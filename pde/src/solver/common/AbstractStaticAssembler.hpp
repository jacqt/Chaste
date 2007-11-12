#ifndef _ABSTRACTSTATICASSEMBLER_HPP_
#define _ABSTRACTSTATICASSEMBLER_HPP_

#include "AbstractAssembler.hpp"
#include "LinearBasisFunction.cpp"
#include "GaussianQuadratureRule.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "LinearSystem.hpp"
#include "GaussianQuadratureRule.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "EventHandler.hpp"
#include <iostream>

/**
 *  AbstractStaticAssembler
 *
 *  Implmentation of main assembler methods so that the virtual base class 
 *  does not need to contain data (which should improve performance).
 * 
 *  Templated over the PROBLEM_DIM so also handles problems with more than one
 *  unknown variable (ie those of the form u_xx + v = 0, v_xx + 2u = 1, where
 *  PROBLEM_DIM is equal to 2)
 *
 *  It defines default code for AssembleSystem, AssembleOnElement and 
 *  AssembleOnSurfaceElement. Each of these work
 *  for any PROBLEM_DIM>=1. Each of these methods work in both the
 *  dynamic case (when there is a current solution available) and the static
 *  case. The same code is used for the nonlinear and linear cases
 *
 *  See also documentation for AbstractAssembler
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractStaticAssembler : virtual public AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>
{
protected:
    /** Mesh to be solved on */
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /** Quadrature rule for use on normal elements */
    GaussianQuadratureRule<ELEMENT_DIM> *mpQuadRule;
    /** Quadrature rule for use on boundary elements */
    GaussianQuadratureRule<ELEMENT_DIM-1> *mpSurfaceQuadRule;
    

    /** Basis function for use with normal elements */
    typedef LinearBasisFunction<ELEMENT_DIM> BasisFunction;
    /** Basis function for use with boundary elements */
    typedef LinearBasisFunction<ELEMENT_DIM-1> SurfaceBasisFunction;

    /**
     *  The CURRENT SOLUTION as a replicated vector for linear dynamic problems. 
     *  (Empty for a static problem).  The CURRENT GUESS for nonlinear problems.
     */
    ReplicatableVector mCurrentSolutionOrGuessReplicated;
    
    /**
     *  The linear system that is assembled in linear pde problems. Not used in
     *  nonlinear problems
     */
    LinearSystem *mpLinearSystem;


    /**
     *  Calculate the contribution, in a PETSc-friendly format, of a single element. It
     *      contains a call to AssembleOnElement 
     * 
     */    
    void AssembleOnElementPetscFormat( Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                       PetscInt* m, PetscInt idxm[], 
                                       PetscInt* n, PetscInt idxn[],
                                       double v[],
                                       PetscInt* m_rhs, PetscInt idx_rhs[], 
                                       double rhs[],
                                       bool assembleVector, bool assembleMatrix,
                                       unsigned num_elem_nodes)
    {
        *m=*n=*m_rhs=0;
        
//        //Assume all elements have the same number of nodes...
//        const unsigned num_elem_nodes = (*iter)->GetNumNodes();
        
        c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> b_elem;
        
        // Call original AssembleOnElement Method        
        AssembleOnElement(rElement, a_elem, b_elem, assembleVector, assembleMatrix);

        PetscInt RangeLo, RangeHi;
        mpLinearSystem->GetOwnershipRange(RangeLo, RangeHi);
                        
        //unsigned n_elem = 0u;
        for (unsigned i=0; i<num_elem_nodes; i++)
        {
            unsigned node1 = rElement.GetNodeGlobalIndex(i);
                    
            if (assembleMatrix)
            {
                for (unsigned j=0; j<num_elem_nodes; j++)
                {
                    unsigned node2 = rElement.GetNodeGlobalIndex(j);
             
                    for (unsigned k=0; k<PROBLEM_DIM; k++)
                    {
                        PetscInt row = PROBLEM_DIM*node1+k;
                        if (row >= RangeLo && row < RangeHi)
                        {
                            if (j==0)
                            {
                                idxm[*m] = row;
                                *m = (*m) + 1;                                                            
                            }
                            for (unsigned l=0; l<PROBLEM_DIM; l++)
                            {
                                if (i==0 && k==0)
                                {
                                    idxn[*n] = PROBLEM_DIM*node2+l;
                                    *n = (*n) + 1;                                                            
                                }

                                unsigned small_mat_row = PROBLEM_DIM*i+k;
                                unsigned small_mat_col = PROBLEM_DIM*j+l;
                                unsigned small_mat_num_cols = PROBLEM_DIM*(ELEMENT_DIM+1);        
                                v[small_mat_row*small_mat_num_cols + small_mat_col] = a_elem(small_mat_row, small_mat_col);
  
                                //std::cout << PROBLEM_DIM*node1+k << " " << PROBLEM_DIM*node2+l << " : " << a_elem(PROBLEM_DIM*i+k,PROBLEM_DIM*j+l) << std::endl;
                            //mpLinearSystem->AddToMatrixElement( PROBLEM_DIM*node1+k,
                            //                                PROBLEM_DIM*node2+l,
                            //                                a_elem(PROBLEM_DIM*i+k,PROBLEM_DIM*j+l) );
                            }
                        }
                    }
                }
            }
        
            if (assembleVector)
            {
                for (unsigned k=0; k<PROBLEM_DIM; k++)
                {
                    PetscInt row = PROBLEM_DIM*node1+k;
                    if (row >= RangeLo && row < RangeHi)
                    {                   
                      idx_rhs[*m_rhs] = PROBLEM_DIM*node1+k;
                      rhs[*m_rhs]       = b_elem(PROBLEM_DIM*i+k);
                      *m_rhs = (*m_rhs) + 1;
                    }
                    //mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node1+k,b_elem(PROBLEM_DIM*i+k));
                }
            }
        }

    }
          
    /**
     *  Calculate the contribution of a single element to the linear system.
     * 
     *  @param rElement The element to assemble on.
     *  @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     *  @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     *  @param currentSolutionOrGuess For the parabolic linear case, the solution at the current 
     *     timestep. NULL for the static linear case. In the nonlinear case, the current
     *     guess.
     *  @param assembleVector a bool stating whether to assemble the load vector (in the 
     *     linear case) or the residual vector (in the nonlinear case)
     *  @param assembleMatrix a bool stating whether to assemble the stiffness matrix (in 
     *     the linear case) or the Jacobian matrix (in the nonlinear case)
     * 
     *  Called by AssembleSystem()
     *  Calls ComputeMatrixTerm() etc
     */
    virtual void AssembleOnElement( Element<ELEMENT_DIM,SPACE_DIM> &rElement,
                                    c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1) > &rAElem,
                                    c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> &rBElem,
                                    bool assembleVector,
                                    bool assembleMatrix)
    {
        GaussianQuadratureRule<ELEMENT_DIM> &quad_rule =
            *(AbstractStaticAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::mpQuadRule);
            
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function.
         */
        const c_matrix<double, SPACE_DIM, SPACE_DIM> *p_inverse_jacobian = NULL;
        double jacobian_determinant = rElement.GetJacobianDeterminant();
        
        // Initialise element contributions to zero
        if ( assembleMatrix || this->ProblemIsNonlinear() ) // don't need to construct grad_phi or grad_u in that case
        {
            p_inverse_jacobian = rElement.GetInverseJacobian();
        }
        
        if (assembleMatrix)
        {
            rAElem.clear();
        }
        
        if (assembleVector)
        {
            rBElem.clear();
        }
        
        const unsigned num_nodes = rElement.GetNumNodes();
        
        // loop over Gauss points
        for (unsigned quad_index=0; quad_index < quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const ChastePoint<ELEMENT_DIM>& quad_point = quad_rule.rGetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM+1> phi = BasisFunction::ComputeBasisFunctions(quad_point);
            c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> grad_phi;
            
            if ( assembleMatrix || this->ProblemIsNonlinear() )
            {
                grad_phi = BasisFunction::ComputeTransformedBasisFunctionDerivatives
                           (quad_point, *p_inverse_jacobian);
            }
            
            // Location of the gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            ChastePoint<SPACE_DIM> x(0,0,0);
            
            c_vector<double,PROBLEM_DIM> u = zero_vector<double>(PROBLEM_DIM);
            c_matrix<double,PROBLEM_DIM,SPACE_DIM> grad_u = zero_matrix<double>(PROBLEM_DIM,SPACE_DIM);
            
            // allow the concrete version of the assembler to interpolate any
            // desired quantities
            this->ResetInterpolatedQuantities();
            
            
            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            for (unsigned i=0; i<num_nodes; i++)
            {
                const Node<SPACE_DIM> *p_node = rElement.GetNode(i);
                const c_vector<double, SPACE_DIM> node_loc = p_node->rGetLocation();
                
                // interpolate x
                x.rGetLocation() += phi(i)*node_loc;
                
                // interpolate u and grad u if a current solution or guess exists
                unsigned node_global_index = rElement.GetNodeGlobalIndex(i);
                if (mCurrentSolutionOrGuessReplicated.size()>0)
                {
                    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
                    {
                        // If we have a current solution (e.g. this is a dynamic problem)
                        // get the value in a usable form.
                        
                        // NOTE - currentSolutionOrGuess input is actually now redundant at this point -
                        
                        // NOTE - following assumes that, if say there are two unknowns u and v, they
                        // are stored in the current solution vector as
                        // [U1 V1 U2 V2 ... U_n V_n]
                        double u_at_node=GetCurrentSolutionOrGuessValue(node_global_index, index_of_unknown);
                        u(index_of_unknown) += phi(i)*u_at_node;
                        
                        if ( this->ProblemIsNonlinear() ) // don't need to construct grad_phi or grad_u in that case
                        {
                            for (unsigned j=0; j<SPACE_DIM; j++)
                            {
                                grad_u(index_of_unknown,j) += grad_phi(j,i)*u_at_node;
             
                            }
                        }
                    }
                }
                
                // allow the concrete version of the assembler to interpolate any
                // desired quantities
                IncrementInterpolatedQuantities(phi(i), p_node);
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            
            ////////////////////////////////////////////////////////////
            // create rAElem and rBElem
            ////////////////////////////////////////////////////////////
            if (assembleMatrix)
            {
                noalias(rAElem) += ComputeMatrixTerm(phi, grad_phi, x, u, grad_u) * wJ;
            }
            
            if (assembleVector)
            {
                noalias(rBElem) += ComputeVectorTerm(phi, grad_phi, x, u, grad_u) * wJ;
            }
        }
    }
    
    
    
    /**
     * Calculate the contribution of a single surface element with Neumann
     * boundary condition to the linear system.
     * 
     * @param rSurfaceElement The element to assemble on.
     * @param rBSurfElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     */
    virtual void AssembleOnSurfaceElement(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
                                          c_vector<double, PROBLEM_DIM*ELEMENT_DIM> &rBSurfElem)
    {
        GaussianQuadratureRule<ELEMENT_DIM-1> &quad_rule =
            *(AbstractStaticAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::mpSurfaceQuadRule);
            
        double jacobian_determinant = rSurfaceElement.GetJacobianDeterminant();
        
        rBSurfElem.clear();
        
        // loop over Gauss points
        for (unsigned quad_index=0; quad_index<quad_rule.GetNumQuadPoints(); quad_index++)
        {
            const ChastePoint<ELEMENT_DIM-1>& quad_point = quad_rule.rGetQuadPoint(quad_index);
            
            c_vector<double, ELEMENT_DIM>  phi = SurfaceBasisFunction::ComputeBasisFunctions(quad_point);
            
            
            /////////////////////////////////////////////////////////////
            // interpolation
            /////////////////////////////////////////////////////////////
            
            // Location of the gauss point in the original element will be
            // stored in x
            ChastePoint<SPACE_DIM> x(0,0,0);
            
            this->ResetInterpolatedQuantities();
            for (unsigned i=0; i<rSurfaceElement.GetNumNodes(); i++)
            {
                const c_vector<double, SPACE_DIM> node_loc = rSurfaceElement.GetNode(i)->rGetLocation();
                x.rGetLocation() += phi(i)*node_loc;
                
                // allow the concrete version of the assembler to interpolate any
                // desired quantities
                IncrementInterpolatedQuantities(phi(i), rSurfaceElement.GetNode(i));
                
                ///\todo: add interpolation of u as well
            }
            
            double wJ = jacobian_determinant * quad_rule.GetWeight(quad_index);
            
            ////////////////////////////////////////////////////////////
            // create rAElem and rBElem
            ////////////////////////////////////////////////////////////
            ///\todo Improve efficiency of Neumann BC implementation.
            noalias(rBSurfElem) += ComputeVectorSurfaceTerm(rSurfaceElement, phi, x) * wJ;
        }
    }
    
    
    
    /**
     *  AssembleSystem - the major method for all assemblers
     * 
     *  Assemble the linear system for a linear PDE, or the residual or Jacobian for
     *  nonlinear PDEs. Loops over each element (and each each surface element if 
     *  there are non-zero Neumann boundary conditions), calls AssembleOnElement() 
     *  and adds the contribution to the linear system.
     * 
     *  Takes in current solution and time if necessary but only used if the problem 
     *  is a dynamic one. This method uses PROBLEM_DIM and can assemble linear systems 
     *  for any number of unknown variables.
     * 
     *  Called by Solve()
     *  Calls AssembleOnElement()
     * 
     *  @param assembleVector  Whether to assemble the RHS vector of the linear system
     *     (i.e. the residual vector for nonlinear problems).
     *  @param assembleMatrix  Whether to assemble the LHS matrix of the linear system
     *     (i.e. the jacobian matrix for nonlinear problems).
     * 
     *  @param currentSolutionOrGuess The current solution in a linear dynamic problem, 
     *     or the current guess in a nonlinear problem. Should be NULL for linear static 
     *     problems.
     * 
     *  @param currentTime The current time for dynamic problems. Not used in static 
     *     problems.
     */
    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix,
                                Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        EventType assemble_event;
        if(assembleMatrix)
        {
            assemble_event = ASSEMBLE_SYSTEM;
            //std::cout << "*";
        }
        else
        {
            assemble_event = ASSEMBLE_RHS;
        }
        
        // Check we've actually been asked to do something!
        assert(assembleVector || assembleMatrix);
        
        // Check the linear system object has been set up correctly
        assert(mpLinearSystem != NULL);
        assert(mpLinearSystem->GetSize() == PROBLEM_DIM * this->mpMesh->GetNumNodes());
        assert(!assembleVector || mpLinearSystem->rGetRhsVector() != NULL);
        assert(!assembleMatrix || mpLinearSystem->rGetLhsMatrix() != NULL);
                   

        // Replicate the current solution and store so can be used in
        // AssembleOnElement
        if (currentSolutionOrGuess != NULL)
        {
            EventHandler::BeginEvent(COMMUNICATION);
            this->mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolutionOrGuess);
            EventHandler::EndEvent(COMMUNICATION);
        }

        EventHandler::BeginEvent(assemble_event);

        // the AssembleOnElement type methods will determine if a current solution or
        // current guess exists by looking at the size of the replicated vector, so
        // check the size is zero if there isn't a current solution
        assert(    ( currentSolutionOrGuess && mCurrentSolutionOrGuessReplicated.size()>0)
                || ( !currentSolutionOrGuess && mCurrentSolutionOrGuessReplicated.size()==0));


        // the concrete class can override this following method if there is
        // work to be done before assembly
        this->PrepareForAssembleSystem(currentSolutionOrGuess, currentTime);

        // Zero the matrix/vector if it is to be assembled
        if (assembleVector)
        {
            mpLinearSystem->ZeroRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->ZeroLhsMatrix();
        }
        
        // Get an iterator over the elements of the mesh
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator
            iter = this->mpMesh->GetElementIteratorBegin();
              
        c_matrix<double, PROBLEM_DIM*(ELEMENT_DIM+1), PROBLEM_DIM*(ELEMENT_DIM+1)> a_elem;
        c_vector<double, PROBLEM_DIM*(ELEMENT_DIM+1)> b_elem;
        
        
        ////////////////////////////////////////////////////////
        // loop over elements
        ////////////////////////////////////////////////////////
        while (iter != this->mpMesh->GetElementIteratorEnd())
        {
            Element<ELEMENT_DIM, SPACE_DIM>& element = **iter;
            
            if (element.GetOwnership() == true)
            {
                AssembleOnElement(element, a_elem, b_elem, assembleVector, assembleMatrix);
                    
                unsigned p_indices[PROBLEM_DIM*(ELEMENT_DIM+1)];                
                element.GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, p_indices);
                    
                if (assembleMatrix)
                {
                    mpLinearSystem->AddLhsMultipleValues(p_indices, a_elem);
                }
                
                if (assembleVector)
                {
                    mpLinearSystem->AddRhsMultipleValues(p_indices, b_elem);
//                    for (unsigned i=0; i<num_elem_nodes; i++)
//                    {
//                        unsigned node1 = element.GetNodeGlobalIndex(i);
//                        for (unsigned k=0; k<PROBLEM_DIM; k++)
//                        {
//                            mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node1+k,b_elem(PROBLEM_DIM*i+k));
//                        }
//                    }
                }
            }
                
            iter++;
        }
        
        // add the integrals associated with Neumann boundary conditions to the linear system
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator
            surf_iter = this->mpMesh->GetBoundaryElementIteratorBegin();
        
        
        ////////////////////////////////////////////////////////
        // loop over surface elements
        ////////////////////////////////////////////////////////
        
        // note, the following condition is not true of Bidomain or Monodomain
        if (this->mpBoundaryConditions->AnyNonZeroNeumannConditions()==true && assembleVector)
        {
            if (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
            {
                const unsigned num_surf_nodes = (*surf_iter)->GetNumNodes();
                c_vector<double, PROBLEM_DIM*ELEMENT_DIM> b_surf_elem;
                
                while (surf_iter != this->mpMesh->GetBoundaryElementIteratorEnd())
                {
                    const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>& surf_element = **surf_iter;
                    
                    ///\todo Check surf_element is in the Neumann surface in an efficient manner
                    /// e.g. by iterating over boundary conditions!
                    if (this->mpBoundaryConditions->HasNeumannBoundaryCondition(&surf_element))
                    {
                        AssembleOnSurfaceElement(surf_element, b_surf_elem);
                        
                        for (unsigned i=0; i<num_surf_nodes; i++)
                        {
                            unsigned node_index = surf_element.GetNodeGlobalIndex(i);
                            
                            for (unsigned k=0; k<PROBLEM_DIM; k++)
                            {
                                mpLinearSystem->AddToRhsVectorElement(PROBLEM_DIM*node_index + k, b_surf_elem(PROBLEM_DIM*i+k));
                            }
                        }
                    }
                    surf_iter++;
                }
            }
        }
        
        if (assembleVector)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->AssembleIntermediateLhsMatrix();
        }
        
        // Apply dirichlet boundary conditions
        this->ApplyDirichletConditions(currentSolutionOrGuess, assembleMatrix);
        
        if (assembleVector)
        {
            mpLinearSystem->AssembleRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->AssembleFinalLhsMatrix();
        }
        
        // overload this method if the assembler has to do anything else
        // required (like setting up a null basis (see BidomainDg0Assembler))
        this->FinaliseAssembleSystem(currentSolutionOrGuess, currentTime);
        
        EventHandler::EndEvent(assemble_event);
    }
    
    
    
    /**
     *  This method is called at the beginning of Solve(). Subclass assemblers can 
     *  use it to check everything has been set up correctly
     */
    virtual void PrepareForSolve()
    {
        assert(mpMesh != NULL);
// commented out because FlaggedMeshAssembler has it's own FlaggedMeshBcc.. - design issue
//        assert(this->mpBoundaryConditions != NULL);
        
        //Set the elements' ownerships according to the node ownership
        DistributedVector::SetProblemSize(this->mpMesh->GetNumNodes());
        this->mpMesh->SetElementOwnerships(DistributedVector::Begin().Global,
                                           DistributedVector::End().Global);
    }
    
    
    /**
     * Accessor method that subclasses of AbstractAssembler (but not us)
     * can use to get to useful data.
     */
    ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rGetMesh()
    {
        assert(mpMesh!=NULL);
        return *mpMesh;
    }

    /**
     * Accessor method that subclasses of AbstractAssembler (but not us)
     * can use to get to useful data.
     */
    LinearSystem** GetLinearSystem()
    {
        return &mpLinearSystem;
    }

    /**
     * Accessor method that subclasses of AbstractAssembler (but not us)
     * can use to get to useful data.
     */
    ReplicatableVector& rGetCurrentSolutionOrGuess()
    {
        return mCurrentSolutionOrGuessReplicated;
    }
    
    
    virtual double GetCurrentSolutionOrGuessValue(unsigned nodeIndex, unsigned indexOfUnknown)
    {
        return mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*nodeIndex + indexOfUnknown];
    }
public:
    /**
     * Default constructor. Uses linear basis functions.
     * 
     * @param numQuadPoints Number of quadrature points to use per dimension.
     */
    AbstractStaticAssembler(unsigned numQuadPoints = 2)
        : AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>()
    {
        // Initialise mesh and bcs to null, so we can check they
        // have been set before attempting to solve
        mpMesh = NULL;
        
        mpQuadRule = NULL;
        mpSurfaceQuadRule = NULL;
        SetNumberOfQuadraturePointsPerDimension(numQuadPoints);
        
        mpLinearSystem = NULL;
    }
    
    /**
     * Set the number of quadrature points to use, per dimension.
     * 
     * This method will throw an exception if the requested number of quadrature
     * points is not supported. (///\todo: There may be a small memory leak if this
     * occurs.)
     * 
     * @param numQuadPoints Number of quadrature points to use per dimension.
     */
    void SetNumberOfQuadraturePointsPerDimension(unsigned numQuadPoints)
    {
        delete mpQuadRule;
        mpQuadRule = new GaussianQuadratureRule<ELEMENT_DIM>(numQuadPoints);
        delete mpSurfaceQuadRule;
        mpSurfaceQuadRule = new GaussianQuadratureRule<ELEMENT_DIM-1>(numQuadPoints);
    }
    
    
    /**
     * Set the mesh.
     */
    void SetMesh(ConformingTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    {
        mpMesh = pMesh;
    }
    
    
    /**
     * Delete any memory allocated by this class.
     */
    virtual ~AbstractStaticAssembler()
    {
        delete mpQuadRule;
        delete mpSurfaceQuadRule;
        delete mpLinearSystem;
    }
};

#endif //_ABSTRACTSTATICASSEMBLER_HPP_
