/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_
#define ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_

#include "NonlinearElasticityAssembler.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "QuadraticBasisFunction.hpp" // not included in NonlinearElasticityAssembler.hpp, just the cpp
#include "LinearBasisFunction.hpp"
#include "AbstractContractionModel.hpp"

/**
 *  AbstractCardiacMechanicsAssembler
 *
 *  Base class to implicit and explicit cardiac mechanics assemblers. Inherits from NonlinearElasticityAssembler
 *  Main method is the overloaded AssembleOnElement which does the extra work needed for cardiac problems. The
 *  child classes hold the contraction models and need to implement a method for getting the active tension from
 *  the model.
 */
template<unsigned DIM>
class AbstractCardiacMechanicsAssembler : public NonlinearElasticityAssembler<DIM>
{
protected:
    static const unsigned STENCIL_SIZE = NonlinearElasticityAssembler<DIM>::STENCIL_SIZE;
    static const unsigned NUM_NODES_PER_ELEMENT = NonlinearElasticityAssembler<DIM>::NUM_NODES_PER_ELEMENT;
    static const unsigned NUM_VERTICES_PER_ELEMENT = NonlinearElasticityAssembler<DIM>::NUM_VERTICES_PER_ELEMENT;

    /**
     *  Vector of contraction model (pointers). One for each quadrature point.
     *  Note the indexing: the i-th entry corresponds to the i-th global quad
     *  point, when looping over elements and then quad points
     */
    std::vector<AbstractContractionModel*> mContractionModelSystems;

    /**
     *  Stored stretches (in fibre direction, at each quadrature point). Should be stored
     *  when GetActiveTensionAndTensionDerivs() is called, and can be used either in that
     *  timestep (implicit solver), or the next timestep (explicit solver)
     */
    std::vector<double> mStretches;

    /** Total number of quad points in the (mechanics) mesh */
    unsigned mTotalQuadPoints;

    /** Whether the material law was passed in or the default used */
    bool mAllocatedMaterialLawMemory;

    /** Current time */
    double mCurrentTime;
    /** Time to which the solver has been asked to solve to */
    double mNextTime;
    /** Time used to integrate the contraction model */
    double mOdeTimestep;

    /** The fibre-sheet-normal directions (in a matrix), if constant (defaults to the identity, ie fibres in the X-direction, sheet in the XY plane) */
    c_matrix<double,DIM,DIM> mConstantFibreSheetDirections;

    /**
     * The fibre-sheet-normal directions (matrices), one for each element. Only non-NULL if SetVariableFibreSheetDirections()
     * is called, if not mConstantFibreSheetDirections is used instead
     */
    std::vector<c_matrix<double,DIM,DIM> >* mpVariableFibreSheetDirections;

    /** (Pointer to) the fibre-sheet matrix for the current element being assembled on */    
    c_matrix<double,DIM,DIM>* mpCurrentElementFibreSheetMatrix;
    /** The fibre direction for the current element being assembled on */    
    c_vector<double,DIM> mCurrentElementFibreDirection;


    /** Check whether the given matrix is orthogonal, by computing A^T A and verifying whether each
      * component matches the identity matrix, to a given tolerance.
      * @param rMatrix reference to the matrix being tested
      */
    void CheckOrthogonality(c_matrix<double,DIM,DIM>& rMatrix);


    /**
     *  Whether the solver is implicit or not (ie whether the contraction model depends on lambda (and depends on
     *  lambda at the current time)). For whether dTa_dLam dependent terms need to be added to the Jacbobian
     */
    virtual bool IsImplicitSolver()=0;

    /**
     * Overloaded AssembleOnElement. Does an bit of initial set up (sets mpCurrentElementFibreSheetMatrix
     * and mCurrentElementFibreDirection), then calls the base class version on NonlinearElasticityAssembler.
     * That then calls the overloaded ComputeStressAndStressDerivative, which
     * computes the passive stress as normal but also computes the extra, active stresses to be
     * added to the total stress.
     * 
     * @param rElement The element to assemble on.
     * @param rAElem The element's contribution to the LHS matrix is returned in this
     *     n by n matrix, where n is the no. of nodes in this element. There is no
     *     need to zero this matrix before calling.
     * @param rAElemPrecond The element's contribution to the matrix passed to PetSC
     *     in creating a preconditioner
     * @param rBElem The element's contribution to the RHS vector is returned in this
     *     vector of length n, the no. of nodes in this element. There is no
     *     need to zero this vector before calling.
     * @param assembleResidual A bool stating whether to assemble the residual vector.
     * @param assembleJacobian A bool stating whether to assemble the Jacobian matrix.
     */
    void AssembleOnElement(Element<DIM, DIM>& rElement,
                           c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElem,
                           c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElemPrecond,
                           c_vector<double,STENCIL_SIZE>& rBElem,
                           bool assembleResidual,
                           bool assembleJacobian);


    /**
     *  Overloaded ComputeStressAndStressDerivative(), which computes the passive part of the
     *  stress as normal but also calls on the contraction model to get the active stress and 
     *  adds it on.
     * 
     *  @param pMaterialLaw The material law for this element
     *  @param rC The Lagrangian deformation tensor (F^T F)
     *  @param rInvC The inverse of C. Should be computed by the user.
     *  @param pressure The current pressure
     *  @param currentQuadPointGlobalIndex The index (assuming an outer loop over elements and an inner 
     *    loop over quadrature points), of the current quadrature point.
     *  @param rT The stress will be returned in this parameter
     *  @param rDTdE the stress derivative will be returned in this parameter, assuming
     *    the final parameter is true
     *  @param computeDTdE A boolean flag saying whether the stress derivative is
     *    required or not.
     */ 
    void ComputeStressAndStressDerivative(AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                          c_matrix<double,DIM,DIM>& rC, 
                                          c_matrix<double,DIM,DIM>& rInvC, 
                                          double pressure, 
                                          unsigned currentQuadPointGlobalIndex,
                                          c_matrix<double,DIM,DIM>& rT,
                                          FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                          bool computeDTdE);
    



    /**
     *  Pure method called in AbstractCardiacMechanicsAssembler::ComputeStressAndStressDerivative(), which needs to provide
     *  the active tension (and other info if implicit (if the contraction model depends on stretch
     *  or stretch rate)) at a particular quadrature point. Takes in the current fibre stretch.
     *
     *  @param currentFibreStretch The stretch in the fibre direction
     *  @param currentQuadPointGlobalIndex quadrature point the integrand is currently being evaluated at in AssembleOnElement
     *  @param assembleJacobian  A bool stating whether to assemble the Jacobian matrix.
     *  @param rActiveTension The returned active tension
     *  @param rDerivActiveTensionWrtLambda The returned dT_dLam, derivative of active tension wrt stretch. Only should be set in implicit assemblers
     *  @param rDerivActiveTensionWrtDLambdaDt The returned dT_dLamDot, derivative of active tension wrt stretch rate.  Only should be set in implicit assemblers
     */
    virtual void GetActiveTensionAndTensionDerivs(double currentFibreStretch,
                                                  unsigned currentQuadPointGlobalIndex,
                                                  bool assembleJacobian,
                                                  double& rActiveTension,
                                                  double& rDerivActiveTensionWrtLambda,
                                                  double& rDerivActiveTensionWrtDLambdaDt)=0;
public:
    /**
     * Constructor
     *
     * @param pQuadMesh A pointer to the mesh.
     * @param outputDirectory The output directory, relative to TEST_OUTPUT
     * @param rFixedNodes The fixed nodes
     * @param pMaterialLaw The material law for the tissue. If NULL the default
     *   (NashHunterPoleZero) law is used.
     */
    AbstractCardiacMechanicsAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                      std::string outputDirectory,
                                      std::vector<unsigned>& rFixedNodes,
                                      AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw);


    /**
     *  Destructor just deletes memory if it was allocated
     */
    ~AbstractCardiacMechanicsAssembler();
    
    
    /** Get the total number of quad points in the mesh. Pure, implemented in concrete assembler */
    unsigned GetTotalNumQuadPoints()
    {
        return mTotalQuadPoints;
    }

    /** Get the quadrature rule used in the elements. */
    virtual GaussianQuadratureRule<DIM>* GetQuadratureRule()
    {
        return this->mpQuadratureRule;
    }


    /**
     *    Set a constant fibre-sheet-normal direction (a matrix) to something other than the default (fibres in X-direction,
     *  sheet in the XY plane)
     *  @param rFibreSheetMatrix The fibre-sheet-normal matrix (fibre dir the first column, normal-to-fibre-in sheet in second
     *  column, sheet-normal in third column).
     */
    void SetConstantFibreSheetDirections(const c_matrix<double,DIM,DIM>& rFibreSheetMatrix)
    {
        mConstantFibreSheetDirections = rFibreSheetMatrix;
        CheckOrthogonality(mConstantFibreSheetDirections);
    }

    /**
     *    Set a variable fibre-sheet-normal direction (matrices), one for each element, from a file.
     *  The file should be a .ortho file (ie each line has the fibre dir, sheet dir, normal dir for that element).
     *  The number of elements must match the number in the MECHANICS mesh!
     *  @param orthoFile the file containing the fibre/sheet directions
     */
    void SetVariableFibreSheetDirections(std::string orthoFile);



    /**
     *  Set the intracellular Calcium concentrations and voltages at each quad point. Pure.
     *
     *  Implicit solvers (for contraction models which are functions of stretch (and maybe
     *  stretch rate) would integrate the contraction model with this Ca/V/t using the current
     *  stretch (ie inside AssembleOnElement, ie inside GetActiveTensionAndTensionDerivs).
     *  Explicit solvers (for contraction models which are NOT functions of stretch can immediately
     *  integrate the contraction models to get the active tension.
     *
     *  @param rCalciumConcentrations Reference to a vector of intracellular calcium concentrations at each quadrature point
     *  @param rVoltages Reference to a vector of voltages at each quadrature point
     */
    void SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations,
                              std::vector<double>& rVoltages);

    /**
     *  Solve for the deformation, integrating the contraction model ODEs.
     *
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    virtual void Solve(double time, double nextTime, double odeTimestep)=0;
    
    
    
    // following method will eventually become 
    // ComputeDeformationGradientsAndStretchesInEachElement(std::vector<c_matrix<double,DIM,DIM> >& rDeformationGradients, 
    //                                                      std::vector<double>& rStretches)

    /**
     *  Compute the stretch in the fibre direction, for each element in the mesh.
     *  Note: using quadratic interpolation for position, the stretches actually vary linearly 
     *  in each element. However, for computational efficiency reasons, when computing the
     *  stretch to pass back to the cell model, we just assume the stretches are constant in
     *  each element (ie ignoring the quadratic correction to the displacement). This means that
     *  the (const) stretch for each element can be computed in advance and stored, and we don't
     *  have to worry about interpolation onto the precise location of the cell-model (electrics-mesh)
     *  node, just which element it is in. 
     *  
     *  To compute this constant stretch we just have to compute F using the deformed positions at the 
     *  vertices only, with linear bases, rather than all the nodes and quadratic bases. 
     * 
     *  @rStretches A reference of a std::vector in which the stretch in each element will be returned.
     *  Must be allocated prior to being passed in
     */
    void ComputeStretchesInEachElement(std::vector<double>& rStretches);
};


template<unsigned DIM>
AbstractCardiacMechanicsAssembler<DIM>::AbstractCardiacMechanicsAssembler(QuadraticMesh<DIM>* pQuadMesh,
                                                                          std::string outputDirectory,
                                                                          std::vector<unsigned>& rFixedNodes,
                                                                          AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw)
   : NonlinearElasticityAssembler<DIM>(pQuadMesh,
                                       pMaterialLaw!=NULL ? pMaterialLaw : new NashHunterPoleZeroLaw<DIM>,
                                       zero_vector<double>(DIM),
                                       DOUBLE_UNSET,
                                       outputDirectory,
                                       rFixedNodes),
                                       mCurrentTime(DBL_MAX),
                                       mNextTime(DBL_MAX),
                                       mOdeTimestep(DBL_MAX)
{
    // compute total num quad points
    mTotalQuadPoints = pQuadMesh->GetNumElements()*this->mpQuadratureRule->GetNumQuadPoints();

    // note that if pMaterialLaw is NULL a new NashHunter law was sent to the
    // NonlinElas constuctor (see above)
    mAllocatedMaterialLawMemory = (pMaterialLaw==NULL);

    // initialise the store of fibre stretches
    mStretches.resize(mTotalQuadPoints, 1.0);

    // initialise fibre/sheet direction matrix to be the identity, fibres in X-direction, and sheet in XY-plane
    mConstantFibreSheetDirections = zero_matrix<double>(DIM,DIM);
    for(unsigned i=0; i<DIM; i++)
    {
        mConstantFibreSheetDirections(i,i) = 1.0;
    }

    mpVariableFibreSheetDirections = NULL;
}


template<unsigned DIM>
AbstractCardiacMechanicsAssembler<DIM>::~AbstractCardiacMechanicsAssembler()
{
    if(mAllocatedMaterialLawMemory)
    {
        assert(this->mMaterialLaws.size()==1); 
        delete this->mMaterialLaws[0];
    }

    if(mpVariableFibreSheetDirections)
    {
        delete mpVariableFibreSheetDirections;
    }
}



template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations,
                                                                  std::vector<double>& rVoltages)

{
    assert(rCalciumConcentrations.size() == this->mTotalQuadPoints);
    assert(rVoltages.size() == this->mTotalQuadPoints);

    ContractionModelInputParameters input_parameters;

    for(unsigned i=0; i<rCalciumConcentrations.size(); i++)
    {
        input_parameters.intracellularCalciumConcentration = rCalciumConcentrations[i];
        input_parameters.voltage = rVoltages[i];

        mContractionModelSystems[i]->SetInputParameters(input_parameters);
    }
}

template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::ComputeStressAndStressDerivative(AbstractIncompressibleMaterialLaw<DIM>* pMaterialLaw,
                                                                         c_matrix<double,DIM,DIM>& rC, 
                                                                         c_matrix<double,DIM,DIM>& rInvC, 
                                                                         double pressure, 
                                                                         unsigned currentQuadPointGlobalIndex,
                                                                         c_matrix<double,DIM,DIM>& rT,
                                                                         FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                                                         bool assembleJacobian)
{
    // 1. Compute T and dTdE for the PASSIVE part of the strain energy.
    pMaterialLaw->SetChangeOfBasisMatrix(*mpCurrentElementFibreSheetMatrix);
    pMaterialLaw->ComputeStressAndStressDerivative(rC,rInvC,pressure,rT,rDTdE,assembleJacobian);

    // 2. Compute the active tension and add to the stress and stress-derivative
    double I4 = inner_prod(mCurrentElementFibreDirection, prod(rC, mCurrentElementFibreDirection));
    double lambda = sqrt(I4);

    double active_tension = 0;
    double d_act_tension_dlam = 0.0;     // Set and used if assembleJacobian==true
    double d_act_tension_d_dlamdt = 0.0; // Set and used if assembleJacobian==true

    GetActiveTensionAndTensionDerivs(lambda, currentQuadPointGlobalIndex, assembleJacobian,
                                     active_tension, d_act_tension_dlam, d_act_tension_d_dlamdt);

    // amend the stress and dTdE using the active tension
    double dTdE_coeff = -2*active_tension/(I4*I4); // note: I4*I4 = lam^4
    if(IsImplicitSolver())
    {
        double dt = mNextTime-mCurrentTime;
        dTdE_coeff += (d_act_tension_dlam + d_act_tension_d_dlamdt/dt)/(lambda*I4); // note: I4*lam = lam^3
    }

    rT += (active_tension/I4)*outer_prod(mCurrentElementFibreDirection,mCurrentElementFibreDirection);

    if(assembleJacobian)
    {
        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned N=0; N<DIM; N++)
            {
                for (unsigned P=0; P<DIM; P++)
                {
                    for (unsigned Q=0; Q<DIM; Q++)
                    {
                        rDTdE(M,N,P,Q) +=  dTdE_coeff * mCurrentElementFibreDirection(M)
                                                      * mCurrentElementFibreDirection(N)
                                                      * mCurrentElementFibreDirection(P)
                                                      * mCurrentElementFibreDirection(Q);
                    }
                }
            }
        }
    }
}




template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::ComputeStretchesInEachElement(std::vector<double>& rStretches)
{
    assert(rStretches.size()==this->mpQuadMesh->GetNumElements());
   
    static c_matrix<double,DIM,NUM_VERTICES_PER_ELEMENT> element_current_displacements;
    static c_matrix<double,DIM,NUM_VERTICES_PER_ELEMENT> grad_lin_phi;
    static c_matrix<double,DIM,DIM> F;      // the deformation gradient, F = dx/dX, F_{iM} = dx_i/dX_M

    static c_matrix<double,DIM,DIM> jacobian;
    static c_matrix<double,DIM,DIM> inverse_jacobian;
    double jacobian_determinant;
    ChastePoint<DIM> quadrature_point; // not needed, but has to be passed in  
    
    // loop over all the elements
    for(unsigned elem_index=0; elem_index<this->mpQuadMesh->GetNumElements(); elem_index++)
    {
        Element<DIM,DIM>* p_elem = this->mpQuadMesh->GetElement(elem_index);

        // get the fibre direction for this element
        mpCurrentElementFibreSheetMatrix = mpVariableFibreSheetDirections ? &(*mpVariableFibreSheetDirections)[elem_index] : &mConstantFibreSheetDirections;
        for(unsigned i=0; i<DIM; i++)
        {
            mCurrentElementFibreDirection(i) = (*mpCurrentElementFibreSheetMatrix)(i,0);
        }

        // get the displacement at this element
        for (unsigned II=0; II<NUM_VERTICES_PER_ELEMENT; II++)
        {
            for (unsigned JJ=0; JJ<DIM; JJ++)
            {
                element_current_displacements(JJ,II) = this->mCurrentSolution[DIM*p_elem->GetNodeGlobalIndex(II) + JJ];
            }
        }

        // set up the linear basis functions
        this->mpQuadMesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_determinant, inverse_jacobian);
        LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quadrature_point, inverse_jacobian, grad_lin_phi);
        
        F = identity_matrix<double>(DIM,DIM); 

        // loop over the vertices and interpolate F, the deformation gradient
        // (note: could use matrix-mult if this becomes inefficient
        for (unsigned node_index=0; node_index<NUM_VERTICES_PER_ELEMENT; node_index++)
        {
            for (unsigned i=0; i<DIM; i++)
            {
                for (unsigned M=0; M<DIM; M++)
                {
                    F(i,M) += grad_lin_phi(M,node_index)*element_current_displacements(i,node_index);
                }
            }
        }
        
        // compute and save the stretch: m^T C m = ||Fm||
        c_vector<double,DIM> deformed_fibre = prod(F, mCurrentElementFibreDirection);
        rStretches[elem_index] = norm_2(deformed_fibre);
    }
}



template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::AssembleOnElement(Element<DIM, DIM>& rElement,
                                                               c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElem,
                                                               c_matrix<double,STENCIL_SIZE,STENCIL_SIZE>& rAElemPrecond,
                                                               c_vector<double,STENCIL_SIZE>& rBElem,
                                                               bool assembleResidual,
                                                               bool assembleJacobian)
{
    // check these have been set
    assert(mCurrentTime != DBL_MAX);
    assert(mNextTime != DBL_MAX);
    assert(mOdeTimestep != DBL_MAX);

    // set up the fibre info for this element
    mpCurrentElementFibreSheetMatrix = mpVariableFibreSheetDirections ? &(*mpVariableFibreSheetDirections)[rElement.GetIndex()] : &mConstantFibreSheetDirections;
    for(unsigned i=0; i<DIM; i++)
    {
        mCurrentElementFibreDirection(i) = (*mpCurrentElementFibreSheetMatrix)(i,0);
    }

    // call base class AssembleOnElement
    NonlinearElasticityAssembler<DIM>::AssembleOnElement(rElement, rAElem, rAElemPrecond, rBElem, assembleResidual, assembleJacobian);
}


template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::SetVariableFibreSheetDirections(std::string orthoFile)
{
    if((orthoFile.length()<7) || orthoFile.substr(orthoFile.length()-6,orthoFile.length()) != ".ortho")
    {
        EXCEPTION("Fibre file must be a .ortho file");
    }

    std::ifstream ifs(orthoFile.c_str());
    if (!ifs.is_open())
    {
        EXCEPTION("Could not open file: " + orthoFile);
    }

    unsigned num_elem_read_from_file;
    ifs >> num_elem_read_from_file;
    assert(num_elem_read_from_file == this->mpQuadMesh->GetNumElements());

    mpVariableFibreSheetDirections = new std::vector<c_matrix<double,DIM,DIM> >(this->mpQuadMesh->GetNumElements(), zero_matrix<double>(DIM,DIM));
    for(unsigned elem_index=0; elem_index<this->mpQuadMesh->GetNumElements(); elem_index++)
    {
        for(unsigned j=0; j<DIM*DIM; j++)
        {
            double data;
            ifs >> data;
            if(ifs.fail())
            {
                std::stringstream error_message;
                error_message << "Error occurred when reading file " << orthoFile
                              << ". Expected " << this->mpQuadMesh->GetNumElements() << " rows and "
                              << "three (not DIM!) columns, ie three components for each fibre, whichever "
                              << "dimension you are in";
                delete mpVariableFibreSheetDirections;
                EXCEPTION(error_message.str());
            }

            (*mpVariableFibreSheetDirections)[elem_index](j/DIM,j%DIM) = data;
        }
    }

    ifs.close();

    for(unsigned elem_index=0; elem_index<this->mpQuadMesh->GetNumElements(); elem_index++)
    {
        CheckOrthogonality((*mpVariableFibreSheetDirections)[elem_index]);
    }
}


template<unsigned DIM>
void AbstractCardiacMechanicsAssembler<DIM>::CheckOrthogonality(c_matrix<double,DIM,DIM>& rMatrix)
{
    c_matrix<double,DIM,DIM>  temp = prod(trans(rMatrix),rMatrix);
    double tol = 1e-4;
    // check temp is equal to the identity
    for(unsigned i=0; i<DIM; i++)
    {
        for(unsigned j=0; j<DIM; j++)
        {
            double val = (i==j ? 1.0 : 0.0);
            if(fabs(temp(i,j)-val)>tol)
            {
                std::stringstream string_stream;
                string_stream << "The given fibre-sheet matrix, " << rMatrix << ", is not orthogonal"
                   << " (A^T A not equal to I to tolerance " << tol << ")";
                EXCEPTION(string_stream.str());
            }
        }
    }
}

#endif /*ABSTRACTCARDIACMECHANICSASSEMBLER_HPP_*/
