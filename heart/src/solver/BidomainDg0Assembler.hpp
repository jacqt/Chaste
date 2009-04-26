/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef _BIDOMAINDG0ASSEMBLER_HPP_
#define _BIDOMAINDG0ASSEMBLER_HPP_

#include "UblasIncludes.hpp"

//#include <iostream>
#include <vector>
#include <petscvec.h>

#include "AbstractDynamicAssemblerMixin.hpp"
#include "AbstractLinearAssembler.hpp"

#include "BidomainPde.hpp"
#include "HeartConfig.hpp"
#include "Element.hpp"
#include "BoundaryElement.hpp"
#include "ChastePoint.hpp"
#include "AbstractMesh.hpp"
#include "BoundaryConditionsContainer.hpp"

/**
 *  BidomainDg0Assembler
 *
 *  The 2 unknowns are voltage and extracellular potential.
 *
 *  This assembler interpolates quantities such as ionic currents and stimuli from
 *  their nodal values (obtained from a BidomainPde) onto a gauss point, and uses
 *  the interpolated values in assembly. The assembler also creates boundary conditions,
 *  which are zero-Neumann boundary conditions on the surface unless
 *  SetFixedExtracellularPotentialNodes() is called.
 *
 *  The user should call Solve() from the superclass AbstractDynamicAssemblerMixin.
 *
 *  NOTE: if any cells have a non-zero extracellular stimulus, phi_e must be fixed at some
 *  nodes (using SetFixedExtracellularPotentialNodes() ), else no solution is possible.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BidomainDg0Assembler
    : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 2, false, BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> >,
      public AbstractDynamicAssemblerMixin<ELEMENT_DIM, SPACE_DIM, 2>
{
public:
    static const unsigned E_DIM = ELEMENT_DIM; /**< The element dimension (to save typing). */
    static const unsigned S_DIM = SPACE_DIM; /**< The space dimension (to save typing). */
    static const unsigned P_DIM = 2u; /**< The problem dimension (to save typing). */

protected:
    // Save typing
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> SelfType; /**< This type (to save typing). */
    typedef AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 2, false, SelfType> BaseClassType; /**< Base class type (to save typing). */

    /// Allow the AbstractStaticAssembler to call our private/protected methods using static polymorphism.
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 2, false, SelfType>;

    /** The PDE to be solved. */
    BidomainPde<SPACE_DIM>* mpBidomainPde;

    HeartConfig* mpConfig;

    // quantities to be interpolated
    double mIionic;
    double mIIntracellularStimulus;
    double mIExtracellularStimulus;

    bool mNullSpaceCreated;

    Vec mExternalVoltageMask;
    std::vector<unsigned> mFixedExtracellularPotentialNodes;

    unsigned mRowForAverageOfPhiZeroed;

    /**
     * Overridden ResetInterpolatedQuantities() method.
     */
    void ResetInterpolatedQuantities();

    /**
     * Create the linear system object if it hasn't been already.
     * Can use an initial solution as PETSc template, or base it on the mesh size.
     * 
     * @param initialSolution an initial guess
     */
    void InitialiseForSolve(Vec initialSolution);

    /**
     * Overridden IncrementInterpolatedQuantities() method.
     *
     * @param phi_i \todo This should be phiI
     * @param pNode
     */
    void IncrementInterpolatedQuantities(double phi_i, const Node<SPACE_DIM>* pNode);

    /**
     * ComputeMatrixTerm()
     *
     * This method is called by AssembleOnElement() and tells the assembler
     * the contribution to add to the element stiffness matrix.
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i \todo should this be rU?
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_matrix<double,2*(ELEMENT_DIM+1),2*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &rU,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     *  ComputeVectorTerm()
     *
     *  This method is called by AssembleOnElement() and tells the assembler
     *  the contribution to add to the element stiffness vector.
     * 
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rGradPhi Basis gradients, rGradPhi(i,j) = d(phi_j)/d(X_i)
     * @param rX The point in space
     * @param u The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     * @param pElement Pointer to the element
     */
    virtual c_vector<double,2*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,2> &u,
        c_matrix<double, 2, SPACE_DIM> &rGradU /* not used */,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement);

    /**
     * ComputeVectorSurfaceTerm()
     *
     * This method is called by AssembleOnSurfaceElement() and tells the
     * assembler what to add to the element stiffness matrix arising
     * from surface element contributions.
     *
     * NOTE: this method has to be implemented but shouldn't ever be called -
     * because all bidomain problems (currently) just have zero Neumann boundary
     * conditions and the AbstractLinearAssmebler::AssembleSystem() method
     * will realise this and not loop over surface elements.
     *
     * @param rSurfaceElement the element which is being considered.
     * @param rPhi The basis functions, rPhi(i) = phi_i, i=1..numBases
     * @param rX The point in space
     */
#define COVERAGE_IGNORE //see NOTE above
    virtual c_vector<double, 2*ELEMENT_DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
        c_vector<double,ELEMENT_DIM> &rPhi,
        ChastePoint<SPACE_DIM> &rX);
#undef COVERAGE_IGNORE

    /**
     *  PrepareForAssembleSystem
     *
     *  Called at the beginning of AbstractLinearAssembler::AssembleSystem()
     *  after the system. Here, used to integrate cell odes.
     */
    virtual void PrepareForAssembleSystem(Vec currentSolution, double time);

    /**
     *  FinaliseAssembleSystem
     *
     *  Called by AbstractLinearAssmebler::AssembleSystem() after the system
     *  has been assembled. Here, used to avoid problems with phi_e drifting
     *  by one of 3 methods: pinning nodes, using a null space, or using an
     *  "average phi_e = 0" row.
     */
    virtual void FinaliseAssembleSystem(Vec currentSolution, double currentTime);

public:

    /**
     * Constructor stores the mesh and pde and sets up boundary conditions.
     * 
     * @param pMesh pointer to the mesh
     * @param pPde pointer to the PDE
     * @param pBoundaryConditions pointer to the boundary conditions
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    BidomainDg0Assembler(AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                         BidomainPde<SPACE_DIM>* pPde,
                         BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 2>* pBcc,
                         unsigned numQuadPoints = 2);

    /**
     * Destructor.
     */
    ~BidomainDg0Assembler();

    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to
     *  zero. This does not necessarily have to be called. If it is not, phi_e
     *  is only defined up to a constant.
     *
     *  @param the nodes to be fixed.
     *
     *  NOTE: currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> fixedExtracellularPotentialNodes);

    void SetRowForAverageOfPhiZeroed(unsigned rowMeanPhiEZero);
};

/**
 * Specialize AssemblerTraits since we define interpolation methods as well as
 * Compute*Term methods.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct AssemblerTraits<BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> >
{
    /** The class in which ComputeVectorTerm is defined. */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CVT_CLS;
    /** The class in which ComputeMatrixTerm is defined. */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> CMT_CLS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined. */
    typedef BidomainDg0Assembler<ELEMENT_DIM, SPACE_DIM> INTERPOLATE_CLS;
};

#endif /*_BIDOMAINDG0ASSEMBLER_HPP_*/
