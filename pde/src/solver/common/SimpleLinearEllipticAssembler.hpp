/*

Copyright (C) University of Oxford, 2008

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
#ifndef _SIMPLELINEARELLIPTICASSEMBLER_HPP_
#define _SIMPLELINEARELLIPTICASSEMBLER_HPP_


#include <vector>
#include <petscvec.h>

#include "LinearSystem.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AbstractLinearAssembler.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "GaussianQuadratureRule.hpp"

#include <boost/mpl/void.hpp>
#include <boost/mpl/if.hpp>

/**
 *  SimpleLinearEllipticAssembler
 *
 *  Assembler for solving AbstractLinearEllipticPdes
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class CONCRETE = boost::mpl::void_>
class SimpleLinearEllipticAssembler : public AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, true, SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE> >
{
public:
    static const unsigned E_DIM = ELEMENT_DIM;
    static const unsigned S_DIM = SPACE_DIM;
    static const unsigned P_DIM = 1u;

    friend class TestSimpleLinearEllipticAssembler;

    typedef SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE> SelfType;
    typedef AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, 1, true, SelfType> BaseClassType;
    /// Allow the AbstractStaticAssembler to call our private/protected methods using static polymorphism.
    friend class AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, 1, true, SelfType>;

protected:
    AbstractLinearEllipticPde<SPACE_DIM>* mpEllipticPde;

protected:
    /**
     *  The term to be added to the element stiffness matrix:
     *
     *   grad_phi[row] \dot ( pde_diffusion_term * grad_phi[col])
     */
    virtual c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> pde_diffusion_term = mpEllipticPde->ComputeDiffusionTerm(rX);

        // if statement just saves computing phi*phi^T if it is to be multiplied by zero
        if(mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)!=0)
        {
            return   prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) )
                   - mpEllipticPde->ComputeLinearInUCoeffInSourceTerm(rX,pElement)*outer_prod(rPhi,rPhi);
        }
        else
        {
            return   prod( trans(rGradPhi), c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1>(prod(pde_diffusion_term, rGradPhi)) );
        }
    }

    /**
     *  The term arising from boundary conditions to be added to the element
     *  stiffness vector
     */
    virtual c_vector<double,1*(ELEMENT_DIM+1)> ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1> &rPhi,
        c_matrix<double, ELEMENT_DIM, ELEMENT_DIM+1> &rGradPhi,
        ChastePoint<SPACE_DIM> &rX,
        c_vector<double,1> &u,
        c_matrix<double,1,SPACE_DIM> &rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
    {
        return mpEllipticPde->ComputeConstantInUSourceTerm(rX) * rPhi;
    }



    virtual c_vector<double, ELEMENT_DIM> ComputeVectorSurfaceTerm(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM> &rSurfaceElement,
            c_vector<double, ELEMENT_DIM> &rPhi,
            ChastePoint<SPACE_DIM> &rX )
    {
        // D_times_gradu_dot_n = [D grad(u)].n, D=diffusion matrix
        double D_times_gradu_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX);
        return rPhi * D_times_gradu_dot_n;
    }



public:
    /**
     * Constructor stores the mesh, pde and boundary conditons, and calls base constructor.
     */
    SimpleLinearEllipticAssembler(AbstractMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  AbstractLinearEllipticPde<SPACE_DIM>* pPde,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1>* pBoundaryConditions,
                                  unsigned numQuadPoints = 2) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,1>(),
            BaseClassType(numQuadPoints)
    {
        // note - we don't check any of these are NULL here (that is done in Solve() instead),
        // to allow the user or a subclass to set any of these later
        mpEllipticPde = pPde;
        this->SetMesh(pMesh);
        this->SetBoundaryConditionsContainer(pBoundaryConditions);
    }

    /**
     * This method is called at the beginning of Solve() in AbstractLinearStaticProblemAssembler
     */
    void PrepareForSolve()
    {
        BaseClassType::PrepareForSolve();
        assert(mpEllipticPde != NULL);
    }
};

/** See the AssemblerTraits for SimpleDg0ParabolicAssembler for why this looks like it does */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, class CONCRETE>
struct AssemblerTraits<SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE> >
{
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CVT_CLS>::type
            CVT_CLS;
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     SimpleLinearEllipticAssembler<ELEMENT_DIM, SPACE_DIM, CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CMT_CLS>::type
            CMT_CLS;
    /**  The class in which IncrementInterpolatedQuantities and ResetInterpolatedQuantities are defined */
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                     AbstractAssembler<ELEMENT_DIM, SPACE_DIM, 1u>,
                     typename AssemblerTraits<CONCRETE>::CMT_CLS>::type
            INTERPOLATE_CLS;
};

#endif //_SIMPLELINEARELLIPTICASSEMBLER_HPP_
