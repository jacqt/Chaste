/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTCABLETESTPROBLEM_HPP_
#define TESTCABLETESTPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include "MixedDimensionMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "ReplicatableVector.hpp"
#include "StiffnessMatrixAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractFeCableIntegralAssembler.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"


template<unsigned DIM>
class CableTestProblemRhsAssembler :  public AbstractFeCableIntegralAssembler<DIM,DIM,1,true,false,NORMAL>
{
private:
    c_vector<double,1*2> ComputeCableVectorTerm(
        c_vector<double, 2>& rPhi,
        c_matrix<double, DIM, 2>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<1,DIM>* pElement)
    {
        return -rPhi;
    }

public:
    CableTestProblemRhsAssembler(MixedDimensionMesh<DIM,DIM>* pMesh)
        : AbstractFeCableIntegralAssembler<DIM,DIM,1,true,false,NORMAL>(pMesh)
    {
    }
};




// solver
template<unsigned DIM>
class CableTestProblemSolver: public AbstractStaticLinearPdeSolver<DIM,DIM,1>
{
private:
	StiffnessMatrixAssembler<DIM,DIM>* mpStiffnessMatrixAssembler;
	CableTestProblemRhsAssembler<DIM>* mpRhsAssembler;

	/** Boundary conditions */
	BoundaryConditionsContainer<DIM,DIM,1>* mpBoundaryConditions;

    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
	{
        // use assemblers to set up KU=b
        mpStiffnessMatrixAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix());
        mpStiffnessMatrixAssembler->Assemble();

        mpRhsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), true);
        mpRhsAssembler->Assemble();

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->SwitchWriteModeLhsMatrix();

        // apply the Dirichlet boundary conditions
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->FinaliseLhsMatrix();
	}


public:
	CableTestProblemSolver(MixedDimensionMesh<DIM,DIM>* pMesh,
                           BoundaryConditionsContainer<DIM,DIM,1>* pBoundaryConditions)
		 : AbstractStaticLinearPdeSolver<DIM,DIM,1>(pMesh),
		   mpBoundaryConditions(pBoundaryConditions)
    {
		// set-up mStiffnessMatrixAssembler and mRhsAssembler
		mpStiffnessMatrixAssembler = new StiffnessMatrixAssembler<DIM,DIM>(pMesh);
		mpRhsAssembler = new CableTestProblemRhsAssembler<DIM>(pMesh);
    }

	~CableTestProblemSolver()
	{
		delete mpStiffnessMatrixAssembler;
		delete mpRhsAssembler;
	}
};

// See #1798
// Solve the problem \nabla^2 u = \delta_{cable}, on a cylindrical geometry r \in [0,1], z\in[0,1]
// with zero-Neumann BCs on the top and bottom surfaces and zero-Dirichlet BCs on the surface r=1.
// The solution is log(r)/2*pi.
//
class TestCableTestProblem : public CxxTest::TestSuite
{
public:
	void TestSolvingTestProblem() throw(Exception)
	{
		std::string mesh_base("mesh/test/data/mixed_dimension_meshes/cylinder");
		TrianglesMeshReader<3,3> reader(mesh_base);
		MixedDimensionMesh<3,3> mesh;
		mesh.ConstructFromMeshReader(reader);

        ConstBoundaryCondition<3>* p_zero_boundary_condition = new ConstBoundaryCondition<3>(0.0);

        // apply boundary conditions u=0 on curved surface of cylinder (zero neumann BCs are
        // applied on the bottom and top surface (z=0,1).
        BoundaryConditionsContainer<3,3,1> bcc;

//// NOTE: this mesh doesn't have a (valid) face file defined, so the face and boundary node
//// data isn't set up in the mesh. Therefore, rather than looping over boundary nodes, we
//// we loop over all nodes, and set the node to be a boundary node before passing it into
//// the BCC.
//        for (MixedDimensionMesh<3,3>::BoundaryNodeIterator iter =
//               mesh.GetBoundaryNodeIteratorBegin();
//             iter != mesh.GetBoundaryNodeIteratorEnd();
//	         iter++)

		for (AbstractTetrahedralMesh<3,3>::NodeIterator current_node = mesh.GetNodeIteratorBegin();
	         current_node != mesh.GetNodeIteratorEnd();
	         ++current_node)
		{
		    Node<3>* p_node = &(*current_node); //Get pointer to the current node from the iterator
            double x = p_node->rGetLocation()[0];
		    double y = p_node->rGetLocation()[1];
            double r = sqrt(x*x+y*y);

		    if(fabs(r-1)<1e-3)
		    {
	            p_node->SetAsBoundaryNode(); // see comment above!
		        bcc.AddDirichletBoundaryCondition(p_node, p_zero_boundary_condition);
		    }
	    }

		// Solver
		CableTestProblemSolver<3> cable_solver(&mesh,&bcc);
		Vec result = cable_solver.Solve();
		double* p_result;
		VecGetArray(result, &p_result);
		// Solution should be u = log(r)/(2*pi)
		for (AbstractTetrahedralMesh<3,3>::NodeIterator current_node = mesh.GetNodeIteratorBegin();
			 current_node != mesh.GetNodeIteratorEnd();
			 ++current_node)
		{
			double x = current_node->GetPoint()[0];
			double y = current_node->GetPoint()[1];
			double r = sqrt(x*x+y*y);


			unsigned local_index = current_node->GetIndex() - mesh.GetDistributedVectorFactory()->GetLow();
			if(r>0.1)
			{
			    // use a tolerance that is weighted by 1-r as accuracy will decrease as
			    // get closer to r=0 for which u=-infty. Visually the solution compared
			    // to the true solution looks pretty good.
                double u = log(r)/(2*M_PI);
			    TS_ASSERT_DELTA(p_result[local_index], u, 0.07*(1-r)+1e-6);
			}
			else
			{
			    // for these nodes r=0 and the true solution is -infinity
			    TS_ASSERT_LESS_THAN(p_result[local_index], -0.4);
			}
		}

		VecRestoreArray(result, &p_result);
		VecDestroy(result);
	}
};

#endif /* TESTCABLETESTPROBLEM_HPP_ */
