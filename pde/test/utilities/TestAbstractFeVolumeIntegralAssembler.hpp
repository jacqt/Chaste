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

#ifndef TESTABSTRACTFEVOLUMEINTEGRALASSEMBLER_HPP_
#define TESTABSTRACTFEVOLUMEINTEGRALASSEMBLER_HPP_

#include "AbstractFeVolumeIntegralAssembler.hpp"
#include "TetrahedralMesh.hpp"
#include "MassMatrixAssembler.hpp"
#include "StiffnessMatrixAssembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscMatTools.hpp"


// Note: PROBLEM_DIM>1 is not tested here, so only in coupled PDE solves


// simple assembler which just returns a c_vector of ones for each quad point
//  => final 'b' vector for each element will be equal to elem_vol*(1,..,1)
template<unsigned DIM>
class BasicVectorAssembler : public AbstractFeVolumeIntegralAssembler<DIM,DIM,1,true,false,NORMAL>
{
private:
    c_vector<double,1*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
        c_vector<double,1*(DIM+1)> vec;
        for (unsigned i=0; i<DIM+1; i++)
        {
            vec(i)=1.0;
        }
        return vec;
    }

public:
    BasicVectorAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<DIM,DIM,1,true,false,NORMAL>(pMesh)
    {
    }
};


// simple assembler which just returns a c_matrix of ones for each quad point
//  => final 'A' matrix for each element will be equal to elem_vol*(1,..,1)
template<unsigned DIM>
class BasicMatrixAssembler : public AbstractFeVolumeIntegralAssembler<DIM,DIM,1,false,true,NORMAL>
{
private:
    c_matrix<double,1*(DIM+1),1*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,1*(DIM+1),1*(DIM+1)> mat;
        for (unsigned i=0; i<DIM+1; i++)
        {
            for (unsigned j=0; j<DIM+1; j++)
            {
                mat(i,j)=1.0;
            }
        }
        return mat;
    }
public:
    BasicMatrixAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<DIM,DIM,1,false,true,NORMAL>(pMesh)
    {
    }
};



// Assembler which does both of the above
template<unsigned DIM>
class BasicVectorAndMatrixAssembler : public AbstractFeVolumeIntegralAssembler<DIM,DIM,1,true,true,NORMAL>
{
private:
    c_matrix<double,1*(DIM+1),1*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,1*(DIM+1),1*(DIM+1)> mat;
        for (unsigned i=0; i<DIM+1; i++)
        {
            for (unsigned j=0; j<DIM+1; j++)
            {
                mat(i,j)=1.0;
            }
        }
        return mat;
    }

    c_vector<double,1*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
        c_vector<double,1*(DIM+1)> vec;
        for (unsigned i=0; i<DIM+1; i++)
        {
            vec(i)=1.0;
        }
        return vec;
    }
public:
    BasicVectorAndMatrixAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<DIM,DIM,1,true,true,NORMAL>(pMesh)
    {
    }
};


// Assembler for checking the current solution is interpolated and passed up into here
class TestingAssembler : public AbstractFeVolumeIntegralAssembler<1,1,1,true,false,NORMAL>
{
private:
    c_vector<double,1*(1+1)> ComputeVectorTerm(
        c_vector<double, 1+1>& rPhi,
        c_matrix<double, 1, 1+1>& rGradPhi,
        ChastePoint<1>& rX,
        c_vector<double,1>& rU,
        c_matrix<double, 1, 1>& rGradU,
        Element<1,1>* pElement)
    {
        TS_ASSERT_LESS_THAN(0.0, rX[0]);             // check x is not 0.0
        TS_ASSERT_DELTA(rU(0), 10 + 10*rX[0], 1e-6); // check u,x interpolated properly
        return zero_vector<double>(1+1);
    }
public:
    TestingAssembler(AbstractTetrahedralMesh<1,1>* pMesh)
        : AbstractFeVolumeIntegralAssembler<1,1,1,true,false,NORMAL>(pMesh)
    {
    }
};


class TestAbstractFeVolumeIntegralAssembler : public CxxTest::TestSuite
{
private:

    template<unsigned DIM>
    void DoTestBasicVectorAssemblers()
    {
        TetrahedralMesh<DIM,DIM> mesh;
        double h = 0.1;
        if (DIM==1)
        {
            mesh.ConstructRegularSlabMesh(h, 1.0);
        }
        else if (DIM==2)
        {
            mesh.ConstructRegularSlabMesh(h, 1.0, 1.0);
        }
        else
        {
            h = 0.5;
            mesh.ConstructRegularSlabMesh(h, 1.0, 1.0, 1.0);
        }

        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes());

        BasicVectorAssembler<DIM> basic_vector_assembler(&mesh);

        basic_vector_assembler.SetVectorToAssemble(vec,true);
        basic_vector_assembler.Assemble();

        PetscVecTools::Finalise(vec);

        ReplicatableVector vec_repl(vec);
        double volume_of_element = DIM==1 ? h : (DIM==2 ? 0.5*h*h : h*h*h/6);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(vec_repl[i], volume_of_element*mesh.GetNode(i)->GetNumContainingElements(), 1e-4);
        }

        VecDestroy(vec);
    }

public:

    // Test vector assembly
    void TestBasicVectorAssemblers() throw(Exception)
    {
        DoTestBasicVectorAssemblers<1>();
        DoTestBasicVectorAssemblers<2>();
        DoTestBasicVectorAssemblers<3>();
    }




    // Test matrix assembly
    // only test the 1d one as here as more difficult to write down correct matrix on paper
    // Tested better with 2d mass matrix below
    void TestBasicMatrixAssemblers() throw(Exception)
    {
        TetrahedralMesh<1,1> mesh;
        double h = 0.1;
        mesh.ConstructRegularSlabMesh(h, h); // must be one element for this test to work

        Mat mat;
        PetscTools::SetupMat(mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 2);

        BasicMatrixAssembler<1> basic_matrix_assembler_1d(&mesh);

        basic_matrix_assembler_1d.SetMatrixToAssemble(mat);
        basic_matrix_assembler_1d.Assemble();

        PetscMatTools::Finalise(mat);

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                TS_ASSERT_DELTA(value, h, 1e-4);
            }
        }

        MatDestroy(mat);
    }

    // Test the ability to assemble both a vector and a matrix
    void TestBasicVectorAndMatrixAssembler() throw(Exception)
    {
        TetrahedralMesh<1,1> mesh;
        double h = 0.1;
        mesh.ConstructRegularSlabMesh(h, h); // must be one element for this test to work

        Mat mat;
        PetscTools::SetupMat(mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 2);

        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes());

        BasicVectorAndMatrixAssembler<1> assembler_1d(&mesh);

        TS_ASSERT_THROWS_THIS(assembler_1d.Assemble(), "Matrix to be assembled has not been set");

        assembler_1d.SetMatrixToAssemble(mat);

        TS_ASSERT_THROWS_THIS(assembler_1d.Assemble(), "Vector to be assembled has not been set");

        assembler_1d.SetVectorToAssemble(vec, true);

        assembler_1d.Assemble();

        PetscVecTools::Finalise(vec);
        PetscMatTools::Finalise(mat);

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);

        ReplicatableVector vec_repl(vec);

        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            TS_ASSERT_DELTA(vec_repl[i], h, 1e-4);
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                TS_ASSERT_DELTA(value, h, 1e-4);
            }
        }

        // call assemble again and check everything is the same (so everything
        // was zeroed, not added to)

        assembler_1d.Assemble();

        PetscVecTools::Finalise(vec);
        PetscMatTools::Finalise(mat);

        ReplicatableVector vec_repl2(vec);

        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            TS_ASSERT_DELTA(vec_repl2[i], h, 1e-4);
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                TS_ASSERT_DELTA(value, h, 1e-4);
            }
        }

        // zero the matrix and call AssembleVector, check the matrix is not
        // assembled.
        PetscVecTools::Zero(vec);
        PetscMatTools::Zero(mat);

        assembler_1d.AssembleVector();

        PetscVecTools::Finalise(vec);
        PetscMatTools::Finalise(mat);

        ReplicatableVector vec_repl3(vec);

        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            TS_ASSERT_DELTA(vec_repl3[i], h, 1e-4);
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                // value is 0.0 now
                TS_ASSERT_DELTA(value, 0.0, 1e-4);
            }
        }

        // now call AssembleMatrix()
        assembler_1d.AssembleMatrix();

        PetscMatTools::Finalise(mat);

        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                TS_ASSERT_DELTA(value, h, 1e-4);
            }
        }

        VecDestroy(vec);
        MatDestroy(mat);

    }

    void TestMassMatrixAssembler1d() throw(Exception)
    {
        TetrahedralMesh<1,1> mesh;
        double h = 0.1;
        mesh.ConstructRegularSlabMesh(h, 0.5);

        Mat mat;
        PetscTools::SetupMat(mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 3);

        MassMatrixAssembler<1,1> assembler(&mesh);

        assembler.SetMatrixToAssemble(mat);
        assembler.Assemble();

        PetscMatTools::Finalise(mat);

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);

        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                if (i>0 && i<mesh.GetNumNodes()-1)
                {
                    // All rows except first and last should look like
                    // [0, .., 0, h/6, 4h/6, h/6, 0, .., 0]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, 4.0*h/6.0, 1e-5);
                    }
                    else if (j==i+1 || j+1==i)
                    {
                        TS_ASSERT_DELTA(value, 1.0*h/6.0, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if (i==0)
                {
                    // top row: [2h/6, h/6, 0, .., 0]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, 2.0*h/6.0, 1e-5);
                    }
                    else if (j==i+1)
                    {
                        TS_ASSERT_DELTA(value, 1.0*h/6.0, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if (i+1==mesh.GetNumNodes())
                {
                    // bottom row: [0, .., 0, h/6, 2h/6]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, 2.0*h/6.0, 1e-5);
                    }
                    else if (j+1==i)
                    {
                        TS_ASSERT_DELTA(value, 1.0*h/6.0, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }

            }
        }

        MatDestroy(mat);
    }


    void TestMassMatrixAssembler2dIncludingScaleFactor() throw(Exception)
    {
        TetrahedralMesh<2,2> mesh;
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements"); // so we know the exact connectivity
        mesh.ConstructFromMeshReader(reader);

        Mat mat;
        PetscTools::SetupMat(mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 4);

        double scale_factor = 2.464525345;
        MassMatrixAssembler<2,2> assembler(&mesh, false, scale_factor);

        assembler.SetMatrixToAssemble(mat);
        assembler.Assemble();

        PetscMatTools::Finalise(mat);

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);

        // check against exact answers

        if (lo==0)
        {
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,0,0)/scale_factor, 1.0/12, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,0,1)/scale_factor, 1.0/24, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,0,2)/scale_factor,    0.0, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,0,3)/scale_factor, 1.0/24, 1e-6);
        }

        if (lo<=1 && 1<hi)
        {
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,1,0)/scale_factor, 1.0/24, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,1,1)/scale_factor, 1.0/6,  1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,1,2)/scale_factor, 1.0/24, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,1,3)/scale_factor, 1.0/12, 1e-6);
        }

        if (lo<=2 && 2<hi)
        {
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,2,0)/scale_factor,    0.0, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,2,1)/scale_factor, 1.0/24, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,2,2)/scale_factor, 1.0/12, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,2,3)/scale_factor, 1.0/24, 1e-6);
        }

        if (lo<=3 && 3<hi)
        {
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,3,0)/scale_factor, 1.0/24, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,3,1)/scale_factor, 1.0/12, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,3,2)/scale_factor, 1.0/24, 1e-6);
            TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,3,3)/scale_factor, 1.0/6 , 1e-6);
        }

        MatDestroy(mat);
    }

    void TestStiffnessMatrixAssembler1d() throw(Exception)
    {
        TetrahedralMesh<1,1> mesh;
        double h = 0.1;
        mesh.ConstructRegularSlabMesh(h, 0.5);

        Mat mat;
        PetscTools::SetupMat(mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 3);

        StiffnessMatrixAssembler<1,1> assembler(&mesh);

        assembler.SetMatrixToAssemble(mat);
        assembler.Assemble();

        PetscMatTools::Finalise(mat);

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);

        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                if (i>0 && i<mesh.GetNumNodes()-1)
                {
                    // All rows except first and last should look like
                    // [0, .., 0, -1/h, 2/h, -1/h, 0, .., 0]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, 2.0/h, 1e-5);
                    }
                    else if (j==i+1 || j+1==i)
                    {
                        TS_ASSERT_DELTA(value, -1.0/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if (i==0)
                {
                    // top row: [1/h, -1/h, 0, .., 0]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, 1.0/h, 1e-5);
                    }
                    else if (j==i+1)
                    {
                        TS_ASSERT_DELTA(value, -1.0/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if (i+1==mesh.GetNumNodes())
                {
                    // bottom row: [0, .., 0, -1/h, 1/h]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, 1.0/h, 1e-5);
                    }
                    else if (j+1==i)
                    {
                        TS_ASSERT_DELTA(value, -1.0/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }

            }
        }

        MatDestroy(mat);
    }

    void TestInterpolationOfPositionAndCurrentSolution() throw(Exception)
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0);

        std::vector<double> u(2);
        u[0] = 10.0;
        u[1] = 20.0;
        Vec current_solution = PetscTools::CreateVec(u);
        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes());

        TestingAssembler assembler(&mesh);

        assembler.SetVectorToAssemble(vec,true);
        assembler.SetCurrentSolution(current_solution);

        // No tests here, there are two TS_ASSERTS insider TestingAssembler::ComputeVectorTerm()
        assembler.Assemble();

        VecDestroy(vec);
        VecDestroy(current_solution);
    }
};
#endif /*TESTABSTRACTFEVOLUMEINTEGRALASSEMBLER_HPP_*/