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

#ifndef TESTBIDOMAINPARALLELMESH_HPP_
#define TESTBIDOMAINPARALLELMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "Hdf5DataReader.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "EventHandler.hpp"
#include "ParallelTetrahedralMesh.hpp"

class TestBidomainParallelMesh : public CxxTest::TestSuite
{
public:

    void TestBidomainProblemWithDistributedMesh1D()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));        
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("DistributedMesh1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("tetrahedral1d");
            
        Vec nondistributed_results;

        PlaneStimulusCellFactory<LuoRudyIModel1991OdeSystem, 1> cell_factory;

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        {
            ///////////////////////////////////////////////////////////////////
            // ParallelTetrahedralMesh
            ///////////////////////////////////////////////////////////////////
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            BidomainProblem<1> nondistributed_problem( &cell_factory );
            nondistributed_problem.SetMesh(&mesh);
            nondistributed_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            // now solve
            nondistributed_problem.Solve();

            VecDuplicate(nondistributed_problem.GetVoltage(), &nondistributed_results);
            VecCopy(nondistributed_problem.GetVoltage(), nondistributed_results);
        }


        ///////////////////////////////////////////////////////////////////
        // ParallelTetrahedralMesh
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetOutputFilenamePrefix("distributed1d");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        ParallelTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        BidomainProblem<1> distributed_problem( &cell_factory );

        distributed_problem.SetMesh(&mesh);

        distributed_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        // now solve
        distributed_problem.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector dist_nondistributed_voltage(nondistributed_results);
        DistributedVector::Stripe nondistributed_voltage(dist_nondistributed_voltage, 0);
        DistributedVector::Stripe nondistributed_potential(dist_nondistributed_voltage, 1);


        DistributedVector dist_distributed_voltage(distributed_problem.GetVoltage());
        DistributedVector::Stripe distributed_voltage(dist_distributed_voltage, 0);
        DistributedVector::Stripe distributed_potential(dist_distributed_voltage, 1);

        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index != DistributedVector::End();
             ++index)
        {
            if (index.Global==0)
            {
                TS_ASSERT_LESS_THAN(0, nondistributed_voltage[index]);
            }
            // the mono and bidomains should agree closely
            TS_ASSERT_DELTA(nondistributed_voltage[index], distributed_voltage[index], 1e-6);
            TS_ASSERT_DELTA(nondistributed_potential[index], distributed_potential[index], 1e-6);
        }

        VecDestroy(nondistributed_results);
    }
};

#endif /*TESTBIDOMAINPARALLELMESH_HPP_*/
