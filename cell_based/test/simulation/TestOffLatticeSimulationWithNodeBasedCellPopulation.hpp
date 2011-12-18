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

#ifndef TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_
#define TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomCellKiller.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "LogFile.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "SmartPointers.hpp"

class TestOffLatticeSimulationWithNodeBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    /**
     * Create a simulation of a NodeBasedCellPopulation with a NodeBasedCellPopulationMechanicsSystem.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayer() throw (Exception)
    {
        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulation");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    void TestSimulationWithBoxes() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulation");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with a NodeBasedCellPopulationMechanicsSystem
     * and a CellKiller. Test that no exceptions are thrown, and write the results to file.
     */
    void TestCellDeath() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationCellPtrDeath");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Add cell killer
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&node_based_cell_population, 0.997877574));
        simulator.AddCellKiller(p_killer);

        // Solve
        simulator.Solve();

        // Check some results
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 17u);

        std::vector<double> node_3_location = simulator.GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 3.4931, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 1.0062, 1e-4);

        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 1.0840, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 1.7208, 1e-4);
    }

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationStandardResult");
        simulator.SetEndTime(2.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) =-1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&node_based_cell_population, zero_vector<double>(2), normal)); // y>0
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        // Solve
        simulator.Solve();

        // Check some results
        std::vector<double> node_3_location = simulator.GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 2.9062, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 0.0172, 1e-4);

        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.8102, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.000, 1e-4);
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationSaveAndLoad");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) =-1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&node_based_cell_population, zero_vector<double>(2), normal)); // y>0
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        // Solve
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);
    }

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 1.0
        OffLatticeSimulation<2>* p_simulator1;
        p_simulator1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestOffLatticeSimulationWithNodeBasedCellPopulationSaveAndLoad", 0.1);

        p_simulator1->SetEndTime(1.0);
        p_simulator1->Solve();

        // Save, then reload and run from 1.0 to 2.5
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(p_simulator1);
        OffLatticeSimulation<2>* p_simulator2
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestOffLatticeSimulationWithNodeBasedCellPopulationSaveAndLoad", 1.0);

        p_simulator2->SetEndTime(2.5);
        p_simulator2->Solve();

        // These results are from time 2.5 in TestStandardResultForArchivingTestBelow()
        std::vector<double> node_3_location = p_simulator2->GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 2.9062, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 0.0172, 1e-4);

        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 3.8102, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.000, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_*/