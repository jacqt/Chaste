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

#ifndef TESTCRYPTSIMULATIONBOUNDARYCONDITION_HPP_
#define TESTCRYPTSIMULATIONBOUNDARYCONDITION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CryptSimulationBoundaryCondition.hpp"
#include "WntConcentration.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class TestCryptSimulationBoundaryCondition : public AbstractCellBasedTestSuite
{
public:

    void TestSetAndGetUseJiggledBottomCells() throw (Exception)
    {
        CryptSimulationBoundaryCondition<1> boundary_condition1d(NULL);
        TS_ASSERT_EQUALS(boundary_condition1d.GetUseJiggledBottomCells(), false);

        boundary_condition1d.SetUseJiggledBottomCells(true);
        TS_ASSERT_EQUALS(boundary_condition1d.GetUseJiggledBottomCells(), true);

        CryptSimulationBoundaryCondition<2> boundary_condition2d(NULL);
        TS_ASSERT_EQUALS(boundary_condition2d.GetUseJiggledBottomCells(), false);

        boundary_condition2d.SetUseJiggledBottomCells(true);
        TS_ASSERT_EQUALS(boundary_condition2d.GetUseJiggledBottomCells(), true);
    }

    void TestConstructorWithCellPopulation1d() throw (Exception)
    {
        // Create 1D cell population
        unsigned num_cells = 5;
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(num_cells-1);

        std::vector<CellPtr> cells1d;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator1d;
        cells_generator1d.GenerateBasic(cells1d, mesh.GetNumNodes());
        MeshBasedCellPopulation<1> crypt1d(mesh, cells1d);

        // Pass this cell population to a boundary condition object
        CryptSimulationBoundaryCondition<1> boundary_condition1d(&crypt1d);

        // Test access to the cell population
        const AbstractCellPopulation<1>* p_population1d = boundary_condition1d.GetCellPopulation();
        TS_ASSERT(p_population1d != NULL);
    }

    void TestConstructorWithCellPopulation2d() throw (Exception)
    {
        // Create a 2d mesh-based cell population
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Pass this cell population to a boundary condition object
        CryptSimulationBoundaryCondition<2> boundary_condition(&crypt);

        // Test access to the cell population
        const AbstractCellPopulation<2>* p_population = boundary_condition.GetCellPopulation();
        TS_ASSERT(p_population != NULL);
    }

    void TestOutputParameters1d() throw(Exception)
    {
        std::string output_directory = "TestOutputParameters1d";
        OutputFileHandler output_file_handler(output_directory, false);

        CryptSimulationBoundaryCondition<1> boundary_condition(NULL);
        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "CryptSimulationBoundaryCondition-1");

        out_stream boundary_condition_parameter_file = output_file_handler.OpenOutputFile("results1d.parameters");
        boundary_condition.OutputCellPopulationBoundaryConditionParameters(boundary_condition_parameter_file);
        boundary_condition_parameter_file->close();

        std::string boundary_condition_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + boundary_condition_results_dir + "results1d.parameters crypt/test/data/TestCryptSimulationBoundaryCondition/results1d.parameters").c_str()), 0);

        // Test OutputCellPopulationBoundaryConditionInfo() method
        out_stream boundary_condition_info_file = output_file_handler.OpenOutputFile("results1d.info");
        boundary_condition.OutputCellPopulationBoundaryConditionInfo(boundary_condition_info_file);
        boundary_condition_info_file->close();

        TS_ASSERT_EQUALS(system(("diff " + boundary_condition_results_dir + "results1d.info crypt/test/data/TestCryptSimulationBoundaryCondition/results1d.info").c_str()), 0);
    }

    void TestOutputParameters2d() throw(Exception)
    {
        std::string output_directory = "TestOutputParameters2d";
        OutputFileHandler output_file_handler(output_directory, false);

        CryptSimulationBoundaryCondition<2> boundary_condition(NULL);
        TS_ASSERT_EQUALS(boundary_condition.GetIdentifier(), "CryptSimulationBoundaryCondition-2");

        out_stream boundary_condition_parameter_file = output_file_handler.OpenOutputFile("results2d.parameters");
        boundary_condition.OutputCellPopulationBoundaryConditionParameters(boundary_condition_parameter_file);
        boundary_condition_parameter_file->close();

        std::string boundary_condition_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + boundary_condition_results_dir + "results2d.parameters crypt/test/data/TestCryptSimulationBoundaryCondition/results2d.parameters").c_str()), 0);

        // Test OutputCellPopulationBoundaryConditionInfo() method
        out_stream boundary_condition_info_file = output_file_handler.OpenOutputFile("results2d.info");
        boundary_condition.OutputCellPopulationBoundaryConditionInfo(boundary_condition_info_file);
        boundary_condition_info_file->close();

        TS_ASSERT_EQUALS(system(("diff " + boundary_condition_results_dir + "results2d.info crypt/test/data/TestCryptSimulationBoundaryCondition/results2d.info").c_str()), 0);
    }

    void TestImposeBoundaryConditionWithNoWnt1d() throw(Exception)
    {
        // Create 1D cell population
        unsigned num_cells = 5;
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(num_cells-1);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Store the node locations
        std::vector<c_vector<double, 1> > node_locations_before(mesh.GetNumNodes());
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Now move the first cell (which should be on x=0, and a stem cell) to the left a bit
        AbstractCellPopulation<1>::Iterator cell_iter = crypt.Begin();

        // Check is a stem cell
        TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), STEM);

        // Check initially at x=0
        TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[0], 0.0, 1e-6);

        // Now move to the left a bit
        crypt.GetNode(0)->rGetModifiableLocation()[0] = -0.1;
        TS_ASSERT_LESS_THAN(crypt.GetLocationOfCellCentre(*cell_iter)[0], 0.0);

        // Pass this cell population to a boundary condition object
        CryptSimulationBoundaryCondition<1> boundary_condition(&crypt);

        // Impose the boundary condition without a Wnt stimulus
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

        /*
         * The first cell should have been pulled back to x=0 since it is a stem cell and
         * there is no Wnt stimulus. It should be unaffected by jiggling, which is not imposed
         * in 1D.
         */
        TS_ASSERT_DELTA(0.0, crypt.GetLocationOfCellCentre(*cell_iter)[0], 1e-3);

        // The nodes should all now be at their original locations
        std::vector<c_vector<double, 1> > node_locations_after(mesh.GetNumNodes());
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            node_locations_after[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(node_locations_before[node_index](0), node_locations_after[node_index](0), 1e-3);
        }
    }

    void TestImposeBoundaryConditionWithNoWntOrJiggling2d() throw(Exception)
    {
        // Create a cell population
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Now move the first cell (which should be on y=0) down a bit
        AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
        TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0, 1e-6);
        crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;
        TS_ASSERT_LESS_THAN(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0);

        // Create a boundary condition object
        CryptSimulationBoundaryCondition<2> boundary_condition(&crypt);
        boundary_condition.SetUseJiggledBottomCells(false);

        // Impose the boundary condition without a Wnt stimulus or jiggled bottom cells
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

        // The first cell should have been moved back by the boundary condition object
        TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0, 1e-4);

        // The nodes should all now be at their original locations
        std::vector<c_vector<double, 2> > node_locations_after(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            node_locations_after[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(node_locations_before[node_index](0), node_locations_after[node_index](0), 1e-3);
            TS_ASSERT_DELTA(node_locations_before[node_index](1), node_locations_after[node_index](1), 1e-3);
        }
    }

    void TestImposeBoundaryConditionWithNoWntButWithJiggling() throw(Exception)
    {
        // Create a cell population
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Move the first cell (which should be on y=0, and is a stem cell) down a bit
        AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
        crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;

        // Also move the second cell (which should be on y=0, and we make a transit cell)
        ++cell_iter;
        TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), STEM);
        cell_iter->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);
        TS_ASSERT_EQUALS(cell_iter->GetCellCycleModel()->GetCellProliferativeType(), TRANSIT);
        TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0, 1e-6);
        crypt.GetNode(1)->rGetModifiableLocation()[1] = -0.1;
        TS_ASSERT_LESS_THAN(crypt.GetLocationOfCellCentre(*cell_iter)[1], 0.0);

        // Create a boundary condition object
        CryptSimulationBoundaryCondition<2> boundary_condition(&crypt);
        boundary_condition.SetUseJiggledBottomCells(true);

        // Impose the boundary condition without a Wnt stimulus but with jiggled bottom cells
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

        /*
         * The first cell should have been pulled up to y=0 since it is a stem cell and
         * there is no Wnt stimulus. It should then be unaffected by the jiggling. The
         * second cell should not have been pulled up since it is a stem cell, but should
         * then have been moved to above y=0 by the jiggling.
         */
        AbstractCellPopulation<2>::Iterator cell_iter2 = crypt.Begin();
        TS_ASSERT_DELTA(0.0, crypt.GetLocationOfCellCentre(*cell_iter2)[1], 1e-3);
        ++cell_iter2;
        TS_ASSERT_LESS_THAN(0.0, crypt.GetLocationOfCellCentre(*cell_iter2)[1]);
    }

    void TestImposeBoundaryConditionWithWnt1d() throw(Exception)
    {
        // Create 1D cell population
        unsigned num_cells = 5;
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(num_cells-1);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up a WntConcentration singleton object
        WntConcentration<1>* p_wnt = WntConcentration<1>::Instance();
        WntConcentration<1>::Instance()->SetType(LINEAR);
        WntConcentration<1>::Instance()->SetCellPopulation(crypt);
        WntConcentration<1>::Instance()->SetCryptLength(crypt.GetWidth(0));
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        // Store the node locations
        std::vector<c_vector<double, 1> > node_locations_before(mesh.GetNumNodes());
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Now move the first cell (which should be on x=0) to the left a bit
        AbstractCellPopulation<1>::Iterator cell_iter = crypt.Begin();
        // Check initially at x=0
        TS_ASSERT_DELTA(crypt.GetLocationOfCellCentre(*cell_iter)[0], 0.0, 1e-6);
        // Now move to the left a bit
        crypt.GetNode(0)->rGetModifiableLocation()[0] = -0.1;
        TS_ASSERT_DELTA(-0.1, crypt.GetLocationOfCellCentre(*cell_iter)[0], 1e-3);

        // Pass this cell population to a boundary condition object
        CryptSimulationBoundaryCondition<1> boundary_condition(&crypt);

        // Impose the boundary condition with a Wnt stimulus
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

        // This node will be moved back to 0.0
        TS_ASSERT_DELTA(0.0, crypt.GetLocationOfCellCentre(*cell_iter)[0], 1e-3);

        // Tidy up
        WntConcentration<1>::Destroy();
    }


    void TestImposeBoundaryConditionWithWntButNoJiggling() throw(Exception)
    {
        // Create a cell population
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Move the first cell (which should be on y=0) down a bit
        AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
        crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;
        c_vector<double, 2> perturbed_location = crypt.GetLocationOfCellCentre(*cell_iter);

        // Create a boundary condition object
        CryptSimulationBoundaryCondition<2> boundary_condition(&crypt);
        boundary_condition.SetUseJiggledBottomCells(false);

        // Set up a WntConcentration singleton object
        WntConcentration<2>* p_wnt = WntConcentration<2>::Instance();
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt.GetWidth(1));
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        // Impose the boundary condition with a Wnt stimulus but without jiggled bottom cells
        boundary_condition.ImposeBoundaryCondition(node_locations_before);

        /*
         * The first cell should not have been moved by the boundary condition object,
         * and hence should remain in its perturbed location.
         */
        c_vector<double, 2> location_after = crypt.GetLocationOfCellCentre(*cell_iter);
        TS_ASSERT_DELTA(location_after[0], perturbed_location[0], 1e-3);
        TS_ASSERT_DELTA(location_after[1], location_after[1], 1e-3);

        // Tidy up
        WntConcentration<2>::Destroy();
    }

    void TestVerifyBoundaryCondition1d() throw(Exception)
    {
        // Create 1D cell population
        unsigned num_cells = 5;
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(num_cells-1);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<1> crypt(mesh, cells);

        // Set up a WntConcentration singleton object
        WntConcentration<1>* p_wnt = WntConcentration<1>::Instance();
        WntConcentration<1>::Instance()->SetType(LINEAR);
        WntConcentration<1>::Instance()->SetCellPopulation(crypt);
        WntConcentration<1>::Instance()->SetCryptLength(crypt.GetWidth(0));
        TS_ASSERT_EQUALS(p_wnt->IsWntSetUp(), true);

        // Store the node locations
        std::vector<c_vector<double, 1> > node_locations_before(mesh.GetNumNodes());
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Now move the first cell (which should be on x=0) to the left
        AbstractCellPopulation<1>::Iterator cell_iter = crypt.Begin();
        crypt.GetNode(0)->rGetModifiableLocation()[0] = -0.1;

        // Create a boundary condition object
        CryptSimulationBoundaryCondition<1> boundary_condition(&crypt);

        // Before imposing the boundary condition, it should not be satisfied
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), false);

        boundary_condition.ImposeBoundaryCondition(node_locations_before);

        // After imposing the boundary condition, it should be satisfied
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);
    }

    void TestVerifyBoundaryCondition2d() throw(Exception)
    {
        // Create a cell population
        CylindricalHoneycombMeshGenerator generator(4, 4, 0, 1.0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true);

        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Store the node locations
        std::vector<c_vector<double, 2> > node_locations_before(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            node_locations_before[node_index] = crypt.GetNode(node_index)->rGetLocation();
        }

        // Now move the first cell (which should be on y=0) down a bit
        AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
        crypt.GetNode(0)->rGetModifiableLocation()[1] = -0.1;

        // Create a boundary condition object
        CryptSimulationBoundaryCondition<2> boundary_condition(&crypt);

        // Before imposing the boundary condition, it should not be satisfied
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), false);

        boundary_condition.ImposeBoundaryCondition(node_locations_before);

        // After imposing the boundary condition, it should be satisfied
        TS_ASSERT_EQUALS(boundary_condition.VerifyBoundaryCondition(), true);
    }

    void TestArchiving1d() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false); // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CryptSimulationBoundaryCondition<1>.arch";

        {
            // Create an output archive
            CryptSimulationBoundaryCondition<1> boundary_condition(NULL);
            boundary_condition.SetUseJiggledBottomCells(true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            AbstractCellPopulationBoundaryCondition<1>* const p_boundary_condition = &boundary_condition;
            output_arch << p_boundary_condition;
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCellPopulationBoundaryCondition<1>* p_boundary_condition;

            // Restore from the archive
            input_arch >> p_boundary_condition;

            // Test we have restored the plane geometry correctly
            TS_ASSERT_EQUALS(static_cast<CryptSimulationBoundaryCondition<1>*>(p_boundary_condition)->GetUseJiggledBottomCells(), true);

            // Tidy up
            delete p_boundary_condition;
       }
    }

    void TestArchiving2d() throw (Exception)
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false); // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CryptSimulationBoundaryCondition<2>.arch";

        {
            // Create an output archive
            CryptSimulationBoundaryCondition<2> boundary_condition(NULL);
            boundary_condition.SetUseJiggledBottomCells(true);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            AbstractCellPopulationBoundaryCondition<2>* const p_boundary_condition = &boundary_condition;
            output_arch << p_boundary_condition;
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCellPopulationBoundaryCondition<2>* p_boundary_condition;

            // Restore from the archive
            input_arch >> p_boundary_condition;

            // Test we have restored the plane geometry correctly
            TS_ASSERT_EQUALS(static_cast<CryptSimulationBoundaryCondition<2>*>(p_boundary_condition)->GetUseJiggledBottomCells(), true);

            // Tidy up
            delete p_boundary_condition;
       }
    }
};

#endif /* TESTCRYPTSIMULATIONBOUNDARYCONDITION_HPP_ */