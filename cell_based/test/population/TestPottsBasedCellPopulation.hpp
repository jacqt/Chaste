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

#ifndef TESTPOTTSBASEDCELLPOPULATION_HPP_
#define TESTPOTTSBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "PottsMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"

class TestPottsBasedCellPopulation : public AbstractCellBasedTestSuite
{
public:

    void TestConstructor() throw(Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test that the mesh and cells are correctly assigned
        TS_ASSERT_EQUALS(&(cell_population.rGetMesh()), p_mesh);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), p_mesh->GetNumNodes());

        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(*cell_iter, cells[i]);
            ++cell_iter;
        }

        // Test that we do not have any update rules present
        TS_ASSERT_EQUALS(cell_population.rGetUpdateRuleCollection().empty(), true);

        // Test that the other member variables of this object are initialised correctly
        TS_ASSERT(cell_population.GetElementTessellation() == NULL);
        TS_ASSERT_DELTA(cell_population.GetTemperature(), 0.1, 1e-6);
        TS_ASSERT_EQUALS(cell_population.GetNumSweepsPerTimestep(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetUpdateNodesInRandomOrder(), true);
        TS_ASSERT_EQUALS(cell_population.GetIterateRandomlyOverUpdateRuleCollection(), false);
    }

    void TestIsCellAssociatedWithADeletedLocation() throw (Exception)
    {
        // Create a Potts-based cell population but do not try to validate
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();
        p_mesh->GetElement(0)->MarkAsDeleted();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells, false, false);

        // Test IsCellAssociatedWithADeletedLocation() method
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            bool is_deleted = cell_population.IsCellAssociatedWithADeletedLocation(*cell_iter);

            if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 0)
            {
                TS_ASSERT_EQUALS(is_deleted, true);
            }
            else
            {
                TS_ASSERT_EQUALS(is_deleted, false);
            }
        }
    }

    void TestValidate() throw (Exception)
    {
        // Create a simple potts-based mesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements()-1);

        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<cells.size(); i++)
        {
            cell_location_indices.push_back(i);
        }

        // This should throw an exception as the number of cells does not equal the number of elements
        std::vector<CellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(PottsBasedCellPopulation<2> cell_population(*p_mesh, cells_copy),
                "Element 3 does not appear to have a cell associated with it");

        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_cell(new Cell(p_state, p_model));

        double birth_time = 0.0 - p_mesh->GetNumElements()-1;
        p_cell->SetBirthTime(birth_time);

        cells.push_back(p_cell);
        cell_location_indices.push_back(p_mesh->GetNumElements()-1);

        // This should pass as the number of cells equals the number of elements
        std::vector<CellPtr> cells_copy2(cells);
        TS_ASSERT_THROWS_NOTHING(PottsBasedCellPopulation<2> cell_population(*p_mesh, cells_copy2));

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Check correspondence between elements and cells
        for (PottsMesh<2>::PottsElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            std::set<unsigned> expected_node_indices;
            unsigned expected_index = iter->GetIndex();

            for (unsigned i=0; i<iter->GetNumNodes(); i++)
            {
                expected_node_indices.insert(iter->GetNodeGlobalIndex(i));
            }

            std::set<unsigned> actual_node_indices;
            unsigned elem_index = iter->GetIndex();
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(elem_index);
            PottsElement<2>* p_actual_element = cell_population.GetElementCorrespondingToCell(p_cell);
            unsigned actual_index = p_actual_element->GetIndex();

            for (unsigned i=0; i<p_actual_element->GetNumNodes(); i++)
            {
                actual_node_indices.insert(p_actual_element->GetNodeGlobalIndex(i));
            }

            TS_ASSERT_EQUALS(actual_index, expected_index);
            TS_ASSERT_EQUALS(actual_node_indices, expected_node_indices);
        }

        // Create another simple potts-based mesh
        PottsMeshGenerator<2> generator2(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh2 = generator2.GetMesh();

        // Create cells
        std::vector<CellPtr> cells2;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, p_mesh2->GetNumElements()+1);

        std::vector<unsigned> cell_location_indices2;
        for (unsigned i=0; i<cells2.size(); i++)
        {
            cell_location_indices2.push_back(i%p_mesh2->GetNumElements()); // Element 0 will have 2 cells
        }

        // This should throw an exception as the number of cells
        // does not equal the number of elements
        TS_ASSERT_THROWS_THIS(PottsBasedCellPopulation<2> cell_population2(*p_mesh2, cells2, false, true, cell_location_indices2),
                "Element 0 appears to have 2 cells associated with it");
    }

    void TestRemoveDeadCellsAndUpdate() throw(Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test RemoveDeadCells() method
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 4u);

        cell_population.Begin()->Kill();

        TS_ASSERT_EQUALS(cell_population.RemoveDeadCells(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 4u);

        // Test that Update() throws no errors
        TS_ASSERT_THROWS_NOTHING(cell_population.Update());
    }

    void TestAddCell() throw(Exception)
    {
        // Create a simple 2D PottsMesh with one cell
        PottsMeshGenerator<2> generator(2, 1, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);

        // Create a new cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_new_cell(new Cell(p_state, p_model));

        // Add new cell to the cell population by dividing the cell
        AbstractCellPopulation<2>::Iterator cell_iter_1 = cell_population.Begin();
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), *cell_iter_1);

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_new_cell), 1u);

        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);

        // Check locations of parent and daughter cell
        AbstractCellPopulation<2>::Iterator cell_iter_2 = cell_population.Begin();
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_2), 0u);
        ++cell_iter_2;
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_2), 1u);

        // Check Elements are as expected
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNodeGlobalIndex(1), 3u);
    }

    void TestUpdateCellLocations()
    {
        // Create a simple 2D PottsMesh with two cells
        PottsMeshGenerator<2> generator(4, 2, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set node selection to non-random lattice sweeping: this will loop over the nodes in index order
        TS_ASSERT_EQUALS(true, cell_population.GetUpdateNodesInRandomOrder());
        cell_population.SetUpdateNodesInRandomOrder(false);

        // Increase temperature: allows swaps to be more likely
        TS_ASSERT_EQUALS(cell_population.GetTemperature(),0.1);
        cell_population.SetTemperature(10.0);

        // Create a volume update rule and pass to the population
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        cell_population.AddUpdateRule(p_volume_constraint_update_rule);

        // Commence lattice sweeping, updating where necessary
        cell_population.UpdateCellLocations(1.0);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNumNodes(), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNumNodes(), 7u);
    }

    void TestUpdateCellLocationsRandomly()
    {
        // Create a simple 2D PottsMesh with two cells
        PottsMeshGenerator<2> generator(4, 2, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set node selection to random lattice sweeping and (for coverage) random iteration over update rules
        cell_population.SetUpdateNodesInRandomOrder(true);
        cell_population.SetIterateRandomlyOverUpdateRuleCollection(true);

        // Increase temperature: allows swaps to be more likely
        TS_ASSERT_EQUALS(cell_population.GetTemperature(),0.1);
        cell_population.SetTemperature(10.0);

        // Create a volume update rule and pass to the population
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        cell_population.AddUpdateRule(p_volume_constraint_update_rule);

        // Commence lattice sweeping, updating where necessary
        cell_population.UpdateCellLocations(1.0);

        // Note that these results differ to those in the above test due to extra random numbers being called
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNumNodes(), 2u);
    }

    ///\todo implement this test (#1666)
//    void TestVoronoiMethods()
//    {
//        // Create a simple 2D PottsMesh
//        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
//        PottsMesh<2>* p_mesh = generator.GetMesh();
//
//        // Create cells
//        std::vector<CellPtr> cells;
//        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());
//
//        // Create cell population
//        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
//
//        // Create element tessellation
//        cell_population.CreateElementTessellation();
//
//        VertexMesh<2,2>* p_tessellation = cell_population.GetElementTessellation();
//
//        TS_ASSERT(p_tessellation != NULL);
//
//        TS_ASSERT_EQUALS(p_tessellation->GetNumNodes(), 25u);
//        TS_ASSERT_EQUALS(p_tessellation->GetNumElements(), 16u);
//    }

    void TestWriteResultsToFileAndOutputCellPopulationParameters()
    {
        // Resetting the maximum cell ID to zero (to account for previous tests)
        Cell::ResetMaxCellId();

        std::string output_directory = "TestPottsBasedCellPopulationWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // For coverage, label one cell
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "PottsBasedCellPopulation-2");

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVariables(true);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetNumSweepsPerTimestep(5);

        TS_ASSERT_EQUALS(cell_population.GetNumSweepsPerTimestep(), 5u);

        // VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        TS_ASSERT_THROWS_NOTHING(cell_population.CreateOutputFiles(output_directory, false));

        cell_population.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(cell_population.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     cell_based/test/data/TestPottsBasedCellPopulationWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     cell_based/test/data/TestPottsBasedCellPopulationWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     cell_based/test/data/TestPottsBasedCellPopulationWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     cell_based/test/data/TestPottsBasedCellPopulationWriters/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     cell_based/test/data/TestPottsBasedCellPopulationWriters/cellages.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat     cell_based/test/data/TestPottsBasedCellPopulationWriters/cellareas.dat").c_str()), 0);

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare output with saved files of what they should look like
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters    cell_based/test/data/TestPottsBasedCellPopulationWriters/results.parameters").c_str()), 0);
    }

    void TestNodeAndMeshMethods() throw(Exception)
    {
        // Create a Potts-based cell population
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test GetNumNodes()
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);

        // Test GetNode()
        for (unsigned index=0; index<cell_population.GetNumNodes(); index++)
        {
            Node<2>* p_node = cell_population.GetNode(index);
            TS_ASSERT_EQUALS(p_node->GetIndex(), index);

            c_vector<double, 2> node_location = p_node->rGetLocation();
            double expected_x = (double)(index%4);
            double expected_y = (double)(index>3) + (double)(index>7) + (double)(index>11);
            TS_ASSERT_DELTA(node_location[0], expected_x, 1e-3);
            TS_ASSERT_DELTA(node_location[1], expected_y, 1e-3);
        }

        // Test GetElement()
        for (unsigned index=0; index<cell_population.GetNumElements(); index++)
        {
            PottsElement<2>* p_element = cell_population.GetElement(index);
            TS_ASSERT_EQUALS(p_element->GetIndex(), index);

            TS_ASSERT_EQUALS(p_element->GetNumNodes(), 4u);
        }

        // Test SetNode()
        ChastePoint<2> unused_point;
        TS_ASSERT_THROWS_THIS(cell_population.SetNode(0, unused_point),
                              "SetNode() cannot be called on a subclass of AbstractOnLatticeCellPopulation.");

        // Test GetWidth() method
        double width_x = cell_population.GetWidth(0);
        TS_ASSERT_DELTA(width_x, 3.0, 1e-4);

        double width_y = cell_population.GetWidth(1);
        TS_ASSERT_DELTA(width_y, 3.0, 1e-4);

        // Test GetNeighbouringNodeIndices() method
        TS_ASSERT_THROWS_THIS(cell_population.GetNeighbouringNodeIndices(10),
            "Cannot call GetNeighbouringNodeIndices() on a PottsBasedCellPopulation, need to go through the PottsMesh instead");
    }

    void TestGetLocationOfCellCentre() throw (Exception)
    {
        // Create a Potts-based cell population
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test GetLocationOfCellCentre()
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();

        for (unsigned i=0; i<4; i++)
        {
            c_vector<double, 2> cell_1_location = cell_population.GetLocationOfCellCentre(*cell_iter);

            double x = 0.5 + 2*(i%2 != 0);
            double y = 0.5 + 2*(i > 1);

            TS_ASSERT_DELTA(cell_1_location[0], x, 1e-6);
            TS_ASSERT_DELTA(cell_1_location[1], y, 1e-6);

            ++cell_iter;
        }
    }

    void TestAddingUpdateRules() throw(Exception)
    {
        // Create a simple 2D PottsMesh with one cell
        PottsMeshGenerator<2> generator(2, 1, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);

        // Create an update rule and pass to the population
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        cell_population.AddUpdateRule(p_volume_constraint_update_rule);

        // Check the update rules are correct
        std::vector<boost::shared_ptr<AbstractPottsUpdateRule<2> > > update_rule_collection = cell_population.rGetUpdateRuleCollection();
        TS_ASSERT_EQUALS(update_rule_collection.size(),1u);
        TS_ASSERT_EQUALS((*update_rule_collection[0]).GetIdentifier(), "VolumeConstraintPottsUpdateRule-2");
    }

    void TestArchiving() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "potts_cell_population_2d.arch";

        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("potts_mesh_2d");

        // Create mesh
        PottsMeshReader<2> mesh_reader("cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d");
        PottsMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Archive cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a Potts-based cell population object
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumElements());

            // Create cell population
            AbstractCellPopulation<2>* const p_cell_population = new PottsBasedCellPopulation<2>(mesh, cells);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create an update rule and pass to the population
            MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
            static_cast<PottsBasedCellPopulation<2>*>(p_cell_population)->AddUpdateRule(p_volume_constraint_update_rule);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Set member variables in order to test that they are archived correctly
            static_cast<PottsBasedCellPopulation<2>*>(p_cell_population)->SetTemperature(0.25);
            static_cast<PottsBasedCellPopulation<2>*>(p_cell_population)->SetNumSweepsPerTimestep(3);
            static_cast<PottsBasedCellPopulation<2>*>(p_cell_population)->SetUpdateNodesInRandomOrder(false);
            static_cast<PottsBasedCellPopulation<2>*>(p_cell_population)->SetIterateRandomlyOverUpdateRuleCollection(true);

            // Archive the cell population
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        // Restore cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            AbstractCellPopulation<2>* p_cell_population;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore the cell population
            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_cell_population;

            // Test that the member variables have been archived correctly
            PottsBasedCellPopulation<2>* p_static_population = static_cast<PottsBasedCellPopulation<2>*>(p_cell_population);
            TS_ASSERT_DELTA(p_static_population->GetTemperature(), 0.25, 1e-6);
            TS_ASSERT_EQUALS(p_static_population->GetNumSweepsPerTimestep(), 3u);
            TS_ASSERT_EQUALS(p_static_population->GetUpdateNodesInRandomOrder(), false);
            TS_ASSERT_EQUALS(p_static_population->GetIterateRandomlyOverUpdateRuleCollection(), true);

            // Test that the update rule has been archived correctly
            std::vector<boost::shared_ptr<AbstractPottsUpdateRule<2> > > update_rule_collection = p_static_population->rGetUpdateRuleCollection();
            TS_ASSERT_EQUALS(update_rule_collection.size(), 1u);
            TS_ASSERT_EQUALS((*update_rule_collection[0]).GetIdentifier(), "VolumeConstraintPottsUpdateRule-2");

            // Tidy up
            delete p_cell_population;
        }
    }
};

#endif /*TESTPOTTSBASEDCELLPOPULATION_HPP_*/