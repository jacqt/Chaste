#ifndef TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_
#define TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CryptSimulation2DPeriodic.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "CryptHoneycombMeshGenerator.hpp"

class TestCryptSimulation2DPeriodicNightly : public CxxTest::TestSuite
{
	void CheckAgainstPreviousRun(std::string resultDirectory, unsigned maxCells, unsigned maxElements)
    {
        std::cout << "Comparing " << resultDirectory << std::endl << std::flush;
        
        ColumnDataReader computed_node_results = ColumnDataReader(resultDirectory+"Results",
                                                             "tabulated_node_results",
                                                             true);
                                                             
        ColumnDataReader expected_node_results = ColumnDataReader("models/test/data/" + resultDirectory+"Results",
                                                             "tabulated_node_results",
                                                             false);
        ColumnDataReader computed_element_results = ColumnDataReader(resultDirectory+"Results",
                                                             "tabulated_element_results",
                                                             true);
                                                             
        ColumnDataReader expected_element_results = ColumnDataReader("models/test/data/" + resultDirectory+"Results",
                                                             "tabulated_element_results",
                                                             false);
                                                                                                                  
        for (unsigned cell=0; cell<maxCells; cell++)
        {
            std::stringstream cell_type_var_name;
            std::stringstream cell_x_position_var_name;
            std::stringstream cell_y_position_var_name;
            cell_type_var_name << "cell_type_" << cell;
            cell_x_position_var_name << "cell_x_position_" << cell;
            cell_y_position_var_name << "cell_y_position_" << cell;
            
            // Vector of Cell Types
            std::vector<double> expected_cell_types = expected_node_results.GetValues(cell_type_var_name.str());
            std::vector<double> computed_cell_types = computed_node_results.GetValues(cell_type_var_name.str());
            
            //Vector of Cell Positions
            std::vector<double> expected_cell_x_positions = expected_node_results.GetValues(cell_x_position_var_name.str());
            std::vector<double> computed_cell_x_positions = computed_node_results.GetValues(cell_x_position_var_name.str());
            
            std::vector<double> expected_cell_y_positions = expected_node_results.GetValues(cell_y_position_var_name.str());
            std::vector<double> computed_cell_y_positions = computed_node_results.GetValues(cell_y_position_var_name.str());
            
            //Comparing expected and computed vector length
            TS_ASSERT_EQUALS(expected_cell_types.size(), computed_cell_types.size());
            TS_ASSERT_EQUALS(expected_cell_x_positions.size(), computed_cell_x_positions.size());
            TS_ASSERT_EQUALS(expected_cell_y_positions.size(), computed_cell_y_positions.size());
            
            //Walkthrough of the expected and computed vectors
            for (unsigned time_step = 0; time_step < expected_cell_types.size(); time_step++)
            {
                TS_ASSERT_EQUALS(expected_cell_types[time_step], computed_cell_types[time_step]);
                TS_ASSERT_DELTA(expected_cell_x_positions[time_step], computed_cell_x_positions[time_step],1e-6);
                TS_ASSERT_DELTA(expected_cell_y_positions[time_step], computed_cell_y_positions[time_step],1e-6);
            }
        }
        
        for (unsigned element=0; element<maxElements; element++)
        {
            std::stringstream nodeA_var_name;
            std::stringstream nodeB_var_name;
            std::stringstream nodeC_var_name;
            nodeA_var_name << "nodeA_" << element;
            nodeB_var_name << "nodeB_" << element;
            nodeC_var_name << "nodeC_" << element;
            
            // Vector of Node A indexes
            std::vector<double> expected_NodeA_numbers = expected_element_results.GetValues(nodeA_var_name.str());
            std::vector<double> computed_NodeA_numbers = computed_element_results.GetValues(nodeA_var_name.str());
            
            // Vector of Node B indexes
            std::vector<double> expected_NodeB_numbers = expected_element_results.GetValues(nodeB_var_name.str());
            std::vector<double> computed_NodeB_numbers = computed_element_results.GetValues(nodeB_var_name.str());
            
            // Vector of Node C indexes
            std::vector<double> expected_NodeC_numbers = expected_element_results.GetValues(nodeC_var_name.str());
            std::vector<double> computed_NodeC_numbers = computed_element_results.GetValues(nodeC_var_name.str());
            
            TS_ASSERT_EQUALS(expected_NodeA_numbers.size(), computed_NodeA_numbers.size());
            TS_ASSERT_EQUALS(expected_NodeB_numbers.size(), computed_NodeB_numbers.size());
            TS_ASSERT_EQUALS(expected_NodeC_numbers.size(), computed_NodeC_numbers.size());
            
            for (unsigned time_step = 0; time_step < expected_NodeA_numbers.size(); time_step++)
            {
                TS_ASSERT_EQUALS(expected_NodeA_numbers[time_step], computed_NodeA_numbers[time_step]);
                TS_ASSERT_EQUALS(expected_NodeB_numbers[time_step], computed_NodeB_numbers[time_step]);
                TS_ASSERT_EQUALS(expected_NodeC_numbers[time_step], computed_NodeC_numbers[time_step]);
            }
        }
    }
    
    
    
public:

    void TestWithBirthOnHoneycombMeshNotPeriodic() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
        int num_cells_depth = 11;
        int num_cells_width = 6;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;

        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth); 
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<int> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        // Any old rubbish here just so the simulation time is set up to set up cell cycle models
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
          
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -random_num_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -1; //hours
            }
            
            MeinekeCryptCell cell(cell_type, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        //CryptSimulation2D simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        simulator.SetPeriodicSides(false);
        simulator.SetGhostNodes(ghost_node_indices);

		SimulationTime::Destroy();
        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DHoneycombMesh", 400u, 800u);
    }
    
    
    //////////////////////////////////////////////////////////////////
    // starting with a small mesh with 1 stem cell and the rest
    // differentiated, check the number of cells at the end of the 
    // simulation is as expected.
    //////////////////////////////////////////////////////////////////
    void Test2DCorrectCellNumbers() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
        // check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellCycleTime(), 24, 1e-12);  
        TS_ASSERT_DELTA(p_params->GetTransitCellCycleTime(), 12, 1e-12);  
        
        int num_cells_width = 6;
        int num_cells_depth = 5;

        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth); 
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::vector<int> ghost_node_indices = generator.GetGhostNodeIndices();

        double crypt_width = num_cells_width - 1.0;
        double crypt_length = num_cells_depth - 1.0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        // Any old rubbish here just so the simulation time is set up to set up cell cycle models
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            
            if (i==14) // middle of bottom row of cells
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -1; 
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -1; //hours
            }
            
            MeinekeCryptCell cell(cell_type, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation2DPeriodic simulator(*p_mesh,cells);
        simulator.SetOutputDirectory("Crypt2DSpringsCorrectCellNumbers");
        simulator.SetEndTime(40); //hours
        simulator.SetMaxCells(800);
        simulator.SetMaxElements(800);

        simulator.SetGhostNodes(ghost_node_indices);
        simulator.SetPeriodicSides(false);
        SimulationTime::Destroy();
        simulator.Solve();
        
        // now count the number of each type of cell
        std::vector<MeinekeCryptCell> cells_after_simulation = simulator.GetCells();
        std::vector<bool> is_ghost_node = simulator.GetGhostNodes();

        int num_stem = 0; 
        int num_transit = 0;
        int num_differentiated = 0;
        
        for (int index = 0; index<p_mesh->GetNumAllNodes(); index++)
        {
            if(!is_ghost_node[index])   //!mesh.GetNode(index)->IsDeleted())
            {
                CryptCellType type = cells_after_simulation[index].GetCellType();
                
                if(type==STEM)
                {
                    num_stem++;
                }
                else if(type==TRANSIT)
                {
                    num_transit++;
                }
                else if(type==DIFFERENTIATED)
                {
                    num_differentiated++;
                }
                else
                {
                    // shouldn't get here
                    TS_ASSERT(false);                
                }
            }
        }
        
        TS_ASSERT_EQUALS(num_stem, 1);
        TS_ASSERT_EQUALS(num_transit, 2);
        
        TS_ASSERT_LESS_THAN(num_differentiated, 23);
        TS_ASSERT_LESS_THAN(18, num_differentiated);
        delete p_mesh;
    } 
    
    void TestWithBirthOnHoneycombMeshPeriodic() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator random_num_gen;
        
        unsigned cells_across = 6;
		unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);  
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh(); 
        std::vector<int> ghost_node_indices = generator.GetGhostNodeIndices(); 

		double crypt_length = (double)cells_up*(sqrt(3)/2)*crypt_width/(double)cells_across;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        // Any old rubbish here just so the simulation time is set up to set up cell cycle models
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(54.0, 9);
        
     // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::cout << "Num Cells = " << num_cells << "\n";
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
          
            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -random_num_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 4)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
             {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -1; //hours
            }
            
            MeinekeCryptCell cell(cell_type, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        

        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicNightly");

		// Set length of simulation here
        simulator.SetEndTime(12.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        
        simulator.SetGhostNodes(ghost_node_indices);
        
        SimulationTime::Destroy();
        
        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DPeriodicNightly", 500u, 1000u);
    }
    
    void TestCalculateCryptBoundaries()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Translate(0.0,-2.0,0.0) ;
        //Create Vector of ghost nodes
        std::vector<int> ghost_node_indices;
        for (int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint().rGetLocation()[0];
            double y = mesh.GetNode(i)->GetPoint().rGetLocation()[1];
            if ((x<2.0)||(x>8.0)||(y>6.0)||(y<0.0))
            {
               ghost_node_indices.push_back(i);
            }
        }
        
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->SetCryptLength(6.0);
        p_params->SetCryptWidth(6.0);
        CryptSimulation2DPeriodic simulator(mesh);
        simulator.SetGhostNodes(ghost_node_indices);

        simulator.CalculateCryptBoundary();
        
        //simulator.DetectNaughtyCellsJoiningPeriodicEdges();
        
        std::vector<unsigned> calculated_boundary_nodes  = simulator.GetCryptBoundary();
        
        std::vector<unsigned> actual_boundary_nodes(24);
      
        actual_boundary_nodes[0] = 24;
        actual_boundary_nodes[1] = 25;
        actual_boundary_nodes[2] = 26;
        actual_boundary_nodes[3] = 27;
        actual_boundary_nodes[4] = 28;
        actual_boundary_nodes[5] = 29;
        actual_boundary_nodes[6] = 30;
        actual_boundary_nodes[7] = 35;
        actual_boundary_nodes[8] = 41;
        actual_boundary_nodes[9] = 46;
        actual_boundary_nodes[10] = 52;
        actual_boundary_nodes[11] = 57;
        actual_boundary_nodes[12] = 63;
        actual_boundary_nodes[13] = 68;
        actual_boundary_nodes[14] = 74;
        actual_boundary_nodes[15] = 79;
        actual_boundary_nodes[16] = 85;
        actual_boundary_nodes[17] = 90;
        actual_boundary_nodes[18] = 91;
        actual_boundary_nodes[19] = 92;
        actual_boundary_nodes[20] = 93;
        actual_boundary_nodes[21] = 94;
        actual_boundary_nodes[22] = 95;
        actual_boundary_nodes[23] = 96;
        
        TS_ASSERT_EQUALS(actual_boundary_nodes.size(),calculated_boundary_nodes.size());
        
        
        for(unsigned i=0; i<calculated_boundary_nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(actual_boundary_nodes[i],calculated_boundary_nodes[i]);
        }
        
        std::vector<unsigned> calculated_left_boundary_nodes = simulator.GetLeftCryptBoundary();
        std::vector<unsigned> calculated_right_boundary_nodes = simulator.GetRightCryptBoundary();
        
        std::vector<unsigned> actual_left_boundary_nodes(7);
        
        actual_left_boundary_nodes[0] = 24;
        actual_left_boundary_nodes[1] = 35;
        actual_left_boundary_nodes[2] = 46;
        actual_left_boundary_nodes[3] = 57;
        actual_left_boundary_nodes[4] = 68;
        actual_left_boundary_nodes[5] = 79;
        actual_left_boundary_nodes[6] = 90;
        
        std::vector<unsigned> actual_right_boundary_nodes(7);
        
        actual_right_boundary_nodes[0] = 24+6;
        actual_right_boundary_nodes[1] = 35+6;
        actual_right_boundary_nodes[2] = 46+6;
        actual_right_boundary_nodes[3] = 57+6;
        actual_right_boundary_nodes[4] = 68+6;
        actual_right_boundary_nodes[5] = 79+6;
        actual_right_boundary_nodes[6] = 90+6;
        
        //* Must uncomment and include this test
        
        TS_ASSERT_EQUALS(actual_left_boundary_nodes.size(),calculated_left_boundary_nodes.size());
        TS_ASSERT_EQUALS(actual_right_boundary_nodes.size(),calculated_right_boundary_nodes.size());
        
        
        for(unsigned i=0; i<calculated_left_boundary_nodes.size(); i++)
        {
            //std::cout<< "calculated_left_boundary_nodes "<< calculated_left_boundary_nodes[i] <<"\n" << std::flush;
            
            TS_ASSERT_EQUALS(actual_left_boundary_nodes[i],calculated_left_boundary_nodes[i]);
            TS_ASSERT_EQUALS(actual_right_boundary_nodes[i],calculated_right_boundary_nodes[i]);
        }
    
        
    }
    
    void TestWithWntDependentCells() throw (Exception)
    {
    	CancerParameters *p_params = CancerParameters::Instance();
		// There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        RandomNumberGenerator random_num_gen;
        
        unsigned cells_across = 6;
		unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;
        
        CryptHoneycombMeshGenerator generator(cells_across, cells_up, crypt_width,thickness_of_ghost_layer);  
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh(); 
        std::vector<int> ghost_node_indices = generator.GetGhostNodeIndices(); 

		double crypt_length = (double)cells_up*(sqrt(3)/2)*crypt_width/(double)cells_across;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
		SimulationTime* p_simulation_time = SimulationTime::Instance();
        // Any old rubbish here just so the simulation time is set up to set up cell cycle models
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1000);
        
     // Set up cells by iterating through the mesh nodes
        unsigned num_cells = p_mesh->GetNumAllNodes();
        std::cout << "Num Cells = " << num_cells << std::endl;
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
			double typical_wnt_cycle = 16.0;
            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours - doesn't matter for stem cell;
            }
            else if (y < 2)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 3)
            {
            	cell_type = TRANSIT;
                generation = 2;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else if (y < 4)
            {
            	cell_type = TRANSIT;
                generation = 3;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            else
             {	// There are no fully differentiated cells in a Wnt simulation!
                cell_type = TRANSIT;
                generation = 4;
                birth_time = -random_num_gen.ranf()*typical_wnt_cycle; //hours
            }
            
            double wnt = 1.0 - y/p_params->GetCryptLength();
            MeinekeCryptCell cell(cell_type, generation, new WntCellCycleModel(wnt));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        std::cout << "Cells Ready." << std::endl;

        CryptSimulation2DPeriodic simulator(*p_mesh, cells);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntNightly");

		// Set length of simulation here
        simulator.SetEndTime(12.0);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Set to re-mesh and birth
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetWntIncluded(true);
        
        simulator.SetGhostNodes(ghost_node_indices);
                
        SimulationTime::Destroy();
        
        simulator.Solve();
        CheckAgainstPreviousRun("Crypt2DPeriodicWntNightly", 500u, 1000u);
    }
    
};


#endif /*TESTCRYPTSIMULATION2DPERIODICNIGHTLY_HPP_*/
