#ifndef TESTCRYPTSIMULATION2DPERIODIC_HPP_
#define TESTCRYPTSIMULATION2DPERIODIC_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "SloughingCellKiller.hpp"

// Possible types of Cell Cycle Model (just for CreateVectorOfCells method)
typedef enum CellCycleType_
{
    FIXED,
    STOCHASTIC,
    WNT,
    TYSONNOVAK
} CellCycleType;


class TestCryptSimulation2DPeriodic : public CxxTest::TestSuite
{
    void CreateVectorOfCells(std::vector<MeinekeCryptCell>& rCells, 
                             ConformingTetrahedralMesh<2,2>& rMesh, 
                             CellCycleType cycleType, 
                             bool randomBirthTimes,
                             double y0 = 0.3,
                             double y1 = 2.0,
                             double y2 = 3.0,
                             double y3 = 4.0)
    {
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        unsigned num_cells = rMesh.GetNumNodes();

        AbstractCellCycleModel* p_cell_cycle_model = NULL;
        double typical_transit_cycle_time;
        double typical_stem_cycle_time;
        
        CancerParameters* p_params = CancerParameters::Instance();
        
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;

            double y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (cycleType==FIXED)
            {
                p_cell_cycle_model = new FixedCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==STOCHASTIC)
            {
                p_cell_cycle_model = new StochasticCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==WNT)
            {
                WntGradient wnt_gradient(LINEAR);
                double wnt = wnt_gradient.GetWntLevel(y);
                p_cell_cycle_model = new WntCellCycleModel(wnt,0);
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==TYSONNOVAK)
            {
                p_cell_cycle_model = new TysonNovakCellCycleModel();
                typical_transit_cycle_time = 1.25;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else
            {
                EXCEPTION("Cell Cycle Type is not recognised");   
            }
            
            
            double birth_time = 0.0;
            
            if (y <= y0)
            {
                cell_type = STEM;
                generation = 0;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_stem_cycle_time; // hours
                }
            }
            else if (y < y1)
            {
                cell_type = TRANSIT;
                generation = 1;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y2)
            {
                cell_type = TRANSIT;
                generation = 2;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y3)
            {
                cell_type = TRANSIT;
                generation = 3;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else
            {
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
                if(cycleType==WNT || cycleType==TYSONNOVAK)
                {
                    // There are no fully differentiated cells!
                    cell_type = TRANSIT;
                    
                }
                else
                {
                    cell_type = DIFFERENTIATED;
                }                
                generation = 4;
            }

            MeinekeCryptCell cell(cell_type, HEALTHY, generation, p_cell_cycle_model);
            
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            rCells.push_back(cell);
        }
    }

    /**
     * Compare 2 meshes to see if they are 'the same'.  Doesn't check everything,
     * but is fairly thorough.  Used for testing serialization.
     */
    template<unsigned DIM>
    void CompareMeshes(ConformingTetrahedralMesh<DIM,DIM>* pMesh1,
                       ConformingTetrahedralMesh<DIM,DIM>* pMesh2)
    {
        TS_ASSERT_EQUALS(pMesh1->GetNumAllNodes(), pMesh2->GetNumAllNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumNodes(), pMesh2->GetNumNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumBoundaryNodes(), pMesh2->GetNumBoundaryNodes());
        TS_ASSERT_EQUALS(pMesh1->GetNumCornerNodes(), pMesh2->GetNumCornerNodes());
        
        for (unsigned i=0; i<pMesh1->GetNumAllNodes(); i++)
        {
            Node<DIM> *p_node = pMesh1->GetNode(i);
            Node<DIM> *p_node2 = pMesh2->GetNode(i);
            TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
            TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
            TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
            for (unsigned j=0; j<DIM; j++)
            {
                TS_ASSERT_DELTA(p_node->rGetLocation()[j], p_node2->rGetLocation()[j], 1e-16);
            }
        }
        
        TS_ASSERT_EQUALS(pMesh1->GetNumElements(), pMesh2->GetNumElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumAllElements(), pMesh2->GetNumAllElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumBoundaryElements(), pMesh2->GetNumBoundaryElements());
        TS_ASSERT_EQUALS(pMesh1->GetNumAllBoundaryElements(), pMesh2->GetNumAllBoundaryElements());
        typename ConformingTetrahedralMesh<DIM,DIM>::ElementIterator it=pMesh1->GetElementIteratorBegin();
        typename ConformingTetrahedralMesh<DIM,DIM>::ElementIterator it2=pMesh2->GetElementIteratorBegin();
        for (;
             it != pMesh1->GetElementIteratorEnd();
             ++it, ++it2)
        {
            Element<DIM,DIM>* p_elt = *it;
            Element<DIM,DIM>* p_elt2 = *it2;
            TS_ASSERT_EQUALS(p_elt->GetNumNodes(), p_elt2->GetNumNodes());
            for (unsigned i=0; i<p_elt->GetNumNodes(); i++)
            {
                TS_ASSERT_EQUALS(p_elt->GetNodeGlobalIndex(i), p_elt2->GetNodeGlobalIndex(i));
            }
        }
    }

public:
    void Test2DCylindrical() throw (Exception)
    {        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);// true = mature cells

        Crypt<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetEndTime(0.1);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxCells(10));
        simulator.SetMaxCells(500);
        TS_ASSERT_THROWS_ANYTHING(simulator.SetMaxElements(10));
        simulator.SetMaxElements(1000);
        
        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(false);
        simulator.SetOutputDirectory("Crypt2DCylindrical");
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        simulator.Solve();

        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_nodes, ghost_cells.size());
        TS_ASSERT_EQUALS(number_of_cells, cells_across*cells_up+1u);  // 6 cells in a row*12 rows + 1 birth
        TS_ASSERT_EQUALS(number_of_nodes, number_of_cells+thickness_of_ghost_layer*2*cells_across); 

        // Coverage of exceptions (after main test to avoid problems with SimulationTime).
        simulator.SetEndTime(10.0);
        simulator.SetOutputDirectory("");
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());

        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    // This is a rubbish test - all cells start at birthTime = 0.
    // So bizarrely the crypt shrinks as the rest lengths are shortened! 
    // But at least it uses Wnt cell cycle and runs reasonably quickly...
    // For a better test with more randomly distributed cell ages see the Nightly test pack.
    void TestWithWntDependentCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 0;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells 
        std::vector<MeinekeCryptCell> cells;        
        CreateVectorOfCells(cells, *p_mesh, WNT, false);

        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicWnt");
        
        // Set length of simulation here
        simulator.SetEndTime(0.3);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        simulator.SetWntGradient(LINEAR);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
        
        simulator.Solve();
        
        std::vector<double> node_35_location = simulator.GetNodeLocation(35);
        //std::vector<double> node_100_location = simulator.GetNodeLocation(100);
        
        TS_ASSERT_DELTA(node_35_location[0], 5.5000 , 1e-4);
        TS_ASSERT_DELTA(node_35_location[1], 4.4000 , 1e-4);
//        TS_ASSERT_DELTA(node_100_location[0], 4.0000 , 1e-4);
//        TS_ASSERT_DELTA(node_100_location[1], 8.0945 , 1e-4);
//          
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    // A better check that the loaded mesh is the same as that saved
    void TestMeshSurvivesSaveLoad() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        
        // Set up a simulation
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, WNT, false);
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DPeriodicWntSaveAndLoad");
        simulator.SetEndTime(0.1);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // Save
        simulator.Save();
        
        // Load
        TissueSimulation<2>* p_simulator;
        p_simulator = TissueSimulation<2>::Load("Crypt2DPeriodicWntSaveAndLoad", 0.0);
        
        // Create an identical mesh for comparison purposes
        HoneycombMeshGenerator generator2(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh2=generator2.GetCylindricalMesh();
        
        // Compare
        CompareMeshes(p_mesh2, &(p_simulator->rGetCrypt().rGetMesh()));
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    // Testing Save (based on previous test)
    void TestSave() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000); 
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, WNT, false);
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory("Crypt2DPeriodicWntSaveAndLoad");
        
        // Our full end time is 0.2, here we run for half the time
        simulator.SetEndTime(0.1);
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        simulator.SetWntGradient(LINEAR);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        simulator.Solve();
        
        // save the results..
        simulator.Save();
        
        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception) 
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 0.2
        TissueSimulation<2>* p_simulator1;
        p_simulator1 = TissueSimulation<2>::Load("Crypt2DPeriodicWntSaveAndLoad", 0.1);
        
        p_simulator1->SetEndTime(0.2);
        p_simulator1->Solve();
        
        // save that then reload
        // and run from 0.2 to 0.3.
        p_simulator1->Save();
        
        TissueSimulation<2>* p_simulator2 = TissueSimulation<2>::Load("Crypt2DPeriodicWntSaveAndLoad", 0.2);
        
        CompareMeshes(&(p_simulator1->rGetCrypt().rGetMesh()),
                      &(p_simulator2->rGetCrypt().rGetMesh()));
        
        p_simulator2->SetEndTime(0.3);
        p_simulator2->Solve();
        
        /* 
         * This checks that these two nodes are in exactly the same location 
         * (after a saved and loaded run) as after a single run
         */
        std::vector<double> node_35_location = p_simulator2->GetNodeLocation(35);
        TS_ASSERT_DELTA(node_35_location[0], 5.5000 , 1e-4);
        TS_ASSERT_DELTA(node_35_location[1], 2.5104 , 1e-4);
        std::vector<double> node_100_location = p_simulator2->GetNodeLocation(100);
        TS_ASSERT_DELTA(node_100_location[0], 4.0000 , 1e-4);
        TS_ASSERT_DELTA(node_100_location[1], 8.0945 , 1e-4);
        
        delete p_simulator1;
        delete p_simulator2;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    /* 
     * Cells are compressed and are trying to spread out. So the cells 
     * at the bottom try to push across y=0 but are prevented from doing so.
     * This is not covered in previous tests because cells are all at age = 0
     * and try to come together to simulate birth.
     * 
     * It is potentially an expensive test computationally, because all 
     * Wnt cells have to run cell cycle models for a large time 
     * to be 'mature' cells which won't shrink together. 
     * Limited this by using only four cells of minimum age.
     */
    void TestWntCellsCannotMoveAcrossYEqualsZero() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in Wnt simulation
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 2;
        unsigned cells_up = 2;
        double crypt_width = 0.5;
        unsigned thickness_of_ghost_layer = 0;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, false, crypt_width/cells_across);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, WNT, true);
        cells[0].SetBirthTime(-1.0);   // Make cell cycle models do minimum work
        cells[1].SetBirthTime(-1.0);
        cells[1].SetMutationState(LABELLED);
        cells[2].SetBirthTime(-1.0);
        cells[2].SetMutationState(APC_ONE_HIT);
        cells[3].SetBirthTime(-1.0);
        cells[3].SetMutationState(BETA_CATENIN_ONE_HIT);

        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        TissueSimulation<2> simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DWntMatureCells");
        // If you want to visualize this use the 'notcylindrical' option
        // (it is too small for it to figure out what's happening on its own).
        
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetWntGradient(LINEAR);
        simulator.SetEndTime(0.01);
                
        simulator.Solve();
        
        // Check that nothing has moved below y=0
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_LESS_THAN(-1e-15,p_mesh->GetNode(cell_iter->GetNodeIndex())->rGetLocation()[1]);
        }
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    
    // This is strange test -- all cells divide within a quick time, it gives
    // good testing of the periodic boundaries though... [comment no longer valid?]
    void TestWithTysonNovakCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        // There is no limit on transit cells in T&N
        p_params->SetMaxTransitGenerations(1000);
        
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        unsigned thickness_of_ghost_layer = 4;
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, TYSONNOVAK, true);
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);
        
        TissueSimulation<2> simulator(crypt); 
               
        AbstractCellKiller<2>* p_cell_killer = new SloughingCellKiller(&crypt);
        simulator.AddCellKiller(p_cell_killer);
                
        simulator.SetOutputDirectory("Crypt2DPeriodicTysonNovak");
        simulator.SetEndTime(0.05);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        simulator.SetDt(0.001);
        
        // Test that labelling a few cells doesn't make any difference to the simulation
        // and therefore log them in the visualizer files for the next test to check.
        simulator.rGetCrypt().rGetCellAtNodeIndex(57).SetMutationState(LABELLED);
        simulator.rGetCrypt().rGetCellAtNodeIndex(56).SetMutationState(APC_ONE_HIT);
        simulator.rGetCrypt().rGetCellAtNodeIndex(51).SetMutationState(APC_TWO_HIT);
        simulator.rGetCrypt().rGetCellAtNodeIndex(63).SetMutationState(BETA_CATENIN_ONE_HIT);
                
        simulator.Solve();
        
        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        std::vector<bool> ghost_cells = crypt.rGetGhostNodes();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
        
        TS_ASSERT_EQUALS(number_of_nodes,ghost_cells.size());
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 81u);
        TS_ASSERT_EQUALS(number_of_nodes, 123u);

        delete p_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    /*
     * This test compares the visualizer output from the previous test with a known file.
     * 
     * Note - if the previous test is changed we need to update the file this test refers to. 
     */
    void TestVisualizerOutput() throw (Exception)
    {
        // work out where the previous test wrote its files
        OutputFileHandler handler("Crypt2DPeriodicTysonNovak",false);
        std::string results_dir = handler.GetTestOutputDirectory() + "results_from_time_0/vis_results";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.vizelements models/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.viznodes models/test/data/Crypt2DPeriodicTysonNovak_vis/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.vizsetup models/test/data/Crypt2DPeriodicTysonNovak_vis/results.vizsetup").c_str()), 0);
    }
   
    
    void TestPrivateFunctionsOf2DCryptSimulation() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator::Instance();
        
        double crypt_length = 9.3;
        double crypt_width = 10.0;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);        
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, mesh, FIXED, false, 0.0, 3.0, 6.5, 8.0);
        
        cells[60].SetBirthTime(-50.0);
        
        Crypt<2> crypt(mesh, cells);
        
        TissueSimulation<2> simulator(crypt);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);
        
        unsigned num_deaths = simulator.DoCellRemoval();
        unsigned num_births = simulator.DoCellBirth();
                                                                
        TS_ASSERT_EQUALS(num_births, 1u);
        TS_ASSERT_EQUALS(num_deaths,11u);
               
        delete p_sloughing_cell_killer;       
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    
    void TestPrivateFunctionsOf2DCryptSimulationOnHoneycombMesh() throw (Exception)
    {
        /*
         ************************************************************************
         ************************************************************************ 
         *     Set up a simulation class to run the individual tests on.
         ************************************************************************
         ************************************************************************ 
         */
        unsigned cells_across = 7;
        unsigned cells_up = 5;
        unsigned thickness_of_ghost_layer = 3;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer,false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        unsigned num_cells = p_mesh->GetNumAllNodes();
        
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator::Instance();
        

        
        // Set up cells by iterating through the mesh nodes
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            double birth_time;
            CryptCellType cell_type;
            unsigned generation;
            double y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -2.0; //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -2.0; //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -2.0;  //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -2.0;  //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = -2.0;  //hours
            }
            
            CryptCellMutationState mutation_state;
            if(i!=60)
            {
                mutation_state = HEALTHY;
            }
            else
            {
                mutation_state = APC_TWO_HIT;
            }
            
            WntGradient wnt_gradient(LINEAR);
            double wnt = wnt_gradient.GetWntLevel(y);
            
            MeinekeCryptCell cell(cell_type, mutation_state, generation, new WntCellCycleModel(wnt,0));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(400);

       /*
        ************************************************************************
        ************************************************************************ 
        *  Test Calculate Velocities on each node
        ************************************************************************
        ************************************************************************ 
        */
                
        std::vector<c_vector<double, 2> > velocities_on_each_node(p_mesh->GetNumAllNodes());
        
        velocities_on_each_node = simulator.CalculateVelocitiesOfEachNode();
 
        for (unsigned i=0; i<p_mesh->GetNumAllNodes(); i++)
        {
            std::set<unsigned>::iterator iter = ghost_node_indices.find(i);
            bool is_a_ghost_node = (iter!=ghost_node_indices.end());

            if (!is_a_ghost_node)
            {
                TS_ASSERT_DELTA(velocities_on_each_node[i][0], 0.0, 1e-4);
                TS_ASSERT_DELTA(velocities_on_each_node[i][1], 0.0, 1e-4);
            }
        }
        
        // Move a node along the x-axis and calculate the force exerted on a neighbour
        c_vector<double,2> old_point = p_mesh->GetNode(59)->rGetLocation();
        Point<2> new_point;
        new_point.rGetLocation()[0] = old_point[0]+0.5;
        new_point.rGetLocation()[1] = old_point[1];
  
        p_mesh->SetNode(59, new_point, false);
        velocities_on_each_node = simulator.CalculateVelocitiesOfEachNode();
        
        TS_ASSERT_DELTA(velocities_on_each_node[60][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantMutant(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[60][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[59][0], (-3+4.0/sqrt(7))*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[59][1], 0.0, 1e-4);

        TS_ASSERT_DELTA(velocities_on_each_node[58][0], 0.5*p_params->GetSpringStiffness()/p_params->GetDampingConstantNormal(), 1e-4);
        TS_ASSERT_DELTA(velocities_on_each_node[58][1], 0.0, 1e-4);
        
        
       /*
        ************************************************************************
        ************************************************************************ 
        *  Test Calculate force on a spring
        ************************************************************************
        ************************************************************************ 
        */
        
        c_vector<double,2> force_on_spring ; // between nodes 59 and 60
        
        // Find one of the elements that nodes 59 and 60 live on
        Point<2> new_point2;
        new_point2.rGetLocation()[0] = new_point[0] + 0.01;
        new_point2.rGetLocation()[1] = new_point[1] + 0.01 ;
        
        unsigned elem_index = p_mesh->GetContainingElementIndex(new_point2,false);
        Element<2,2>* p_element = p_mesh->GetElement(elem_index);
        
        force_on_spring = simulator.CalculateForceBetweenNodes(p_element->GetNodeGlobalIndex(1),p_element->GetNodeGlobalIndex(0));
        
        TS_ASSERT_DELTA(force_on_spring[0], 0.5*p_params->GetSpringStiffness(), 1e-4);
        TS_ASSERT_DELTA(force_on_spring[1], 0.0, 1e-4);
        
        
       /*
        ************************************************************************
        ************************************************************************ 
        *  Test UpdateNodePositions
        ************************************************************************
        ************************************************************************ 
        */
        
        Point<2> point_of_node60 = p_mesh->GetNode(60)->rGetLocation();
        
        simulator.SetDt(0.01);
        simulator.UpdateNodePositions(velocities_on_each_node);
        
        TS_ASSERT_DELTA(p_mesh->GetNode(60)->rGetLocation()[0],point_of_node60.rGetLocation()[0]+force_on_spring[0]/p_params->GetDampingConstantMutant() *0.01, 1e-4);
        TS_ASSERT_DELTA(p_mesh->GetNode(60)->rGetLocation()[1],point_of_node60.rGetLocation()[1], 1e-4);
        
       /*
        ************************************************************************
        ************************************************************************ 
        * Test UpdateCellTypes 
        ************************************************************************
        ************************************************************************ 
        */
        simulator.SetWntGradient(LINEAR);
        simulator.UpdateCellTypes();
        
        
        for (Crypt<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            CryptCellType cell_type;
            cell_type = cell_iter->GetCellType();
            if (!cell_type==STEM)
            {
                //std::cout << "Cell type = " << cell_type << std::endl;
                WntCellCycleModel *p_this_model = static_cast<WntCellCycleModel*>(cell_iter->GetCellCycleModel());
                double beta_cat_level = p_this_model->GetProteinConcentrations()[6]+ p_this_model->GetProteinConcentrations()[7];
                //std::cout << "Cell " << i << ", beta-cat = " << beta_cat_level << std::endl;
                if (beta_cat_level > 0.4127)
                {
                    TS_ASSERT_EQUALS(cell_type,TRANSIT);
                }
                else
                {
                    TS_ASSERT_EQUALS(cell_type,DIFFERENTIATED);
                }
            }
        
        }
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestCalculateDividingCellCentreLocationsConfMesh() throw (Exception)
    {
        double separation = 0.1;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a parent node
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);

        ConformingTetrahedralMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> conf_cells;
        CreateVectorOfCells(conf_cells, conf_mesh, TYSONNOVAK, true);

        Crypt<2> conf_crypt(conf_mesh, conf_cells);

        Crypt<2>::Iterator conf_iter = conf_crypt.Begin();

        TissueSimulation<2> simulator(conf_crypt);
                     
        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(conf_iter);
        c_vector<double, 2> new_parent_location = conf_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = conf_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), separation, 1e-7);

        SimulationTime::Destroy();
    }
    

    void TestCalculateDividingCellCentreLocationsConfMeshStemCell() throw (Exception)
    {
        double separation = 0.1;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a parent node
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=0.0;                                    // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        ConformingTetrahedralMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> conf_cells;
        CreateVectorOfCells(conf_cells, conf_mesh, TYSONNOVAK, true);        
        Crypt<2> conf_crypt(conf_mesh, conf_cells);

        Crypt<2>::Iterator conf_iter = conf_crypt.Begin();

        TissueSimulation<2> simulator(conf_crypt);
        
        // repeat two times for coverage
        // need vector from parent to daughter to have both +ve and -ve y component
        // different branches will execute to make sure daughter stays in crypt ie. +ve y component
        for (unsigned repetitions=0; repetitions<=1; repetitions++)
        {                
            c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(conf_iter);
            c_vector<double, 2> new_parent_location = conf_mesh.GetNode(0)->rGetLocation();
            c_vector<double, 2> parent_to_daughter = conf_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
    
            // The parent stem cell should stay where it is and the daughter be introduced at positive y.
            TS_ASSERT_DELTA(new_parent_location[0], location[0], 1e-7);
            TS_ASSERT_DELTA(new_parent_location[1], location[1], 1e-7);
            TS_ASSERT(daughter_location[1]>=location[1]);
            TS_ASSERT_DELTA(norm_2(parent_to_daughter), 0.5*separation, 1e-7);
        }

        SimulationTime::Destroy();
    }

    void TestCalculateDividingCellCentreLocationsCylindricalMesh() throw (Exception)
    {
        double separation = 0.1;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a mesh
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cyl_cells;
        CreateVectorOfCells(cyl_cells, cyl_mesh, TYSONNOVAK, true);        
        Crypt<2> cyl_crypt(cyl_mesh, cyl_cells);

        Crypt<2>::Iterator cyl_iter = cyl_crypt.Begin();

        TissueSimulation<2> simulator(cyl_crypt);                
        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(cyl_iter);
        c_vector<double, 2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), separation, 1e-7);

        SimulationTime::Destroy();
    }

    void TestCalculateDividingCellCentreLocationsCylindricalMeshStemCell() throw (Exception)
    {
        double separation = 0.1;

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Make a mesh
        c_vector<double ,2> location;
        location[0]=1.0;
        location[1]=0.0;                                    // <- y=0
        Node<2>* p_node = new Node<2>(0u,location, false);
        Cylindrical2dMesh cyl_mesh(6.0);
        cyl_mesh.AddNode(p_node);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cyl_cells;
        CreateVectorOfCells(cyl_cells, cyl_mesh, TYSONNOVAK, true);        
        Crypt<2> cyl_crypt(cyl_mesh, cyl_cells);

        Crypt<2>::Iterator cyl_iter = cyl_crypt.Begin();

        TissueSimulation<2> simulator(cyl_crypt);                
        c_vector<double, 2> daughter_location = simulator.CalculateDividingCellCentreLocations(cyl_iter);
        c_vector<double, 2> new_parent_location = cyl_mesh.GetNode(0)->rGetLocation();
        c_vector<double, 2> parent_to_daughter = cyl_mesh.GetVectorFromAtoB(new_parent_location, daughter_location);
        
        // The parent stem cell should stay where it is and the daughter be introduced at positive y.
        TS_ASSERT_DELTA(new_parent_location[0], location[0], 1e-7);
        TS_ASSERT_DELTA(new_parent_location[1], location[1], 1e-7);
        TS_ASSERT(daughter_location[1]>=location[1]);
        TS_ASSERT_DELTA(norm_2(parent_to_daughter), 0.5*separation, 1e-7);

        SimulationTime::Destroy();
    }
    
    // short test which sets mNoBirth for coverage
    void TestNoBirth() throw (Exception)
    {
        std::string output_directory = "Crypt2DCylindricalNoBirth";        
        unsigned cells_across = 2;
        unsigned cells_up = 3;
        double crypt_width = 2.0;
        unsigned thickness_of_ghost_layer = 0;
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, true, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);// true = mature cells

        Crypt<2> crypt(*p_mesh, cells);               
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(2.0); // long enough for a cell to be born were SetNoBirth not called
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);
        
        // These are for coverage and use the defaults
        simulator.SetDt(1.0/120.0);
        simulator.SetReMeshRule(true);
        simulator.SetNoBirth(true);
        
        AbstractCellKiller<2>* p_sloughing_cell_killer = new SloughingCellKiller(&crypt, true);
        simulator.AddCellKiller(p_sloughing_cell_killer);

        simulator.Solve();

        // test we have the same number of cells and nodes at the end of each time
        // (if we do then the boundaries are probably working!)
        unsigned number_of_cells = crypt.GetNumRealCells();
        unsigned number_of_nodes = crypt.rGetMesh().GetNumNodes();
                
        TS_ASSERT_EQUALS(number_of_cells, cells_across*cells_up); 
        TS_ASSERT_EQUALS(number_of_nodes, number_of_cells+thickness_of_ghost_layer*2*cells_across); 

        delete p_sloughing_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    

    // Death on a non-periodic mesh
    // Massive amount of random death eventually leading to every cell being killed off..
    // Note that birth does occur too.
    void TestRandomDeathOnNonPeriodicCrypt() throw (Exception)
    {
        unsigned cells_across = 2;
        unsigned cells_up = 1;
        unsigned thickness_of_ghost_layer = 1;
        
        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        CreateVectorOfCells(cells, *p_mesh, FIXED, true);
              
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);
        
        simulator.SetOutputDirectory("Crypt2DRandomDeathNonPeriodic");
        simulator.SetEndTime(0.5);
        simulator.SetMaxCells(500);
        simulator.SetMaxElements(1000);

        AbstractCellKiller<2>* p_random_cell_killer = new RandomCellKiller<2>(&crypt, 0.1);
        simulator.AddCellKiller(p_random_cell_killer);

        simulator.Solve();
        
        // there should be no cells left after this amount of time
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 0u);
    
        delete p_random_cell_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCRYPTSIMULATION2DPERIODIC_HPP_*/
