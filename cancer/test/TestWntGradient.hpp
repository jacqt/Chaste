#ifndef TESTWNTGRADIENT_HPP_
#define TESTWNTGRADIENT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"
#include "WntConcentration.hpp"
#include "MeshBasedTissue.cpp"
#include "TissueCell.hpp"
#include "WntCellCycleModel.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestWntConcentration : public CxxTest::TestSuite
{
private:
    
    void setUp()
    {
        // Initialise singleton classes
        SimulationTime::Instance()->SetStartTime(0.0);
        CancerParameters::Instance()->Reset();
    }
    void tearDown()
    {        
        // Clear up singleton classes
        SimulationTime::Destroy();
        WntConcentration::Destroy();
    }
    
public:    
    void TestNoWntConcentration() throw(Exception)
    {
        WntConcentration* p_wnt = WntConcentration::Instance();
        p_wnt->SetType(NONE);
        
        TS_ASSERT_EQUALS(p_wnt->GetType(), NONE);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
    }
    
    void TestLinearWntConcentration() throw(Exception)
    {
        WntConcentration* p_wnt = WntConcentration::Instance();
        p_wnt->SetType(LINEAR);
        
        CancerParameters *params = CancerParameters::Instance();
                
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 1.0-height/params->GetCryptLength(), 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
    }
    
    
    void TestOffsetLinearWntConcentration() throw(Exception)
    {
        WntConcentration* p_wnt = WntConcentration::Instance();
        
        p_wnt->SetType(LINEAR);
        CancerParameters *params = CancerParameters::Instance();
        params->SetTopOfLinearWntConcentration(1.0/3.0);
        
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0 , 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        // under a third of the way up the crypt.
        params->SetCryptLength(22.0);
        height = 7.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 1.0 - height/((1.0/3.0)*params->GetCryptLength()) , 1e-9);
        // more than a third of the way up the crypt.
        height = 10.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
    }
    
    
    void TestRadialWntConcentration() throw(Exception)
    {
        WntConcentration* p_wnt = WntConcentration::Instance();
        p_wnt->SetType(RADIAL);

        CancerParameters *params = CancerParameters::Instance();
        
        // Test GetWntLevel(double) method
        double height = 100;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-4);
        
        height = -1e-12;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-4);        
        
        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0454 , 1e-4);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        
        params->SetCryptLength(22.0);
        height = 7.0;
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 0.6818, 1e-4);
                
        // Test GetWntLevel(TissueCell*) method
                
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // translate mesh so that its centre is at (0,0)
        mesh.Translate(-0.5,-0.5); 
        
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_model);
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);        
        }
        
        // create a crypt
        MeshBasedTissue<2> crypt(mesh,cells);        
        CancerParameters::Instance()->SetCryptLength(1.0);
        p_wnt->SetTissue(crypt);
        
        MeshBasedTissue<2>::Iterator cell_iter = crypt.Begin();
        
        double wnt_gradient_at_cell0 = p_wnt->GetWntLevel(&(*cell_iter));
        
        while(cell_iter!=crypt.End())
        {
            TS_ASSERT_DELTA(p_wnt->GetWntLevel(&(*cell_iter)), wnt_gradient_at_cell0, 1e-12);
            ++cell_iter;
        }                
    }
    
    void TestArchiveWntConcentration()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "wnt_grad.arch";
        
        // Create an ouput archive
        {            
            WntConcentration::Instance()->SetType(LINEAR);
            
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            
            output_arch << static_cast<const WntConcentration&>(*WntConcentration::Instance());
            
            WntConcentration::Destroy();
        }
        
        {
            WntConcentration* p_wnt = WntConcentration::Instance();
            
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            
            // Restore from the archive
            input_arch >> *p_wnt;
           
            double height = 21.0;
            double wnt_level = p_wnt->GetWntLevel(height);
            
            TS_ASSERT_DELTA(wnt_level, 1.0-height/CancerParameters::Instance()->GetCryptLength(), 1e-9);
        }
    }
    
    
    void TestSingletonnessOfWntConcentration()
    {
        CancerParameters *params = CancerParameters::Instance();
        
        WntConcentration* p_wnt = WntConcentration::Instance();
        p_wnt->SetType(NONE);
        
        double height = 5;
        double wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        TS_ASSERT_THROWS_ANYTHING(p_wnt->SetType(NONE));
        
        WntConcentration::Destroy();   
 
        p_wnt = WntConcentration::Instance();
        p_wnt->SetType(LINEAR);
        
        height = 100;
        wnt_level = 0.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 0.0, 1e-9);
        
        height = -1e-12;    // for cells very close to 0 on negative side.
        wnt_level = p_wnt->GetWntLevel(height);
        TS_ASSERT_DELTA(wnt_level, 1.0, 1e-9);        
        
        height = 21.0;
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level, 1.0-height/params->GetCryptLength(), 1e-9);
        
        params->SetCryptLength(10.0);
        wnt_level = p_wnt->GetWntLevel(height);
        
        TS_ASSERT_DELTA(wnt_level , 0.0 , 1e-9);
        
        TS_ASSERT_THROWS_ANYTHING(p_wnt->SetConstantWntValueForTesting(-10));      
    }
    
    
    void TestWntInitialisationSetup() throw(Exception)
    {                
        // create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_2_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        std::vector<WntCellCycleModel*> models;
        
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            WntCellCycleModel* p_model = new WntCellCycleModel();
            TissueCell cell(STEM, HEALTHY, p_model);
            double birth_time = 0.0-i;
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);

            cells.push_back(cell);
            models.push_back(p_model);
        }
        
        // Create the crypt
        MeshBasedTissue<2> crypt(mesh,cells);
        
        CancerParameters::Instance()->SetCryptLength(1.0);

        WntConcentration::Instance()->SetType(LINEAR);
        WntConcentration::Instance()->SetTissue(crypt);
        
        // As there is no tissue simulation we must explicitly initialise the cells
        crypt.InitialiseCells();
        
        MeshBasedTissue<2>::Iterator iter = crypt.Begin();
        
        while(iter!=crypt.End())
        {
            const WntCellCycleModel* p_model = (WntCellCycleModel*) iter->GetCellCycleModel();
            std::vector<double> proteins = p_model->GetProteinConcentrations();
        
            if(iter.GetNode()->rGetLocation()[1]==0.0)
            {
                TS_ASSERT_DELTA(proteins[5],4.975124378109454e-03, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7],6.002649406788524e-01, 1e-3);
                TS_ASSERT_DELTA(proteins[8],1.00, 1e-3);
            }
            else
            {
                TS_ASSERT_DELTA(iter.GetNode()->rGetLocation()[1],1.0,1e-12);
                TS_ASSERT_DELTA(proteins[5],1.000, 1e-3);
                TS_ASSERT_DELTA(proteins[6]+proteins[7],0.0074, 1e-3);
                TS_ASSERT_DELTA(proteins[8],0.00, 1e-3);
            }

            ++iter;
        }
    }
};

#endif /*TESTWNTGRADIENT_HPP_*/
