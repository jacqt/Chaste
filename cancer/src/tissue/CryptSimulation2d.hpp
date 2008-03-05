#ifndef CRYPTSIMULATION2D_HPP_
#define CRYPTSIMULATION2D_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "TissueSimulation.cpp"
#include "SimpleDataWriter.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

class CryptSimulation2d : public TissueSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestCryptSimulation2d;

private :
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {   
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>      
        archive & boost::serialization::base_object<TissueSimulation<2> >(*this);
        archive & mUseJiggledBottomCells;
        
        assert((mpMechanicsSystem == NULL) ||
               (&(mpMechanicsSystem->rGetTissue()) == &mrTissue));
    }
    
    /** Whether to use a flat bottom surface or to jiggle the cells on the bottom surface */
    bool mUseJiggledBottomCells;
    
    /** The file that the values of beta catenin is written out to. */
    out_stream mBetaCatResultsFile;    
    
    MeshBasedTissue<2>* mpStaticCastTissue;
    
    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent 
     * and daughter in opposing directions on this axis.
     * 
     * @param node_index The parent node index
     * 
     * @return daughter_coords The coordinates for the daughter cell.
     * 
     */
    c_vector<double, 2> CalculateDividingCellCentreLocations(AbstractTissue<2>::Iterator parentCell)
    {     
        double separation = CancerParameters::Instance()->GetDivisionSeparation();
        c_vector<double, 2> parent_coords = parentCell.rGetLocation();
        c_vector<double, 2> daughter_coords;
            
        // pick a random direction and move the parent cell backwards by 0.5*sep in that
        // direction and return the position of the daughter cell (0.5*sep forwards in the
        // random vector direction
    
        // Make a random direction vector of the required length
        c_vector<double, 2> random_vector;        
        
        double random_angle = RandomNumberGenerator::Instance()->ranf();
        random_angle *= 2.0*M_PI;
        
        random_vector(0) = 0.5*separation*cos(random_angle);
        random_vector(1) = 0.5*separation*sin(random_angle);
        
        c_vector<double, 2> proposed_new_parent_coords = parent_coords-random_vector;
        c_vector<double, 2> proposed_new_daughter_coords = parent_coords+random_vector;
        
        if (   (proposed_new_parent_coords(1) >= 0.0)
            && (proposed_new_daughter_coords(1) >= 0.0))
        {
            // We are not too close to the bottom of the tissue
            // move parent
            parent_coords = proposed_new_parent_coords;
            daughter_coords = proposed_new_daughter_coords;
        }
        else
        {   
            proposed_new_daughter_coords = parent_coords+2.0*random_vector;
            while (proposed_new_daughter_coords(1) < 0.0)
            {
                random_angle = RandomNumberGenerator::Instance()->ranf();
                random_angle *= 2.0*M_PI;
                
                random_vector(0) = separation*cos(random_angle);
                random_vector(1) = separation*sin(random_angle);
                proposed_new_daughter_coords = parent_coords+random_vector;
            }
            daughter_coords = proposed_new_daughter_coords;
        }
        
        assert(daughter_coords(1)>=0.0); // to make sure dividing cells stay in the tissue
        assert(parent_coords(1)>=0.0);   // to make sure dividing cells stay in the tissue
                
        // set the parent to use this location
        ChastePoint<2> parent_coords_point(parent_coords);
        mrTissue.MoveCell(parentCell, parent_coords_point);
        return daughter_coords;           
    }
    
    
    /**
     * Moves each node to a new position for this timestep
     *
     * @param rDrDt the x and y force components on each node.
     */
    void UpdateNodePositions(const std::vector< c_vector<double, 2> >& rDrDt)
    {
        // update ghost positions first because they do not affect the real cells
        mpStaticCastTissue->UpdateGhostPositions(mDt);
        
        // Iterate over all cells to update their positions.
        for (AbstractTissue<2>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            TissueCell& cell = *cell_iter;
            unsigned index = cell.GetNodeIndex();
            
            ChastePoint<2> new_point(mrTissue.GetNode(index)->rGetLocation() + mDt*rDrDt[index]);
                        
            bool is_wnt_included = WntConcentration::Instance()->IsWntSetUp();
            if (!is_wnt_included) WntConcentration::Destroy();
            
            // stem cells are fixed if no wnt, so reset the x-value to the old x-value           
            if ((cell.GetCellType()==STEM) && (!is_wnt_included))
            {
                new_point.rGetLocation()[0] = mrTissue.GetNode(index)->rGetLocation()[0];
                new_point.rGetLocation()[1] = mrTissue.GetNode(index)->rGetLocation()[1];
            }
            
            // for all cells - move up if below the bottom surface
            if (new_point.rGetLocation()[1] < 0.0) 
            {
                new_point.rGetLocation()[1] = 0.0; 
                if (mUseJiggledBottomCells)
                {
                   /*  
                    * Here we give the cell a push upwards so that it doesn't  
                    * get stuck on y=0 for ever (ticket:422). 
                    *  
                    * Note that all stem cells may get moved to same height and  
                    * random numbers try to ensure we aren't left with the same  
                    * problem at a different height! 
                    */ 
                    new_point.rGetLocation()[1] = 0.05*mpRandomGenerator->ranf();
                } 
            } 
            
            // move the cell
            assert(new_point[1]>=0.0);
            mrTissue.MoveCell(cell_iter, new_point); 
                    
        }
    }
    
    
    void WriteVisualizerSetupFile()
    {
        *mpSetupFile << "MeshWidth\t" << mpStaticCastTissue->rGetMesh().GetWidth(0u);// get furthest distance between nodes in the x-direction
    }
    
    
    void SetupWriteBetaCatenin()
    {
        OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/",false);
        mBetaCatResultsFile = output_file_handler.OpenOutputFile("results.vizbCat");
        *mpSetupFile << "BetaCatenin\n";
    }
    
    
    void WriteBetaCatenin(double time)
    {
        *mBetaCatResultsFile <<  time << "\t";
        
        double global_index;
        double x;
        double y;
        double b_cat_membrane;
        double b_cat_cytoplasm;
        double b_cat_nuclear;
        
        for (AbstractTissue<2>::Iterator cell_iter = mrTissue.Begin();
             cell_iter != mrTissue.End();
             ++cell_iter)
        {
            global_index = (double) cell_iter.GetNode()->GetIndex();
            x = cell_iter.rGetLocation()[0];
            y = cell_iter.rGetLocation()[1];
            b_cat_membrane = cell_iter->GetCellCycleModel()->GetMembraneBoundBetaCateninLevel();
            b_cat_cytoplasm = cell_iter->GetCellCycleModel()->GetCytoplasmicBetaCateninLevel();
            b_cat_nuclear = cell_iter->GetCellCycleModel()->GetNuclearBetaCateninLevel();
            
            *mBetaCatResultsFile << global_index << " " << x << " " << y << " " << b_cat_membrane << " " << b_cat_cytoplasm << " " << b_cat_nuclear << " ";
        }

        *mBetaCatResultsFile << "\n";
    }    
    
    
    void SetupSolve()
    {
        if (   ( mrTissue.Begin() != mrTissue.End() )  // there are any cells
            && ( mrTissue.Begin()->GetCellCycleModel()->UsesBetaCat()) ) // assume all the cells are the same
        {
            SetupWriteBetaCatenin();            
            double current_time = SimulationTime::Instance()->GetDimensionalisedTime();            
            WriteBetaCatenin(current_time);
        }
    }
    
    
    void PostSolve()
    {
        SimulationTime *p_time = SimulationTime::Instance();
                
        if ((p_time->GetTimeStepsElapsed()+1)%mSamplingTimestepMultiple==0)
        {
            if (   ( mrTissue.Begin() != mrTissue.End() )  // there are any cells
                && ( mrTissue.Begin()->GetCellCycleModel()->UsesBetaCat()) ) // assume all the cells are the same
            {
                double time_next_step = p_time->GetDimensionalisedTime() + p_time->GetTimeStep();
                WriteBetaCatenin(time_next_step);
            }
        }        
    }

    void AfterSolve()
    {
        if (   ( mrTissue.Begin() != mrTissue.End() )  // there are any cells
            && ( mrTissue.Begin()->GetCellCycleModel()->UsesBetaCat()) ) // assume all the cells are the same
        {
            mBetaCatResultsFile->close();
        }
    }
    
public :            

    /** 
     *  Constructor
     * 
     *  @param rTissue A tissue facade class (contains a mesh and cells)
     *  @param deleteTissue whether to delete the tissue on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     */
    CryptSimulation2d(AbstractTissue<2>& rTissue, 
                      AbstractDiscreteTissueMechanicsSystem<2>* pMechanicsSystem=NULL, 
                      bool deleteTissueAndMechanicsSystem=false,
                      bool initialiseCells=true)
        : TissueSimulation<2>(rTissue, pMechanicsSystem, deleteTissueAndMechanicsSystem, initialiseCells),
          mUseJiggledBottomCells(false)
    {
        mpStaticCastTissue = static_cast<MeshBasedTissue<2>*>(&mrTissue);
    }
    
    
    void UseJiggledBottomCells()
    {            
        mUseJiggledBottomCells = true;                
    }
    
    /**
     * Saves the whole tissue simulation for restarting later.
     *
     * Puts it in the folder mOutputDirectory/archive/
     * and the file "tissue_sim_at_time_<SIMULATION TIME>.arch"
     *
     * First archives simulation time then the simulation itself.
     * 
     * Note that this method has to be implemented in this class,
     * so you save the right sort of simulation to the archive.
     * Not really sure why this is needed, but...
     */
    void Save()
    {
        CommonSave(this);
    }

    /**
     * Loads a saved tissue simulation to run further.
     *
     * @param rArchiveDirectory the name of the simulation to load
     * (specified originally by simulation.SetOutputDirectory("wherever"); )
     * @param rTimeStamp the time at which to load the simulation (this must
     * be one of the times at which simulation.Save() was called) 
     * 
     * Note that this method has to be implemented in this class, since it's a static method.
     */
    static CryptSimulation2d* Load(const std::string& rArchiveDirectory, const double& rTimeStamp)
    {
	    std::string archive_filename = TissueSimulation<2>::GetArchivePathname(rArchiveDirectory, rTimeStamp);

        // Create an input archive
        std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
        boost::archive::text_iarchive input_arch(ifs);

	    TissueSimulation<2>::CommonLoad(input_arch);

        CryptSimulation2d* p_sim; 
        input_arch >> p_sim;
         
        if (p_sim->rGetTissue().GetNumNodes()!=p_sim->rGetTissue().rGetCells().size()) 
        { 
            #define COVERAGE_IGNORE 
            std::stringstream string_stream; 
            string_stream << "Error in Load(), number of nodes (" << p_sim->rGetTissue().GetNumNodes() 
                          << ") is not equal to the number of cells (" << p_sim->rGetTissue().rGetCells().size()  
                          << ")"; 
            EXCEPTION(string_stream.str()); 
            #undef COVERAGE_IGNORE 
        } 
          
        return p_sim;         
    }       
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(CryptSimulation2d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2d * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const AbstractTissue<2> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    
    const AbstractDiscreteTissueMechanicsSystem<2> * p_spring_system = &(t->rGetMechanicsSystem());
    ar & p_spring_system;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation2d * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    AbstractTissue<2>* p_tissue;
    ar >> p_tissue;
    
    AbstractDiscreteTissueMechanicsSystem<2>* p_spring_system;
    ar >> p_spring_system;
    
    
    // invoke inplace constructor to initialize instance
    ::new(t)CryptSimulation2d(*p_tissue, p_spring_system, true, false);
}
}
} // namespace ...

#endif /*CRYPTSIMULATION2D_HPP_*/

