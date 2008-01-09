#ifndef _ABSTRACTSIMPLECELLCYCLEMODEL_HPP_
#define _ABSTRACTSIMPLECELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class contains all the things common to simple cell cycle models
 * 
 * i.e. models where the length of cell cycle phases are determined when 
 * the cell cycle model is created, 
 * rather than evaluated 'on the fly' by ODEs and suchlike.
 * 
 * N.B. Whether or not the cell should actually divide may depend on 
 * Wnt / Oxygen etc. in subclasses...
 */
class AbstractSimpleCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
    }
    
protected:

    /** 
     * Protected constructor for creating an identical daughter cell 
     * (with the same G1 duration).
     */
    AbstractSimpleCellCycleModel(double g1Duration, unsigned generation)
    {
        mG1Duration = g1Duration;
        mGeneration = generation;
    }
    
    /** 
     * Subclasses can override this function if they wish,
     * this just allocates the cancer parameter default values for each
     * of the different cell types' G1 durations.
     */
    virtual void SetG1Duration();
	
	
public:
    /**
     * Default constructor - creates an AbstractSimpleCellCycleModel
     */
    AbstractSimpleCellCycleModel()
    {
    }
        
    /**
     * Default destructor
     */
    virtual ~AbstractSimpleCellCycleModel()
    {}
    
    virtual void ResetModel();
    
    /**
     * Default UpdateCellCyclePhase function for a simple cell cycle model.
     * 
     * Can be overridden if they should do something more subtle. 
     */
    virtual void UpdateCellCyclePhase();
    
    /** 
     * Set the new cell's G1 duration once it has been created after division. 
     * The duration will be based on cell type. 
     */ 
    void InitialiseDaughterCell();     
    virtual void Initialise();
};

BOOST_IS_ABSTRACT(AbstractSimpleCellCycleModel)

#endif //_ABSTRACTSIMPLECELLCYCLEMODEL_HPP_
