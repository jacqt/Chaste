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
#ifndef FIXEDCELLCYCLEMODEL_HPP_
#define FIXEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleMeinekeCellCycleModel.hpp"

/**
 *  Fixed cell cycle model
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  CancerParameters::StemCellG1Duration + SG2MDuration
 *  and CancerParameters::TransitCellG1Duration + SG2MDuration)
 */
class FixedCellCycleModel : public AbstractSimpleMeinekeCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleMeinekeCellCycleModel>(*this);
    }

    /**
     * Private constructor for identical cells.
     * 
     * @param g1Duration 
     * @param generation
     */
    FixedCellCycleModel(double g1Duration, unsigned generation):
        AbstractSimpleMeinekeCellCycleModel(g1Duration, generation){};

public:

    /**
     * Default constructor. Note that mBirthTime is set in 
     * AbstractCellCycleModel() and mG1Duration is set in 
     * AbstractSimpleCellCycleModel().
     */
    FixedCellCycleModel();

    /** 
     * Overridden builder method to create new instances of 
     * the cell cycle model.
     */
    AbstractCellCycleModel* CreateDaughterCellCycleModel();

};

// Declare identifier for the serializer
BOOST_CLASS_EXPORT(FixedCellCycleModel)


#endif /*FIXEDCELLCYCLEMODEL_HPP_*/
