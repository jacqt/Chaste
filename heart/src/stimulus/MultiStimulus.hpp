/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef MULTISTIMULUS_HPP_
#define MULTISTIMULUS_HPP_

#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class provides a stimulus function which is the
 * sum of an arbitrary number of stimuli.
 *
 * After creation it behaves like a ZeroStimulus until
 * any number of stimuli are added.
 */

class MultiStimulus : public AbstractStimulusFunction
{
private:

    /** Vector of stimuli. */
    std::vector<AbstractStimulusFunction*> mStimuli;

public:
    /**
     * Combine a stimulus with the existing ones.
     *
     * @param pStimulus pointer to the stimulus to be added.
     */
     void AddStimulus(AbstractStimulusFunction* pStimulus);

    /**
     * Get the magnitude of the multiple stimuli at time 'time'
     *
     * @param time  time at which to return the stimulus
     * @return  Magnitude of stimulus at time 'time'.
     */
     double GetStimulus(double time);
};

#endif /*MULTISTIMULUS_HPP_*/
