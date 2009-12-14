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

#ifndef _FoxModel2002Modified_
#define _FoxModel2002Modified_

/** @file
 *
 * Model: fox_model_2002
 *
 * Processed by pycml - CellML Tools in Python
 *     (translate: $Revision: 698 $, pycml: $Revision: 693 $)
 * on Wed May  9 07:45:54 2007
 *
 * <autogenerated>
 */

#include <cmath>
#include <cassert>
#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"

/**
 * Generated from CellML, and P_Ca parameter modified.
 */
class FoxModel2002Modified : public AbstractCardiacCell
{
public:

    /**
     * Constructor
     * 
     * @param pSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    FoxModel2002Modified(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                         boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~FoxModel2002Modified(void);

    /**
     * Calculates the ionic current
     * 
     * @returns the total ionic current
     */
    double GetIIonic();

    /**
     * Evaluate the derivatives of the state variables
     * 
     * @param var_Environment__time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives (
            double var_environment__time,
            const std::vector<double> &rY,
            std::vector<double> &rDY);
};

#endif
