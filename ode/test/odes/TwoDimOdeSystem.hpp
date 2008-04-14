/*
Copyright (C) University of Oxford, 2008

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TWODIMODESYSTEM_HPP_
#define TWODIMODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"

class TwoDimOdeSystem : public AbstractOdeSystem
{
public :
    TwoDimOdeSystem()
    {
        mNumberOfStateVariables=2;
        
        mInitialConditions.push_back(1);
        mInitialConditions.push_back(2);
        
        mStateVariables.push_back(3);
        mStateVariables.push_back(4);
    }
    
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
    {
        rDY.assign(rY.begin(), rY.end());
    }
};


#endif /*TWODIMODESYSTEM_HPP_*/
