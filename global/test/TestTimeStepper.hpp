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


#ifndef TESTTIMESTEPPER_HPP_
#define TESTTIMESTEPPER_HPP_

#include <cxxtest/TestSuite.h>

#include "TimeStepper.hpp"
#include <assert.h>
#include <math.h>

#include <cfloat>

class TestTimeStepper : public CxxTest::TestSuite
{
public:
    void TestAdvance()
    {
        const double smidge=1e-10;
        
        double startTime=0.0;
        double endTime=2.0;
        double timeStep=3.7e-05;

        TS_ASSERT_THROWS_ANYTHING(TimeStepper(endTime, startTime, timeStep));

        TimeStepper stepper(startTime, endTime, timeStep);
        
        TS_ASSERT_EQUALS( stepper.EstimateTimeSteps(),
                          (unsigned) ceil((endTime - startTime)/timeStep) );
        
        
        double real_time_step = timeStep;
        unsigned time_step_number = 0;
        double current_time = startTime;
        
        /* We'll trap for stopping times that are close to the end time
         * in order to avoid having a timestep of 1e-14 (or whatever) at
         * the end in the case of rounding errors.
         */
        double close_to_end_time = endTime - smidge*timeStep;
        
        while (current_time < endTime)
        {
            TS_ASSERT(!stepper.IsTimeAtEnd());
            
            time_step_number++;
            // Determine what the value time step should really be like
            double to_time = startTime+time_step_number*timeStep;
    
            if (to_time >= close_to_end_time)
            {
                real_time_step = endTime - current_time;
                // std::cout<<"InternalSolve "<<timeStep<<" "<<real_time_step<<"\n";
                to_time = endTime;
            }
            
            //std::cout << stepper.GetNextTimeStep()-real_time_step << std::endl;
    
            TS_ASSERT_EQUALS(stepper.GetNextTimeStep(), real_time_step);
            TS_ASSERT_EQUALS(stepper.GetTime(), current_time);
            TS_ASSERT_EQUALS(stepper.GetNextTime(), to_time);
                                
            // Determine the new current time
            current_time = to_time;
            stepper.AdvanceOneTimeStep();
    
            TS_ASSERT_EQUALS(current_time, stepper.GetTime());
        }
        
        TS_ASSERT(stepper.IsTimeAtEnd());
        TS_ASSERT(stepper.GetTimeStepsElapsed()==time_step_number);
    }
};

#endif /*TESTTIMESTEPPER_HPP_*/
