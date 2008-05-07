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


#ifndef TESTCARDIACELECTROMECHANICSPROBLEMEXPERIMENTAL_HPP_
#define TESTCARDIACELECTROMECHANICSPROBLEMEXPERIMENTAL_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem1d.hpp"
#include "CardiacElectroMechanicsProblem.hpp"

class TestCardiacElectroMechanicsProblemExperimental : public CxxTest::TestSuite
{
public:
    void xTest1dCompareExplicitVsImplicit() throw(Exception)
    {
        double time_step = 0.01;
        PlaneStimulusCellFactory<1> cell_factory(time_step, -1000*1000);

        // put these inside braces so that there destructor gets called and all output
        // files get closed
        {
            // instabilities appear at about 6.8
            CardiacElectroMechanicsProblem1d explicit_problem(&cell_factory, 5, true, 1,  "ExplicitCardiacElectroMech");
            explicit_problem.Solve();

            CardiacElectroMechanicsProblem1d implicit_problem(&cell_factory, 5, false, 1, "ImplicitCardiacElectroMech");
            implicit_problem.Solve();
        }

        std::string file1 = OutputFileHandler::GetChasteTestOutputDirectory() + 
                            "ExplicitCardiacElectroMech/length.txt";

        std::string file2 = OutputFileHandler::GetChasteTestOutputDirectory() + 
                            "ImplicitCardiacElectroMech/length.txt";
        
        std::ifstream ifs1(file1.c_str());
        std::ifstream ifs2(file2.c_str());
        
        double length1;
        double length2;
        
        ifs1 >> length1;
        ifs2 >> length2;

        double counter = 0;
        
        TS_ASSERT(!ifs1.eof());
        TS_ASSERT(!ifs2.eof());
        
        while(!ifs1.eof())
        {
            std::cout << length1 << " " << length2 << "\n";
            TS_ASSERT_DELTA(length1, length2,  fabs(length1*1e-5));
            ifs1 >> length1;
            ifs2 >> length2;
            
            // hardcoded test
            if(counter==450)
            {
                TS_ASSERT_DELTA(length2, 0.999378, 1e-5);
            }
            counter++;
        }
    }
    
    void Test2dCompareExplicitVsImplicit() throw(Exception)
    {
        PlaneStimulusCellFactory<2> cell_factory(0.01, -1000*1000);

        unsigned dt=128;
        for(unsigned i=0;i<4;i++)
        {
            std::stringstream name;
            name << "CardiacElectroMech_Time_MORE_" << dt;
            
            CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 200, 10, false, dt, 0.01, name.str());
            implicit_problem.SetNoElectricsOutput();
            implicit_problem.Solve();
            
            dt *= 2;
        }


//        unsigned num_nodes_in_each_dir=8;
//        for(unsigned i=0;i<1;i++)
//        {
//            std::stringstream name;
//            name << "CardiacElectroMech_Space_NEW_" << num_nodes_in_each_dir;
//            
//            CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 200, num_nodes_in_each_dir, false, 100, 1.0, name.str());
//            implicit_problem.SetNoElectricsOutput();
//            implicit_problem.Solve();
//            
//            num_nodes_in_each_dir *= 2;
//        }

//        double nhs_ode_time_step = 0.01;
//
//        for(unsigned i=0; i<8; i++)
//        {
//            std::stringstream name;
//            name << "CardiacElectroMech_OdeTimeStepNew_" << i;
//
//            std::cout << nhs_ode_time_step << "\n";
//            
//            CardiacElectroMechanicsProblem<2> implicit_problem(&cell_factory, 200, 10, false, 128, nhs_ode_time_step, name.str());
//            implicit_problem.SetNoElectricsOutput();
//            implicit_problem.Solve();
//            
//            nhs_ode_time_step *= 2;
//        }
    }
};
#endif /*TESTCARDIACELECTROMECHANICSPROBLEMEXPERIMENTAL_HPP_*/
