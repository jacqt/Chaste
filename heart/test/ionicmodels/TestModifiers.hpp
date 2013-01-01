/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef _TESTMODIFIERS_HPP_
#define _TESTMODIFIERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/shared_ptr.hpp>
#include "Exception.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ZeroStimulus.hpp"
#include "Modifiers.hpp"
#include "Shannon2004.hpp"

class TestModifiers : public CxxTest::TestSuite
{
public:
    void TestAccessingParametersWithoutModifiers() throw(Exception)
    {
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus());
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellShannon2004FromCellML* p_shannon = new CellShannon2004FromCellML(p_solver, p_stimulus);

        // We should now have all of the following methods available as an alternative to using 'modifiers'
        TS_ASSERT_DELTA(p_shannon->GetParameter("membrane_fast_sodium_current_conductance"),16.0,1e-5);
        TS_ASSERT_DELTA(p_shannon->GetParameter("membrane_L_type_calcium_current_conductance"),5.4e-4,1e-5);
        TS_ASSERT_DELTA(p_shannon->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance"),0.03,1e-5);
        TS_ASSERT_DELTA(p_shannon->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance"),0.07,1e-5);

        delete p_shannon;
     }

    void TestAssigningModifiersToACellModel() throw(Exception)
    {
        boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus());
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellShannon2004FromCellML* p_shannon = new CellShannon2004FromCellML(p_solver, p_stimulus);

        TS_ASSERT_EQUALS(p_shannon->HasModifier("Alan"), false);
        TS_ASSERT_THROWS_THIS(p_shannon->GetModifier("Alan"), "There is no modifier called Alan in this model.");


        // Default modifier shouldn't do anything to the value inputted to calc()
        TS_ASSERT_EQUALS(p_shannon->HasModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance"), true);
        TS_ASSERT_DELTA(p_shannon->GetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance")->Calc(123,0),123,1e-9);

        // Make a new modifier
        boost::shared_ptr<AbstractModifier> p_new_modifier(new FixedModifier(-90.0));

        TS_ASSERT_THROWS_THIS(p_shannon->SetModifier("Alan",p_new_modifier), "There is no modifier called Alan in this model.");

        // Assign it to the Shannon model
        p_shannon->SetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance",p_new_modifier);

        // We should now get a new answer to this.
        TS_ASSERT_DELTA(p_shannon->GetModifier("membrane_rapid_delayed_rectifier_potassium_current_conductance")->Calc(0,0),-90,1e-9);

        delete p_shannon;
    }

    void TestDummyModifiers(void) throw(Exception)
    {
        DummyModifier dummymod;

        double parameter = 2;
        double returned = dummymod.Calc(parameter, 0.0);

        TS_ASSERT_DELTA(parameter, returned, 1e-9);
    }

    void TestFactorModifiers(void) throw(Exception)
    {
        double factor = 2;
        FactorModifier mod(factor);

        double parameter = 2;
        double returned = mod.Calc(parameter, 0.0);

        TS_ASSERT_DELTA(parameter*factor, returned, 1e-9);
    }

    void TestFixedModifiers(void) throw(Exception)
    {
        double fixed = 32;
        FixedModifier mod(fixed);

        double parameter = 2;
        double returned = mod.Calc(parameter, 1.0);

        TS_ASSERT_DELTA(fixed, returned, 1e-9);
    }

    void TestTimeModifiers(void) throw(Exception)
    {
        // This class just provides an example of how to make a time modifier you might want.
        TimeModifier mod;
        double parameter = 2;

        for (unsigned i=0; i<7; i++)
        {
            double time = (double)i;
            double returned = mod.Calc(parameter, time);
            TS_ASSERT_DELTA(parameter*sin(time), returned, 1e-9);
        }
    }

};


#endif //_TESTMODIFIERS_HPP_
