#ifndef TESTCONVERGENCE_HPP_
#define TESTCONVERGENCE_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "BackwardEulerLuoRudyIModel1991.hpp"
#include "ConvergenceTester.hpp"

class TestConvergence : public CxxTest::TestSuite
{   
public:

    
    void Test1D() throw(Exception)
    {
        ConvergenceTester<BackwardEulerLuoRudyIModel1991, BidomainProblem<1>, 1> tester;
    }

//    void Test2D() throw(Exception)
//    {
//        ConvergenceTester<LuoRudyIModel1991OdeSystem, BidomainProblem<2>, 2> tester;
//    }
//    
//    void Test3D() throw(Exception)
//    {
//        ConvergenceTester<LuoRudyIModel1991OdeSystem, BidomainProblem<3>, 3> tester;
//    }    
};

#endif /*TESTCONVERGENCE_HPP_*/
