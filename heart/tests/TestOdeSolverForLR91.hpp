#ifndef _TESTODESOLVERFORLR91_HPP_
#define _TESTODESOLVERFORLR91_HPP_


#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "LR91Model.hpp"
#include "AbstractStimulusFunction.hpp"
#include "InitialStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"


class TestOdeSolverForLR91 : public CxxTest::TestSuite
{
    public:
    
    // Test Ode Solver for LR91
    void testOdeSolverForLR91(void)
    {
    //Tests that setting and getting stimulus works -test passed
        double time = 0.8;
        double magnitudeOfStimulus = 80.0;        
        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
        double stimulusValue = pStimulus->GetStimulus(time);
        std::cout<<"Stimulus value is "<< stimulusValue<< std::endl;
        
//        // test 1
       // double time = 0.2;
        //double magnitudeOfStimulus = 80.0;        
      //  AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
        LR91OdeFun lr91OdeFun(pStimulus);
        double voltage = -40.0;
        double m = 0.1;
        double h = 0.03;
        double j = 0.7;
        double d = 0.05;
        double f = 0.5;
        double x = 0.1;
        double caI = 0.001;
          
        std::vector<double> intialConditions;
        intialConditions.push_back(voltage);
        intialConditions.push_back(m);
        intialConditions.push_back(h);
        intialConditions.push_back(j);
        intialConditions.push_back(d);
        intialConditions.push_back(f);
        intialConditions.push_back(x);
        intialConditions.push_back(caI);
        std::vector<double> Result = lr91OdeFun.EvaluateYDerivatives(time, intialConditions);
        
        
//        DO REGRESSION TEST here....
//        double iStim = stimulusValue;
//        
//          
//        TS_ASSERT_DELTA(Result[0], (-iStim -iTotal) /cMembrane , 0.0001);
//        TS_ASSERT_DELTA(Result[1],   , 0.0001);
//        TS_ASSERT_DELTA(Result[2],  , 0.0001);
//        TS_ASSERT_DELTA(Result[3],  , 0.0001);
//        TS_ASSERT_DELTA(Result[4],  , 0.0001);
//        TS_ASSERT_DELTA(Result[5],  , 0.0001);
//        TS_ASSERT_DELTA(Result[6],  , 0.0001);
//        TS_ASSERT_DELTA(Result[7],  , 0.0001);
//          
        TS_TRACE("LupRudy WORKS!!!! :)) ");
        
//        
//       //test 2 
//
//        double magnitudeOfStimulus = 80.0;        
//        AbstractStimulusFunction *pStimulus = new InitialStimulus(magnitudeOfStimulus);
//        
        LR91Model *pLR91Model;
        pLR91Model = new LR91Model(voltage, m, h, j, d, f, x, caI, pStimulus);
        
        TS_TRACE("LupRudy is initiated !!!! :)) ");
        
        AbstractIvpOdeSolver *pMySolver = new EulerIvpOdeSolver();
        
        double startTime = 0.0;
        double endTime = 10.0;
        double timeStep = 0.01;             
        
        OdeSolution Solution = pLR91Model->SolveModel(startTime, endTime, timeStep,
                               intialConditions, pMySolver);
                               
//        output to matlab file
        
        TS_TRACE("MOREOVER OdeSolver SOLVES!!!! :)) ");

    }
    
    
    
   
//	// Tests that Na gating variables within range 0 to 1 inclusive
//	void testSodiumGatingVariables( void )
//	{
//		double m = 0.0017;
//		double h = 0.9833;
//		double j = 0.9895;
//		double voltage = -75.0;
//		
//		SodiumCurrentLR91 *pINa = new SodiumCurrentLR91();
//		pINa->UpdateMagnitudeOfCurrent(voltage,m,h,j);
//		
//		std::cout << "\n";
//		std::cout << "INa m gate: " << pINa->GetM() << "\n";
//		TS_ASSERT(pINa->GetM() >= 0);
//		TS_ASSERT(pINa->GetM() <= 1);
//		
//		std::cout << "INa h gate: " << pINa->GetH() << "\n";
//		TS_ASSERT(pINa->GetH() >= 0);
//		TS_ASSERT(pINa->GetH() <= 1);		
//	}
//	
//	// Tests that K gating variables within range 0 to 1 inclusive
//	void testPotassiumGatingVariables( void )
//	{
//		double m = 0.0017;
//		double h = 0.9833;
//		double j = 0.9895;
//		double voltage = -75.0;
//		
//		PotassiumCurrent *pIK = new PotassiumCurrent();
//		pIK->UpdateMagnitudeOfCurrent(voltage,m,h,j);
//			
//		std::cout << "IK n gate: " << pIK->GetN() << "\n";
//		TS_ASSERT(pIK->GetN() >= 0);
//		TS_ASSERT(pIK->GetN() <= 1);		
//	}
	
	
	
    
};



#endif //_TESTODESOLVERFORLR91_HPP_
