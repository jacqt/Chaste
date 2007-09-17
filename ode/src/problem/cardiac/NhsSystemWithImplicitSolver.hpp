#ifndef NHSSYSTEMWITHIMPLICITSOLVER_HPP_
#define NHSSYSTEMWITHIMPLICITSOLVER_HPP_

#include "NHSCellularMechanicsOdeSystem.hpp"

// todo: doxygen

/** 
 *  NHS system with build in implicit solver. Jon Whiteley's method, which breaks down
 *  the multivariable implicit solve into a sequence of 1d implicit solves
 */
class NhsSystemWithImplicitSolver : public NHSCellularMechanicsOdeSystem
{
private:
    const static double mTolerance = 1e-6; 
    double mDt;
    bool mUseImplicitExplicitSolveForZ;

    std::vector<double> mTempStoredStateVariables;
    std::vector<double> mCurrentStateVars;
    
    double mActiveTensionInitialGuess;
    double mActiveTensionSolution;

    void ImplicitSolveForActiveTension();
    double CalcActiveTensionResidual(double activeTensionGuess);

    double ImplicitSolveForCaTrop(double newActiveTension);

    double ImplicitSolveForZ(double newCaTrop);
    double CalcZResidual(double z, double newCaTrop);
    double ImplicitExplicitSolveForZ(double newCaTrop);

    double ImplicitSolveForQ();

public :  
    NhsSystemWithImplicitSolver();

    void SetActiveTensionInitialGuess(double activeTensionInitialGuess);
    
    void SolveDoNotUpdate(double startTime, double endTime, double timestep);
    
    void UpdateStateVariables();
    
    void UseImplicitExplicitSolveForZ(bool useImplicitExplicitSolveForZ = true);
};

#endif /*NHSSYSTEMWITHIMPLICITSOLVER_HPP_*/
