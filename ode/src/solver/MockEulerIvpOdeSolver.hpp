#ifndef _MOCKEULERIVPODESOLVER_HPP_
#define _MOCKEULERIVPODESOLVER_HPP_

#include "EulerIvpOdeSolver.hpp"

/**
 * This 'mock' class is only used in testing. It is the same
 * as the EulerIvpOdeSolver, but also keeps a count of how many
 * times it has been called. This is useful to check ode solving
 * has been parallelised.
 */


class MockEulerIvpOdeSolver : public EulerIvpOdeSolver
{
    private:
    int mCallCount;
    
    public:
    MockEulerIvpOdeSolver();
                      
    int GetCallCount();
    
    virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem, 
                          double startTime,
                          double endTime,
                          double timeStep,
                          std::vector<double> initialConditions);
                          
    virtual ~MockEulerIvpOdeSolver()
    {
    }                        

};

#endif //_MOCKEULERIVPODESOLVER_HPP_
