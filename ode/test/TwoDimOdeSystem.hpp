#ifndef TWODIMODESYSTEM_HPP_
#define TWODIMODESYSTEM_HPP_

#include "AbstractOdeSystem.hpp"

class TwoDimOdeSystem : public AbstractOdeSystem
{
public :
    TwoDimOdeSystem()
    {
        mInitialConditions.push_back(1);
        mInitialConditions.push_back(2);
        
        mStateVariables.push_back(3);
        mStateVariables.push_back(4);
    } 

    std::vector<double> TwoDimOdeSystem::EvaluateYDerivatives (double time, const std::vector<double> &rY)
    {   
        // do nothing
    }
};


#endif /*TWODIMODESYSTEM_HPP_*/
