#include "AbstractOdeBasedCellCycleModel.hpp"

AbstractOdeBasedCellCycleModel::AbstractOdeBasedCellCycleModel(double lastTime)
        : mpOdeSystem(NULL),
          mLastTime(lastTime),
          mDivideTime(lastTime),
          mReadyToDivide(false),
          mFinishedRunningOdes(false),
          mG2PhaseStartTime(DBL_MAX)
{
    AbstractCellCycleModel::SetBirthTime(lastTime);
}


AbstractOdeBasedCellCycleModel::~AbstractOdeBasedCellCycleModel()
{
    if (mpOdeSystem!=NULL)
    {
        delete mpOdeSystem;   
    }
}


void AbstractOdeBasedCellCycleModel::SetBirthTime(double birthTime)
{
    AbstractCellCycleModel::SetBirthTime(birthTime);
    mLastTime = birthTime;
    mDivideTime = birthTime;
}


std::vector<double> AbstractOdeBasedCellCycleModel::GetProteinConcentrations() const
{
    assert(mpOdeSystem!=NULL);
    return mpOdeSystem->rGetStateVariables();
}


void AbstractOdeBasedCellCycleModel::SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations)
{
    assert(mpOdeSystem!=NULL);
    assert(proteinConcentrations.size()==mpOdeSystem->rGetStateVariables().size());
    mLastTime = lastTime;
    mpOdeSystem->SetStateVariables(proteinConcentrations);
}


bool AbstractOdeBasedCellCycleModel::ReadyToDivide()
{
    assert(mpOdeSystem!=NULL);
    assert(mpCell!=NULL);
    
    double current_time = SimulationTime::Instance()->GetDimensionalisedTime();
    
//    we want to start solving the ODEs at the end of M phase - could possibly 
//    hijack mLastTime to do this
    //std::cout << "mBirthTime = " << mBirthTime << " current time = "<< current_time << "\n" << std::flush; 
    if (mCurrentCellCyclePhase == M_PHASE)
    {
        double m_duration = GetMDuration();
        if (current_time - mBirthTime >= m_duration)
        {
            mCurrentCellCyclePhase = G_ONE_PHASE;
            mLastTime = m_duration + mBirthTime;
        } 
        else 
        {
            return false;
        }
        //return false;
    }
    
    if (current_time>mLastTime)
    {
        if (!mFinishedRunningOdes)
        {   
            //std::cout << "Running ODEs to time " << current_time << "\n" << std::flush;
            mFinishedRunningOdes = SolveOdeToTime(current_time);
            
            for (unsigned i=0 ; i<mpOdeSystem->GetNumberOfStateVariables() ; i++)
            {
                if (mpOdeSystem->rGetStateVariables()[i]< 0)
                {
                    #define COVERAGE_IGNORE
                    std::cout << "Protein["<< i <<"] = "<< mpOdeSystem->rGetStateVariables()[i] << "\n";
                    EXCEPTION("A protein concentration has gone negative\nCHASTE predicts that the CellCycleModel numerical method is probably unstable.");
                    #undef COVERAGE_IGNORE
                }
            }
                        
            if (mFinishedRunningOdes)
            {
                mCurrentCellCyclePhase = S_PHASE;  
                mG2PhaseStartTime = GetOdeStopTime() + GetSDuration();
                mDivideTime = mG2PhaseStartTime + GetG2Duration();
                //std::cout << "Ode Stop time = " << GetOdeStopTime() << "\n" << std::flush;
                //std::cout << "mDivideTime = " << mDivideTime << "\n" << std::flush;
//              need to do some clever business here - instead of a divide time, we should
//              get back the separate pahse durations for S and G2
//              (need to do this in each concrete instance)
                if (current_time >= mDivideTime)
                {
                    mReadyToDivide = true;
                }
                if (current_time >= mG2PhaseStartTime)
                {
                    mCurrentCellCyclePhase = G_TWO_PHASE;
                }
            }
            mLastTime = current_time;   // This is the last time the ODEs were evaluated.
        }
        else
        {   // ODE model finished, just increasing time until division...
            if (current_time >= mDivideTime)
            {
                mReadyToDivide = true;
            }
            if (current_time >= mG2PhaseStartTime)
            {
                mCurrentCellCyclePhase = G_TWO_PHASE;
            }
        }
    }
    return mReadyToDivide;
}


void AbstractOdeBasedCellCycleModel::ResetModel()
{
    assert(mReadyToDivide);
    mBirthTime = mDivideTime;
    mLastTime = mDivideTime;
    mFinishedRunningOdes = false;
    mReadyToDivide = false;
    mCurrentCellCyclePhase = M_PHASE;
}
