#ifndef TYSONNOVAKCELLCYCLEMODEL_HPP_
#define TYSONNOVAKCELLCYCLEMODEL_HPP_
#include "AbstractCellCycleModel.hpp"
#include "TysonNovak2001OdeSystem.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
/**
 *  Tyson Novak cell cycle model
 *
 *  Time taken to progress through the cycle is actually deterministic as ODE system
 *  independent of external factors.
 */
class TysonNovakCellCycleModel : public AbstractCellCycleModel
{
private:    
    TysonNovak2001OdeSystem mOdeSystem;
    BackwardEulerIvpOdeSolver* mpSolver;
    double mLastTime;
    double mDivideTime;
    std::vector <double> mProteinConcentrations;
	bool mReadyToDivide;
    
    TysonNovakCellCycleModel(std::vector<double> parentProteinConcentrations, double divideTime);
    
    
public:
    
    TysonNovakCellCycleModel();
    
    ~TysonNovakCellCycleModel();
        
    virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    virtual void ResetModel();
    
    virtual void SetBirthTime(double birthTime);
    
    std::vector< double > GetProteinConcentrations();
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
};





#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
