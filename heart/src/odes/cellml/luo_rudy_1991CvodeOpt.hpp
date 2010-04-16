#ifdef CHASTE_CVODE
#ifndef CELLLUO_RUDY_1991FROMCELLMLCVODEOPT_HPP_
#define CELLLUO_RUDY_1991FROMCELLMLCVODEOPT_HPP_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: luo_rudy_1991
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 8678, pycml: 8676)
//! on Fri Apr 16 16:49:37 2010
//! 
//! <autogenerated>

#include "AbstractCvodeCell.hpp"
#include "AbstractStimulusFunction.hpp"

class Cellluo_rudy_1991FromCellMLCvodeOpt : public AbstractCvodeCell
{
    // 
    // Settable parameters and readable variables
    // 
    double var_membrane__I_stim;
    double var_membrane__i_Na;
    double var_membrane__i_si;
    double var_membrane__i_K;
    double var_membrane__i_K1;
    double var_membrane__i_Kp;
    double var_membrane__i_b;
    
public:
    double Get_membrane__I_stim();
    double Get_membrane__i_Na();
    double Get_membrane__i_si();
    double Get_membrane__i_K();
    double Get_membrane__i_K1();
    double Get_membrane__i_Kp();
    double Get_membrane__i_b();
    Cellluo_rudy_1991FromCellMLCvodeOpt(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~Cellluo_rudy_1991FromCellMLCvodeOpt();
    void VerifyStateVariables();
    double GetIIonic();
    void EvaluateRhs(double var_environment__time, const N_Vector rY, N_Vector rDY);
};


#endif // CELLLUO_RUDY_1991FROMCELLMLCVODEOPT_HPP_
#endif // CHASTE_CVODE
