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
#include "StochasticDivisionRuleCellCycleModel.hpp"


StochasticDivisionRuleCellCycleModel::StochasticDivisionRuleCellCycleModel(double g1Duration,
                                                                           unsigned generation,
                                                                           bool dividedSymmetrically)
    : AbstractSimpleCellCycleModel(g1Duration, generation),
      mDividedSymmetrically(dividedSymmetrically)
{
}


StochasticDivisionRuleCellCycleModel::StochasticDivisionRuleCellCycleModel(bool dividedSymmetrically)
    : mDividedSymmetrically(dividedSymmetrically)
{
}

       
void StochasticDivisionRuleCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);

    CancerParameters* p_params = CancerParameters::Instance();
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mpCell->GetCellType())
    {
        case STEM:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetStemCellG1Duration(), 1.0);
            break;
        case TRANSIT:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetTransitCellG1Duration(), 1.0);
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < p_params->GetMinimumGapDuration())
    {
        #define COVERAGE_IGNORE
        mG1Duration = p_params->GetMinimumGapDuration();
        #undef COVERAGE_IGNORE
    }
}


void StochasticDivisionRuleCellCycleModel::ResetForDivision()
{
    if (mGeneration+1u > CancerParameters::Instance()->GetMaxTransitGenerations())
    {
        mpCell->SetCellType(DIFFERENTIATED);
    }

    /**
     * If dealing with a stem cell, we may have symmetric division.
     * We therefore neglect the possibility of de-differentiation.
     * 
     * NB. This code must be implemented before the call to
     * AbstractSimpleCellCycleModel::ResetForDivision(), because 
     * that method sets the G1 duration based on the cell type.
     */
    if (mpCell->GetCellType() == STEM)
    {
        double test_number = RandomNumberGenerator::Instance()->ranf(); // U(0,1)
        double sym_div_prob = CancerParameters::Instance()->GetSymmetricDivisionProbability();

        // If undergoing symmetric division...
        if (test_number < sym_div_prob)
        {
            mDividedSymmetrically = true;

            // Check if the daughter cells are both STEM or TRANSIT.
            // We assign an equal probability to each of these events.
            if (test_number < 0.5*sym_div_prob)
            {
                mpCell->SetCellType(STEM);
            }
            else
            {
                mpCell->SetCellType(TRANSIT);
            }
        }
        else
        {
            mDividedSymmetrically = false;
        }
    }

    AbstractSimpleCellCycleModel::ResetForDivision();

    if (mpCell->GetCellType() == STEM)
    {
        mGeneration = 0;
    }
}


void StochasticDivisionRuleCellCycleModel::InitialiseDaughterCell()
{
    // If the cell was born out of symmetric division,
    // then do not alter generation or cell type
    if (mDividedSymmetrically == false)
    {
        if (mGeneration == 0)
        {
            mGeneration = 1;
        }
        // Daughter cell is always a TRANSIT or DIFFERENTIATED
        mpCell->SetCellType(TRANSIT);
        if (mGeneration > CancerParameters::Instance()->GetMaxTransitGenerations())
        {
            mpCell->SetCellType(DIFFERENTIATED);
        }
    }
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}


AbstractCellCycleModel* StochasticDivisionRuleCellCycleModel::CreateDaughterCellCycleModel()
{
    // Use a private constructor that doesn't reset mG1Duration
    return new StochasticDivisionRuleCellCycleModel(mG1Duration, mGeneration, mDividedSymmetrically);
}


bool StochasticDivisionRuleCellCycleModel::DividedSymmetrically()
{
    return mDividedSymmetrically;
}
