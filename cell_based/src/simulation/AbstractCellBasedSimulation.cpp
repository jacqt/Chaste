/*

Copyright (C) University of Oxford, 2005-2011

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

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"

#include <typeinfo>

template<unsigned DIM>
AbstractCellBasedSimulation<DIM>::AbstractCellBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                              bool deleteCellPopulationAndCellKillersInDestructor,
                                              bool initialiseCells)
    : mEndTime(0.0),  // hours - this is set later on
      mrCellPopulation(rCellPopulation),
      mDeleteCellPopulationAndCellKillersInDestructor(deleteCellPopulationAndCellKillersInDestructor),
      mInitialiseCells(initialiseCells),
      mNoBirth(false),
      mUpdateCellPopulation(true),
      mOutputDirectory(""),
      mSimulationOutputDirectory(mOutputDirectory),
      mNumBirths(0),
      mNumDeaths(0),
      mSamplingTimestepMultiple(1)
{
    // This line sets a random seed of 0 if it wasn't specified earlier.
    mpRandomGenerator = RandomNumberGenerator::Instance();

    // Different time steps are used for cell-centre and vertex-based simulations
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        mDt = 1.0/120.0; // 30 seconds
    }
    else
    {
        mDt = 0.002; // smaller time step required for convergence/stability
    }

    if (mInitialiseCells)
    {
        mrCellPopulation.InitialiseCells();
    }
}

template<unsigned DIM>
AbstractCellBasedSimulation<DIM>::~AbstractCellBasedSimulation()
{
    if (mDeleteCellPopulationAndCellKillersInDestructor)
    {
        for (typename std::vector<AbstractCellKiller<DIM>*>::iterator it=mCellKillers.begin();
             it != mCellKillers.end();
             ++it)
        {
            delete *it;
        }

        delete &mrCellPopulation;
    }
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::DoCellBirth()
{
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // Check if this cell is ready to divide - if so create a new cell etc.
        if (cell_iter->GetAge() > 0.0)
        {
            if (cell_iter->ReadyToDivide())
            {
                // Create a new cell
                CellPtr p_new_cell = cell_iter->Divide();

                // Call method that determines how cell division occurs and returns a vector
                c_vector<double, DIM> new_location = CalculateCellDivisionVector(*cell_iter);

                // Add new cell to the cell population
                mrCellPopulation.AddCell(p_new_cell, new_location, *cell_iter);

                // Update counter
                num_births_this_step++;
            }
        }
    }
    return num_births_this_step;
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step = 0;

    /*
     * This labels cells as dead or apoptosing. It does not actually remove the cells,
     * mrCellPopulation.RemoveDeadCells() needs to be called for this.
     */
    for (typename std::vector<AbstractCellKiller<DIM>*>::iterator killer_iter = mCellKillers.begin();
         killer_iter != mCellKillers.end();
         ++killer_iter)
    {
        (*killer_iter)->TestAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrCellPopulation.RemoveDeadCells();

    return num_deaths_this_step;
}

template<unsigned DIM>
c_vector<double, DIM> AbstractCellBasedSimulation<DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    /**
     * \todo Could remove this dynamic_cast by moving the code block below into
     * AbstractCentreBasedCellPopulation::AddCell(), allowing it to be overruled by
     * this method when overridden in subclasses. See also comment on #1093.
     */
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&mrCellPopulation))
    {
        // Location of parent and daughter cells
        c_vector<double, DIM> parent_coords = mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, DIM> daughter_coords;

        // Get separation parameter
        double separation = static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&mrCellPopulation)->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, DIM> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        switch (DIM)
        {
            case 1:
            {
                double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

                random_vector(0) = 0.5*separation*random_direction;
                break;
            }
            case 2:
            {
                double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

                random_vector(0) = 0.5*separation*cos(random_angle);
                random_vector(1) = 0.5*separation*sin(random_angle);
                break;
            }
            case 3:
            {
                double random_zenith_angle = M_PI*RandomNumberGenerator::Instance()->ranf(); // phi
                double random_azimuth_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf(); // theta

                random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(2) = 0.5*separation*cos(random_zenith_angle);
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }

        parent_coords = parent_coords - random_vector;
        daughter_coords = parent_coords + random_vector;

        // Set the parent to use this location
        ChastePoint<DIM> parent_coords_point(parent_coords);
        unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else
    {
        // Do something for Vertex / Potts Models here
        return zero_vector<double>(DIM);
    }
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
}

template<unsigned DIM>
double AbstractCellBasedSimulation<DIM>::GetDt()
{
    return mDt;
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::GetNumBirths()
{
    return mNumBirths;
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::GetNumDeaths()
{
    return mNumDeaths;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetEndTime(double endTime)
{
    assert(endTime > 0);
    mEndTime = endTime;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}

template<unsigned DIM>
std::string AbstractCellBasedSimulation<DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
}

template<unsigned DIM>
AbstractCellPopulation<DIM>& AbstractCellBasedSimulation<DIM>::rGetCellPopulation()
{
    return mrCellPopulation;
}

template<unsigned DIM>
const AbstractCellPopulation<DIM>& AbstractCellBasedSimulation<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetUpdateCellPopulationRule(bool updateCellPopulation)
{
    mUpdateCellPopulation = updateCellPopulation;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetNoBirth(bool noBirth)
{
    mNoBirth = noBirth;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::AddCellKiller(AbstractCellKiller<DIM>* pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}

template<unsigned DIM>
std::vector<double> AbstractCellBasedSimulation<DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for (unsigned i=0; i<DIM; i++)
    {
        location.push_back(mrCellPopulation.GetNode(rNodeIndex)->rGetLocation()[i]);
    }
    return location;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::Solve()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);

    if (current_time > 0) // use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }

    if (mOutputDirectory=="")
    {
        EXCEPTION("OutputDirectory not set");
    }

    double time_now = p_simulation_time->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    mSimulationOutputDirectory = results_directory;

    ///////////////////////////////////////////////////////////
    // Set up Simulation
    ///////////////////////////////////////////////////////////

    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", true);

    mrCellPopulation.CreateOutputFiles(results_directory+"/", false);

    mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

    SetupSolve();

    // Age the cells to the correct time. Note that cells are created with
    // negative birth times so that some are initially almost ready to divide.
    LOG(1, "Setting up cells...");
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // We don't use the result; this call is just to force the cells to age
        // to the current time running their cell-cycle models to get there
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");

    // Write initial conditions to file for the visualizer
    WriteVisualizerSetupFile();

    *mpVizSetupFile << std::flush;

    mrCellPopulation.WriteResultsToFiles();

    OutputSimulationSetup();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);

    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////

    while ((p_simulation_time->GetTimeStepsElapsed() < num_time_steps) && !(StoppingEventHasOccurred()) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

        /**
         * This function calls:
         * DoCellRemoval()
         * DoCellBirth()
         * CellPopulation::Update()
         */
        UpdateCellPopulation();

        /////////////////////////////////////
        // Update Cell Locations and topology
        /////////////////////////////////////
        //CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
        UpdateCellLocationsAndTopology();
        //CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);

        //////////////////////////////////////////
        // PostSolve, which may be implemented by
        // child classes (e.g. to solve PDEs)
        //////////////////////////////////////////
        PostSolve();

        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();

        ////////////////////////////
        // Output current results
        ////////////////////////////
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

        // Write results to file
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple == 0)
        {
            mrCellPopulation.WriteResultsToFiles();
        }

        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    }

    LOG(1, "--END TIME = " << SimulationTime::Instance()->GetTime() << "\n");
    /* Carry out a final update so that cell population is coherent with new cell positions.
     * NB cell birth/death still need to be checked because they may be spatially-dependent.*/
    UpdateCellPopulation();

    AfterSolve();

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

    mrCellPopulation.CloseOutputFiles();

    *mpVizSetupFile << "Complete\n";
    mpVizSetupFile->close();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
}

template<unsigned DIM>
bool AbstractCellBasedSimulation<DIM>::StoppingEventHasOccurred()
{
    return false;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::UpdateCellPopulation()
{
    /////////////////////////
    // Remove dead cells
    /////////////////////////
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::DEATH);
    unsigned deaths_this_step = DoCellRemoval();
    mNumDeaths += deaths_this_step;
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::DEATH);

    /////////////////////////
    // Divide cells
    /////////////////////////
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::BIRTH);
    unsigned births_this_step = DoCellBirth();
    mNumBirths += births_this_step;
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::BIRTH);

    // This allows NodeBasedCellPopulation::Update() to do the minimum amount of work
    bool births_or_death_occurred = ((births_this_step>0) || (deaths_this_step>0));

    ////////////////////////////
    // Update topology of cell population
    ////////////////////////////
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
    if (mUpdateCellPopulation)
    {
        LOG(1, "\tUpdating cell population...");
        mrCellPopulation.Update(births_or_death_occurred);
        LOG(1, "\tdone.\n");
    }
    else if (births_or_death_occurred)
    {
        EXCEPTION("CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::OutputSimulationSetup()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);

    // Output machine information
    ExecutableSupport::WriteMachineInfoFile(this->mSimulationOutputDirectory + "/system_info");

    // Output Chaste provenance information
    out_stream build_info_file = output_file_handler.OpenOutputFile("build.info");
    ExecutableSupport::WriteLibraryInfo(build_info_file);
    build_info_file->close();

    // Output simulation parameter and setup details
    out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

    // Output simulation details
    std::string simulation_type = GetIdentifier();

    *parameter_file << "<Chaste>\n";
    *parameter_file << "\n\t<" << simulation_type << ">\n";
    OutputSimulationParameters(parameter_file);
    *parameter_file << "\t</" << simulation_type << ">\n";
    *parameter_file << "\n";

    // Output cell population details (includes cell-cycle model details)
    mrCellPopulation.OutputCellPopulationInfo(parameter_file);

    // Loop over cell killers
    *parameter_file << "\n\t<CellKillers>\n";
    for (typename std::vector<AbstractCellKiller<DIM>*>::iterator iter = mCellKillers.begin();
         iter != mCellKillers.end();
         ++iter)
    {
        // Output cell killer details
        (*iter)->OutputCellKillerInfo(parameter_file);
    }
    *parameter_file << "\t</CellKillers>\n";

    // This is used to output information about subclasses.
    OutputAdditionalSimulationSetup(parameter_file);

    *parameter_file << "\n</Chaste>\n";
    parameter_file->close();
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<Dt>" << mDt << "</Dt>\n";
    *rParamsFile << "\t\t<EndTime>" << mEndTime << "</EndTime>\n";
    *rParamsFile << "\t\t<SamplingTimestepMultiple>" << mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
}

////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////

template class AbstractCellBasedSimulation<1>;
template class AbstractCellBasedSimulation<2>;
template class AbstractCellBasedSimulation<3>;