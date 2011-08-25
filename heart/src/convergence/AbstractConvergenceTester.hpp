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


#ifndef ABSTRACTCONVERGENCETESTER_HPP_
#define ABSTRACTCONVERGENCETESTER_HPP_

#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"

#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "AbstractTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "Hdf5DataReader.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "CuboidMeshConstructor.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "ConstBoundaryCondition.hpp"
#include "StimulusBoundaryCondition.hpp"

#include "Warnings.hpp"

typedef enum StimulusType_
{
    PLANE=0,
    QUARTER,
    NEUMANN
} StimulusType;

/**
 * RampedQuarterStimulusCellFactory stimulates a quarter of a mesh of width mMeshWidth
 * ie all the cells in 0 < x <= mMeshWidth/4
 */
template <class CELL, unsigned DIM>
class RampedQuarterStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    /** define a new vector of stimuli - one for each step in x-direction*/
    std::vector< boost::shared_ptr<SimpleStimulus> > mpStimuli;
    /** Width (x-width) of mesh*/
    double mMeshWidth;
    /** Step size of mesh is derived from the width and the number of elements across*/
    double mStepSize;
    /** The number of stimulated levels covering 0 <= x < mMeshWidth/4.
     * Note that the nodes on the quarter level are not included and are unstimulated
     */
    unsigned mLevels;
public:

    /**
     * Constructor.
     * @param meshWidth x-width of mesh
     * @param numElemAcross this allows us to deduce the mesh step size.
     * @param fullStim  the maximum stimulus level
     */
    RampedQuarterStimulusCellFactory(double meshWidth, unsigned numElemAcross, double fullStim=-1e7)
        : AbstractCardiacCellFactory<DIM>(),
          mMeshWidth(meshWidth),
          mStepSize(meshWidth/numElemAcross),
          mLevels(numElemAcross/4)
    {
        assert(numElemAcross%4 == 0); //numElemAcross is supposed to be a multiple of 4

        for (unsigned level=0; level<mLevels; level++)
        {
            double this_stim = fullStim - (level*fullStim)/mLevels;
            //this_stim is full_stim at the zero level and would be zero at level=mLevels
            mpStimuli.push_back((boost::shared_ptr<SimpleStimulus>)new SimpleStimulus(this_stim, 0.5));
        }
    }


    /**
     * Create cell model
     *
     * @param node Global node index
     */
    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        double x = this->GetMesh()->GetNode(node)->GetPoint()[0];
        double d_level = x/mStepSize;
        unsigned level = (unsigned) d_level;
        assert(fabs(level-d_level) < DBL_MAX); //x ought to really be a multiple of the step size

        if (level < mLevels)
        {
            return new CELL(this->mpSolver, this->mpStimuli[level]);
        }
        else
        {
            return new CELL(this->mpSolver, this->mpZeroStimulus);
        }
    }
};


/**
 * AbstractUntemplatedConvergenceTester
 * contains core functionality used in more convergence testers
 */
class AbstractUntemplatedConvergenceTester
{
protected:
    /** Mesh width (for cuboid mesh)*/
    double mMeshWidth;
public:
    /** OdeTimeStep to be varied in OdeConvergenceTester etc*/
    double OdeTimeStep;
    /** PdeTimeStep to be varied in PdeConvergenceTester etc*/
    double PdeTimeStep;
    /** Mesh number - mesh 0 has 4 elements in each space dimension  0.05cm on a 0.2cm mesh
     *              - mesh 1 has 8 (0.025cm)
     *              - mesh 2 has 16 (0.0125cm)
     */
    unsigned MeshNum;
    double RelativeConvergenceCriterion; /**< Main convergence test is LastDifference < RelativeConvergenceCriterion */
    double LastDifference; /**< Used to store and retrieve the difference between success runs (relative 2-norm of Vm at 3rd quarter node).*/
    double Apd90FirstQn; /**< Used to store and retrieve the APD90 of a node in the first quarter x-value*/
    double Apd90ThirdQn; /**< Used to store and retrieve the APD90 of a node in the third quarter x-value*/
    double ConductionVelocity;  /**< Used to store and retrieve the conduction velocity between the first & third quarter nodes*/
    /** Set to true once the result of a test is known (either read from file
     * or produced by the coarsest run of the tester).
     */
    bool PopulatedResult;
    /** true if converging to a known standard result
     *  used in StimulusConverger in projects/jmpf
     */
    bool FixedResult;
    /** true if the plane stimulus should applied with an exact value (not scaled by space-step)
     *  used in StimulusConverger in projects/jmpf
     */
    bool UseAbsoluteStimulus;
    /**
     * A value to be used with a plane stimulus (not scaled by space-step)
     *  used in StimulusConverger in projects/jmpf
     */
    double AbsoluteStimulus;
    bool SimulateFullActionPotential; /**< Set it true in order to run simulations for long enough to get a whole AP (and thus get convergence history for APD90).  Note that this slackens the relative L2 norm criterion, since we are in plateau phase for longer*/
    bool Converged; /**< Set to true when convergence has been reached */
    StimulusType Stimulus; /**< The type of stimulus: PLANE (x=0), REGION (first quarter in x) or NEUMANN (monodomain only)*/
    double NeumannStimulus; /**< Quantity of face stimulus to use in the Neumann case */

    AbstractUntemplatedConvergenceTester();

    /**
     * Run the same test at different levels of refinement until
     * some convergence criterion is met.
     * @param nameOfTest The name of the convergence test (typically the name in the suite) for use in naming files.
     */
    virtual void Converge(std::string nameOfTest)=0;

    virtual ~AbstractUntemplatedConvergenceTester();
};

/**
 * Use template specialization to set the appropriate conductivities for the problem.
 */
template<class CARDIAC_PROBLEM, unsigned DIM>
void SetConductivities(CARDIAC_PROBLEM& rCardiacProblem);

template<unsigned DIM>
void SetConductivities(BidomainProblem<DIM>& rProblem)
{
    c_vector<double, DIM> conductivities;
    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 1.75;
    }
    HeartConfig::Instance()->SetIntracellularConductivities(conductivities);

    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 7.0;
    }
    HeartConfig::Instance()->SetExtracellularConductivities(conductivities);
}

template<unsigned DIM>
void SetConductivities(MonodomainProblem<DIM>& rProblem)
{
    c_vector<double, DIM> conductivities;
    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 1.75;
    }
    HeartConfig::Instance()->SetIntracellularConductivities(conductivities);
}

/**
 * AbstractConvergenceTester
 * Run convergence for a particular cell type, mono/bidomain and dimension
 */
template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class AbstractConvergenceTester : public AbstractUntemplatedConvergenceTester
{
public:
    /**
     * \todo This is a scarily long method; could do with some parts extracted?
     * @param nameOfTest The name of the convergence test (typically the name in the suite) for use in naming files.
     */
    void Converge(std::string nameOfTest)
    {
        std::cout << "=========================== Beginning Test...==================================\n";
        // Create the meshes on which the test will be based
        const std::string mesh_dir = "ConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);
        ReplicatableVector voltage_replicated;

        unsigned file_num=0;

        // Create a file for storing conduction velocity and AP data and write the header
        OutputFileHandler conv_info_handler("ConvergencePlots", false);
        out_stream p_conv_info_file;
        if (PetscTools::AmMaster())
        {
            p_conv_info_file = conv_info_handler.OpenOutputFile(nameOfTest+"_info.csv");
            (*p_conv_info_file) << "#Abcisa\t"
                                << "l2-norm-full\t"
                                << "l2-norm-onset\t"
                                << "Max absolute err\t"
                                << "APD90_1st_quad\t"
                                << "APD90_3rd_quad\t"
                                << "Conduction velocity (relative diffs)" << std::endl;
        }
        SetInitialConvergenceParameters();

        double prev_apd90_first_qn=0.0;
        double prev_apd90_third_qn=0.0;
        double prev_cond_velocity=0.0;
        std::vector<double> prev_voltage;
        std::vector<double> prev_times;
        PopulateStandardResult(prev_voltage, prev_times);

        do
        {
            CuboidMeshConstructor<DIM> constructor;

            //If the printing time step is too fine, then simulations become I/O bound without much improvement in accuracy
            double printing_step = this->PdeTimeStep;
#define COVERAGE_IGNORE
            while (printing_step < 1.0e-4)
            {
                printing_step *= 2.0;
                std::cout<<"Warning: PrintingTimeStep increased\n";
            }
#undef COVERAGE_IGNORE
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(this->OdeTimeStep, this->PdeTimeStep, printing_step);
#define COVERAGE_IGNORE
            if (SimulateFullActionPotential)
            {
                HeartConfig::Instance()->SetSimulationDuration(400.0);
            }
            else
            {
                HeartConfig::Instance()->SetSimulationDuration(8.0);
            }
#undef COVERAGE_IGNORE
            HeartConfig::Instance()->SetOutputDirectory("Convergence");
            HeartConfig::Instance()->SetOutputFilenamePrefix("Results");

            DistributedTetrahedralMesh<DIM, DIM> mesh;
            constructor.Construct(mesh, this->MeshNum, mMeshWidth);

            unsigned num_ele_across = (unsigned) SmallPow(2, this->MeshNum+2); // number of elements in each dimension

            AbstractCardiacCellFactory<DIM>* p_cell_factory=NULL;

            switch (this->Stimulus)
            {
                case NEUMANN:
                {
                    p_cell_factory = new ZeroStimulusCellFactory<CELL, DIM>();
                    break;
                }
                case PLANE:
                {
                    if (this->UseAbsoluteStimulus)
                    {
                        #define COVERAGE_IGNORE
                        p_cell_factory = new GeneralPlaneStimulusCellFactory<CELL, DIM>(0, this->AbsoluteStimulus, true);
                        #undef COVERAGE_IGNORE
                    }
                    else
                    {
                        p_cell_factory = new GeneralPlaneStimulusCellFactory<CELL, DIM>(num_ele_across, constructor.GetWidth(), false, this->AbsoluteStimulus);
                    }
                    break;
                }
                case QUARTER:
                {
                    ///\todo consider reducing all stimuli to match this one.
                    p_cell_factory = new RampedQuarterStimulusCellFactory<CELL, DIM>(constructor.GetWidth(), num_ele_across, this->AbsoluteStimulus/10.0);
                    break;
                }
            }


            CARDIAC_PROBLEM cardiac_problem(p_cell_factory);
            cardiac_problem.SetMesh(&mesh);

            // Calculate positions of nodes 1/4 and 3/4 through the mesh
            unsigned third_quadrant_node;
            unsigned first_quadrant_node;
            switch(DIM)
            {
                case 1:
                {
                    first_quadrant_node = (unsigned) (0.25*constructor.GetNumElements());
                    third_quadrant_node = (unsigned) (0.75*constructor.GetNumElements());
                    break;
                }
                case 2:
                {
                    unsigned n= (unsigned) SmallPow (2, this->MeshNum+2);
                    first_quadrant_node =   (n+1)*(n/2)+  n/4 ;
                    third_quadrant_node =   (n+1)*(n/2)+3*n/4 ;
                    break;
                }
                case 3:
                {
                    const unsigned first_quadrant_nodes_3d[5]={61, 362, 2452, 17960, 137296};
                    const unsigned third_quadrant_nodes_3d[5]={63, 366, 2460, 17976, 137328};
                    assert(this->PdeTimeStep<5);
                    first_quadrant_node = first_quadrant_nodes_3d[this->MeshNum];
                    third_quadrant_node = third_quadrant_nodes_3d[this->MeshNum];
                    break;
                }

                default:
                    NEVER_REACHED;
            }

            double mesh_width=constructor.GetWidth();

            // We only need the output of these two nodes
            std::vector<unsigned> nodes_to_be_output;
            nodes_to_be_output.push_back(first_quadrant_node);
            nodes_to_be_output.push_back(third_quadrant_node);
            cardiac_problem.SetOutputNodes(nodes_to_be_output);

            // The results of the tests were originally obtained with the following conductivity
            // values. After implementing fibre orientation the defaults changed. Here we set
            // the former ones to be used.
            SetConductivities(cardiac_problem);

            cardiac_problem.Initialise();

            boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM>);
            SimpleStimulus stim(NeumannStimulus, 0.5);
            if (Stimulus==NEUMANN)
            {

                StimulusBoundaryCondition<DIM>* p_bc_stim = new StimulusBoundaryCondition<DIM>(&stim);

                // get mesh
                AbstractTetrahedralMesh<DIM, DIM> &r_mesh = cardiac_problem.rGetMesh();
                // loop over boundary elements
                typename AbstractTetrahedralMesh<DIM, DIM>::BoundaryElementIterator iter;
                iter = r_mesh.GetBoundaryElementIteratorBegin();
                while (iter != r_mesh.GetBoundaryElementIteratorEnd())
                {
                    double x = ((*iter)->CalculateCentroid())[0];
                    ///\todo remove magic number? (#1884)
                    if (x*x<=1e-10)
                    {
                        p_bcc->AddNeumannBoundaryCondition(*iter, p_bc_stim);
                    }
                    iter++;
                }
                // pass the bcc to the problem
                cardiac_problem.SetBoundaryConditionsContainer(p_bcc);
            }

            DisplayRun();
            double time_before=MPI_Wtime();
            //// use this to get some info printed out
            //cardiac_problem.SetWriteInfo();

            try
            {
                cardiac_problem.Solve();
            }
            catch (Exception e)
            {
                WARNING("This run threw an exception.  Check convergence results\n");
                std::cout << e.GetMessage() << std::endl;
            }

            std::cout << "Time to solve = " << MPI_Wtime()-time_before << " seconds\n";

            OutputFileHandler results_handler("Convergence", false);
            Hdf5DataReader results_reader = cardiac_problem.GetDataReader();

            {
                std::vector<double> transmembrane_potential = results_reader.GetVariableOverTime("V", third_quadrant_node);
                std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();

                OutputFileHandler plot_file_handler("ConvergencePlots", false);
                if (PetscTools::AmMaster())
                {
                    // Write out the time series for the node at third quadrant
                    {
                        std::stringstream plot_file_name_stream;
                        plot_file_name_stream<< nameOfTest << "_Third_quadrant_node_run_"<< file_num << ".csv";
                        out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                        for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                        {
                            (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";
                        }
                        p_plot_file->close();
                    }
                    // Write time series for first quadrant node
                    {
                        std::vector<double> transmembrane_potential_1qd=results_reader.GetVariableOverTime("V", first_quadrant_node);
                        std::vector<double> time_series_1qd = results_reader.GetUnlimitedDimensionValues();
                        std::stringstream plot_file_name_stream;
                        plot_file_name_stream<< nameOfTest << "_First_quadrant_node_run_"<< file_num << ".csv";
                        out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                        for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                        {
                            (*p_plot_file) << time_series_1qd[data_point] << "\t" << transmembrane_potential_1qd[data_point] << "\n";
                        }
                        p_plot_file->close();
                    }
                }
                // calculate conduction velocity and APD90 error
                PropagationPropertiesCalculator ppc(&results_reader);

                try
                {
                    #define COVERAGE_IGNORE
                    if (SimulateFullActionPotential)
                    {
                        Apd90FirstQn = ppc.CalculateActionPotentialDuration(90.0, first_quadrant_node);
                        Apd90ThirdQn = ppc.CalculateActionPotentialDuration(90.0, third_quadrant_node);
                    }
                    #undef COVERAGE_IGNORE
                    ConductionVelocity  = ppc.CalculateConductionVelocity(first_quadrant_node,third_quadrant_node,0.5*mesh_width);
                }
                catch (Exception e)
                {
                    #define COVERAGE_IGNORE
                    std::cout << "Warning - this run threw an exception in calculating propagation.  Check convergence results\n";
                    std::cout << e.GetMessage() << std::endl;
                    #undef COVERAGE_IGNORE
                }
                double cond_velocity_error = 1e10;
                double apd90_first_qn_error = 1e10;
                double apd90_third_qn_error = 1e10;

                if (this->PopulatedResult)
                {
                    if (prev_cond_velocity != 0.0)
                    {
                        cond_velocity_error = fabs(ConductionVelocity - prev_cond_velocity) / prev_cond_velocity;
                    }
#define COVERAGE_IGNORE
                    if (prev_apd90_first_qn != 0.0)
                    {
                        apd90_first_qn_error = fabs(Apd90FirstQn - prev_apd90_first_qn) / prev_apd90_first_qn;
                    }
                    if (prev_apd90_third_qn != 0.0)
                    {
                        apd90_third_qn_error = fabs(Apd90ThirdQn - prev_apd90_third_qn) / prev_apd90_third_qn;
                    }
                    if (apd90_first_qn_error == 0.0)
                    {
                        apd90_first_qn_error = DBL_EPSILON; //Avoid log zero on plot
                    }
                    if (apd90_third_qn_error == 0.0)
                    {
                        apd90_third_qn_error = DBL_EPSILON; //Avoid log zero on plot
                    }
#undef COVERAGE_IGNORE
                    if (cond_velocity_error == 0.0)
                    {
                        cond_velocity_error = DBL_EPSILON; //Avoid log zero on plot
                    }
                }

                prev_cond_velocity = ConductionVelocity;
                prev_apd90_first_qn = Apd90FirstQn;
                prev_apd90_third_qn = Apd90ThirdQn;

                // calculate l2norm
                double max_abs_error = 0;
                double sum_sq_abs_error =0;
                double sum_sq_prev_voltage = 0;
                double sum_sq_abs_error_full =0;
                double sum_sq_prev_voltage_full = 0;
                if (this->PopulatedResult)
                {
                    //If the PDE step is varying then we'll have twice as much data now as we use to have

                    unsigned time_factor=(time_series.size()-1) / (prev_times.size()-1);
                    assert (time_factor == 1 || time_factor == 2 || time_factor == 8);
                    //Iterate over the shorter time series data
                    for (unsigned data_point = 0; data_point<prev_times.size(); data_point++)
                    {
                        unsigned this_data_point=time_factor*data_point;

                        assert(time_series[this_data_point] == prev_times[data_point]);
                        double abs_error = fabs(transmembrane_potential[this_data_point]-prev_voltage[data_point]);
                        max_abs_error = (abs_error > max_abs_error) ? abs_error : max_abs_error;
                        //Only do resolve the upstroke...
                        sum_sq_abs_error_full += abs_error*abs_error;
                        sum_sq_prev_voltage_full += prev_voltage[data_point] * prev_voltage[data_point];
                        if (time_series[this_data_point] <= 8.0)
                        {
                            //In most regular cases we always do this, since the simulation stops at ms
                            sum_sq_abs_error += abs_error*abs_error;
                            sum_sq_prev_voltage += prev_voltage[data_point] * prev_voltage[data_point];
                        }
                    }

                }
                if (!this->PopulatedResult || !FixedResult)
                {
                    prev_voltage = transmembrane_potential;
                    prev_times = time_series;
                }

                if (this->PopulatedResult)
                {

                    if (PetscTools::AmMaster())
                    {
                        (*p_conv_info_file) << std::setprecision(8)
                                            << Abscissa() << "\t"
                                            << sqrt(sum_sq_abs_error_full/sum_sq_prev_voltage_full) << "\t"
                                            << sqrt(sum_sq_abs_error/sum_sq_prev_voltage) << "\t"
                                            << max_abs_error << "\t"
                                            << Apd90FirstQn <<" "<< apd90_first_qn_error <<""<< "\t"
                                            << Apd90ThirdQn <<" "<< apd90_third_qn_error <<""<< "\t"
                                            << ConductionVelocity <<" "<< cond_velocity_error  <<""<< std::endl;
                    }
                    // convergence criterion
                    this->Converged = sum_sq_abs_error/sum_sq_prev_voltage<this->RelativeConvergenceCriterion;
                    this->LastDifference=sum_sq_abs_error/sum_sq_prev_voltage;
#define COVERAGE_IGNORE
                    if (time_series.size() == 1u)
                    {
                        std::cout << "Failed after successful convergence - give up this convergence test\n";
                        break;
                    }
#undef COVERAGE_IGNORE

                }

                if (!this->PopulatedResult)
                {
                    this->PopulatedResult=true;

                }
            }

            // Get ready for the next test by halving the time step
            if (!this->Converged)
            {
                UpdateConvergenceParameters();
                file_num++;
            }
            delete p_cell_factory;
        }
        while (!GiveUpConvergence() && !this->Converged);


        if (PetscTools::AmMaster())
        {
            p_conv_info_file->close();

            std::cout << "Results: " << std::endl;
            MPIABORTIFNON0(system, "cat " + conv_info_handler.GetOutputDirectoryFullPath() + nameOfTest + "_info.csv");
        }

    }

    void DisplayRun()
    {
        unsigned num_ele_across = (unsigned) SmallPow(2, this->MeshNum+2);// number of elements in each dimension
        double scaling = mMeshWidth/(double) num_ele_across;

        std::cout<<"================================================================================"<<std::endl;
        std::cout<<"Solving in "<<DIM<<"D\t";
        std::cout<<"Space step "<< scaling << " cm (mesh " << this->MeshNum << ")" << "\n";
        std::cout<<"PDE step "<<this->PdeTimeStep<<" ms"<<"\t";
        std::cout<<"ODE step "<<this->OdeTimeStep<<" ms"<<"\t";
        if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
        {
            std::cout<<"KSP absolute "<<HeartConfig::Instance()->GetAbsoluteTolerance()<<"\t";
        }
        else
        {
            std::cout<<"KSP relative "<<HeartConfig::Instance()->GetRelativeTolerance()<<"\t";
        }
        switch (this->Stimulus)
        {
            case PLANE:
            std::cout<<"Stimulus = Plane\n";
            break;

            case QUARTER:
            std::cout<<"Stimulus = Ramped first quarter\n";
            break;

            case NEUMANN:
            std::cout<<"Stimulus = Neumann\n";
            break;

        }
        EXPECT0(system, "date");//To keep track of what Nightly things are doing
        ///\todo The UseAbsoluteStimulus is temporary, while we are sorting out
        ///3D stimulus.  It is to be removed later (along with StimulusConvergenceTester)
        if (this->UseAbsoluteStimulus)
        {
            #define COVERAGE_IGNORE
            std::cout<<"Using absolute stimulus of "<<this->AbsoluteStimulus<<std::endl;
            #undef COVERAGE_IGNORE
        }
        std::cout << std::flush;
        //HeartEventHandler::Headings();
        //HeartEventHandler::Report();

    }

public:
    virtual ~AbstractConvergenceTester() {}

    /** Initial values of parameters at the beginning of the convergence test (the parameter to be varied will be larger than the expected value at convergence)*/
    virtual void SetInitialConvergenceParameters()=0;
    /** Update the parameter which is being varied*/
    virtual void UpdateConvergenceParameters()=0;
    /** Assess whether to abort the convergence test (convergence is unlikely to happen).*/
    virtual bool GiveUpConvergence()=0;
    /** The value of the parameter which is being varied*/
    virtual double Abscissa()=0;
    /** This is currently used as stub for convergence testers which need to converge
     * to a known standardised result (the StimulusConvergence tester in projects/jmpf).
     *
     * @param result a standard vector to be sized and filled with V_m values by this method (in subclass)
     * @param times a standard vector to be sized and filled with times values by this method (in subclass)
     */
    virtual void PopulateStandardResult(std::vector<double> &result, std::vector<double> &times)
    {
        assert(this->PopulatedResult==false);
        assert(result.size()==0);
        assert(times.size()==0);
    }

    /**
     * @return when the convergence criterion is met
     */
    bool IsConverged()
    {
        return Converged;
    }

    /**
     * @param meshWidth set the dimension of the cuboid mesh (default value is 0.2cm)
     */
    void SetMeshWidth(double meshWidth)
    {
        mMeshWidth=meshWidth;
    }
};

#endif /*ABSTRACTCONVERGENCETESTER_HPP_*/
