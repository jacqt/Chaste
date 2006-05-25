#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_


#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "BoundaryConditionsContainer.hpp"
#include "BidomainDg0Assembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "BidomainPde.hpp"
#include "AbstractCardiacCellFactory.hpp"


/**
 * Class which specifies and solves a bidomain problem.
 * 
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N), 
 * where V_j is the voltage at node j and phi_j is the 
 * extracellular potential at node j.
 */
template<int SPACE_DIM>
class BidomainProblem
{
private:
    std::string mMeshFilename;

    /** 
     *  Start time defaults to 0, pde timestep defaults to 0.01 (ms), the
     *  end time is not defaulted and must be set
     */
    double mStartTime;
    double mEndTime;
    double mPdeTimeStep;

    /** data is not written if output directory or output file prefix are not set*/ 
    std::string  mOutputDirectory, mOutputFilenamePrefix;  

    BidomainPde<SPACE_DIM>* mpBidomainPde;

    AbstractCardiacCellFactory<SPACE_DIM>* mpCellFactory;
    
    Vec mVoltage; // Current solution
    int mLo, mHi; // Ownership range for the double size potential vector mVoltage

    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> mMesh;
    
public:
    
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should 
     * create cells.
     */
    BidomainProblem(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory)
    : mMeshFilename(""),     // i.e. undefined
      mOutputDirectory(""),  // i.e. undefined
      mOutputFilenamePrefix(""),   // i.e. undefined
      mpBidomainPde(NULL)
    {
        mpCellFactory = pCellFactory;
                
        mStartTime   = 0.0;  // ms
        mPdeTimeStep = 0.01; // ms
        mEndTime     = -1;   // so can check has been set
        
        //initialise these to -1
        mLo = -1;
        mHi = -1; 
    }

    /**
     * Destructor
     */
    ~BidomainProblem()
    { 
        delete mpBidomainPde;
    }
    
    
    /** Initialise the system. Must be called before Solve() */
    void Initialise()
    {
        assert( mMeshFilename!="" );

        TrianglesMeshReader mesh_reader(mMeshFilename);
        mMesh.ConstructFromMeshReader(mesh_reader);
        
        mpCellFactory->SetMesh( &mMesh );
        
        mpBidomainPde = new BidomainPde<SPACE_DIM>( mpCellFactory, mStartTime, mPdeTimeStep );
    }
     

    /**
     * Solve the problem
     */
    void Solve()
    {
        assert( mpBidomainPde != NULL );  // if pde is NULL, Initialise() probably hasn't been called
        assert( mStartTime < mEndTime );
        
        // Linear solver
        SimpleLinearSolver linear_solver;
    
        // Assembler
        BidomainDg0Assembler<SPACE_DIM,SPACE_DIM> bidomain_assembler(mpBidomainPde, &mMesh, &linear_solver);
        
        // initial condition;   
        Vec initial_condition;
        int lo, hi;
        mpBidomainPde->GetOwnershipRange(lo, hi);
        VecCreateMPI(PETSC_COMM_WORLD, 2*(hi-lo) , 2*mMesh.GetNumNodes(), &initial_condition);  
        double* p_initial_condition;
        VecGetArray(initial_condition, &p_initial_condition); 

        VecGetOwnershipRange(initial_condition, &mLo, &mHi);
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_initial_condition[2*local_index  ] = mpBidomainPde->GetCardiacCell(global_index)->GetVoltage();
            p_initial_condition[2*local_index+1] = 0;
        }


        VecRestoreArray(initial_condition, &p_initial_condition);      
        VecAssemblyBegin(initial_condition);
        VecAssemblyEnd(initial_condition);
        //VecView(initial_condition, PETSC_VIEWER_STDOUT_WORLD);
    
        //  Write data to a file <mOutputFilenamePrefix>_xx.dat, 'xx' refers to 
        //  'xx'th time step using ColumnDataWriter         
        ///\todo: get writer to write V_m and \phi_e seperately?
        ParallelColumnDataWriter *p_test_writer = NULL;
       
        int time_var_id = 0;
        int voltage_var_id = 0;
        bool write_files = false;

        if (mOutputFilenamePrefix.length() > 0)
        {
            write_files = true;

            p_test_writer = new ParallelColumnDataWriter(mOutputDirectory,mOutputFilenamePrefix);

            p_test_writer->DefineFixedDimension("Node", "dimensionless", 2*mMesh.GetNumNodes() );
            time_var_id = p_test_writer->DefineUnlimitedDimension("Time","msecs");
        
            voltage_var_id = p_test_writer->DefineVariable("Vm_And_Phi_e","mV");
            p_test_writer->EndDefineMode();
        }

        double current_time = mStartTime;        
        int big_steps = 0; 

        if (write_files)
        {        
            p_test_writer->PutVariable(time_var_id, current_time); 
            p_test_writer->PutVector(voltage_var_id, initial_condition);
        }                        
                        
        while( current_time < mEndTime )
        {
            bidomain_assembler.SetTimes(current_time, current_time+mPdeTimeStep, mPdeTimeStep);
            bidomain_assembler.SetInitialCondition( initial_condition );
            
            mVoltage = bidomain_assembler.Solve(); //(mMesh, mpBidomainPde, bcc);
                                    
            // Free old initial condition
            VecDestroy(initial_condition);

            // Initial condition for next loop is current solution
            initial_condition = mVoltage;
            
            // Writing data out to the file <mOutputFilenamePrefix>.dat                 
            if (write_files)
            {
                p_test_writer->AdvanceAlongUnlimitedDimension(); //creates a new file
                p_test_writer->PutVariable(time_var_id, current_time); 
                p_test_writer->PutVector(voltage_var_id, mVoltage);
            }

            mpBidomainPde->ResetAsUnsolvedOdeSystem();

            current_time += mPdeTimeStep;
                
            big_steps++;
        }

        // close the file that stores voltage values            
        if (write_files)
        {
            p_test_writer->Close();
            delete p_test_writer;
        }
    }
    
    
    void SetStartTime(const double &rStartTime)
    {
        mStartTime = rStartTime;
    }

    void SetEndTime(const double &rEndTime)
    {
        mEndTime = rEndTime;
    }

    void SetPdeTimeStep(double pdeTimeStep)
    {
        assert(0.0 < pdeTimeStep);
        mPdeTimeStep=pdeTimeStep;
    }
    
    double GetPdeTimeStep()
    {
        return mPdeTimeStep;   
    }

    void SetMeshFilename(const std::string &rMeshFilename)
    {
        mMeshFilename = rMeshFilename;
    }
        
    void SetOutputDirectory(const std::string &rOutputDirectory)
    {
        mOutputDirectory = rOutputDirectory;
    }
        
    void SetOutputFilenamePrefix(const std::string &rOutputFilenamePrefix)
    {
        mOutputFilenamePrefix = rOutputFilenamePrefix;
    }
    
 
    /** Get the final solution vector. This is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N). 
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     * 
     *  This vector is distributed over all processes,
     *  with the current process owning the [lo, ..., hi-1] components of the vector.
     */
    void GetVoltageArray(double **pVoltageArray, int &lo, int &hi)
    {
        //check these have been set
        assert(mLo > -1);
        assert(mHi > -1);        
        VecGetArray(mVoltage, pVoltageArray);
        lo=mLo;
        hi=mHi; 
    }
    
    /** call this after GetVoltageArray to avoid memory leaks*/
    void RestoreVoltageArray(double **pVoltageArray)
    {
       VecRestoreArray(mVoltage, pVoltageArray);      
       VecAssemblyBegin(mVoltage);
       VecAssemblyEnd(mVoltage);
       VecDestroy(mVoltage);
    }
    
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM> & rGetMesh()
    {
        return mMesh;   
    }
    
    /** 
     *  Get the pde. Can only be called after Initialise()
     */
    BidomainPde<SPACE_DIM>* GetBidomainPde() 
    {
        assert(mpBidomainPde!=NULL);
        return mpBidomainPde;  
    }
};


#endif /*BIDOMAINPROBLEM_HPP_*/
