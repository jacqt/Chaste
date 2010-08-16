/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef PETSCTOOLS_HPP_
#define PETSCTOOLS_HPP_

/**
 * @file
 * Contains the PetscTools class, and the TERMINATE and NEVER_REACHED
 * macros.  The macros are defined here (rather than Exception.hpp)
 * because they rely on MPI functionality to operate correctly when
 * running in parallel.
 */

#include <string>
#include <vector>
#include <cstdlib> // For EXIT_FAILURE

#ifndef SPECIAL_SERIAL
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#endif //SPECIAL_SERIAL

#define EXIT_IF_PARALLEL if(!PetscTools::IsSequential()){TS_TRACE("This test does not pass in parallel yet.");return;}
#define EXIT_IF_SEQUENTIAL if(PetscTools::IsSequential()){TS_TRACE("This test is not meant to be executed in sequential.");return;}


#ifndef NDEBUG
// Uncomment this to trace calls to PetscTools::Barrier
//#define DEBUG_BARRIERS
#endif


/**
 * A helper class of static methods.
 */
class PetscTools
{
private:
    /** Whether PETSc has been initialised. */
    static bool mPetscIsInitialised;

    /** The total number of processors. */
    static unsigned mNumProcessors;

#ifdef DEBUG_BARRIERS
    /** Used to debug number of barriers */
    static unsigned mNumBarriers;
#endif

    /** Which processors we are. */
    static unsigned mRank;

    /** Private method makes sure that (if this is the first use within a test) then PETSc has been probed */
    static inline void CheckCache()
    {
        if (mNumProcessors == 0)
        {
            ResetCache();
        }
    }

public:

    /**
     * As a convention, we consider processor 0 the master process
     */
    static const unsigned MASTER_RANK=0;

    /**
     * Reset our cached values: whether PETSc is initialised,
     * how many processors there are, and which one we are.
     */
    static void ResetCache();

    /**
     * Just returns whether there is one process or not.
     */
    static bool IsSequential();

    /**
     *  Returns total number of processors
     */
    static unsigned GetNumProcs();

    /**
     * Return our rank.
     *
     * If PETSc has not been initialized, returns 0.
     */
    static unsigned GetMyRank();

    /**
     * Just returns whether it is the master process or not.
     *
     * If not running in parallel, always returns true.
     */
    static bool AmMaster();

    /**
     * Just returns whether it is the right-most process or not.
     *
     * If not running in parallel, always returns true.
     */
    static bool AmTopMost();

    /**
     * If MPI is set up, perform a barrier synchronisation.
     * If not, it's a noop.
     *
     * @param callerId  only used in debug mode; printed before & after the barrier call
     */
    static void Barrier(const std::string callerId="");

#ifndef SPECIAL_SERIAL
    /**
     * Create a vector of the specified size. SetFromOptions is called.
     *
     * @param size  the size of the vector
     * @param localSize  the local number of items owned by this process
     * 
     */
    static Vec CreateVec(int size, int localSize=PETSC_DECIDE);

    /**
     * Create a Vec from the given data.
     *
     * @param data  some data
     */
    static Vec CreateVec(std::vector<double> data);

    /**
     * Create a vector of the specified size with all values set to be the given
     * constant. SetFromOptions is called.
     *
     * @param size  the size of the vector
     * @param value  the value to set each entry
     */
    static Vec CreateAndSetVec(int size, double value);

    /**
     * Set up a matrix - set the size using the given parameters. The number of local rows 
     * and columns is by default PETSC_DECIDE. SetFromOptions is called.
     *
     * @param rMat the matrix
     * @param numRows the number of rows in the matrix
     * @param numColumns the number of columns in the matrix
     * @param rowPreallocation the max number of nonzero entries expected on a row
     *   A value of 0 is allowed: no preallocation is then done and the user must 
     *   preallocate the memory for the matrix themselves.
     * @param numLocalRows the number of local rows (defaults to PETSC_DECIDE)
     * @param numLocalColumns the number of local columns (defaults to PETSC_DECIDE)
     */
    static void SetupMat(Mat& rMat, int numRows, int numColumns,
                         unsigned rowPreallocation,
                         int numLocalRows=PETSC_DECIDE,
                         int numLocalColumns=PETSC_DECIDE);

    /**
     * Boolean AND of a flags between processes
     *
     * @param flag is set to true on this process.
     * @return true if any process' flag is set to true.
     */
    static bool ReplicateBool(bool flag);

    /**
     * Ensure exceptions are handled cleanly in parallel code, by causing all processes to
     * throw if any one does.
     *
     * @param flag is set to true if this process has thrown.
     */
    static void ReplicateException(bool flag);

    /**
     * Dumps a given Petsc object to disk.
     *
     * @param rMat a matrix
     * @param rOutputFileFullPath where to dump the matrix to disk
     */
    static void DumpPetscObject(const Mat& rMat, const std::string& rOutputFileFullPath);

    /**
     * Dumps a given Petsc object to disk.
     *
     * @param rVec a vector
     * @param rOutputFileFullPath where to dump the vector to disk
     */
    static void DumpPetscObject(const Vec& rVec, const std::string& rOutputFileFullPath);

    /**
     * Read a previously dumped Petsc object from disk.
     *
     * @param rMat a matrix
     * @param rOutputFileFullPath where to read the matrix from
     * @param rParallelLayout If provided, rMat will have the same parallel layout. Its content is irrelevant. 
     */
    static void ReadPetscObject(Mat& rMat, const std::string& rOutputFileFullPath, Vec rParallelLayout=NULL);

    /**
     * Read a previously dumped Petsc object from disk.
     *
     * @param rVec a vector
     * @param rOutputFileFullPath where to read the matrix from
     * @param rParallelLayout If provided, rMat will have the same parallel layout. Its content is irrelevant. 
     */
    static void ReadPetscObject(Vec& rVec, const std::string& rOutputFileFullPath, Vec rParallelLayout=NULL);
    
#endif //SPECIAL_SERIAL (ifndef)

    /**
     * Level 4 error (Termination).  Execution cannot continue from this point and hence
     * should be terminated (even when running with NDEBUG).
     * 
     * @param rMessage An error message to appear on the screen
     * @param rFilename  which source file produced the termination error
     * @param lineNumber  which line number of the source file produced the termination error
     */
    static void Terminate(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber);


};

#define TERMINATE(message) PetscTools::Terminate(message, __FILE__, __LINE__)

/*
 *  The exit statement at the end of NEVER_REACHED is not really needed but prevents g++ from complaining about 
 * uninitialised variables when you have code that looks like:
 * 
 *   RelativeTo::Value relative_to;
 *   switch (rPath.relative_to())
 *   {
 *       case cp::relative_to_type::cwd:
 *           relative_to = RelativeTo::CWD;
 *           break;
 *       case cp::relative_to_type::chaste_test_output:
 *           relative_to = RelativeTo::ChasteTestOutput;
 *           break;
 *       case cp::relative_to_type::chaste_source_root:
 *           relative_to = RelativeTo::ChasteSourceRoot;
 *           break;
 *       case cp::relative_to_type::absolute:
 *           relative_to = RelativeTo::Absolute;
 *           break;
 *       default:
 *           NEVER_REACHED;
 *           break;
 *   }
 * 
 *  relative_to is considered potentially uninitialised in the default branch unless the compiler finds a exit, 
 * assert or throw statement. 
 */
#define NEVER_REACHED  PetscTools::Terminate("Should have been impossible to reach this line of code", __FILE__, __LINE__); exit(EXIT_FAILURE)

#endif /*PETSCTOOLS_HPP_*/
