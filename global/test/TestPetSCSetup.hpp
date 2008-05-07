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


#ifndef _TESTPETSCSETUP_HPP_
#define _TESTPETSCSETUP_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

/**
 * This tests that the initialisation of PETSc does something.
 */


class TestPetSCSetup : public CxxTest::TestSuite
{
public:
    void TestPetscIsThere()
    {
        PetscTruth is_there;
        PetscInitialized(&is_there);
        TS_ASSERT( is_there == PETSC_TRUE );
    }
    
    void TestPetscExceptions()
    {
        int err=0;
        TS_ASSERT_THROWS_NOTHING(PETSCEXCEPT(err));
        
        Vec v;
        err = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, -1, &v);
        VecDestroy(v);
        //#define PETSC_ERR_ARG_WRONGSTATE   73   /* object in argument is in wrong */
        TS_ASSERT_THROWS_ANYTHING(PETSCEXCEPT(err));
        
        err=PETSC_ERR_FILE_OPEN;
        //#define PETSC_ERR_FILE_OPEN        65   /* unable to open file */
        TS_ASSERT_THROWS_ANYTHING(PETSCEXCEPT(err));
        
        //See if we can do it without a temporary
        TS_ASSERT_THROWS_ANYTHING(
            PETSCEXCEPT(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, -1, &v)));
        VecDestroy(v);
            
        //This test give back an "unknown error" message
        TS_ASSERT_THROWS_ANYTHING( PETSCEXCEPT(-3));
    }
    
    
    void TestKspExceptionsForCoverage()
    {
        TS_ASSERT_THROWS_NOTHING(  KSPEXCEPT(2) );
        //These next few lines are designed to force the coverage test to pass.
        //Some are hard to throw in normal circumstances --
        //"Unknown KSP error code" ought never to be thrown.
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_ITS) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_DTOL) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_BREAKDOWN_BICG) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_NONSYMMETRIC) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(KSP_DIVERGED_INDEFINITE_PC) );
        TS_ASSERT_THROWS_ANYTHING( KSPEXCEPT(-735827) );
    }
};



#endif // _TESTPETSCSETUP_HPP_
