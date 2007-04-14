#ifndef TESTPARALLELVECTOR_HPP_
#define TESTPARALLELVECTOR_HPP_

#include <cxxtest/TestSuite.h>
#include <petscvec.h>
#include <petscmat.h>

#include "PetscSetupAndFinalize.hpp"
#include "PetscException.hpp"

#include "ParallelVector.hpp"

class TestParallelVector : public CxxTest::TestSuite
{
public:
    
    void TestReadAndRestore()
    {
        // create a 10 element petsc vector
        Vec vec;

        VecCreate(PETSC_COMM_WORLD, &vec);
        VecSetSizes(vec, PETSC_DECIDE, 10);
        VecSetFromOptions(vec);
        
        // calculate my range
        PetscInt petsc_lo, petsc_hi;
        VecGetOwnershipRange(vec,&petsc_lo,&petsc_hi);
        unsigned lo=(unsigned)petsc_lo;
        unsigned hi=(unsigned)petsc_hi;   
        
        // stuff the vector with values: global index * local index
        double* p_vec;
        VecGetArray(vec, &p_vec);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            p_vec[local_index] = local_index*global_index;
        }
        VecRestoreArray(vec, &p_vec);
        VecAssemblyBegin(vec);
        VecAssemblyEnd(vec);
        
        // create vec portion
        ParallelVector parallel_vector(vec);
        
        // check that the portion has the correct range
        TS_ASSERT_EQUALS(ParallelVector::Begin().Global,lo);
        TS_ASSERT_EQUALS(ParallelVector::End().Global,hi);
        
        for (ParallelVector::Iterator index = ParallelVector::Begin();
             index!= ParallelVector::End();
             ++index)
        {
            TS_ASSERT_EQUALS(parallel_vector[index], index.Local*index.Global);
        }
        
        // now write something to the vector - local_index * global_index
        for (ParallelVector::Iterator index = ParallelVector::Begin();
             index!= ParallelVector::End();
             ++index)
        {
            parallel_vector[index] =  -(double)(index.Local*index.Global);
        }
        
        TS_ASSERT( ! (ParallelVector::Begin() == ParallelVector::End()) );        
        
        // ask the portion to restore the main vector
        parallel_vector.Restore();
        
        // check that the vector was restored
        VecGetArray(vec, &p_vec);
        
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            unsigned local_index = global_index - lo;
            TS_ASSERT_EQUALS(p_vec[local_index], -(double)local_index*global_index);
        }

    }
           
        
};

#endif /*TESTPARALLELVECTOR_HPP_*/
