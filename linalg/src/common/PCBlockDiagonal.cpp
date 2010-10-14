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

#include "PCBlockDiagonal.hpp"
#include "Exception.hpp"

		//#include "Timer.hpp"

PCBlockDiagonal::PCBlockDiagonal(KSP& rKspObject)
{
    PCBlockDiagonalCreate(rKspObject);
    PCBlockDiagonalSetUp();
}

PCBlockDiagonal::~PCBlockDiagonal()
{
    MatDestroy(mPCContext.A11_matrix_subblock);
    MatDestroy(mPCContext.A22_matrix_subblock);

    PCDestroy(mPCContext.PC_amg_A11);
    PCDestroy(mPCContext.PC_amg_A22);

    VecDestroy(mPCContext.x1_subvector);
    VecDestroy(mPCContext.y1_subvector);

    VecDestroy(mPCContext.x2_subvector);
    VecDestroy(mPCContext.y2_subvector);
    
    VecScatterDestroy(mPCContext.A11_scatter_ctx);
    VecScatterDestroy(mPCContext.A22_scatter_ctx);    
}

void PCBlockDiagonal::PCBlockDiagonalCreate(KSP& rKspObject)
{
    KSPGetPC(rKspObject, &mPetscPCObject);

    Mat system_matrix, dummy;
    MatStructure flag;
    KSPGetOperators(rKspObject, &system_matrix, &dummy, &flag);

    PetscInt num_rows, num_columns;
    MatGetSize(system_matrix, &num_rows, &num_columns);
    assert(num_rows==num_columns);

    PetscInt num_local_rows, num_local_columns;
    MatGetLocalSize(system_matrix, &num_local_rows, &num_local_columns);

    // odd number of rows: impossible in Bidomain.
    // odd number of local rows: impossible if V_m and phi_e for each node are stored in the same processor.
    if ((num_rows%2 != 0) || (num_local_rows%2 != 0))
    {
        TERMINATE("Wrong matrix parallel layout detected in PCLDUFactorisation.");
    }

    // Allocate memory     
    unsigned subvector_num_rows = num_rows/2;    
    unsigned subvector_local_rows = num_local_rows/2;
    mPCContext.x1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.x2_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.y1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.y2_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);

    // Create scatter contexts
    {
        // Needed by VecScatterCreate in order to find out parallel layout.
        Vec dummy_vec = PetscTools::CreateVec(num_rows, num_local_rows);

        IS A11_rows, A22_rows;        
        ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 2, &A11_rows);
        ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 1, 2, &A22_rows);   
        
        IS all_vector;    
        ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 1, &all_vector);
        
        VecScatterCreate(dummy_vec, A11_rows, mPCContext.x1_subvector, all_vector, &mPCContext.A11_scatter_ctx);    
        VecScatterCreate(dummy_vec, A22_rows, mPCContext.x2_subvector, all_vector, &mPCContext.A22_scatter_ctx);
        
        ISDestroy(A11_rows);
        ISDestroy(A22_rows);
        ISDestroy(all_vector);
        
        VecDestroy(dummy_vec);        
    }
            
    // Get matrix sublock A11        
    {
        // Work out local row range for subblock A11 (same as x1 or y1) 
        PetscInt low, high, global_size;
        VecGetOwnershipRange(mPCContext.x1_subvector, &low, &high);
        VecGetSize(mPCContext.x1_subvector, &global_size);        
        assert(global_size == num_rows/2);       
        
        IS A11_local_rows;
        IS A11_columns;
        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low, 2, &A11_local_rows);
        ISCreateStride(PETSC_COMM_WORLD, global_size, 0, 2, &A11_columns);
    
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1) //PETSc 3.1
        MatGetSubMatrix(system_matrix, A11_local_rows, A11_columns,
			MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A11_local_rows, A11_columns, PETSC_DECIDE, 
			MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#endif

    
        ISDestroy(A11_local_rows);
        ISDestroy(A11_columns);
    }

    // Get matrix sublock A22
    {
        // Work out local row range for subblock A22 (same as x2 or y2) 
        PetscInt low, high, global_size;
        VecGetOwnershipRange(mPCContext.x2_subvector, &low, &high);
        VecGetSize(mPCContext.x2_subvector, &global_size);        
        assert(global_size == num_rows/2);       
        
        IS A22_local_rows;
        IS A22_columns;
        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low+1, 2, &A22_local_rows);
        ISCreateStride(PETSC_COMM_WORLD, global_size, 1, 2, &A22_columns);
    
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1) //PETSc 3.1
        MatGetSubMatrix(system_matrix, A22_local_rows, A22_columns,
			MAT_INITIAL_MATRIX, &mPCContext.A22_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A22_local_rows, A22_columns, PETSC_DECIDE, 
			MAT_INITIAL_MATRIX, &mPCContext.A22_matrix_subblock);
#endif
    
        ISDestroy(A22_local_rows);
        ISDestroy(A22_columns);
    }

    // Register call-back function and its context
    PCSetType(mPetscPCObject, PCSHELL);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PCShellSetApply(mPetscPCObject, PCBlockDiagonalApply, (void*) &mPCContext);
#else
    // Register PC context so it gets passed to PCBlockDiagonalApply
    PCShellSetContext(mPetscPCObject, &mPCContext);
    // Register call-back function
    PCShellSetApply(mPetscPCObject, PCBlockDiagonalApply);
#endif

}

void PCBlockDiagonal::PCBlockDiagonalSetUp()
{
    // These options will get read by PCSetFromOptions
    PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
    PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");
    PetscOptionsSetValue("-pc_hypre_type", "boomeramg");

    // Set up amg preconditioner for block A11
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A11));
    PCSetType(mPCContext.PC_amg_A11, PCBJACOBI);
//     PCSetType(mPCContext.PC_amg_A11, PCHYPRE);
//     PCHYPRESetType(mPCContext.PC_amg_A11, "euclid");

    PCSetOperators(mPCContext.PC_amg_A11, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
    PCSetFromOptions(mPCContext.PC_amg_A11);
    PCSetUp(mPCContext.PC_amg_A11);

    // Set up amg preconditioner for block A22
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A22));

    /* Full AMG in the block */
    PCSetType(mPCContext.PC_amg_A22, PCHYPRE);
    PCHYPRESetType(mPCContext.PC_amg_A22, "boomeramg");

    //    PetscOptionsSetValue("-pc_hypre_type", "boomeramg");
    PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
    PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");
    PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type", "HMIS");
    //    PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all","Jacobi");
    //PetscOptionsSetValue("-pc_hypre_boomeramg_max_levels","10");
    //PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl", "1");
    //    PetscOptionsSetValue("-pc_hypre_boomeramg_print_statistics","");
    //    PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type","standard");

    //    PCHYPRESetType(mPCContext.PC_amg_A22, "euclid");



    /* Block Jacobi with AMG at each block */
    //     PCSetType(mPCContext.PC_amg_A22, PCBJACOBI);
    
    //     PetscOptionsSetValue("-sub_pc_type", "hypre");
    
    //     PetscOptionsSetValue("-sub_pc_hypre_type", "boomeramg");
    //     PetscOptionsSetValue("-sub_pc_hypre_boomeramg_max_iter", "1");
    //     PetscOptionsSetValue("-sub_pc_hypre_boomeramg_strong_threshold", "0.0");
    
    //     PetscOptionsSetValue("-pc_hypre_type", "boomeramg");
    //     PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
    //     PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");

    /* Additive Schwarz with AMG at each block */
//     PCSetType(mPCContext.PC_amg_A22, PCASM);
    
//     PetscOptionsSetValue("-pc_asm_type", "basic");
//     PetscOptionsSetValue("-pc_asm_overlap", "1");

//     PetscOptionsSetValue("-sub_ksp_type", "preonly");    

//     PetscOptionsSetValue("-sub_pc_type", "hypre");
    
//     PetscOptionsSetValue("-sub_pc_hypre_type", "boomeramg");
//     PetscOptionsSetValue("-sub_pc_hypre_boomeramg_max_iter", "1");
//     PetscOptionsSetValue("-sub_pc_hypre_boomeramg_strong_threshold", "0.0");
 

    PCSetOperators(mPCContext.PC_amg_A22, mPCContext.A22_matrix_subblock, mPCContext.A22_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
    PCSetFromOptions(mPCContext.PC_amg_A22);
    PCSetUp(mPCContext.PC_amg_A22);
}

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1) //PETSc 3.1
PetscErrorCode PCBlockDiagonalApply(PC pc_object, Vec x, Vec y)
{
  void* pc_context;

  PCShellGetContext(pc_object, &pc_context);   
#else
PetscErrorCode PCBlockDiagonalApply(void* pc_context, Vec x, Vec y)
{
#endif

    // Cast the context pointer to PCBlockDiagonalContext
    PCBlockDiagonal::PCBlockDiagonalContext* block_diag_context = (PCBlockDiagonal::PCBlockDiagonalContext*) pc_context;
    assert(block_diag_context!=NULL);

    /*
     * Scatter x = [x1 x2]'
     */
//PETSc-3.x.x or PETSc-2.3.3
//    Timer::Reset();
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A11_scatter_ctx, x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A11_scatter_ctx, x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A11_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A11_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_scatter_ctx, x, block_diag_context->x2_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A22_scatter_ctx, x, block_diag_context->x2_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x2_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x2_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_scatter_ctx);
#endif
    //    Timer::PrintAndReset("scatter");

    /*
     *  y1 = AMG(A11)*x1
     *  y2 = AMG(A22)*x2
     */
    PCApply(block_diag_context->PC_amg_A11, block_diag_context->x1_subvector, block_diag_context->y1_subvector);
    //    Timer::PrintAndReset("first prec");    
    PCApply(block_diag_context->PC_amg_A22, block_diag_context->x2_subvector, block_diag_context->y2_subvector);
    //    Timer::PrintAndReset("second prec");        

    /*
     * Gather y = [y1 y2]'
     */
//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A11_scatter_ctx, block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A11_scatter_ctx, block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A11_scatter_ctx);
    VecScatterEnd(block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A11_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_scatter_ctx, block_diag_context->y2_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A22_scatter_ctx, block_diag_context->y2_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y2_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_scatter_ctx);
    VecScatterEnd(block_diag_context->y2_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_scatter_ctx);
#endif
    //    Timer::PrintAndReset("gather");

    return 0;
}
