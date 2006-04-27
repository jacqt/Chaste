#ifndef REPLICATABLEVECTOR_HPP_
#define REPLICATABLEVECTOR_HPP_

#include <vector>
#include <petscvec.h>

class ReplicatableVector
{
private:
    /**
     * The wrapped vector.
     */
    std::vector<double> mData;

public:
    /**
     * Default constructor.
     * Note that the vector will need to be resized before it can be used.
     */
    ReplicatableVector()
    {
    }

    /**
     * Constructor to make a vector of given size.
     */
    ReplicatableVector(unsigned size)
    {
        resize(size);
    }

    /**
     * Return the size of the vector.
     */
    unsigned size(void)
    {
        return mData.size();
    }
    
    /**
     * Resize the vector.
     * 
     * @param size  The number of elements to allocate memory for.
     */
    void resize(unsigned size)
    {
        mData.resize(size);
    }
    
    /**
     * Access the vector.
     */
    double& operator[](unsigned index)
    {
        return mData[index];
    }
    
    /**
     * Replicate the given vector over all processes.
     * 
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes, storing it in this object.
     * 
     * This object must be resized to the size of the global vector before this method
     * is called.
     * 
     * @param lo  The start of our ownership range
     * @param hi  One past the end of our ownership range
     * @param input_array  The local portion of the array to be replicated.  Should
     *    contain hi-lo entries.  (If your input vector is larger, just pass the address
     *    of the first entry in your ownership range.)
     */
    void ReplicateVector(unsigned lo, unsigned hi, double *input_array)
    {
        unsigned size = mData.size();
        // Set up an array for MPI replication to use
        double input_vector_local_array[size];
        for (unsigned global_index=0; global_index<size; global_index++)
        {
            if (lo <= global_index && global_index < hi)
            { 
                unsigned local_index = global_index - lo;
                input_vector_local_array[global_index] = input_array[local_index]; 
            } 
            else 
            {
                input_vector_local_array[global_index] = 0.0;
            }
        }
        
        // Replicate
        MPI_Allreduce(input_vector_local_array, &(mData[0]), size,
                      MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }

    /**
     * Replicate this vector over all processes.
     * 
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes.
     * 
     * @param lo  The start of our ownership range
     * @param hi  One past the end of our ownership range
     */
    void Replicate(unsigned lo, unsigned hi)
    {
        ReplicateVector(lo, hi, &(mData[lo]));
    }

    /**
     * Replicate the given PETSc vector over all processes.
     * 
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes, storing it in this object.
     * 
     * Our data vector will automatically be resized to fit the whole PETSc vector.
     * 
     * @param vec  The PETSc vector to replicate.
     */
    void ReplicatePetscVector(Vec vec)
    {
        int lo, hi, size;  
        VecGetOwnershipRange(vec, &lo, &hi);
        VecGetSize(vec, &size);       
        double *input_array;
        VecGetArray(vec, &input_array);
        resize(size);
        ReplicateVector(lo, hi, input_array);
        VecRestoreArray(vec, &input_array);
    }
};

#endif /*REPLICATABLEVECTOR_HPP_*/
