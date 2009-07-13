/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef PARALLELTETRAHEDRALMESH_HPP_
#define PARALLELTETRAHEDRALMESH_HPP_

#include <map>
#include <vector>
#include <set>

#include "AbstractTetrahedralMesh.hpp"
#include "Node.hpp"
#include "AbstractTetrahedralMeshReader.hpp"

/*
 *  The following definition fixes an odd incompatibility of METIS 4.0 and Chaste. Since
 * the library was compiled with a plain-C compiler, it fails to link using a C++ compiler.
 * Note that METIS 4.0 fails to compile with g++ or icpc, so a C compiler should be used.
 *
 * Somebody had this problem before: http://www-users.cs.umn.edu/~karypis/.discus/messages/15/113.html?1119486445
 *
 * Note that it is necessary to define the function header before the #include statement.
*/
extern "C" {
extern void METIS_PartMeshNodal(int*, int*, int*, int*, int*, int*, int*, int*, int*);
};
#include "metis.h"

/**
 * Parallel implementation of a mesh
 * Nodes are distributed such that each process has
 * A set of nodes (possibly reordered) with contiguous global indices
 * A local copy of all the elements supporting those nodes
 * A local copy of ghost/halo nodes which are all the nodes used in the supporting elements, but not owned outright.
 */ 
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ParallelTetrahedralMesh : public AbstractTetrahedralMesh< ELEMENT_DIM, SPACE_DIM>
{
    friend class TestParallelTetrahedralMesh;

public:

     /** Definition of partition types. 
      * "DUMB" is using naturally mesh ordering with PETSC_DECIDE.
      * "METIS_BINARY" via METIS file dump and a call to the Metis binary partdmesh.
      * "METIS_LIBRARY" is a call to the sequential METIS library 
      * */
    typedef enum
    {
        DUMB=0,
        METIS_BINARY,
        METIS_LIBRARY
    } PartitionType;

private:

    /** The total number of elements in the mesh. */
    unsigned mTotalNumElements;

    /** The total number of boundary elements in the mesh. */
    unsigned mTotalNumBoundaryElements;

    /** The total number of nodes in the mesh. */
    unsigned mTotalNumNodes;

    /** Vector of pointer to halo nodes used by this process. */
    std::vector<Node<SPACE_DIM>* > mHaloNodes;

    /** A map from node global index to local index used by this process. */
    std::map<unsigned, unsigned> mNodesMapping;

    /** A map from halo node global index to local index used by this process. */
    std::map<unsigned, unsigned> mHaloNodesMapping;

    /** A map from element global index to local index used by this process. */
    std::map<unsigned, unsigned> mElementsMapping;

    /** A map from boundary element global index to local index used by this process. */
    std::map<unsigned, unsigned> mBoundaryElementsMapping;

    /** Partition type (given by enum PartitionType). */
    PartitionType mMetisPartitioning;

public:

    /**
     * Constructor.
     * 
     * @param metisPartitioning defaults to METIS_LIBRARY
     */
    ParallelTetrahedralMesh(PartitionType metisPartitioning=METIS_LIBRARY);

    /**
     * Destructor.
     */
    virtual ~ParallelTetrahedralMesh();

    /**
     * Construct the mesh using a MeshReader.
     * This method must be overridden in concrete classes.
     *
     * @param rMeshReader the mesh reader
     * @param cullInternalFaces whether to cull internal faces (defaults to false)
     */
    void ConstructFromMeshReader(AbstractTetrahedralMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader,
                                 bool cullInternalFaces=false);

    /**
     * Get the number of nodes that are entirely owned by the local process.
     * (Does not include halo nodes).
     */
    unsigned GetNumLocalNodes() const;
    
    /**
     * Get the number of Elements which are owned by this process (have at least one entirely
     * locally-owned node).
     */
    unsigned GetNumLocalElements() const;

    /**
     * Get the total number of nodes that are actually in use (globally).
     */
    unsigned GetNumNodes() const;

    /**
     * Get the total number of elements that are actually in use (globally).
     */
    unsigned GetNumElements() const;
    
    /**
     * Get the total number of boundary elements that are actually in use (globally).
     */
    unsigned GetNumBoundaryElements() const;

    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     *
     * @param lo is the lowest node number owned by the process
     * @param hi is one higher than the highest node number owned by the process
     * ie. this process owns nodes [lo..hi)
     * and element is "owned" if one or more of its nodes are owned
     */
    void SetElementOwnerships(unsigned lo, unsigned hi);

private:

    /**
     * Add the most recently constructed node to the global->local node mapping
     * 
     * @param index is the global index of node to be registered
     */
    void RegisterNode(unsigned index);

    /**
     * Add the most recently constructed halo node to the global->local halo node mapping
     * 
     * @param index is the global index of halo node to be registered
     */
    void RegisterHaloNode(unsigned index);
    
    /**
     * Add the most recently constructed element to the global->local element mapping
     * 
     * @param index is the global index of element to be registered
     */
    void RegisterElement(unsigned index);

     /**
     * Add the most recently constructed boundary element to the global->local boundary element mapping
     * 
     * @param index is the global index of boundary element to be registered
     */
    void RegisterBoundaryElement(unsigned index);

    /**
     * Overridden solve node mapping method.
     *
     * @param index the global index of the node
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Overridden solve halo node mapping method.
     *
     * @param index the global index of the node
     */
    unsigned SolveHaloNodeMapping(unsigned index);

    /**
     * Overridden solve element mapping method.
     *
     * @param index the global index of the element
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Overridden solve boundary element mapping method.
     *
     * @param index the global index of the boundary element
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /**
     * Compute a parallel partitioning of a given mesh
     * using specialised methods below based on the value
     * of mMetisPartitioning
     * 
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is a set to be filled with the indices of nodes owned by this process
     * @param rHaloNodesOwned is a set to be filled with the indices of halo nodes owned by this process
     * @param rElementsOwned is a set to be filled with the indices of elements owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     * @param rNodePermutation a vector to be filled with the permutation applied to the node numbering by the partitioning method
     * \todo Make it clear which way the permutation applies
     *  
     */
    void ComputeMeshPartitioning(AbstractTetrahedralMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                 std::set<unsigned>& rNodesOwned,
                                 std::set<unsigned>& rHaloNodesOwned,
                                 std::set<unsigned>& rElementsOwned,
                                 std::vector<unsigned>& rProcessorsOffset,
                                 std::vector<unsigned>& rNodePermutation);

    /**
     * Specialised method to compute a parallel partitioning of a given mesh
     * (called by ComputeMeshPartitioning, based on the value of mMetisPartitioning
     * 
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     */
    void DumbNodePartitioning(AbstractTetrahedralMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                              std::set<unsigned>& rNodesOwned);

    /**
     * Specialised method to compute a parallel partitioning of a given mesh
     * (called by ComputeMeshPartitioning, based on the value of mMetisPartitioning
     * 
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     * @param rNodePermutation a vector to be filled with the permutation applied to the node numberig by the partitioning method
     * 
     */
    void MetisBinaryNodePartitioning(AbstractTetrahedralMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                     std::set<unsigned>& rNodesOwned, 
                                     std::vector<unsigned>& rProcessorsOffset,
                                     std::vector<unsigned>& rNodePermutation);

    /**
     * Specialised method to compute a parallel partitioning of a given mesh
     * (called by ComputeMeshPartitioning, based on the value of mMetisPartitioning
     * 
     * @param rMeshReader is the reader pointing to the mesh to be read in and partitioned
     * @param rNodesOwned is an empty set to be filled with the indices of nodes owned by this process
     * @param rProcessorsOffset a vector of length NumProcs to be filled with the index of the lowest indexed node owned by each process
     * @param rNodePermutation a vector to be filled with the permutation applied to the node numberig by the partitioning method
     * 
     */
    void MetisLibraryNodePartitioning(AbstractTetrahedralMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader,
                                      std::set<unsigned>& rNodesOwned, 
                                      std::vector<unsigned>& rProcessorsOffset,
                                      std::vector<unsigned>& rNodePermutation);

    /**
     * Reorder the node indices in this mesh by applying a given permutation
     * 
     * @param rNodePermutation the index permutation to apply
     */
    void ReorderNodes(std::vector<unsigned>& rNodePermutation);
};

#endif /*PARALLELTETRAHEDRALMESH_HPP_*/
