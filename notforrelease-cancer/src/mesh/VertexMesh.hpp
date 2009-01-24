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
#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

#include <iostream>
#include <map>
#include <algorithm>

#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>

#include "AbstractMesh.hpp"
#include "VertexMeshReader2d.hpp"
#include "VertexElement.hpp"
#include "NodeMap.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based tissue simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh
{
    friend class TestVertexMesh;
    
private:

    /** Vector of pointers to nodes. */
    std::vector<Node<SPACE_DIM>*> mNodes;

    /** Vector of pointers to VertexElements. */    
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;
   
    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangment. */
    double mThresholdDistance;
    
    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;
    
     /** Indices of elements that have been deleted. These indices can be reused when adding new elements. */
    std::vector<unsigned> mDeletedElementIndices;
    
    /** Whether nodes have been added to the mesh. */
    bool mAddedNodes;
    
    /** Whether Elements have been added to the mesh. */
    bool mAddedElements;

    /** Create correspondences between VertexElements and Nodes in the mesh. */
    void SetupVertexElementsOwnedByNodes();
    
    /**
     * Helper method for ReMesh to Identify the type of swap
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     * 
     * @param rNodeA one of the nodes to perform the swap with 
     * @param rNodeB the other node to perform the swap
     */  
    void IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);
    
    /**
     * Helper method for ReMesh to merge nodes when needed 
     * Move node with smallest global index to center and remove other node
     * 
     * @param rNodeA one of the nodes to perform the swap with 
     * @param rNodeB the other node to perform the swap
     */  
    void PerformNodeMerge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB,
                          std::set<unsigned> ElementsContainingNodes);
    
    /**
     * Helper method for ReMesh to perform the T1 Swap
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     * 
     * @param rNodeA one of the nodes to perform the swap with 
     * @param rNodeB the other node to perform the swap
     */  
    void PerformT1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, 
                       std::set<unsigned> ElementsContainingNodes);
    
    /**
     * Method to divide an element in half 
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     * 
     * @param pElement the element to divide
     */  
    void DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement);
    
    /**
     * Method to divide an element given 2 nodes in which to divide the element with 
     * 
     * \todo This method currently assumes SPACE_DIM = 2 (see #866)
     * 
     * @param rElement the element to divide
     * @param nodeAIndex the local index of node where to divide
     * @param nodeBindex the local index of node where to divide
     * 
     * @return the index of the new element
     */  
    unsigned DivideElement(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned nodeAIndex, unsigned nodeBIndex);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mThresholdDistance;      
    }
    
public:

    /**
     * Default constructor.
     * 
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param thresholdDistance the minimum threshold distance for element rearrangment (defaults to 0.01)
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes, 
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
               double thresholdDistance=0.01);

    /**
     * Constructor for use by serializer.
     *
     * @param thresholdDistance the minimum threshold distance for element rearrangment (defaults to 0.01)
     */
    VertexMesh(double thresholdDistance=0.01);
    
    /**
     * Destructor.
     */
    ~VertexMesh();
    
    /**
     * @return mThresholdDistance
     */
    double GetThresholdDistance() const;
    
    /**
     * Set method for mThresholdDistance.
     * 
     * @param thresholdDistance
     */
    void SetThresholdDistance(double thresholdDistance);

    /**
     * Helper constructor, creates a rectangular vertex-based mesh.
     * 
     * @param numAcross number of VertexElements across
     * @param numUp number of VertexElements up
     * @param thresholdDistance the minimum threshold distance for element rearrangment (defaults to 0.01)
     */
    VertexMesh(unsigned numAcross, unsigned numUp, double thresholdDistance=0.01);

    /**
     * @return the number of Nodes in the mesh.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    unsigned GetNumElements() const;

    /**
     * @param index the global index of a specified node
     * 
     * @return a pointer to the node
     */
    Node<SPACE_DIM>* GetNode(unsigned index) const;
    
    /**
     * @param index  the global index of a specified vertex element
     * 
     * @return a pointer to the vertex element
     */    
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Add a node to the mesh.
     *
     * Note: After calling this one or more times, you must then call ReMesh.
     * 
     * @param pNewNode pointer to the new node
     * @return the global index of the new node in the mesh.
     */
    unsigned AddNode(Node<SPACE_DIM> *pNewNode);


    /**
     * Add an element to the mesh.
     *
     */
    unsigned AddElement(VertexElement<ELEMENT_DIM, SPACE_DIM> *pNewElement);


    /**
     *  Move the node with a particular index to a new point in space.
     * 
      * @param nodeIndex the index of the node to be moved
      * @param point the new target location of the node
      */
    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);

    /**
     * Delete mNodes and mElements.
     */
    void Clear();

    /** 
     * Add a node on the edge between two nodes.
     * 
     * @param pNodeA a pointer to one node
     * @param pNodeB a pointer to the other nodes
     */
    void DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);
        
    /**
     * Re-mesh the mesh.
     * 
     * @param elementMap a NodeMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created 
     *                   with the correct size, GetNumElements()
     */
    void ReMesh(NodeMap& elementMap);

    /**
     * Alternative version of remesh which takes no parameters does not require a NodeMap. 
     * Note: inherited classes should overload ReMesh(NodeMap&).
     */
    void ReMesh();

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     * 
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Given a node and one of its containing elements, find a set containing 
     * the indices of those neighbouring node(s) that are NOT also in the element.
     * 
     * Note that we allow for more than one such index, since there is no reason 
     * a priori to assume that each node is contained by exactly three elements.
     * 
     * @param nodeIndex global index of the node
     * @param elemIndex global index of the element
     * 
     * @return its neighbouring nodes that are not in the element
     */
    std::set<unsigned> GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex);

    /** 
     * Construct the mesh using a mesh reader.
     * 
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(VertexMeshReader2d& rMeshReader);

};

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VertexMesh
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const VertexMesh<ELEMENT_DIM, SPACE_DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const double threshold_distance = t->GetThresholdDistance();
    ar << threshold_distance;
}

/**
 * De-serialize constructor parameters and initialise VertexMesh.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, VertexMesh<ELEMENT_DIM, SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double threshold_distance;
    ar >> threshold_distance;

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexMesh<ELEMENT_DIM, SPACE_DIM>(threshold_distance);
}
}
} // namespace ...

#endif /*VERTEXMESH_HPP_*/
