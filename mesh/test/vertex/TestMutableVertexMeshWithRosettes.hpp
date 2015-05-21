/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTMUTABLEVERTEXMESHWITHROSETTES_HPP_
#define TESTMUTABLEVERTEXMESHWITHROSETTES_HPP_

#include <cxxtest/TestSuite.h>

#include "MutableVertexMeshWithRosettes.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "ArchiveOpener.hpp"
#include "Debug.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestMutableVertexMeshWithRosettes : public CxxTest::TestSuite
{
private:

    MutableVertexMeshWithRosettes<2,2>* ConstructFiveCellRosette()
    {
        //Make 11 nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.00000, 0.00000));
        nodes.push_back(new Node<2>(2, true, 0.80902, 0.58779));
        nodes.push_back(new Node<2>(3, true, 0.30902, 0.95106));
        nodes.push_back(new Node<2>(4, true, -0.30902, 0.95106));
        nodes.push_back(new Node<2>(5, true, -0.80902, 0.58779));
        nodes.push_back(new Node<2>(6, true, -1.00000, 0.00000));
        nodes.push_back(new Node<2>(7, true, -0.80902, -0.58779));
        nodes.push_back(new Node<2>(8, true, -0.30902, -0.95106));
        nodes.push_back(new Node<2>(9, true, 0.30902, -0.95106));
        nodes.push_back(new Node<2>(10, true, 0.80902, -0.58779));

        //Make 5 quadrangular cells
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);
        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[5]);
        nodes_elem_3.push_back(nodes[6]);
        nodes_elem_3.push_back(nodes[7]);
        std::vector<Node<2>*> nodes_elem_4;
        nodes_elem_4.push_back(nodes[0]);
        nodes_elem_4.push_back(nodes[7]);
        nodes_elem_4.push_back(nodes[8]);
        nodes_elem_4.push_back(nodes[9]);
        std::vector<Node<2>*> nodes_elem_5;
        nodes_elem_5.push_back(nodes[0]);
        nodes_elem_5.push_back(nodes[9]);
        nodes_elem_5.push_back(nodes[10]);
        nodes_elem_5.push_back(nodes[1]);

        //Make 5 vertex elements
        std::vector<VertexElement<2,2>* > vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_4));
        vertex_elements.push_back(new VertexElement<2,2>(4, nodes_elem_5));

        return new MutableVertexMeshWithRosettes<2,2>(nodes, vertex_elements);
    }

    MutableVertexMeshWithRosettes<2,2>* ConstructProtorosette()
    {
        //Make 9 nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.00000, 0.00000));
        nodes.push_back(new Node<2>(2, true, 0.70711, 0.70711));
        nodes.push_back(new Node<2>(3, true, 0.00000, 1.00000));
        nodes.push_back(new Node<2>(4, true, -0.70711, 0.70711));
        nodes.push_back(new Node<2>(5, true, -1.00000, 0.00000));
        nodes.push_back(new Node<2>(6, true, -0.70711, -0.70711));
        nodes.push_back(new Node<2>(7, true, 0.00000, -1.00000));
        nodes.push_back(new Node<2>(8, true, 0.70711, -0.70711));

        //Make 4 quadrangular cells
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        std::vector<Node<2>*> nodes_elem_2;
        nodes_elem_2.push_back(nodes[0]);
        nodes_elem_2.push_back(nodes[3]);
        nodes_elem_2.push_back(nodes[4]);
        nodes_elem_2.push_back(nodes[5]);
        std::vector<Node<2>*> nodes_elem_3;
        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[5]);
        nodes_elem_3.push_back(nodes[6]);
        nodes_elem_3.push_back(nodes[7]);
        std::vector<Node<2>*> nodes_elem_4;
        nodes_elem_4.push_back(nodes[0]);
        nodes_elem_4.push_back(nodes[7]);
        nodes_elem_4.push_back(nodes[8]);
        nodes_elem_4.push_back(nodes[1]);

        //Make 4 vertex elements
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_4));

        return new MutableVertexMeshWithRosettes<2,2>(nodes, vertex_elements);
    }

public:

    void TestSetAndGetMethods() throw (Exception)
    {
        MutableVertexMeshWithRosettes<2,2>* p_mesh = ConstructFiveCellRosette();

        // Use all three set methods
        p_mesh->SetProtorosetteFormationProbability(0.123);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.234);
        p_mesh->SetRosetteResolutionProbabilityPerTimestep(0.345);

        // Check that we get the correct values
        TS_ASSERT_DELTA(p_mesh->GetProtorosetteFormationProbability(), 0.123, 1e-10);
        TS_ASSERT_DELTA(p_mesh->GetProtorosetteResolutionProbabilityPerTimestep(), 0.234, 1e-10);
        TS_ASSERT_DELTA(p_mesh->GetRosetteResolutionProbabilityPerTimestep(), 0.345, 1e-10);

        // Test that the probability exceptions are working as expected
        TS_ASSERT_THROWS_THIS(p_mesh->SetProtorosetteFormationProbability(-1.234), "Attempting to assign a negative probability.");
        TS_ASSERT_THROWS_THIS(p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(-2.345), "Attempting to assign a negative probability.");
        TS_ASSERT_THROWS_THIS(p_mesh->SetRosetteResolutionProbabilityPerTimestep(-3.456), "Attempting to assign a negative probability.");

        TS_ASSERT_THROWS_THIS(p_mesh->SetProtorosetteFormationProbability(1.0123), "Attempting to assign a probability greater than one.");
        TS_ASSERT_THROWS_THIS(p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(1.0234), "Attempting to assign a probability greater than one.");
        TS_ASSERT_THROWS_THIS(p_mesh->SetRosetteResolutionProbabilityPerTimestep(1.0345), "Attempting to assign a probability greater than one.");
    }

    void TestArchiving() throw (Exception)
    {
        // Set archiving location
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mvmwr.arch";
        ArchiveLocationInfo::SetMeshFilename("mvmwr");

        // Create mesh
        MutableVertexMeshWithRosettes<2,2>* p_mesh = ConstructFiveCellRosette();

        // Set member variables
        p_mesh->SetProtorosetteFormationProbability(0.123);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.234);
        p_mesh->SetRosetteResolutionProbabilityPerTimestep(0.345);

        AbstractMesh<2,2>* const p_abstract_mesh = p_mesh;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects
         * changing during the save, and so object tracking leading to wrong results.
         *
         * E.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.
         */

        // Create an output archive
        {
            TS_ASSERT_EQUALS((static_cast<MutableVertexMesh<2,2>*>(p_abstract_mesh))->GetNumNodes(), 11u);
            TS_ASSERT_EQUALS((static_cast<MutableVertexMesh<2,2>*>(p_abstract_mesh))->GetNumElements(), 5u);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost
            (*p_arch) << p_abstract_mesh;
        }

        {
            // De-serialize and compare
            AbstractMesh<2,2>* p_abstract_mesh_2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_abstract_mesh_2;

            MutableVertexMeshWithRosettes<2,2>* p_mesh_original = static_cast<MutableVertexMeshWithRosettes<2,2>*>(p_abstract_mesh);
            MutableVertexMeshWithRosettes<2,2>* p_mesh_loaded = static_cast<MutableVertexMeshWithRosettes<2,2>*>(p_abstract_mesh_2);

            // Test member variables were archived correctly
            TS_ASSERT_DELTA(p_mesh_original->GetProtorosetteFormationProbability(), 0.123, 1e-10);
            TS_ASSERT_DELTA(p_mesh_loaded->GetProtorosetteFormationProbability(), 0.123, 1e-10);

            TS_ASSERT_DELTA(p_mesh_original->GetProtorosetteResolutionProbabilityPerTimestep(), 0.234, 1e-10);
            TS_ASSERT_DELTA(p_mesh_loaded->GetProtorosetteResolutionProbabilityPerTimestep(), 0.234, 1e-10);

            TS_ASSERT_DELTA(p_mesh_original->GetRosetteResolutionProbabilityPerTimestep(), 0.345, 1e-10);
            TS_ASSERT_DELTA(p_mesh_loaded->GetRosetteResolutionProbabilityPerTimestep(), 0.345, 1e-10);

            // Compare the loaded mesh against the original
            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), 11u);
            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), p_mesh_loaded->GetNumNodes());

            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), 5u);
            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), p_mesh_loaded->GetNumElements());

            for (unsigned node_idx = 0 ; node_idx < p_mesh_original->GetNumNodes() ; node_idx++)
            {
                Node<2>* p_node = p_mesh_original->GetNode(node_idx);
                Node<2>* p_node2 = p_mesh_loaded->GetNode(node_idx);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());
                TS_ASSERT_DELTA(p_node->rGetLocation()[0], p_node2->rGetLocation()[0], 1e-10);
                TS_ASSERT_DELTA(p_node->rGetLocation()[1], p_node2->rGetLocation()[1], 1e-10);
            }

            for (unsigned elem_idx = 0 ; elem_idx < p_mesh_original->GetNumElements() ; elem_idx++)
            {
                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_idx)->GetNumNodes(),
                                 p_mesh_loaded->GetElement(elem_idx)->GetNumNodes());

                for (unsigned local_index=0; local_index<p_mesh_original->GetElement(elem_idx)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_idx)->GetNodeGlobalIndex(local_index),
                                     p_mesh_loaded->GetElement(elem_idx)->GetNodeGlobalIndex(local_index));
                }
            }

            // Tidy up
            delete p_mesh_original;
            delete p_mesh_loaded;
        }
    }

    void TestHandleHighOrderJunctions() throw (Exception)
    {
        // Make an empty mesh
        MutableVertexMeshWithRosettes<2,2> mesh;

        // Add 6 nodes, three each for two triangular elements
        mesh.AddNode(new Node<2>(0, true, 0.0, 0.0));
        mesh.AddNode(new Node<2>(1, true, 1.0, 0.0));
        mesh.AddNode(new Node<2>(2, true, 1.0, 1.0));
        mesh.AddNode(new Node<2>(3, true, 2.0, 0.0));
        mesh.AddNode(new Node<2>(4, true, 3.0, 0.0));
        mesh.AddNode(new Node<2>(5, true, 2.0, 1.0));

        // Create two vectors of nodes, one representing each of two triangular elements
        std::vector<Node<2>* > nodes_elem_1;
        nodes_elem_1.push_back(mesh.GetNode(0));
        nodes_elem_1.push_back(mesh.GetNode(1));
        nodes_elem_1.push_back(mesh.GetNode(2));

        std::vector<Node<2>* > nodes_elem_2;
        nodes_elem_2.push_back(mesh.GetNode(3));
        nodes_elem_2.push_back(mesh.GetNode(4));
        nodes_elem_2.push_back(mesh.GetNode(5));

        // Add three copies of each element to the mesh
        mesh.AddElement(new VertexElement<2,2>(0, nodes_elem_1));
        mesh.AddElement(new VertexElement<2,2>(1, nodes_elem_1));
        mesh.AddElement(new VertexElement<2,2>(2, nodes_elem_1));

        mesh.AddElement(new VertexElement<2,2>(3, nodes_elem_2));
        mesh.AddElement(new VertexElement<2,2>(4, nodes_elem_2));
        mesh.AddElement(new VertexElement<2,2>(5, nodes_elem_2));

        // Get pointers to our two test-nodes
        Node<2>* p_node_0 = mesh.GetNode(0);
        Node<2>* p_node_5 = mesh.GetNode(5);

        /**
         * At this point, nodes 0 and 5, say, are in three elements and HandleHighOrderJunctions() should
         * throw an exception
         */
        TS_ASSERT_THROWS_THIS(mesh.HandleHighOrderJunctions(p_node_0, p_node_5), "Neither node is high order");

        /**
         * If we now make a further copy of each element containing, both node 0 and node 5 will be contained in four
         * elements, and we expect an exception to be thrown
         */
        mesh.AddElement(new VertexElement<2,2>(6, nodes_elem_1));
        mesh.AddElement(new VertexElement<2,2>(7, nodes_elem_2));
        TS_ASSERT_THROWS_THIS(mesh.HandleHighOrderJunctions(p_node_0, p_node_5), "Both nodes involved in a swap event are contained in more than three elements");

        /**
         * If we now delete element 7, node 0 will still be contained in four elements while node 5 will only be
         * contained in three.  This is expected behaviour and so no exception should be thrown.
         */
        mesh.DeleteElementPriorToReMesh(7);
        TS_ASSERT_THROWS_NOTHING(mesh.HandleHighOrderJunctions(p_node_0, p_node_5));
    }

    void TestHandleAdditionalRemodellingBehaviour() throw (Exception)
    {
        // Make an empty mesh
        MutableVertexMeshWithRosettes<2,2> mesh;

        // Add 2 nodes
        mesh.AddNode(new Node<2>(0, true, 0.0, 0.0));
        mesh.AddNode(new Node<2>(1, true, 1.0, 0.0));

        // Get pointers to the two nodes
        Node<2>* p_node_a = mesh.GetNode(0);
        Node<2>* p_node_b = mesh.GetNode(1);

        // Create an empty set
        std::set<unsigned> empty_set;

        // Check that, for all case numbers not equal to four, we trigger the exception
        for( unsigned case_num = 0 ; case_num < 10 ; case_num++ )
        {
            if(case_num == 4)
            {
                continue;
            }

            TS_ASSERT_THROWS_THIS(mesh.HandleAdditionalRemodellingBehaviour(p_node_a, p_node_b, empty_set, case_num),
                                  "No functionality for this case yet");
        }
    }

    void TestPerformRosetteRankIncrease() throw (Exception)
    {
        // Create the standard five-cell rosette
        MutableVertexMeshWithRosettes<2,2>* p_mesh = ConstructFiveCellRosette();

        /**
         * Modify the mesh to incorporate an additional element which will go on to increase the rosette rank
         */

        // One new node will be needed
        unsigned new_node_idx = p_mesh->AddNode(new Node<2>(11, false, 0.5, 0.0));

        // One new element will be needed, consisting of four nodes
        std::vector<Node<2>* > nodes_new_elem;
        nodes_new_elem.push_back(p_mesh->GetNode(new_node_idx));
        nodes_new_elem.push_back(p_mesh->GetNode(10));
        nodes_new_elem.push_back(p_mesh->GetNode(1));
        nodes_new_elem.push_back(p_mesh->GetNode(2));

        unsigned new_elem_idx = p_mesh->AddElement(new VertexElement<2,2>(5, nodes_new_elem));
        p_mesh->GetElement(new_elem_idx)->RegisterWithNodes();

        // Add new node in to elements 0 and 4, and remove node with global index 1 from elements 0 and 4
        VertexElement<2,2>* p_elem_0 = p_mesh->GetElement(0);
        VertexElement<2,2>* p_elem_4 = p_mesh->GetElement(4);
        VertexElement<2,2>* p_elem_n = p_mesh->GetElement(new_elem_idx);

        p_elem_0->AddNode(p_mesh->GetNode(new_node_idx), p_elem_0->GetNodeLocalIndex(0));
        p_elem_4->AddNode(p_mesh->GetNode(new_node_idx), p_elem_4->GetNodeLocalIndex(1));

        p_elem_0->DeleteNode(p_elem_0->GetNodeLocalIndex(1));
        p_elem_4->DeleteNode(p_elem_4->GetNodeLocalIndex(1));

        /**
         * Now the mesh is as desired, we can get the necessary numbers, perform the rosette rank increase, and assert
         * that everything is as intended after the operation
         */

        unsigned num_nodes_before = p_mesh->GetNumNodes();
        assert(num_nodes_before == 12);

        unsigned num_nodes_elem_0_before = p_elem_0->GetNumNodes();
        assert(num_nodes_elem_0_before == 4);

        unsigned num_nodes_elem_4_before = p_elem_4->GetNumNodes();
        assert(num_nodes_elem_4_before == 4);

        unsigned num_nodes_new_elem_before = p_elem_n->GetNumNodes();
        assert(num_nodes_new_elem_before == 4);

        c_vector<double, 2> node_0_location_before = p_mesh->GetNode(0)->rGetLocation();

        // Perform the rosette rank increase
        p_mesh->PerformRosetteRankIncrease(p_mesh->GetNode(0), p_mesh->GetNode(new_node_idx));

        // The mesh should have lost one node
        TS_ASSERT_EQUALS(num_nodes_before - 1, p_mesh->GetNumNodes());

        // Elements 0 and 4 should have lost a node, but the new element should not have lost any
        TS_ASSERT_EQUALS(num_nodes_elem_0_before - 1, p_elem_0->GetNumNodes());
        TS_ASSERT_EQUALS(num_nodes_elem_4_before - 1, p_elem_4->GetNumNodes());
        TS_ASSERT_EQUALS(num_nodes_new_elem_before, p_elem_n->GetNumNodes());

        // The node with global index 0 should remain in the same location
        TS_ASSERT_DELTA(node_0_location_before[0], p_mesh->GetNode(0)->rGetLocation()[0], 1e-10);
        TS_ASSERT_DELTA(node_0_location_before[1], p_mesh->GetNode(0)->rGetLocation()[1], 1e-10);

        // The node with global index 0 should now be included in the new element
        TS_ASSERT_LESS_THAN(p_elem_n->GetNodeLocalIndex(0), UINT_MAX);
    }

    void TestPerformProtorosetteResolution() throw (Exception)
    {
        // Let us first create a protorosette
        MutableVertexMeshWithRosettes<2,2>* p_mesh = ConstructProtorosette();

        VertexElement<2,2>* p_elem_0 = p_mesh->GetElement(0);
        VertexElement<2,2>* p_elem_1 = p_mesh->GetElement(1);
        VertexElement<2,2>* p_elem_2 = p_mesh->GetElement(2);
        VertexElement<2,2>* p_elem_3 = p_mesh->GetElement(3);

        Node<2>* p_node_0 = p_mesh->GetNode(0);

        // Perform the protorosette resolution
        p_mesh->PerformProtorosetteResolution(p_node_0);

        /**
         * We now need to check that this has been done correctly
         */

        // The number of nodes in the mesh should have increase by one (from 9 to 10)
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 10);
        Node<2>* p_node_9 = p_mesh->GetNode(9);

        // Two elements (either 0 & 2, or 1 & 3) should now have five nodes, while the other two should still have four
        TS_ASSERT( ( p_elem_0->GetNumNodes() == 4 ) || ( p_elem_0->GetNumNodes() == 5 ) );

        // The central node is still one of the two in the centre which are now joined by an edge.  Due to the initial
        // symmetry, this edge will have length CellRearrangementThreshold * CellRearrangementRatio
        double node_spacing = p_mesh->GetCellRearrangementThreshold() * p_mesh->GetCellRearrangementRatio();
        c_vector<double, 2> node_0_pos = p_node_0->rGetLocation();
        c_vector<double, 2> node_9_pos = p_node_9->rGetLocation();

        TS_ASSERT_DELTA(norm_2(node_0_pos), 0.5 * node_spacing, 1e-10);
        TS_ASSERT_DELTA(norm_2(node_9_pos), 0.5 * node_spacing, 1e-10);
        TS_ASSERT_DELTA(norm_2(node_9_pos - node_0_pos), node_spacing, 1e-10);

        // Because the division axis is random, we are left with two possible cases
        if (p_elem_0->GetNumNodes() == 4)
        {
            TS_ASSERT_EQUALS(p_elem_0->GetNumNodes(), 4);
            TS_ASSERT_EQUALS(p_elem_2->GetNumNodes(), 4);
            TS_ASSERT_EQUALS(p_elem_1->GetNumNodes(), 5);
            TS_ASSERT_EQUALS(p_elem_3->GetNumNodes(), 5);

            if (node_0_pos[0] > 0)
            {
                TS_ASSERT_DELTA(node_0_pos[0], node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_0_pos[1], node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_9_pos[0], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_9_pos[1], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
            }
            else // (node_0_pos[0] < 0)
            {
                TS_ASSERT_DELTA(node_9_pos[0], node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_9_pos[1], node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_0_pos[0], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_0_pos[1], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
            }
        }
        else //p_elem_0->GetNumNodes() == 5
        {
            TS_ASSERT_EQUALS(p_elem_0->GetNumNodes(), 5);
            TS_ASSERT_EQUALS(p_elem_2->GetNumNodes(), 5);
            TS_ASSERT_EQUALS(p_elem_1->GetNumNodes(), 4);
            TS_ASSERT_EQUALS(p_elem_3->GetNumNodes(), 4);

            if (node_0_pos[0] > 0)
            {
                TS_ASSERT_DELTA(node_0_pos[0], node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_0_pos[1], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_9_pos[0], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_9_pos[1], node_spacing * sqrt(2.0) * 0.25, 1e-10);
            }
            else // (node_0_pos[0] < 0)
            {
                TS_ASSERT_DELTA(node_9_pos[0], node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_9_pos[1], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_0_pos[0], -1.0 * node_spacing * sqrt(2.0) * 0.25, 1e-10);
                TS_ASSERT_DELTA(node_0_pos[1], node_spacing * sqrt(2.0) * 0.25, 1e-10);
            }
        }
    }

    void TestPerformRosetteRankDecrease() throw (Exception)
    {
        // Let us first create a protorosette
        MutableVertexMeshWithRosettes<2,2>* p_mesh = ConstructFiveCellRosette();

        VertexElement<2,2>* p_elem_0 = p_mesh->GetElement(0);
        VertexElement<2,2>* p_elem_1 = p_mesh->GetElement(1);
        VertexElement<2,2>* p_elem_2 = p_mesh->GetElement(2);
        VertexElement<2,2>* p_elem_3 = p_mesh->GetElement(3);
        VertexElement<2,2>* p_elem_4 = p_mesh->GetElement(4);

        Node<2>* p_node_0 = p_mesh->GetNode(0);
        Node<2>* p_node_1 = p_mesh->GetNode(1);
        Node<2>* p_node_2 = p_mesh->GetNode(2);
        Node<2>* p_node_3 = p_mesh->GetNode(3);
        Node<2>* p_node_4 = p_mesh->GetNode(4);
        Node<2>* p_node_5 = p_mesh->GetNode(5);
        Node<2>* p_node_6 = p_mesh->GetNode(6);
        Node<2>* p_node_7 = p_mesh->GetNode(7);
        Node<2>* p_node_8 = p_mesh->GetNode(8);
        Node<2>* p_node_9 = p_mesh->GetNode(9);
        Node<2>* p_node_10 = p_mesh->GetNode(10);

        // Perform the protorosette resolution
        p_mesh->PerformRosetteRankDecrease(p_node_0);

        /**
         * We now need to check that this has been done correctly
         */

        // The number of nodes in the mesh should have increase by one (from 11 to 12)
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 12);
        Node<2>* p_node_11 = p_mesh->GetNode(11);

        // Two elements should now have five nodes, and the other three should have four
        std::vector<VertexElement<2,2>* > five_node_elems;
        std::vector<VertexElement<2,2>* > four_node_elems;

        for (unsigned elem_idx = 0 ; elem_idx < 5 ; elem_idx++)
        {
            VertexElement<2,2>* p_current_element = p_mesh->GetElement(elem_idx);

            if (p_current_element->GetNumNodes() == 4)
            {
                four_node_elems.push_back(p_current_element);
            }
            else if (p_current_element->GetNumNodes() == 5)
            {
                five_node_elems.push_back(p_current_element);
            }
            else
            {
                TS_FAIL("No element should contain this number of nodes");
            }
        }

        TS_ASSERT_EQUALS(five_node_elems.size(), 2);
        TS_ASSERT_EQUALS(four_node_elems.size(), 3);

        // The five-node elements should be separated by one (0,2 or 1,3 or 2,4 or 3,0)
        TS_ASSERT( ( (five_node_elems[0]->GetIndex() + 2) % 5 == five_node_elems[1]->GetIndex() ) ||
                   ( (five_node_elems[1]->GetIndex() + 2) % 5 == five_node_elems[0]->GetIndex() ) );

        // We can also check the location of the new node, which should be a distance from the origin given by
        // CellRearrangementThreshold * CellRearrangementRatio.  Node 0 should not have moved
        double node_spacing = p_mesh->GetCellRearrangementThreshold() * p_mesh->GetCellRearrangementRatio();
        c_vector<double, 2> node_0_pos = p_node_0->rGetLocation();
        c_vector<double, 2> node_11_pos = p_node_11->rGetLocation();

        TS_ASSERT_DELTA(norm_2(node_0_pos), 0.0, 1e-10);
        TS_ASSERT_DELTA(norm_2(node_11_pos), node_spacing, 1e-10);
    }

    void TestCheckForRosettes() throw (Exception)
    {
        // Let us first create reference meshes for four and five cell rosettes
        MutableVertexMeshWithRosettes<2,2>* p_ref_rosette = ConstructFiveCellRosette();
        MutableVertexMeshWithRosettes<2,2>* p_ref_protorosette = ConstructProtorosette();

        // We also need to create meshes which we will CheckForRosettes
        MutableVertexMeshWithRosettes<2,2>* p_rosette = ConstructFiveCellRosette();
        MutableVertexMeshWithRosettes<2,2>* p_protorosette = ConstructProtorosette();

        /**
         * Initially, set the parameters such that nothing should happen when we check for rosettes
         */
        p_rosette->SetRosetteResolutionProbabilityPerTimestep(0.0);
        p_protorosette->SetProtorosetteResolutionProbabilityPerTimestep(0.0);

        // Check for rosettes, and then verify that nothing has yet changed
        p_rosette->CheckForRosettes();
        p_protorosette->CheckForRosettes();

        TS_ASSERT_EQUALS(p_ref_rosette->GetNumNodes(), p_rosette->GetNumNodes());
        for( unsigned node_idx = 0 ; node_idx < p_ref_rosette->GetNumNodes() ; node_idx++ )
        {
            Node<2>* current_ref_node = p_ref_rosette->GetNode(node_idx);
            Node<2>* current_node = p_rosette->GetNode(node_idx);

            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[0], current_node->rGetLocation()[0], 1e-10);
            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[1], current_node->rGetLocation()[1], 1e-10);
        }

        TS_ASSERT_EQUALS(p_ref_protorosette->GetNumNodes(), p_protorosette->GetNumNodes());
        for( unsigned node_idx = 0 ; node_idx < p_ref_protorosette->GetNumNodes() ; node_idx++ )
        {
            Node<2>* current_ref_node = p_ref_protorosette->GetNode(node_idx);
            Node<2>* current_node = p_protorosette->GetNode(node_idx);

            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[0], current_node->rGetLocation()[0], 1e-10);
            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[1], current_node->rGetLocation()[1], 1e-10);
        }

        /**
         * Next, set the parameters such that resolution and rank-decrease will happen when we check for rosettes
         */
        p_rosette->SetRosetteResolutionProbabilityPerTimestep(1.0);
        p_protorosette->SetProtorosetteResolutionProbabilityPerTimestep(1.0);

        // Check for rosettes, and then verify that the number of nodes in each mesh has increased by one
        p_rosette->CheckForRosettes();
        p_protorosette->CheckForRosettes();

        // In the rosette mesh, all nodes should remain in the same place, except the new one
        TS_ASSERT_EQUALS(p_ref_rosette->GetNumNodes() + 1, p_rosette->GetNumNodes());
        for( unsigned node_idx = 0 ; node_idx < p_ref_rosette->GetNumNodes() ; node_idx++ )
        {
            Node<2>* current_ref_node = p_ref_rosette->GetNode(node_idx);
            Node<2>* current_node = p_rosette->GetNode(node_idx);

            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[0], current_node->rGetLocation()[0], 1e-10);
            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[1], current_node->rGetLocation()[1], 1e-10);
        }

        // In the protorosette mesh, all nodes should remain in the same place, except node 0 and the new one
        TS_ASSERT_EQUALS(p_ref_protorosette->GetNumNodes() + 1, p_protorosette->GetNumNodes());
        for( unsigned node_idx = 1 ; node_idx < p_ref_protorosette->GetNumNodes() ; node_idx++ )
        {
            Node<2>* current_ref_node = p_ref_protorosette->GetNode(node_idx);
            Node<2>* current_node = p_protorosette->GetNode(node_idx);

            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[0], current_node->rGetLocation()[0], 1e-10);
            TS_ASSERT_DELTA(current_ref_node->rGetLocation()[1], current_node->rGetLocation()[1], 1e-10);
        }
    }
};

#endif /*TESTMUTABLEVERTEXMESHWITHROSETTES_HPP_*/