/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef TESTNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_
#define TESTNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "ArchiveOpener.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TrianglesMeshReader.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestNodeBasedCellPopulationParallelMethods : public AbstractCellBasedTestSuite
{
private:

    NodesOnlyMesh<3>* mpNodesOnlyMesh;

    NodeBasedCellPopulation<3>* mpNodeBasedCellPopulation;

    void setUp()
    {
        AbstractCellBasedTestSuite::setUp();

        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0 + (double)PetscTools::GetMyRank()));

        mpNodesOnlyMesh = new NodesOnlyMesh<3>;
        mpNodesOnlyMesh->ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mpNodesOnlyMesh->GetNumNodes());

        mpNodeBasedCellPopulation = new NodeBasedCellPopulation<3>(*mpNodesOnlyMesh, cells);

        delete nodes[0];
    }

    void tearDown()
    {
        delete mpNodeBasedCellPopulation;
        delete mpNodesOnlyMesh;

        AbstractCellBasedTestSuite::tearDown();
    }
public:

    void TestGetCellAndNodePair() throw (Exception)
    {
        std::pair<CellPtr, Node<3>* > pair = mpNodeBasedCellPopulation->GetCellNodePair(0);

        TS_ASSERT_EQUALS(pair.second->GetIndex(), 0u);

        CellPtr p_returned_cell = pair.first;
        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->GetLocationIndexUsingCell(p_returned_cell), 0u);
    }

    void TestAddNodeAndCellsToSend() throw (Exception)
    {
        unsigned index_of_node_to_send = 0;
        mpNodeBasedCellPopulation->AddNodeAndCellToSendRight(index_of_node_to_send);
        mpNodeBasedCellPopulation->AddNodeAndCellToSendLeft(index_of_node_to_send);

        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellsToSendRight.size(), 1u);
        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellsToSendLeft.size(), 1u);

        unsigned node_right_index = (*mpNodeBasedCellPopulation->mCellsToSendRight.begin()).second->GetIndex();
        TS_ASSERT_EQUALS(node_right_index, 0u);

        unsigned node_left_index = (*mpNodeBasedCellPopulation->mCellsToSendLeft.begin()).second->GetIndex();
        TS_ASSERT_EQUALS(node_left_index, 0u);

    }
///\todo #1902 fix archive include headers for boost 1-34
//    void TestSendAndRecieveCells() throw (Exception)
//    {
//        unsigned index_of_node_to_send = 0;
//        mpNodeBasedCellPopulation->AddNodeAndCellToSendRight(index_of_node_to_send);
//        mpNodeBasedCellPopulation->AddNodeAndCellToSendLeft(index_of_node_to_send);
//
//        TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mCellCommunicationTag, 123u);
//
//        TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvRight));
//        TS_ASSERT(!(mpNodeBasedCellPopulation->mpCellsRecvLeft));
//
//        mpNodeBasedCellPopulation->SendCellsToNeighbourProcesses();
//
//        if (!PetscTools::AmTopMost())
//        {
//            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvRight->size(), 1u);
//
//            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvRight->begin()).second->GetIndex();
//            TS_ASSERT_EQUALS(index, 0u);
//        }
//        if (!PetscTools::AmMaster())
//        {
//            TS_ASSERT_EQUALS(mpNodeBasedCellPopulation->mpCellsRecvLeft->size(), 1u);
//
//            unsigned index = (*mpNodeBasedCellPopulation->mpCellsRecvLeft->begin()).second->GetIndex();
//            TS_ASSERT_EQUALS(index, 0u);
//        }
//    }
};

#endif /*TESTNODEBASEDCELLPOPULATIONPARALLELMETHODS_HPP_*/