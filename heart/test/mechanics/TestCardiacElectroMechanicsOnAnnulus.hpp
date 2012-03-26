/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTCARDIACELECTROMECHANICSONANNULUS_HPP_
#define TESTCARDIACELECTROMECHANICSONANNULUS_HPP_


#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "NumericFileComparison.hpp"
#include "Hdf5DataReader.hpp"

// small region stimulus, essentially irrelevant
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-5e5, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];
        if (fabs(x)<0.02+1e-6 && y<-0.4+1e-6) // stimulating small region
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

};

class TestCardiacElectroMechanicsOnAnnulus : public CxxTest::TestSuite
{
public:
    void TestOnAnnulus() throw (Exception)
    {
        TetrahedralMesh<2,2> electrics_mesh;
        QuadraticMesh<2> mechanics_mesh;

        // could (should?) use finer electrics mesh, but keeping electrics simulation time down
        TrianglesMeshReader<2,2> reader1("mesh/test/data/annuli/circular_annulus_960_elements");
        electrics_mesh.ConstructFromMeshReader(reader1);

        TrianglesMeshReader<2,2> reader2("mesh/test/data/annuli/circular_annulus_960_elements_quad",2 /*quadratic elements*/);
        mechanics_mesh.ConstructFromMeshReader(reader2);

        PointStimulus2dCellFactory cell_factory;

        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double x = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double y = mechanics_mesh.GetNode(i)->rGetLocation()[1];

            if (fabs(x)<1e-6 && fabs(y+0.5)<1e-6)  // fixed point (0.0,-0.5)
            {
                fixed_nodes.push_back(i);
            }
        }

//// ///\todo #2027 investigate SNES by introducing longer end time
       HeartConfig::Instance()->SetSimulationDuration(20.0);

       ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);

       problem_defn.SetContractionModel(KERCHOFFS2003,0.1);
       problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
       problem_defn.SetZeroDisplacementNodes(fixed_nodes);
       problem_defn.SetMechanicsSolveTimestep(1.0);

       problem_defn.SetVariableFibreSheetDirectionsFile("heart/test/data/fibre_tests/circular_annulus_960_elements.ortho", false);

//// ///\todo #2027 introduce internal pressure
//       std::vector<BoundaryElement<1,2>*> boundary_elems;
//       for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
//               = mechanics_mesh.GetBoundaryElementIteratorBegin();
//             iter != mechanics_mesh.GetBoundaryElementIteratorEnd();
//             ++iter)
//       {
//            ChastePoint<2> centroid = (*iter)->CalculateCentroid();
//            double r = sqrt( centroid[0]*centroid[0] + centroid[1]*centroid[1] );
//
//            if (r < 0.4)
//            {
//                BoundaryElement<1,2>* p_element = *iter;
//                boundary_elems.push_back(p_element);
//            }
//       }
//       problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, -1);

       CardiacElectroMechanicsProblem<2> problem(COMPRESSIBLE,
                                                 &electrics_mesh,
                                                 &mechanics_mesh,
                                                 &cell_factory,
                                                 &problem_defn,
                                                 "TestEmOnAnnulus");

       problem.Solve();

       // we don't test anything, we mainly just want to verify it solves OK


       MechanicsEventHandler::Headings();
       MechanicsEventHandler::Report();
    }
};

#endif // TESTCARDIACELECTROMECHANICSONANNULUS_HPP_

/*

Exact fibre directions HACKED IN!

chaste, no LU
       Assemble            Solve           Update          AllMech          NonMech           Output            Total
   3.510 (  3%)   112.848 ( 91%)     1.738 (  1%)   118.772 ( 96%)     5.172 (  4%)     0.156 (  0%)   124.228 (100%)  (seconds)

chaste, LU
       Assemble            Solve           Update          AllMech          NonMech           Output            Total
   3.531 (  3%)   101.908 ( 90%)     1.771 (  2%)   107.875 ( 95%)     5.474 (  5%)     0.206 (  0%)   113.691 (100%)  (seconds)

snes, LU
       Assemble            Solve           Update          AllMech          NonMech           Output            Total
   0.000 (  0%)     0.000 (  0%)     0.000 (  0%)     5.194 ( 48%)     5.316 ( 49%)     0.171 (  2%)    10.896 (100%)

 */
