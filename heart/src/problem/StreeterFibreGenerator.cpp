/*

Copyright (C) University of Oxford, 2005-2011

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

#include "UblasCustomFunctions.hpp"

#include "StreeterFibreGenerator.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
//#include "HeartRegionCodes.hpp"

template<unsigned SPACE_DIM>
double StreeterFibreGenerator<SPACE_DIM>::GetAveragedThicknessLocalNode(
        const unsigned nodeIndex, const std::vector<double>& wallThickness) const
{
    if (nodeIndex < mrMesh.GetDistributedVectorFactory()->GetLow() ||
        nodeIndex >= mrMesh.GetDistributedVectorFactory()->GetHigh() )
    {
        return 0.0;  //Don't calculate this for nodes which aren't local
    }

    // Initialise the average with the value corresponding to the current node
    double average = wallThickness[nodeIndex];
    unsigned nodes_visited = 1;

    // Use a set to store visited nodes
    std::set<unsigned> visited_nodes;
    visited_nodes.insert(nodeIndex);

    Node<SPACE_DIM>* p_current_node = mrMesh.GetNode(nodeIndex);

    // Loop over the elements containing the given node
    for(typename Node<SPACE_DIM>::ContainingElementIterator element_iterator = p_current_node->ContainingElementsBegin();
        element_iterator != p_current_node->ContainingElementsEnd();
        ++element_iterator)
    {
        // Get a pointer to the container element
        Element<SPACE_DIM,SPACE_DIM>* p_containing_element = mrMesh.GetElement(*element_iterator);

       // Loop over the nodes of the element
       for(unsigned node_local_index=0;
           node_local_index<p_containing_element->GetNumNodes();
           node_local_index++)
       {
            Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
            unsigned neighbour_node_index = p_neighbour_node->GetIndex();

            // Check if the neighbour node has already been visited
            if (visited_nodes.find(neighbour_node_index) == visited_nodes.end())
            {
                average += wallThickness[neighbour_node_index];
                visited_nodes.insert(neighbour_node_index);
                nodes_visited++;
            }
       }
    }

    return average/nodes_visited;
}




template<unsigned SPACE_DIM>
double StreeterFibreGenerator<SPACE_DIM>::GetFibreMaxAngle(
        const c_vector<HeartRegionType, SPACE_DIM+1>& nodesRegionsForElement) const
{
    unsigned lv=0, rv=0;

    for (unsigned index=0; index<SPACE_DIM+1; index++)
    {
        switch (nodesRegionsForElement[index])
        {
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_WALL:
            case HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM:
                lv++;
                break;

            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_WALL:
            case HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM:
                rv++;
                break;

            case HeartGeometryInformation<SPACE_DIM>::UNKNOWN:
            default:
                NEVER_REACHED;
        }
    }

    // If most of the nodes are in the right ventricle
    if (rv>lv)
    {
        return M_PI/4.0;
    }

    // Anywhere else
    return M_PI/3.0;
}

template<unsigned SPACE_DIM>
StreeterFibreGenerator<SPACE_DIM>::StreeterFibreGenerator(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh)
    : mrMesh(rMesh),
      mpGeometryInfo(NULL),
      mApexToBase(zero_vector<double>(SPACE_DIM))
{
}

template<unsigned SPACE_DIM>
StreeterFibreGenerator<SPACE_DIM>::~StreeterFibreGenerator()
{
    delete mpGeometryInfo;
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::SetSurfaceFiles(
            const std::string &rEpicardiumFile,
            const std::string &rRightVentricleFile,
            const std::string &rLeftVentricleFile,
            bool indexFromZero)
{
    // Compute the distance map of each surface
     mpGeometryInfo = new HeartGeometryInformation<SPACE_DIM>(mrMesh, rEpicardiumFile, rLeftVentricleFile, rRightVentricleFile, indexFromZero);
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::GenerateOrthotropicFibreOrientation(
            std::string outputDirectory,
            std::string fibreOrientationFile,
            bool logInfo)
{
    assert(SPACE_DIM == 3);
    if (mpGeometryInfo == NULL)
    {
        EXCEPTION("Files defining the heart surfaces not set");
    }
    unsigned num_elements = mrMesh.GetNumElements();

    // Open files
    OutputFileHandler handler(outputDirectory, false);
    out_stream p_fibre_file, p_regions_file, p_thickness_file, p_ave_thickness_file;

    //Make sure that only the master process makes an empty output file
    //and writes the log files
    if (PetscTools::AmMaster())
    {
        //Open file and close it, but don't write to it yet
        p_fibre_file=handler.OpenOutputFile(fibreOrientationFile);
        // First line of the fibre file: number of elements of the mesh
        *p_fibre_file << num_elements << std::endl;
        p_fibre_file->close();
    }
    else
    {
        logInfo=false;
    }

    if (logInfo)
    {
        p_regions_file  = handler.OpenOutputFile("node_regions.data");
        p_thickness_file = handler.OpenOutputFile("wall_thickness.data");
        p_ave_thickness_file = handler.OpenOutputFile("averaged_thickness.data");
    }


    //We expect that the apex to base has been set
    if (fabs(norm_2(mApexToBase)) < DBL_EPSILON)
    {
        EXCEPTION("Apex to base vector has not been set");
    }

    // Compute wall thickness parameter
    unsigned num_nodes = mrMesh.GetNumNodes();
    std::vector<double> wall_thickness(num_nodes);
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        double dist_epi, dist_endo;

        HeartRegionType node_region = mpGeometryInfo->GetHeartRegion(node_index);

        switch(node_region)
        {
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_WALL:
                dist_epi = mpGeometryInfo->rGetDistanceMapEpicardium()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_WALL:
                dist_epi = mpGeometryInfo->rGetDistanceMapEpicardium()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM:
                dist_epi = mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM:
                dist_epi = mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::UNKNOWN:
                #define COVERAGE_IGNORE
                std::cerr << "Wrong distances node: " << node_index << "\t"
                          << "Epi " << mpGeometryInfo->rGetDistanceMapEpicardium()[node_index] << "\t"
                          << "RV " << mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index] << "\t"
                          << "LV " << mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index]
                          << std::endl;

                // Make wall_thickness=0 as in Martin's code
                dist_epi = 1;
                dist_endo = 0;
                break;
                #undef COVERAGE_IGNORE

            default:
                NEVER_REACHED;
        }

        wall_thickness[node_index] = dist_endo / (dist_endo + dist_epi);

        if (std::isnan(wall_thickness[node_index]))
        {
            #define COVERAGE_IGNORE
            /*
             *  A node contained on both epicardium and lv (or rv) surfaces has wall thickness 0/0.
             *  By setting its value to 0 we consider it contained only on the lv (or rv) surface.
             */
            wall_thickness[node_index] = 0;
            #undef COVERAGE_IGNORE
        }

        if (logInfo)
        {
            *p_regions_file << node_region*100 << "\n";
            *p_thickness_file << wall_thickness[node_index] << "\n";
        }
    }

    /*
     *  For each node, average its value of e with the values of all the neighbours
     */
    std::vector<double> my_averaged_wall_thickness(num_nodes);
    std::vector<double> averaged_wall_thickness(num_nodes);
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        my_averaged_wall_thickness[node_index] = GetAveragedThicknessLocalNode(node_index, wall_thickness);
    }

    //Non-local information appear as zeros in the vector
    MPI_Allreduce(&my_averaged_wall_thickness[0], &averaged_wall_thickness[0], num_nodes,
                  MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    if (logInfo)
    {
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
             *p_ave_thickness_file << averaged_wall_thickness[node_index] << "\n";
        }

    }

    /*
     *  Compute the gradient of the averaged wall thickness at the centroid of each tetrahedron
     */
    c_vector<double,SPACE_DIM> grad_ave_wall_thickness;

    bool fibre_file_is_open= false;
    for (unsigned element_index=0; element_index<num_elements; element_index++)
    {
        //PRINT_VARIABLE(element_index);
        if (mrMesh.CalculateDesignatedOwnershipOfElement(element_index))
        {
            //PRINT_VARIABLE(element_index);
            //Locally own this element
            //Make sure that writing will proceed
            //The same file is written to by multiple processes, therefore we need
            //this barrier - Barrier happens BEFORE opening and AFTER closing

            PetscTools::Barrier("GenerateOrthotropicFibreOrientation");
            if (fibre_file_is_open != true)
            {
                //Open and append
                p_fibre_file = handler.OpenOutputFile(fibreOrientationFile, std::ios::app);
                fibre_file_is_open=true;
            }


            Element<SPACE_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(element_index);

            /*
             *  The gradient of the averaged thickness at the element is:
             *
             *     grad_ave_wall_thickness[element_index] = ave' * BF * inv(J)
             *
             *  being : ave, averaged thickness values of the nodes defining the element
             *          J,   the Jacobian of the element as defined in class Element.
             *                               (-1 -1 -1)
             *          BF,  basis functions ( 1  0  0)
             *                               ( 0  1  0)
             *                               ( 0  0  1)
             *
             *  Defined as u in Streeter paper
             */

            c_vector<double, SPACE_DIM+1> elem_nodes_ave_thickness;
            double element_averaged_thickness = 0.0;
            c_vector<HeartRegionType, SPACE_DIM+1> elem_nodes_region;

            for (unsigned local_node_index=0; local_node_index<SPACE_DIM+1; local_node_index++)
            {
                // Get node's global index
                unsigned global_node_index = p_element->GetNode(local_node_index)->GetIndex();

                elem_nodes_ave_thickness[local_node_index] = averaged_wall_thickness[global_node_index];
                elem_nodes_region[local_node_index] = mpGeometryInfo->GetHeartRegion(global_node_index);

                // Calculate wall thickness averaged value for the element
                element_averaged_thickness +=  wall_thickness[global_node_index];
            }

            element_averaged_thickness /= SPACE_DIM+1;

            c_matrix<double, SPACE_DIM+1, SPACE_DIM> basis_functions( zero_matrix<double>(4u,3u) );
            basis_functions(0,0) = basis_functions(0,1) = basis_functions(0,2) = -1.0;
            basis_functions(1,0) = basis_functions(2,1) = basis_functions(3,2) =  1.0;

            c_matrix<double, SPACE_DIM+1, SPACE_DIM> temp;
            c_matrix<double, SPACE_DIM, SPACE_DIM> jacobian, inverse_jacobian;
            double jacobian_det;
            mrMesh.GetInverseJacobianForElement(element_index, jacobian, jacobian_det, inverse_jacobian);
            noalias(temp) = prod (basis_functions, inverse_jacobian);
            noalias(grad_ave_wall_thickness) = prod(elem_nodes_ave_thickness, temp);

            grad_ave_wall_thickness /= norm_2(grad_ave_wall_thickness);

            /*
             *   Normal to the gradient (v in Streeter paper) which is then the circumferential direction
             * (it will be the fibre direction after rotation)
             *
             *  Computed as the cross product with the base-apex direction (originally assumed base-apex axis is x). The output vector is not normal,
             * since the angle between them may be != 90, normalise it.
             */
            c_vector<double, SPACE_DIM> fibre_direction = VectorProduct(grad_ave_wall_thickness, mApexToBase);
            fibre_direction /= norm_2(fibre_direction);

            /*
             *  Longitude direction (w in Streeter paper)
             */
            c_vector<double, SPACE_DIM> longitude_direction = VectorProduct(grad_ave_wall_thickness, fibre_direction);

            /*
             *  Compute fibre to v angle: alpha = R*(1-2e)^3
             *
             *    R is the maximum angle between the fibre and the v axis (heart region dependant)
             *    (1 - 2e)^3 scales it by a value in [-1, 1] defining the rotation of the fibre based
             *       on the position in the wall
             */
            double alpha = GetFibreMaxAngle(elem_nodes_region) * SmallPow( (1 - 2*element_averaged_thickness), 3 );

            /*
             *  Apply alpha rotation about the u axis to the orthonormal basis
             *
             *               ( u(1) v(1) w(1) )
             *   (u, v, w) = ( u(2) v(2) w(2) )
             *               ( u(3) v(3) w(3) )
             *
             *  The following matrix defines a rotation about the u axis
             *
             *                 ( 1        0           0      ) (u')
             *   R = (u, v, w) ( 0    cos(alpha) -sin(alpha) ) (v')
             *                 ( 0    sin(alpha)  cos(alpha) ) (w')
             *
             *  The rotated basis is computed like:
             *
             *                                             ( 1        0           0      )
             *  (u, v_r, w_r ) = R * (u, v, w) = (u, v, w) ( 0    cos(alpha) -sin(alpha) )
             *                                             ( 0    sin(alpha)  cos(alpha) )
             *
             *  Which simplifies to:
             *
             *   v_r =  v*cos(alpha) + w*sin(alpha)
             *   w_r = -v*sin(alpha) + w*cos(alpha)
             */
            c_vector<double, SPACE_DIM> rotated_fibre_direction = fibre_direction*cos(alpha) + longitude_direction*sin(alpha);
            c_vector<double, SPACE_DIM> rotated_longitude_direction = -fibre_direction*sin(alpha) + longitude_direction*cos(alpha);


            /*
             * Test the orthonormality of the basis
             */
            assert( fabs(norm_2(rotated_fibre_direction) - 1) < 100*DBL_EPSILON );
            assert( fabs(norm_2(grad_ave_wall_thickness) - 1) < 100*DBL_EPSILON );
            assert( fabs(norm_2(rotated_longitude_direction) - 1) < 100*DBL_EPSILON );

            assert( fabs(inner_prod(rotated_fibre_direction, grad_ave_wall_thickness)) < 100*DBL_EPSILON );
            assert( fabs(inner_prod(rotated_fibre_direction, rotated_longitude_direction)) < 100*DBL_EPSILON);
            assert( fabs(inner_prod(grad_ave_wall_thickness, rotated_longitude_direction)) < 100*DBL_EPSILON);

            /*
             *  Output the direction of the myofibre, the transverse to it in the plane of the myocite laminae and the normal to this laminae (in that order)
             *
             *  See Fig. 1 "Laminar Structure of the Heart: a mathematical model" LeGrice et al. 97
             *
             */
            *p_fibre_file << rotated_fibre_direction[0]     << " " << rotated_fibre_direction[1]     << " "  << rotated_fibre_direction[2]     << " "
                          << grad_ave_wall_thickness[0]     << " " << grad_ave_wall_thickness[1]     << " "  << grad_ave_wall_thickness[2]     << " "
                          << rotated_longitude_direction[0] << " " << rotated_longitude_direction[1] << " "  << rotated_longitude_direction[2] << std::endl;
        }
        else
        {
            //We don't locally own this element
            //Give up the file pointer
            if (fibre_file_is_open)
            {
                p_fibre_file->close();
                fibre_file_is_open=false;
            }
            PetscTools::Barrier("GenerateOrthotropicFibreOrientation");

        }
    }

    if (fibre_file_is_open)
    {
        p_fibre_file->close();
    }

    if (logInfo)
    {
        p_regions_file->close();
        p_thickness_file->close();
        p_ave_thickness_file->close();
    }

}


template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::SetApexToBase(const c_vector<double, SPACE_DIM>& apexToBase)
{
    double norm = norm_2(apexToBase);
    if (norm < DBL_EPSILON)
    {
        EXCEPTION("Apex to base vector should be non-zero");
    }
    mApexToBase = apexToBase / norm;
}
template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::SetApexToBase(unsigned axis)
{
    if (axis >= SPACE_DIM)
    {
        EXCEPTION("Apex to base coordinate axis was out of range");
    }
    mApexToBase = zero_vector<double>(SPACE_DIM);
    mApexToBase[axis] = 1.0;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////
#define COVERAGE_IGNORE
template class StreeterFibreGenerator<3>;
#undef COVERAGE_IGNORE