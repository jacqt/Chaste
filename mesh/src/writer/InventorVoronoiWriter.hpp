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


#ifndef INVENTORVORONOIWRITER_HPP_
#define INVENTORVORONOIWRITER_HPP_

#include "VoronoiTessellation.hpp"
#include "OutputFileHandler.hpp"


const std::string INVENTOR_HEADER="#Inventor V2.0 ascii \n\
\n\
Separator { \n\
  Info { \n\
    string \"tetrahedron.iv generated by IVREAD.\" \n\
    string \"Original data in file tetrahedron.wrl.\" \n\
  } \n\
  Separator { \n\
    LightModel { \n\
      model BASE_COLOR \n\
    } \n\
\n\
    Material { \n\
      ambientColor  0.2 0.2 0.2 \n\
\n\
      emissiveColor 0.0 0.0 0.0 \n\
      specularColor 1.0 1.0 1.0 \n\
      shininess     0.2 \n\
\n\
    } \n\
    ShapeHints { \n\
      vertexOrdering COUNTERCLOCKWISE \n\
      shapeType UNKNOWN_SHAPE_TYPE \n\
      faceType CONVEX \n\
      creaseAngle 6.28319 \n\
    } \n\
\n\
    DrawStyle { \n\
        style           LINES \n\
        lineWidth           3 \n\
        linePattern       255 \n\
    } \n\
\n\
    Coordinate3 { \n\
      point [ \n\
";


////  This is what we will need to change the header to changing \n to be \n with a backslash on it:
//const std::string INVENTOR_HEADER="#Inventor V2.1 ascii \n
//\n
//Separator { \n
//    IndexedFaceSet { \n
//    vertexProperty      VertexProperty { \n
//        vertex  [ \n
//";
//
//const std::string INVENTOR_MID="                ] \n
//                texCoord    [  ] \n
//                orderedRGBA [  ] \n
//                materialBinding PER_VERTEX \n
//                normalBinding   PER_PART \n
//    }  \n
//      coordIndex [ \n
//";
//
// const std::string INVENTOR_FOOTER="      ]   \n
//      \n
//  }   \n
// Translation {  \n
//    translation 3 0 0  \n
//    }  \n
//}   \n
//";
////  NOTE: You also need to change INVENTOR_MID used between cells to "\n\"
// 
// 
 
      
const std::string INVENTOR_MID="      ] \n\
    } \n\
    IndexedFaceSet { \n\
      coordIndex [ \n\
";

const std::string INVENTOR_FOOTER="      ] \n\
    } \n\
  } \n\
} \n\
";



class InventorVoronoiWriter
{
protected:
    OutputFileHandler *mpOutputFileHandler; /**< Output file handler */
    std::string mBaseName; /**< Base name for the input files */
    
public:
    /** Constructor */
    InventorVoronoiWriter(const std::string &rDirectory,
                          const std::string &rBaseName,
                          const bool clearOutputDir=true)
            : mBaseName(rBaseName)
    {
        mpOutputFileHandler = new OutputFileHandler(rDirectory, clearOutputDir);
    }
    
    /** Destructor */
    ~InventorVoronoiWriter()
    {
        delete mpOutputFileHandler;
    }
    
    /**
     *  Write the voronoi tessellation in Inventor format
     */
    void Write(const VoronoiTessellation<3>& rTessellation)
    {
        // open inventor file
        std::string file_name = this->mBaseName+".iv";
        out_stream p_file = this->mpOutputFileHandler->OpenOutputFile(file_name);
        
        // write out header part of file      

        *p_file << INVENTOR_HEADER;
        
        // write out vertices
        // and construct map from pointer to vertex to vertex number
        
        std::map< c_vector<double, 3>*, unsigned> vertex_number_map;
        for ( unsigned vertex_number=0;
              vertex_number<rTessellation.mVertices.size();
              vertex_number++ )
        {
            c_vector<double ,3>& vertex=*(rTessellation.mVertices[vertex_number]);
            *p_file << "        " << vertex(0) << " " << vertex(1) << " " << vertex(2) << ",\n";
            
            vertex_number_map[rTessellation.mVertices[vertex_number]]=vertex_number;
        }
        
        *p_file << INVENTOR_MID; //  CHANGE this to: *p_file << "\n";
        
        // write out faces;
        
        for (unsigned face_number=0;
             face_number < rTessellation.mFaces.size();
             face_number++)
        {
            *p_file << "        ";
            Face<3>& face=*(rTessellation.mFaces[face_number]);
            for (unsigned vertex_local_number=0;
                 vertex_local_number < face.mVertices.size();
                 vertex_local_number++)
            {
                // note this assumes we can definitely find the vertex in the map
                unsigned vertex_number=vertex_number_map[face.mVertices[vertex_local_number]];
                *p_file << vertex_number << ", ";
            }
            *p_file << "\n";
        }
        *p_file << INVENTOR_FOOTER;
    }
    
    /**
     *  Scale the vertex of each cell toward the centre of that cell by the given scaleFactor
     *  and write.
     */
    void ScaleAndWrite(VoronoiTessellation<3>& rTessellation, double scaleFactor)
    {
        if ((scaleFactor <= 0.0) || (scaleFactor > 1.0))
        {
            EXCEPTION("scaleFactor should be between 0 and 1");
        }
        
        // open inventor file
        std::string file_name = this->mBaseName+".iv";
        out_stream p_file = this->mpOutputFileHandler->OpenOutputFile(file_name);
        
        // write out header part of file      
        *p_file << INVENTOR_HEADER;
        
        unsigned global_vertex_number = 0;

        // the face data which will be written to file afterwards
        std::vector<std::vector<unsigned> > new_faces_data;
        std::vector<unsigned> number_faces_per_cell;
        
        // loop over cells and write out scaled vertices for each one, storing face info as we go
        for (unsigned cell_index = 0; cell_index<rTessellation.mVoronoiCells.size(); cell_index++)
        {
            c_vector<double, 3>& r_cell_centre = rTessellation.mVoronoiCells[cell_index].rGetVoronoiCellCentre();

            // map from position to (new) global vertex number, for this cell only
            std::map< c_vector<double, 3>*, unsigned> vertex_number_map;

            const VoronoiCell& r_cell = rTessellation.mVoronoiCells[cell_index];
            
            for (unsigned face_number=0; face_number<r_cell.mFaces.size(); face_number++)
            {
                std::vector<unsigned> face_vertex_data;

                Face<3>& r_face = *(r_cell.mFaces[face_number]);
                for (unsigned face_vertex_number=0; face_vertex_number<r_face.mVertices.size(); face_vertex_number++)
                {
                    unsigned global_number_for_this_vertex;
                    
                    // see if vertex is in the map
                    std::map< c_vector<double,3>*,unsigned>::iterator iter = vertex_number_map.find(r_face.mVertices[face_vertex_number]);
                    if(iter!=vertex_number_map.end())
                    {
                        global_number_for_this_vertex = iter->second;
                    }
                    else
                    {
                        global_number_for_this_vertex = global_vertex_number;

                        // not in the map, so add it to map
                        vertex_number_map[r_face.mVertices[face_vertex_number]] = global_number_for_this_vertex;
                        global_vertex_number++;
                        
                        // scale the vertex and print out the new position
                        c_vector<double,3> new_vertex = *(r_face.mVertices[face_vertex_number]);
                        new_vertex = scaleFactor*(new_vertex - r_cell_centre) + r_cell_centre;

                        *p_file << "        " << new_vertex(0) << " " << new_vertex(1) << " " << new_vertex(2) << ",\n";
                    }
                    
                    // store this vertex's global number as a vertex for this face
                    face_vertex_data.push_back(global_number_for_this_vertex);
                }

                // add the vertex data for this face to the store
                new_faces_data.push_back( face_vertex_data );
            }
            
            // store how many faces were in this cell
            number_faces_per_cell.push_back(r_cell.mFaces.size());
        }

        unsigned index=0;

        // write the face info
        // Loop over cells
        for(unsigned i=0; i<number_faces_per_cell.size(); i++)
        {
            if(number_faces_per_cell[i]>0)
            {
                *p_file << INVENTOR_MID;
                // Loop over faces
                for(unsigned j=0; j<number_faces_per_cell[i]; j++)
                {
                    *p_file << "        ";
                    assert(index<new_faces_data.size());
                    if ((rTessellation.rGetCell(i).mOrientations)[j])
                    {
                        for(unsigned k=0; k<new_faces_data[index].size(); k++)
                        {
                            *p_file << new_faces_data[index][k] << ", ";
                        }
                     }
                     else
                     {
                        for(unsigned k=0; k<new_faces_data[index].size(); k++)
                        {
                            unsigned l = new_faces_data[index].size() - k - 1;
                            *p_file << new_faces_data[index][l] << ", ";
                        }
                     }
                     *p_file << "\n";
                     index++;
                }
            }
        }

        *p_file << INVENTOR_FOOTER;
    }
};

#endif /*INVENTORVORONOIWRITER_HPP_*/
