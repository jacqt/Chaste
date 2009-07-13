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


#ifndef _TRIANGLESMESHREADER_HPP_
#define _TRIANGLESMESHREADER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include "AbstractTetrahedralMeshReader.hpp"

/**
 * Concrete version of the AbstractCachedMeshReader class.
 * Once constructed the public methods of the AbstractCachedMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshReader : public AbstractTetrahedralMeshReader<ELEMENT_DIM,SPACE_DIM>
{
private:

    bool mIndexFromZero;            /**< True if input data is numbered from zero, false otherwise */

    std::string mFilesBaseName;     /**< The base name for mesh files. */

    std::ifstream mNodesFile;       /**< The nodes file for the mesh. */
    std::ifstream mElementsFile;    /**< The elements file for the mesh. */
    std::ifstream mFacesFile;       /**< The faces (edges) file for the mesh. */

    unsigned mNumNodes;             /**< Number of nodes in the mesh. */
    unsigned mNumElements;          /**< Number of elements in the mesh. */
    unsigned mNumFaces;             /**< Number of faces in the mesh. */

    unsigned mNodesRead;            /**< Number of nodes read in. */
    unsigned mElementsRead;         /**< Number of elements read in. */
    unsigned mFacesRead;            /**< Number of faces read in. */
    unsigned mBoundaryFacesRead;    /**< Number of boundary faces read in. */

    unsigned mNumNodeAttributes;    /**< Is the number of attributes stored at each node. */
    unsigned mMaxNodeBdyMarker;     /**< Is the maximum node boundary marker. */
    unsigned mNumElementNodes;      /**< Is the number of nodes per element. */
    unsigned mNumElementAttributes; /**< Is the number of attributes stored for each element. */
    unsigned mNumFaceAttributes;    /**< Is the number of attributes stored for each face. */

    unsigned mOrderOfElements;      /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mOrderOfBoundaryElements; /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mNodesPerElement;      /**< The number of nodes contained in each element. */
    unsigned mNodesPerBoundaryElement; /**< The number of nodes in each boundary element. */

public:

    /**
     * Constructor.
     *
     * @param pathBaseName  the base name of the files from which to read the mesh data
     * @param orderOfElements  the order of each element: 1 for linear, 2 for quadratic (defaults to 1)
     * @param orderOfBoundaryElements the order of each boundary element: 1 for linear, 2 for quadratic (defaults to 1. May
     *  or may not be different to orderOfElements (Note tetgen with the -o2 flag creates quadratic elements but doesn't 
     *  create quadratic faces, hence the need for this third parameter)
     *   
     */
    TrianglesMeshReader(std::string pathBaseName, unsigned orderOfElements=1, unsigned orderOfBoundaryElements=1);

    /** Returns the number of elements in the mesh */
    unsigned GetNumElements() const;

    /** Returns the number of nodes in the mesh */
    unsigned GetNumNodes() const;

    /** Returns the number of faces in the mesh (synonym of GetNumEdges()) */
    unsigned GetNumFaces() const;

    /** Returns the number of edges in the mesh (synonym of GetNumFaces()) */
    unsigned GetNumEdges() const;

    /** Returns the number of attributes in the mesh */
    unsigned GetNumElementAttributes() const;

    /** Returns the number of attributes in the mesh */
    unsigned GetNumFaceAttributes() const;

    /** Resets pointers to beginning*/
    void Reset();

    /** Returns a vector of the coordinates of each node in turn */
    std::vector<double> GetNextNode();

    /** Returns a vector of the nodes of each element (and any attribute infomation, if there is any) in turn */
    ElementData GetNextElementData();

    /** Returns a vector of the nodes of each face in turn (synonym of GetNextEdgeData()) */
    ElementData GetNextFaceData();

    /** Returns a vector of the nodes of each edge in turn (synonym of GetNextFace()). */
    ElementData GetNextEdgeData();

private:

    /** Open mesh files. */
    void OpenFiles();

    /** Open node file. \todo Change name to OpenNodesFile for consistency with OpenElementsFile and OpenFacesFile? (#991) */
    void OpenNodeFile();

    /** Open elements file. */
    void OpenElementsFile();

    /** Open faces file. */
    void OpenFacesFile();

    /** Read the header from each mesh file. */
    void ReadHeaders();

    /** Close mesh files. */
    void CloseFiles();

    /**
     * Read in the next line.
     *
     * @param fileStream
     * @param rRawLine
     */
    void GetNextLineFromStream(std::ifstream& fileStream, std::string& rRawLine);

    /** Get method for mFilesBaseName. */
    std::string GetMeshFileBaseName();

};

#endif //_TRIANGLESMESHREADER_HPP_
