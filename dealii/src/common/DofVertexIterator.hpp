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


#ifndef DOFVERTEXITERATOR_HPP_
#define DOFVERTEXITERATOR_HPP_

#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <cassert>

/**
 *  DofVertexIterator
 *
 *  See also TriangulationVertexIterator
 *
 *  Looping over nodes (vertices) in deal.II meshes is a pain, you have to loop over
 *  elements (cells), then loop over each vertex of that cell, checking you haven't
 *  already done that vertex. This class takes care of this for you.
 *
 *  The difference between this class and TriangulationVertexIterator is that the cell
 *  is obtained through a DoFHandler<DIM>::active_cell_iterator rather than a
 *  Triangulation<DIM>::active_cell_iterator, which means the dofs for the vertex
 *  can be obtained.
 *
 *  Note: vertices are reached by looping over elements (internally, in this class), 
 *  therefore this class only reaches ACTIVE vertices (ie ones that are part of an
 *  element). (Unactive vertices arise if the mesh is coarsened)
 * 
 *  Note: this class doesn't act like a typical std::iterator class, in that the
 *  result doesn't point to a vertex, you don't call ++ etc. (although in the
 *  future ++ might be implemented to call Next()
 *
 *  Usage:
 *
 *  DofVertexIterator<2> iter(&mesh,&dof);
 *  while(!iter.ReachedEnd())
 *  {
 *      Point<2> vertex = iter.GetVertex(); //etc
 *
 *      double value_at_node = solutions( iter.GetDof(0) ); // here solutions is a Vector<double>
 *
 *      iter.Next();
 *  }
 *
 *  The vertex number, current cell, local number of the vertex in that cell can
 *  also be obtained.
 *
 *  TODO: perhaps extract commonality between this and TriangulationVertexIterator. 
 *  (although the code is almost identical mpCell are different types in the two classes)
 */
template<unsigned DIM>
class DofVertexIterator
{
private :
    Triangulation<DIM>* mpMesh;
    DoFHandler<DIM>* mpDofHandler;
    std::vector<bool> mVertexTouched;
    typename DoFHandler<DIM>::active_cell_iterator mpCurrentCell;
    
    unsigned mCurrentVertexIndex;
    bool mReachedEnd;
    
    void NextNode()
    {
        if (mCurrentVertexIndex < GeometryInfo<DIM>::vertices_per_cell-1)
        {
            mCurrentVertexIndex++;
        }
        else
        {
            mCurrentVertexIndex = 0;
            mpCurrentCell++;
        }
    }
    
public :
    DofVertexIterator(Triangulation<DIM>* pMesh, DoFHandler<DIM>* pDofHandler) : 
              mpMesh(pMesh),
              mpDofHandler(pDofHandler),
              mVertexTouched(mpMesh->n_vertices(),false),
              mpCurrentCell(mpDofHandler->begin_active())
    {
        assert(pMesh && pDofHandler);
        
        mCurrentVertexIndex = 0;
        mReachedEnd = (mpCurrentCell==mpDofHandler->end());
        
        // set the current node as having been touched
        if(!mReachedEnd)
        {
            mVertexTouched[GetVertexGlobalIndex()] = true;
        }
    }
    
    /**
     *  Move to the next vertex
     */
    void Next()
    {
        bool found=false;
        while ( (found==false) && (!mReachedEnd) )
        {
            NextNode();
            
            mReachedEnd = (mpCurrentCell==mpDofHandler->end());
            if ( !mReachedEnd && !mVertexTouched[mpCurrentCell->vertex_index(mCurrentVertexIndex)] )
            {
                found = true;
                mVertexTouched[GetVertexGlobalIndex()] = true;
            }
        }
    }
    
    /**
     *  The method returns true if the last vertex has been passed
     */
    bool ReachedEnd()
    {
        return mReachedEnd;
    }
    
    /**
     *  Get the position of the vertex
     */
    Point<DIM>& GetVertex()
    {
        assert(!mReachedEnd);
        return mpCurrentCell->vertex(mCurrentVertexIndex);
    }
    
    /**
     *  Get the global vertex index
     */
    unsigned GetVertexGlobalIndex()
    {
        assert(!mReachedEnd);
        return mpCurrentCell->vertex_index(mCurrentVertexIndex);
    }
    
    /**
     *  Get the index of the current vertex in the current cell
     *  To be used with GetCell()
     */
    unsigned GetLocalVertexIndexForCell()
    {
        assert(!mReachedEnd);
        return mCurrentVertexIndex;
    }
    
    /**
     *  Get the current cell. Together with GetLocalVertexIndexForCell() this
     *  specifies the vertex. Needed for calling other methods on the cell. 
     */
    typename DoFHandler<DIM>::active_cell_iterator GetCell()
    {
        assert(!mReachedEnd);
        return mpCurrentCell;
    }
    
    /**
     * Get the Dof associated with the ith unknown at the current vertex
     */
    unsigned GetDof(unsigned i)
    {
        return mpCurrentCell->vertex_dof_index(mCurrentVertexIndex,i);
    }
    
    /**
     *  Reset to the first vertex so the iterator can be used again
     */
    void Reset()
    {
        // resize mVertexTouched in case the mesh has been refined..
        if(mVertexTouched.size()!=mpMesh->n_vertices())
        {
            mVertexTouched.resize(mpMesh->n_vertices());
        }
        
        for (unsigned i=0; i<mpMesh->n_vertices(); i++)
        {
            mVertexTouched[i] = false;
        }
        mpCurrentCell = mpDofHandler->begin_active();
        mCurrentVertexIndex = 0;
        mReachedEnd = (mpCurrentCell==mpDofHandler->end());
        
        // set the current node as having been touched
        if(!mReachedEnd)
        {
            mVertexTouched[GetVertexGlobalIndex()] = true;
        }
    }
};

#endif /*DOFVERTEXITERATOR_HPP_*/
