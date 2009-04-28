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
#ifndef ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_
#define ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_

#include <map>
#include <set>
#include "AbstractBoundaryCondition.hpp"
#include "Node.hpp"
//#include "ConstBoundaryCondition.hpp"
//#include "TetrahedralMesh.hpp"
//#include "LinearSystem.hpp"
//#include "PetscException.hpp"

/**
 * Helper struct storing an operator for computing whether one node 
 * has a lower index than another.
 */
template<unsigned SPACE_DIM>
struct LessThanNode
{
    /**
     * Less-then node index comparison operator.
     * 
     * @param n1 pointer to a node
     * @param n2 pointer to a node
     */
    bool operator()(const Node<SPACE_DIM> * const &n1, const Node<SPACE_DIM> * const &n2)
    {
        return (n1->GetIndex() < n2->GetIndex() );
    }
};

/**
 * Abstract boundary conditions container.
 */
template<unsigned ELEM_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractBoundaryConditionsContainer
{
protected:
    std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >
    *mpDirichletMap[PROBLEM_DIM]; /**< List (map) of Dirichlet boundary conditions */

    typename std::map< const Node<SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*, LessThanNode<SPACE_DIM> >::const_iterator
    mDirichIterator; /**< Internal iterator over Dirichlet boundary conditions */

public:
    /**
     * Constructor allocates memory for the Dirichlet boundary conditions lists.
     */
    AbstractBoundaryConditionsContainer();

    /**
     * Destructor.
     */
    ~AbstractBoundaryConditionsContainer();

    /** {
        assert(indexOfUnknown < PROBLEM_DIM);

        this->mDirichIterator = this->mpDirichletMap[indexOfUnknown]->find(pNode);

        return (this->mDirichIterator != this->mpDirichletMap[indexOfUnknown]->end());
    }
     * Return whether any Dirichlet conditions are defined.
     */
    bool HasDirichletBoundaryConditions();

    /**
     * Delete list of Dirichlet boundary conditions.
     * 
     * @param deletedConditions (optional)
     */
    void DeleteDirichletBoundaryConditions(std::set<const AbstractBoundaryCondition<SPACE_DIM>*> deletedConditions = std::set<const AbstractBoundaryCondition<SPACE_DIM>*>());

    /**
     * Obtain value of Dirichlet boundary condition at specified node.
     *
     * This is unlikely to be needed by the user, the methods ApplyDirichletToLinearProblem or
     * ApplyDirichletToNonlinearProblem can be called instead to apply all Dirichlet boundary conditions
     * at the same time.
     * 
     * @param pBoundaryNode pointer to a boundary node
     * @param indexOfUnknown index of the unknown for which to obtain the value of the boundary condition (defaults to 0)
     */
    double GetDirichletBCValue(const Node<SPACE_DIM>* pBoundaryNode, unsigned indexOfUnknown = 0);

    /**
     * Test if there is a Dirichlet boundary condition defined on the given node.
     *
     * \todo Perhaps have flag in node object for efficiency?
     * 
     * @param pNode pointer to a node
     * @param indexOfUnknown index of the unknown for which to obtain the value of the boundary condition (defaults to 0)
     */
    bool HasDirichletBoundaryCondition(const Node<SPACE_DIM>* pNode, unsigned indexOfUnknown = 0);
};

#endif /*ABSTRACTBOUNDARYCONDITIONSCONTAINER_HPP_*/
