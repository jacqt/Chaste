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

#ifndef CONTINUUMMECHANICSPROBLEMDEFINITION_HPP_
#define CONTINUUMMECHANICSPROBLEMDEFINITION_HPP_

#include "QuadraticMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractIncompressibleMaterialLaw.hpp"
#include "AbstractCompressibleMaterialLaw.hpp"
#include "CompressibilityType.hpp"


/**
 *  Simple enumeration for denoting the type of body force, constant or defined via a function
 */
typedef enum BodyForceType_
{
    CONSTANT_BODY_FORCE,
    FUNCTIONAL_BODY_FORCE
} BodyForceType;

/**
 *  Simple enumeration for denoting the type of traction (Neumann) boundary condition:
 */
typedef enum TractionBoundaryConditionType_
{
    NO_TRACTIONS,
    ELEMENTWISE_TRACTION,
    FUNCTIONAL_TRACTION,
    PRESSURE_ON_DEFORMED
} TractionBoundaryConditionType;


/**
 *  A class for specifying various parts of a continuum mechanics problem, Dirichlet node information (which
 *  nodes are in space in a solid problem, which nodes have fixed flow in a fluids problem), the body force (per unit mass)
 *  (usually acceleration due to gravity or zero), the traction boundary conditions, and the density.
 */
template<unsigned DIM>
class ContinuumMechanicsProblemDefinition
{
public:
    /** Special value for Dirichlet nodes, indicating that a Dirichlet boundary condition
     *  in a particular dimension is not specified */
    static const double FREE; // set to DBL_MAX

protected:
    /** The mesh being solved on */
    QuadraticMesh<DIM>& mrMesh;

    /** Density of the body (constant throughout body) */
    double mDensity;

    //////////////////////////////
    // body force
    //////////////////////////////

    /** The body force type */
    BodyForceType mBodyForceType;

    /** The constant body force, only used if mBodyForceType is set appropriately */
    c_vector<double,DIM> mConstantBodyForce;

    /** The body force as a function of space and time, only used if mBodyForceType is set appropriately */
    c_vector<double,DIM> (*mpBodyForceFunction)(c_vector<double,DIM>& rX, double t);

    //////////////////////////////
    // tractions
    //////////////////////////////

    /** The traction (Neumann) boundary condition type */
    TractionBoundaryConditionType mTractionBoundaryConditionType;

    /** The surface elements on which tractions are applied */
    std::vector<BoundaryElement<DIM-1,DIM>*> mTractionBoundaryElements;

    /** The tractions on each surface element (only used if mTractionBoundaryConditionType is set appropriately) */
    std::vector<c_vector<double,DIM> > mElementwiseTractions;

    /** If the tractions are specified to correspond to a pressure acting on the surface: the pressure for the given
     *  boundary elements (only used if mTractionBoundaryConditionType is set appropriately) */
    double mNormalPressure;

    /** The tractions as a function of space and time (only used if mTractionBoundaryConditionType is set appropriately) */
    c_vector<double,DIM> (*mpTractionBoundaryConditionFunction)(c_vector<double,DIM>& rX, double t);

    ///////////////////////////////////////////
    // Dirichlet boundary conditions
    ///////////////////////////////////////////

    /**
     * All nodes (including non-vertices) which have a dirichlet boundary condition (ie position
     * prescribed in solid mechanics problems, flow prescribed in fluids problems). */
    std::vector<unsigned> mDirichletNodes;

    /** The values at the nodes with Dirichlet boundary conditions (displacement  */
    std::vector<c_vector<double,DIM> > mDirichletNodeValues;


public:
    /**
     * Constructor initialises the body force to zero and density to 1.0
     * @param rMesh  is the mesh being solved on
     */
    ContinuumMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh);

    /** Destructor */
    virtual ~ContinuumMechanicsProblemDefinition()
    {
    }


    /**
     *  Set the density
     *  @param density
     */
    void SetDensity(double density);

    /**
     * Get the density
     */
    double GetDensity();

    /**
     * Set a list of nodes (indices) to be given zero Dirichlet boundary condition
     * @param rZeroDirichletNodes the nodes at which the value (displacement/flow) is zero
     */
    void SetZeroDirichletNodes(std::vector<unsigned>& rZeroDirichletNodes);

    // No method here for setting non-zero Dirichlet BCs - these are in the subclasses..

    /**
     *  Get the Dirichlet nodes
     */
    std::vector<unsigned>& rGetDirichletNodes();

    /**
     * Get the Dirichlet node values.
     */
    std::vector<c_vector<double,DIM> >& rGetDirichletNodeValues();


    /**
     * Set the body force to be used - this version sets a constant body force
     * @param bodyForce the constant body force
     */
    void SetBodyForce(c_vector<double,DIM> bodyForce);

    /**
     * Set the body force to be used - this version sets a functional body force
     * @param pFunction a function of space and time returning a vector
     */
    void SetBodyForce(c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t));

    /**
     * Get the body force at a particular point and time.
     * Note: The user can either call this, or check what type of body force has been set using
     * GetBodyForceType() and then call GetConstantBodyForce() or EvaluateBodyForceFunction(X,t).
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> GetBodyForce(c_vector<double,DIM>& rX, double t = 0.0);

    /**
     * Get the body force type
     */
    BodyForceType GetBodyForceType();

    /**
     * Get the constant body force (error if GetBodyForceType()!=CONSTANT_BODY_FORCE)
     */
    c_vector<double,DIM> GetConstantBodyForce();

    /**
     * Evaluate the body force function (error if GetBodyForceType()!=FUNCTIONAL_BODY_FORCE)
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> EvaluateBodyForceFunction(c_vector<double,DIM>& rX, double t);

    /**
     * Get the traction (Neumann) boundary condition type
     */
    TractionBoundaryConditionType GetTractionBoundaryConditionType();

    /**
     * Set traction (Neumann) boundary conditions. This version takes in a vector of
     *  boundary elements, and corresponding tractions for each one.
     * @param rTractionBoundaryElements the boundary elements
     * @param rElementwiseTractions corresponding tractions
     */
    void SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*>& rTractionBoundaryElements,
                                       std::vector<c_vector<double,DIM> >& rElementwiseTractions);

    /**
     * Set traction (Neumann) boundary conditions. This version takes in a vector of
     *  boundary elements, and a function to be evaluated at points in these boundary
     *  elements
     * @param rTractionBoundaryElements the boundary elements
     * @param pFunction the traction function (a function of space and time, returning a vector)
     */
    void SetTractionBoundaryConditions(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                       c_vector<double,DIM> (*pFunction)(c_vector<double,DIM>& rX, double t));

    /**
     * Set traction (Neumann) boundary conditions. This version says that pressures should be applied
     * on surfaces in the DEFORMED body in the outward normal direction.
     *
     * @param rTractionBoundaryElements The boundary elements
     * @param normalPressure the corresponding pressure
     */
    void SetApplyNormalPressureOnDeformedSurface(std::vector<BoundaryElement<DIM-1,DIM>*> rTractionBoundaryElements,
                                                 double normalPressure);

    /**
     * Get the vector of traction boundary elements
     */
    std::vector<BoundaryElement<DIM-1,DIM>*>& rGetTractionBoundaryElements();

    /**
     *  Get the element-wise tractions vector (corresponding to
     *  vector returned by rGetTractionBoundaryElements())
     *  (error if GetTractionBoundaryConditionType()!=ELEMENTWISE_TRACTION)
     */
    std::vector<c_vector<double,DIM> >& rGetElementwiseTractions();

    /**
     *  Get the pressure for the boundary elements (corresponding to
     *  vector returned by rGetTractionBoundaryElements())
     *  (error if GetTractionBoundaryConditionType()!=PRESSURE_ON_DEFORMED)
     */
    double GetNormalPressure();


    /**
     * Evaluate the traction boundary condition function (error if GetTractionBoundaryConditionType()!=FUNCTIONAL_TRACTION)
     *
     * @param rX spatial location
     * @param t current time
     */
    c_vector<double,DIM> EvaluateTractionFunction(c_vector<double,DIM>& rX, double t);

    /**
     * Check all variables are set appropriately. Exceptions are thrown if any are not.
     * Derived classes can override but should call this version as well.
     */
    virtual void Validate();
};


#endif // CONTINUUMMECHANICSPROBLEMDEFINITION_HPP_