#ifndef _ABSTRACTCOUPLEDPDE_HPP_
#define _ABSTRACTCOUPLEDPDE_HPP_

#include "AbstractLinearParabolicPde.hpp"
#include <vector>
#include "petscvec.h"

typedef std::vector<double> odeVariablesType;

template <int SPACE_DIM>
class AbstractCoupledPde : public AbstractLinearParabolicPde<SPACE_DIM>
{

public:
    // timestep used in the ode solvers        
    double mSmallTimeStep;

    // timestep used by the pde solver
    double mBigTimeStep;

	// simulation time
    double mTime;   

    AbstractIvpOdeSolver *mpOdeSolver;

        // number of nodes in the mesh 
    int mNumNodes;
        
    // Lowest value of index that this part of the global object stores
    int mOwnershipRangeLo;
        
    // One more than the local highest index
    int mOwnershipRangeHi;
    
    /** mOdeVarsAtNode[i] is a vector of the current values of the
     *  voltage, gating variables, intracellular Calcium concentration at node 
     *  i. The voltage returned by the ode solver is not used later since the pde
     *  solves for the voltage.
     */
     // Distributed
    std::vector<odeVariablesType>            mOdeVarsAtNode;

    
public:
    // Replicated
  	std::vector<double>	solutionCache;
    
        //Constructor
    AbstractCoupledPde(int numNodes, AbstractIvpOdeSolver *pOdeSolver, double tStart, double bigTimeStep, double smallTimeStep)
    {
        assert(smallTimeStep < bigTimeStep + 1e-10);
        assert(numNodes > 0);
        
        mNumNodes=numNodes;
        mBigTimeStep=bigTimeStep;
        mpOdeSolver=pOdeSolver;
        mSmallTimeStep=smallTimeStep;
     
        mTime = tStart;
        
        // Resize vectors to the appropriate size for each process
        // Create a PETSc vector and use the ownership range of the PETSc vector to size our C++ vectors
        Vec tempVec;
        VecCreate(PETSC_COMM_WORLD, &tempVec);
        VecSetSizes(tempVec, PETSC_DECIDE, numNodes);
        VecSetFromOptions(tempVec);
        VecGetOwnershipRange(tempVec,&mOwnershipRangeLo,&mOwnershipRangeHi);
        VecDestroy(tempVec); // no longer needed
        
        mOdeVarsAtNode.resize(mOwnershipRangeHi-mOwnershipRangeLo);
      
        solutionCache.resize(mNumNodes);


     }
     
     virtual void PrepareForAssembleSystem(Vec currentSolution)
     {	
     	std::cout<<"AbstractCoupledPde::PrepareForAssembleSystem\n";
     }

};        
        
#endif //_ABSTRACTCOUPLEDPDE_HPP_
