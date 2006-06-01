#ifndef _VARYINGDIFFUSIONANDSOURCETERMPDE_
#define _VARYINGDIFFUSIONANDSOURCETERMPDE_

#include "AbstractLinearEllipticPde.hpp"
#include "Point.hpp"
#include <cmath>

/**
 * A more complex linear elliptic PDE used in tests. The source and diffusion terms
 * depend on x.
 */

template <int SPACE_DIM>
class VaryingDiffusionAndSourceTermPde : public AbstractLinearEllipticPde<SPACE_DIM>
{
private:
	double DistanceFromOrigin(Point<SPACE_DIM> x)
	{
		double sum=0;
		for (int i=0; i<SPACE_DIM; i++)
		{
			sum += x[i]*x[i];
		}
		return sqrt(sum);
	}
	
public:
	double ComputeLinearSourceTerm(Point<SPACE_DIM> x)
	{
		return pow(DistanceFromOrigin(x),3);
	}
    
    double ComputeNonlinearSourceTerm(Point<SPACE_DIM> , double )
    {
    	return 0.0;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(Point<SPACE_DIM> x)
    {
    	return pow(DistanceFromOrigin(x),2)*identity_matrix<double>(SPACE_DIM);
    }
};

#endif //_VARYINGDIFFUSIONANDSOURCETERMPDE_
