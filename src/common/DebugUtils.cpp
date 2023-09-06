#include "DebugUtils.h"

#include <CGAL/point_generators_3.h>

namespace SPGMT
{
	namespace Debug
	{
		PointSet3 Uniform3DCubeSampling(const double anHalfSide, const int aSampleCount)
		{
			PointSet3 samples;
			samples.reserve(aSampleCount);
			CGAL::Random_points_in_cube_3<Point3> generator{ anHalfSide };
			std::copy_n(generator, aSampleCount, samples.point_back_inserter());
			// (NOT MANDATORY) Use a random permutation to hide the creation history of the point set.
			// CGAL::cpp98::random_shuffle(samples.begin(), samples.end()); 
			return samples;
		}
	}
}