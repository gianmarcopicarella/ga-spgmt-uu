#pragma once

#include "Types.h"
#include <variant>

namespace SPGMT
{
	using Range = std::pair<int, int>;

	struct Result
	{
        struct RangeWrapper
        {
            Range myRange;
            int myRefIndex;
        };

        std::unordered_map<int, std::vector<int>> mySortedPlanesIndices;
        std::vector<RangeWrapper> myRangeWrappers;
	};

	struct PlaneWrapper
	{
		PlaneWrapper(const Plane& aPlane) : myPlane(aPlane)
		{
			const Vec3 up{ 0,0,1 };
			myHasNegativeNormal = CGAL::scalar_product(up, aPlane.orthogonal_vector()) < 0.f;
		}
		bool has_on_positive_side(const Point3& aPoint) const
		{
			return myHasNegativeNormal ?
				myPlane.has_on_negative_side(aPoint): 
				myPlane.has_on_positive_side(aPoint);
		}
		bool has_on_negative_side(const Point3& aPoint) const
		{
			return myHasNegativeNormal ?
				myPlane.has_on_positive_side(aPoint) :
				myPlane.has_on_negative_side(aPoint);
		}
		bool has_on(const Point3& aPoint) const
		{
			return myPlane.has_on(aPoint);
		}
		const Plane myPlane;
	private:
		bool myHasNegativeNormal;
	};

	std::vector<PlaneWrapper> GetPlaneWrappers(const std::vector<Plane>& somePlanes);
    Result BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints);
}