#pragma once

#include "Types.h"

#include <vector>

namespace SPGMT
{
	namespace Utils
	{
		template<typename T, typename F>
		bool AreItemsParallel(const std::vector<T>& someItems, F&& aGetter)
		{
			const auto& refDir = aGetter(someItems[0]);
			return std::all_of(std::next(someItems.begin()), someItems.end(), [&](const auto& anItem) {
				const auto& dir = aGetter(anItem);
				return refDir == dir || refDir == -dir;
			});
		}

		template<typename T>
		bool AreItemsUnique(const std::vector<T>& someItems)
		{
			auto areItemsUnique{ true };
			for (int i = 0; i < someItems.size() && areItemsUnique; ++i)
			{
				for (int k = i + 1; k < someItems.size(); ++k)
				{
					if (someItems[i] == someItems[k])
					{
						areItemsUnique = false;
						break;
					}
				}
			}
			return areItemsUnique;
		}

		bool ArePlanesNonVertical(const std::vector<Plane>& somePlanes);
		bool ArePlanesUniformlyOriented(const std::vector<Plane>& somePlanes);

		void FlipPlaneNormalsIfFacingDownwards(std::vector<Plane>& someOutPlanes);
		bool IsPlaneFacingUp(const Plane& aPlane);

		Point3 SingleDualMapping(const Plane& aPlane);
		Plane SingleDualMapping(const Point3& aPoint);

		std::vector<Point3> DualMapping(const std::vector<Plane>& somePlanes);
		std::vector<Plane> DualMapping(const std::vector<Point3>& somePoints);
	}
}