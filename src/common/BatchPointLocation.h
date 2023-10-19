#pragma once

#include "Types.h"
#include <variant>

namespace SPGMT
{
	using Range = std::pair<int, int>;

	struct BatchPointResult
	{
        struct RangeWrapper
        {
            Range myRange;
            int myRefIndex;
        };

        std::unordered_map<int, std::vector<int>> mySortedPlanesIndices;
        std::vector<RangeWrapper> myRangeWrappers;
	};

	BatchPointResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints);

    BatchPointResult ParallelBatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints);
}