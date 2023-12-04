#pragma once

#include "Types.h"

namespace SPGMT
{
	using Range = std::pair<size_t, size_t>;

	struct BatchPointResult
	{
        using SlabCache = std::unordered_map<size_t, std::vector<size_t>>;
        std::unordered_map<size_t, SlabCache> mySortedPlanesCache;

        struct RangeWrapper
        {
            Range myRange;
            size_t myCacheIndex;
            size_t myIndex;
        };

        std::vector<RangeWrapper> myRangeWrappers;
	};

    BatchPointResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints);
}