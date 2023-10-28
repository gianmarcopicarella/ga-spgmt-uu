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

        std::unordered_map<int, std::vector<size_t>> mySortedPlanesIndices;
        std::vector<RangeWrapper> myRangeWrappers;



        // New storage ---------------
        using SlabCache = std::unordered_map<size_t, std::vector<size_t>>;
        std::vector<SlabCache> mySortedPlanesCache;

        struct NewRangeWrapper
        {
            Range myRange;
            size_t myCacheIndex;
            size_t myIndex;
        };

        std::vector<NewRangeWrapper> myNewRangeWrappers;

	};

	BatchPointResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints);

    BatchPointResult ParallelBatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints);
}