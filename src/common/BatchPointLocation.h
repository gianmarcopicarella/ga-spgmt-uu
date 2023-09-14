#pragma once

#include "Types.h"
#include <variant>

namespace SPGMT
{
	using Range = std::pair<int, int>;

	struct OnePlaneResult
	{
		// True -> Plane above point, False -> Plane below or point contained in plane
		std::vector<bool> myIsPointCovered;
	};

	struct ParallelPlanesResult
	{

		std::vector<int> mySortedPlanesIndices;
		std::vector<Range> myRanges;
	};

	struct BaseResult
	{
		struct ZoneRange
		{
			Range myRange;
			int myZoneIndex;
			bool myIsFaceZone;
		};

		std::unordered_map<int, std::vector<int>> mySortedPlanesPerFaceZone;
		std::unordered_map<int, std::vector<int>> mySortedPlanesPerLineZone;
		// Each pair contains the zone index to be used to retrieve the sorted planes list
		// and the range of planes above the point in that list
		std::vector<ZoneRange> myZoneRangesPairs;
	};

	using LocationResult = std::variant<std::monostate, OnePlaneResult, ParallelPlanesResult, BaseResult>;

	namespace Debug
	{
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
	}

	LocationResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints);
}