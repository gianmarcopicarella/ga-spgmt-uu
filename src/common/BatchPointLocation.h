#pragma once

#include "Types.h"

namespace SPGMT
{
	enum class ResultType
	{
		NONE = 0,
		PARALLEL_PLANES = 1,
		SINGLE_PLANE = 2
	};

	struct LocationResult
	{
		ResultType myType{ ResultType::NONE };
		std::vector<std::pair<int, int>> myRanges;
		std::vector<int> mySortedPlanesIndices;
	};

	namespace Debug
	{
		template<typename T>
		bool AreItemsUnique(const std::vector<T>& someItems)
		{
			auto areItemsUnique { true };
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

	void BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints, LocationResult& anOutResult);
}