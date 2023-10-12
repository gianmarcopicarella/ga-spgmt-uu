#pragma once

#include <vector>
#include <CGAL/Random.h>
#include <map>
#include "Types.h"
#include "BruteForce.h"

namespace SPGMT
{
	template<typename T>
	std::vector<int> SampleWithProbability(const std::vector<T>& someItems, const double aProbability)
	{
		std::vector<int> samplesIndices;
		CGAL::Random random{};

		for (int i = 0; i < someItems.size(); ++i)
		{
			if (random.get_double() <= aProbability)
			{
				samplesIndices.push_back(i);
			}
		}

		return samplesIndices;
	}

	LowerEnvelope3d ComputeLowerEnvelopeSmart(const std::vector<Plane>& somePlanes);

}