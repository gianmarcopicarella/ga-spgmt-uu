#pragma once

#include <vector>
#include <CGAL/Random.h>

#include "Types.h"

namespace SPGMT
{
	struct Vertex;

	template<typename T>
	std::vector<int> SampleWithProbability(const std::vector<T>& someItems, const FT aProbability)
	{
		std::vector<int> samplesIndices;
		CGAL::Random random{};

		for (int i = 0; i < someItems.size(); ++i)
		{
			if (random.uniform_real() <= aProbability)
			{
				samplesIndices.push_back(i);
			}
		}

		return samplesIndices;
	}

	struct TriangulationData
	{
		std::vector<int> myVerticesIndices;
	};

	std::vector<Vertex> ComputeLowerEnvelopeSmart(const std::vector<Plane>& somePlanes);

	struct Face;
	std::vector<TriangulationData> TriangulateLowerEnvelopeFaces(const std::vector<Face>& someFaces, const bool aComputeAllEdges = true);

}