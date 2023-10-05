#pragma once

#include <vector>
#include <CGAL/Random.h>
#include <map>
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
		std::map<std::pair<int, int>, std::vector<int>> myTrianglesIndicesAtEdge;
	};

	std::vector<Vertex> ComputeLowerEnvelopeSmart(const std::vector<Plane>& somePlanes);

	struct Face;
	TriangulationData TriangulateLowerEnvelopeFaces(const std::vector<Face>& someFaces, const bool aComputeAllEdges = true);

}