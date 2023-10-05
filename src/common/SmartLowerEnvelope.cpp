#include "SmartLowerEnvelope.h"

#include "BruteForce.h"

namespace SPGMT
{
	namespace
	{
		void locUpdateTrianglesLocalityMap(const int aStartingIndex, const int anEndIndex, const int aTriangleIndex,
			std::map<std::pair<int, int>, std::vector<int>>& someOutTrianglesIndicesAtEdge)
		{
			const auto& edgeKey = std::make_pair(std::min(aStartingIndex, anEndIndex), std::max(aStartingIndex, anEndIndex));
			auto& firstEdgeEntry = someOutTrianglesIndicesAtEdge.find(edgeKey);
			if (firstEdgeEntry == someOutTrianglesIndicesAtEdge.end())
			{
				firstEdgeEntry = someOutTrianglesIndicesAtEdge.insert(std::make_pair(edgeKey, std::vector<int>{})).first;
			}
			firstEdgeEntry->second.push_back(aTriangleIndex);
		}

		// No reference passed
		void locRecursiveLowerEnvelope(const std::vector<Plane> someBasePlanes, const std::vector<int> someSampledPlanesIndices, const int aSampleThreshold)
		{
			if (someSampledPlanesIndices.size() <= aSampleThreshold)
			{
				// Skip sampling and compute LE and Triangles
			}
			else
			{
				// Sample someSampledPlanesIndices and compute LE and triangles
			}


		}
	}

	std::vector<Vertex> ComputeLowerEnvelopeSmart(const std::vector<Plane>& somePlanes)
	{
		// If the input size is within a maximum threshold 
		// then use the brute force algorithm
		constexpr auto threshold = 10;
		if (somePlanes.size() <= threshold)
		{
			return ComputeLowerEnvelope(somePlanes);
		}

		// Otherwise use the smart algorithm

		std::vector<Vertex> result;




		return result;
	}

	TriangulationData TriangulateLowerEnvelopeFaces(const std::vector<Face>& someFaces, const bool aComputeAllEdges)
	{
		CGAL_precondition(someFaces.size() > 0);
		TriangulationData result;
		for (int i = 0; i < someFaces.size(); ++i)
		{
			const auto& vertexIndices = someFaces[i].myVertexIndices;
			CGAL_precondition(vertexIndices.size() > 0);
			const auto sourceVertexIdx = vertexIndices.front();
			// The triangulation is guaranteed to be CCW because the first and second vertices always bound the face to their left
			// So every triangulation will pick the vertices in CCW order by design
			for(int k = 1; k < vertexIndices.size() - 1; ++k)
			{
				result.myVerticesIndices.push_back(sourceVertexIdx);
				result.myVerticesIndices.push_back(vertexIndices[k]);
				result.myVerticesIndices.push_back(vertexIndices[k + 1]);
				const auto triangleIndex = result.myVerticesIndices.size() / 3;
				locUpdateTrianglesLocalityMap(sourceVertexIdx, vertexIndices[k], triangleIndex, result.myTrianglesIndicesAtEdge);
				locUpdateTrianglesLocalityMap(vertexIndices[k], vertexIndices[k + 1], triangleIndex, result.myTrianglesIndicesAtEdge);
				locUpdateTrianglesLocalityMap(vertexIndices[k + 1], sourceVertexIdx, triangleIndex, result.myTrianglesIndicesAtEdge);
			}
		}
		return result;
	}
}