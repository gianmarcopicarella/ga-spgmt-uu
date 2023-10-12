#include "SmartLowerEnvelope.h"

#include "BatchPointLocation.h"
#include "Utils.h"

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
			// Only two triangles can share the same edge
			CGAL_postcondition(firstEdgeEntry->second.size() < 3);
		}

		LowerEnvelope3d locRecursiveLowerEnvelope(
			const std::vector<Plane>& somePlanes, 
			const std::vector<int>& somePlanesFromBase,
			const std::vector<int>& someSampledPlanesFromBase)
		{
			constexpr auto exponent = 1.f / 8.f;
			const int samplingSetSize = std::pow(somePlanesFromBase.size(), exponent);
			
			std::vector<int> currentSampledPlanesFromBase;
			if (someSampledPlanesFromBase.size() <= samplingSetSize)
			{
				// Skip sampling and compute LE and Triangles
				currentSampledPlanesFromBase = someSampledPlanesFromBase;
			}
			else
			{
				// Sample someSampledPlanesIndices and compute LE and triangles
				const auto samplingProbability = (float)samplingSetSize / (float)somePlanesFromBase.size();
				currentSampledPlanesFromBase = SampleWithProbability(someSampledPlanesFromBase, samplingProbability);
			}

			// Compute lower envelope and triangles
			/*std::vector<Plane> planes;
			std::transform(currentSampledPlanesFromBase.begin(), currentSampledPlanesFromBase.end(),
				std::back_inserter(planes), [&somePlanes](const auto& aPlaneIdx) { return somePlanes[aPlaneIdx]; });
			const auto& lowerEnvelope = ComputeLowerEnvelope(planes);
			
			if (lowerEnvelope.size() == 1)
			{
				CGAL_precondition(lowerEnvelope.front().myType == VertexType::INFINITE);
				return lowerEnvelope;
			}

			//const auto& lowerEnvelopeFaces = ExtractLowerEnvelopeFaces(lowerEnvelope);
			
			// Prepare query points for batch point location 
			{
				const auto& pointToPlaneDuality = [](const auto& aVertex) { return Utils::SingleDualMapping(aVertex.myPoint); };
				const auto& planeToPointDuality = [somePlanes](const auto& aPlaneIdx) { 
					return Utils::SingleDualMapping(somePlanes[aPlaneIdx]); };
				std::vector<Plane> dualizedPoints;
				std::vector<Point3> dualizedPlanes;

				std::transform(lowerEnvelope.begin(), lowerEnvelope.end(), std::back_inserter(dualizedPoints), pointToPlaneDuality);
				std::transform(somePlanesFromBase.begin(), somePlanesFromBase.end(), std::back_inserter(dualizedPlanes), planeToPointDuality);

				const auto& locationResult = BatchPointLocation(dualizedPoints, dualizedPlanes);
			}*/

			// TEMPORARY RETURN
			return LowerEnvelope3d{};


			// 1. Edge case) If there is just one line with both ends at infinity
			// else triangulation can proceed

			//const auto& lowerEnvelopeTriangulation = TriangulateLowerEnvelopeFaces(lowerEnvelopeFaces);

			// Prepare input for BatchPointLocation
			// Execute BatchPointLocation
			// Collect conflict list per vertex

			//for (int i = 0; i < lowerEnvelopeTriangulation.myVerticesIndices.size(); i += 3)
			//{
			//	// Collect recursive sets H' and S'
			//	// If S' is not empty then call recursion
			//}

			// Assemble each triangle result together and return the result
		}
	}

	LowerEnvelope3d ComputeLowerEnvelopeSmart(const std::vector<Plane>& somePlanes)
	{
		// If the input size is within a maximum threshold 
		// then use the brute force algorithm
		constexpr auto threshold = 10;
		if (somePlanes.size() <= threshold)
		{
			return ComputeLowerEnvelope(somePlanes);
		}

		// Otherwise use the smart algorithm
		std::vector<int> allPlanesIndices(somePlanes.size());
		std::iota(allPlanesIndices.begin(), allPlanesIndices.end(), 0);
		
		const auto samplingProbability = std::pow(somePlanes.size(), 1.f / 9.f) / (float)somePlanes.size();
		const auto& sampledPlanesIndices = SampleWithProbability(allPlanesIndices, samplingProbability);

		return locRecursiveLowerEnvelope(somePlanes, allPlanesIndices, sampledPlanesIndices);
	}
}