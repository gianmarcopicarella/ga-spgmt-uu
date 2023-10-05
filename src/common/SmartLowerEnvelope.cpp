#include "SmartLowerEnvelope.h"

#include "BruteForce.h"

namespace SPGMT
{
	namespace
	{
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

	std::vector<TriangulationData> TriangulateLowerEnvelopeFaces(const std::vector<Face>& someFaces, const bool aComputeAllEdges)
	{
		CGAL_precondition(someFaces.size() > 0);

		std::vector<TriangulationData> result;

		for (int i = 0; i < someFaces.size(); ++i)
		{
			const auto& vertexIndices = someFaces[i].myVertexIndices;
			CGAL_precondition(vertexIndices.size() > 0);
			const auto sourceVertexIdx = vertexIndices.front();

			// The triangulation is guaranteed to be CCW because the first and second vertices always bound the face to their left
			// So every triangulation will pick the vertices in CCW order by design
			TriangulationData td;
				
			for(int k = 1; k < vertexIndices.size() - 1; ++k)
			{
				td.myVerticesIndices.push_back(sourceVertexIdx);
				
				/*if (aComputeAllEdges)
				{*/
					td.myVerticesIndices.push_back(vertexIndices[k]);
					td.myVerticesIndices.push_back(vertexIndices[k + 1]);
				/*}
				else
				{
					td.myVerticesIndices.push_back(vertexIndices[k + 1]);
				}*/
			}

			/*if (!aComputeAllEdges && someFaces[i].myType == FaceType::BOUNDED)
			{
				td.myVerticesIndices.pop_back();
				td.myVerticesIndices.pop_back();
			}*/

			result.push_back(td);
		}

		return result;
	}
}