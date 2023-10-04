#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "../common/BatchPointLocation.h"
#include "../common/BruteForce.h"
#include "../common/DebugUtils.h"

#include "../common/Serialize.h"
#include "../common/Visualization.h"
/*
Currently this function must behave correctly for the following cases:
1) No Plane intersections found -> Parallel planes (vertical/horizontal cases too) DONE
2) No Plane intersections found -> Only one plane. DONE
*/

TEST_CASE("BatchPointLocation with no planes returns an empty list of ranges", "[BatchPointLocation]")
{
	SECTION("100 Random Points")
	{
		constexpr auto sampleCount = 100;
		std::vector<SPGMT::Plane> planes;
		std::vector<SPGMT::Point3> points = SPGMT::Debug::Uniform3DCubeSampling(100.f, sampleCount);

		const auto result = SPGMT::BatchPointLocation(planes, points);

		REQUIRE((result.myRangeWrappers.empty() && result.mySortedPlanesIndices.empty()));
	}
}

TEST_CASE("BatchPointLocation with one plane returns a list with pair <0, -1> if point is above plane, else <-1, -1>", "[BatchPointLocation]")
{
	// Parallel Random planes
	SECTION("Multiple Parallel Planes and Random Points around them")
	{
		constexpr auto planeSamplesCount = 200;
		constexpr auto pointSamplesCount = 300;
		constexpr auto minPlaneDistance = 20.f;

		auto planes = SPGMT::Debug::RandomParallelPlanesSampling(planeSamplesCount, minPlaneDistance);
		std::vector<SPGMT::Point3> points;

		for (auto i = 0; i < planes.size(); ++i)
		{
			const auto& planePoints = SPGMT::Debug::RandomPointsPartitionedByPlane(pointSamplesCount, planes[i], minPlaneDistance - 1.f, false);
			std::copy(planePoints.mySamples.begin(), planePoints.mySamples.end(), std::back_inserter(points));
		}

		// run function
		const auto result = SPGMT::BatchPointLocation(planes, points);

		REQUIRE(result.myRangeWrappers.size() == (planeSamplesCount * pointSamplesCount));
		REQUIRE((result.mySortedPlanesIndices.size() == 1 && result.mySortedPlanesIndices.at(0).size() == planeSamplesCount));

		// Check ranges
		for (auto i = 0; i < result.myRangeWrappers.size(); ++i)
		{
			// This is the index related to the sorted list of planes in unpackedResult, NOT the list "planes"
			const int firstPlaneAboveIdx = result.myRangeWrappers[i].myRange.first;
			const int rangeEnd = result.myRangeWrappers[i].myRange.second;

			if (firstPlaneAboveIdx != -1)
			{
				// Check that all planes before are below the point and all planes in the range are above the point
				for (int k = 0; k < firstPlaneAboveIdx; ++k)
				{
					const int planeIdx = result.mySortedPlanesIndices.at(0)[k];
					const auto requirement = planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}

				for (int k = firstPlaneAboveIdx; k < rangeEnd; ++k)
				{
					const int planeIdx = result.mySortedPlanesIndices.at(0)[k];
					const auto requirement = !planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}
			}
			else
			{
				// Check that all planes are below the point (here order doesn't matter)
				for (int k = 0; k < planes.size(); ++k)
				{
					const auto requirement = planes[k].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}
			}
		}
	}

	// Random planes
	SECTION("Single Random Plane and Random Points around it")
	{
		constexpr auto planeSamplesCount = 1;
		constexpr auto pointSamplesCount = 10000;

		for (const auto& plane : SPGMT::Debug::RandomPlaneSampling(planeSamplesCount))
		{
			// Generate 1000 points partitioned by the plane
			const auto& points = SPGMT::Debug::RandomPointsPartitionedByPlane(pointSamplesCount, plane);

			// run function
			const auto result = SPGMT::BatchPointLocation({ plane }, points.mySamples);
			REQUIRE(result.myRangeWrappers.size() == pointSamplesCount);
			REQUIRE((result.mySortedPlanesIndices.size() == 1 && result.mySortedPlanesIndices.at(0).size() == planeSamplesCount));

			// Check result
			for (auto i = 0; i < pointSamplesCount; ++i)
			{
				if (result.myRangeWrappers[i].myRange.first == 0)
				{
					const auto requirement = !plane.has_on_positive_side(points.mySamples[i]);
					REQUIRE(requirement);
				}
			}
		}
	}
}

TEST_CASE("BatchPointLocation with multiple parallel 2D lines when projecting plane intersections", "[BatchPointLocation]")
{
	SECTION("One horizontal plane and some random parallel planes intersecting it")
	{
		constexpr auto minPlaneDistance = 20.f;
		constexpr auto planesCount = 30;
		constexpr auto pointSamplesCount = 0;
		constexpr auto specialPointSamplesCount = 10000;
		constexpr auto allowSamplesOverPlane = true;

		auto planes = SPGMT::Debug::RandomParallelPlanesSampling(planesCount, minPlaneDistance);
		{
			SPGMT::Plane horizontalPlane{ SPGMT::Point3{0,0,0}, SPGMT::Dir3{0,0,1} };
			planes.push_back(horizontalPlane);
		}

		std::vector<SPGMT::Point3> points;

		const auto specialPoints = SPGMT::Debug::SamplePointsAlongPlaneIntersections(planes, specialPointSamplesCount * planesCount);
		std::copy(specialPoints.begin(), specialPoints.end(), std::back_inserter(points));

		for (auto i = 0; i < planes.size(); ++i)
		{
			const auto& planePoints = SPGMT::Debug::RandomPointsPartitionedByPlane(
				pointSamplesCount, planes[i], minPlaneDistance - 1.f, allowSamplesOverPlane);
			std::copy(planePoints.mySamples.begin(), planePoints.mySamples.end(), std::back_inserter(points));
		}

		// run function
		const auto result = SPGMT::BatchPointLocation(planes, points);

		REQUIRE(result.myRangeWrappers.size() == points.size());

		// Check ranges
		for (auto i = 0; i < result.myRangeWrappers.size(); ++i)
		{
			const auto& rangeWrapper = result.myRangeWrappers[i];
			const auto& sortedPlanesIndices = result.mySortedPlanesIndices.at(rangeWrapper.myRefIndex);

			// This is the index related to the sorted list of planes in unpackedResult, NOT the list "planes"
			const int firstPlaneAboveIdx = rangeWrapper.myRange.first;
			const int rangeEnd = rangeWrapper.myRange.second;

			if (firstPlaneAboveIdx != planes.size())
			{
				// Check that all planes before are below the point and all planes in the range are above the point
				for (int k = 0; k < firstPlaneAboveIdx; ++k)
				{
					const int planeIdx = sortedPlanesIndices[k];
					const auto requirement = planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}

				for (int k = firstPlaneAboveIdx; k < rangeEnd; ++k)
				{
					const int planeIdx = sortedPlanesIndices[k];
					const auto requirement = !planes[planeIdx].has_on_positive_side(points[i]); // negative or on the plane
					REQUIRE(requirement);
				}
			}
			else
			{
				// Check that all planes are below the point (here order and zone dont matter)
				for (int k = 0; k < planes.size(); ++k)
				{
					const auto requirement = planes[k].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}
			}
		}
	}
}

TEST_CASE("BatchPointLocation with some random planes", "[BatchPointLocation]")
{
	SECTION("random planes")
	{
		constexpr auto minPlaneDistance = 20.f;
		constexpr auto planesCount = 30;
		constexpr auto pointSamplesCount = 10000;
		constexpr auto allowSamplesOverPlane = true;

		auto planes = SPGMT::Debug::RandomPlaneSampling(planesCount);

		std::vector<SPGMT::Point3> points;

		for (auto i = 0; i < planes.size(); ++i)
		{
			const auto& planePoints = SPGMT::Debug::RandomPointsPartitionedByPlane(
				pointSamplesCount, planes[i], minPlaneDistance - 1.f, allowSamplesOverPlane);
			std::copy(planePoints.mySamples.begin(), planePoints.mySamples.end(), std::back_inserter(points));
		}

		const auto specialPoints = SPGMT::Debug::SamplePointsAlongPlaneIntersections(planes, 1000);
		std::copy(specialPoints.begin(), specialPoints.end(), std::back_inserter(points));

		const auto specialPointsVertices = SPGMT::Debug::SampleTriplePlaneIntersectionPoints(planes, 1000);
		std::copy(specialPointsVertices.begin(), specialPointsVertices.end(), std::back_inserter(points));

		// run function
		const auto result = SPGMT::BatchPointLocation(planes, points);

		REQUIRE(result.myRangeWrappers.size() == points.size());

		// Check ranges
		for (auto i = 0; i < result.myRangeWrappers.size(); ++i)
		{
			const auto& zoneRange = result.myRangeWrappers[i];
			const auto& sortedPlanesIndices =
				result.mySortedPlanesIndices.at(zoneRange.myRefIndex);

			// This is the index related to the sorted list of planes in unpackedResult, NOT the list "planes"
			const int firstPlaneAboveIdx = zoneRange.myRange.first;
			const int rangeEnd = zoneRange.myRange.second;

			if (firstPlaneAboveIdx != planes.size())
			{
				// Check that all planes before are below the point and all planes in the range are above the point
				for (int k = 0; k < firstPlaneAboveIdx; ++k)
				{
					const int planeIdx = sortedPlanesIndices[k];
					const auto requirement = planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}

				for (int k = firstPlaneAboveIdx; k < rangeEnd; ++k)
				{
					const int planeIdx = sortedPlanesIndices[k];
					const auto requirement = !planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}
			}
			else
			{
				// Check that all planes are below the point (here order and zone dont matter)
				for (int k = 0; k < planes.size(); ++k)
				{
					const auto requirement = planes[k].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}
			}
		}

	}
}

TEST_CASE("ComputeLowerEnvelope with some parallel planes", "[ComputeLowerEnvelope]")
{
	SECTION("100 random, non-vertical parallel planes")
	{
		constexpr auto planesCount = 100;
		const auto& planes = SPGMT::Debug::RandomParallelPlanesSampling(planesCount);

		const auto result = SPGMT::ComputeLowerEnvelope(planes);

		REQUIRE((
			result.size() == 1 &&
			result[0].myType == SPGMT::VertexType::INFINITE &&
			result[0].myLowestLeftPlanes.size() == 1));

		const auto resultPlaneIdx = result[0].myLowestLeftPlanes[0];
		for (int i = 0; i < planes.size(); ++i)
		{
			const auto requirement = !planes[i].has_on_positive_side(planes[resultPlaneIdx].point());
			REQUIRE(requirement);
		}
	}

	SECTION("100 random, non-vertical parallel planes and one horizontal plane")
	{
		using namespace SPGMT;

		constexpr auto planesCount = 100;
		auto planes = SPGMT::Debug::RandomParallelPlanesSampling(planesCount);
		{
			SPGMT::Plane horizontalPlane{ SPGMT::Point3{0,0,0}, SPGMT::Dir3{0,0,1} };
			planes.push_back(horizontalPlane);
		}

		const auto result = SPGMT::ComputeLowerEnvelope(planes);

		REQUIRE((result.size() == 2 &&
			result[0].myType == SPGMT::VertexType::INFINITE &&
			result[1].myType == SPGMT::VertexType::INFINITE &&
			result[0].mySortedNeighboursIndices.size() == 1 &&
			result[1].mySortedNeighboursIndices.size() == 1));

		const Vec3 up{ 0,0,1 };

		for (int i = 0; i < result.size(); ++i)
		{
			const Vec3 lineSegment{ result[i].myPoint, result[result[i].mySortedNeighboursIndices[0]].myPoint };
			const Vec3 left = -CGAL::cross_product(lineSegment, up);
			const Point3 sample = result[i].myPoint + left * 100.f;
			const Line3 rayLine{ sample, up };

			std::vector<std::pair<FT, int>> zPlanes;

			for (int j = 0; j < planes.size(); ++j)
			{
				const auto intersection = CGAL::intersection(rayLine, planes[j]);
				CGAL_precondition(intersection.has_value());
				zPlanes.push_back(std::make_pair(boost::get<Point3>(&*intersection)->z(), j));
			}

			const auto lowestPlaneIdx = result[i].myLowestLeftPlanes[0];

			std::sort(zPlanes.begin(), zPlanes.end(), [](auto& a, auto& b) {return a.first < b.first; });
			REQUIRE(zPlanes[0].second == lowestPlaneIdx);
		}

		//Serialization::SerializeLowerEnvelope("single_line_LE.ply", result);
	}
}

SPGMT::Point2 locProjectXY(const SPGMT::Point3& aPoint)
{
	return SPGMT::Point2{ aPoint.x(), aPoint.y() };
}
int locFindVertex(const std::vector<SPGMT::Vertex>& someVertices, const SPGMT::Point2& aPoint, const SPGMT::VertexType aVertexType)
{
	const auto matchIt = std::find_if(someVertices.begin(), someVertices.end(), [&aPoint, aVertexType](const auto& aVertex) {
		return aVertex.myType == aVertexType && locProjectXY(aVertex.myPoint) == aPoint;
		});

	if (matchIt != someVertices.end())
	{
		return std::distance(someVertices.begin(), matchIt);
	}

	return -1;
}
int locFindVertex(const std::vector<SPGMT::Vertex>& someVertices, const std::vector<int>& someIndices, const SPGMT::Point2& aPoint, const SPGMT::VertexType aVertexType)
{
	for (int i = 0; i < someIndices.size(); ++i)
	{
		if (someVertices[someIndices[i]].myType == aVertexType)
		{
			const auto& projection = locProjectXY(someVertices[someIndices[i]].myPoint);
			if (aPoint == projection)
			{
				return someIndices[i];
			}
		}
	}

	return -1;
}
int locFindVertexToInfinity(const std::vector<SPGMT::Vertex>& someVertices, const std::vector<int>& someIndices, const SPGMT::Ray2& aRay)
{
	for (int i = 0; i < someIndices.size(); ++i)
	{
		if (someVertices[someIndices[i]].myType == SPGMT::VertexType::INFINITE)
		{
			const auto& projection = locProjectXY(someVertices[someIndices[i]].myPoint);
			if (aRay.has_on(projection))
			{
				return someIndices[i];
			}
		}
	}

	return -1;
}
void locIsLowerEnvelopeProjectionCorrect(const std::vector<SPGMT::Vertex>& someVertices, const std::vector<SPGMT::Plane>& somePlanes)
{
	using VertexIter = Envelope_diagram_2::Vertex_const_iterator;
	using HalfEdgeCwIter = typename Envelope_diagram_2::Halfedge_around_vertex_const_circulator;
	const auto& expectedResult = SPGMT::Debug::GetLowerEnvelopeOfPlanes(somePlanes);
	CGAL_precondition(expectedResult.is_valid());
	std::set<int> usedVerticesIndices;

	for (VertexIter it = expectedResult.vertices_begin(); it != expectedResult.vertices_end(); ++it)
	{
		const auto currentIndex = locFindVertex(someVertices, it->point(), SPGMT::VertexType::FINITE);
		REQUIRE(currentIndex != -1);
		REQUIRE(usedVerticesIndices.find(currentIndex) == usedVerticesIndices.end());
		usedVerticesIndices.insert(currentIndex);

		const auto& currentVertexXY = locProjectXY(someVertices[currentIndex].myPoint);
		HalfEdgeCwIter neighbourIt = it->incident_halfedges();
		std::set<int> usedNeighbourVerticesIndices;

		do
		{
			const auto& halfEdgeCurve = neighbourIt->ccb()->curve();
			if (halfEdgeCurve.is_segment())
			{
				const auto& segment = halfEdgeCurve.segment();
				const auto& pointToFind = currentVertexXY == segment.source() ? segment.target() : segment.source();
				const auto targetIdx = locFindVertex(someVertices, someVertices[currentIndex].mySortedNeighboursIndices,
					pointToFind, SPGMT::VertexType::FINITE);
				REQUIRE(targetIdx != -1);
				REQUIRE(usedNeighbourVerticesIndices.find(targetIdx) == usedNeighbourVerticesIndices.end());
				usedNeighbourVerticesIndices.insert(targetIdx);

			}
			else if (halfEdgeCurve.is_ray())
			{
				const auto& ray = halfEdgeCurve.ray();
				CGAL_precondition(ray.source() == currentVertexXY);
				const int infinityIdx = locFindVertexToInfinity(someVertices, someVertices[currentIndex].mySortedNeighboursIndices, ray);
				REQUIRE(infinityIdx != -1);
				REQUIRE(usedNeighbourVerticesIndices.find(infinityIdx) == usedNeighbourVerticesIndices.end());
				usedNeighbourVerticesIndices.insert(infinityIdx);
			}
			else
			{
				CGAL_precondition(false);
			}
		} while (++neighbourIt != it->incident_halfedges());

		// All vertex neighbours must be explored
		REQUIRE(usedNeighbourVerticesIndices.size() == someVertices[currentIndex].mySortedNeighboursIndices.size());
	}

	// Check that every finite vertex in my result has been visited
	const auto finiteVerticesCounter = std::count_if(someVertices.begin(), someVertices.end(),
		[](const auto& aVertex) {return aVertex.myType == SPGMT::VertexType::FINITE; });
	REQUIRE(finiteVerticesCounter == usedVerticesIndices.size());

	// Check that the number of faces is correct
	const auto& faces = SPGMT::ExtractLowerEnvelopeFaces(someVertices);
	REQUIRE(faces.size() == expectedResult.number_of_faces());

	// Check that each face index is unique
	{
		std::vector<int> faceIndices;
		std::transform(faces.begin(), faces.end(), std::back_inserter(faceIndices), [](const auto& aFace) {return aFace.myPlaneIndex; });
		const auto initialSize = faceIndices.size();
		faceIndices.erase(std::unique(faceIndices.begin(), faceIndices.end()), faceIndices.end());
		const auto finalSize = faceIndices.size();
		REQUIRE(initialSize == finalSize);
	}
}


TEST_CASE("ComputeLowerEnvelope with some random planes", "[ComputeLowerEnvelope]")
{
	SECTION("10 random dual planes")
	{
		using namespace SPGMT;
		constexpr auto planesCount = 30;
		constexpr auto halfSide = 50.f;
		const auto& points = SPGMT::Debug::Uniform3DCubeSampling(halfSide, planesCount);
		const auto& planes = SPGMT::Debug::GetDualPlanes(points);
		const auto result = SPGMT::ComputeLowerEnvelope(planes);
		locIsLowerEnvelopeProjectionCorrect(result, planes);

		//SPGMT::Visualization::VisualizeLowerEnvelope(result);
		//Serialization::SerializeLowerEnvelope("random_planes_LE.ply", result);
	}
}