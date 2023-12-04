#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>

#include "../common/BatchPointLocation.h"
#include "../common/BruteForce.h"
#include "../common/SmartLowerEnvelope.h"
#include "../common/DebugUtils.h"
#include "../common/Visualization.h"
#include "../common/Utils.h"

#include <iostream>
#include <chrono>

// mem test
#include <CGAL/Compact_container.h>

// thread test
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>

//#include <oneapi/tbb/tbbmalloc_proxy.h> // It improves deallocation


#define THREADS_COUNT 14

TEST_CASE("BatchPointLocation with some random planes", "[BatchPointLocation]")
{
	//std::cout << tbb::TBB_malloc_replacement_log(nullptr) << std::endl;

	// Limit the number of threads to two for all oneTBB parallel interfaces
	oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, THREADS_COUNT);

	/*
	SECTION("100 Random Points")
	{
		constexpr auto sampleCount = 100;
		std::vector<SPGMT::Plane> planes;
		std::vector<SPGMT::Point3> points = SPGMT::Debug::Uniform3DCubeSampling(100.f, sampleCount);

		const auto result = SPGMT::BatchPointLocation(planes, points);

		REQUIRE((result.myRangeWrappers.empty() && result.mySortedPlanesCache.empty()));
	}

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
		REQUIRE((result.mySortedPlanesCache.size() == 1 && result.mySortedPlanesCache.at(0).at(0).size() == planeSamplesCount));

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
					const int planeIdx = result.mySortedPlanesCache.at(0).at(0)[k];
					const auto requirement = planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}

				for (int k = firstPlaneAboveIdx; k < rangeEnd; ++k)
				{
					const int planeIdx = result.mySortedPlanesCache.at(0).at(0)[k];
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
			REQUIRE((result.mySortedPlanesCache.size() == 1 && result.mySortedPlanesCache.at(0).size() == planeSamplesCount));

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

	SECTION("One horizontal plane and some random parallel planes intersecting it")
	{
		constexpr auto minPlaneDistance = 20.f;
		constexpr auto planesCount = 30;
		constexpr auto pointSamplesCount = 10;
		constexpr auto specialPointSamplesCount = 10;
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
			const auto& sortedPlanesIndices = result.mySortedPlanesCache[rangeWrapper.myCacheIndex].at(rangeWrapper.myIndex);

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

	
	SECTION("One horizontal plane and some random parallel planes intersecting it")
	{
		constexpr auto minPlaneDistance = 20.f;
		constexpr auto planesCount = 30;
		constexpr auto pointSamplesCount = 0;
		constexpr auto specialPointSamplesCount = 100;
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

		// NEW CHECK
		using namespace SPGMT;
		REQUIRE(result.myRangeWrappers.size() == points.size());

		// Check ranges
		for (auto i = 0; i < result.myRangeWrappers.size(); ++i)
		{
			const auto& zoneRange = result.myRangeWrappers[i];
			const auto& sortedPlanesIndices =
				result.mySortedPlanesCache[zoneRange.myCacheIndex].at(zoneRange.myIndex);

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
	*/

	SECTION("random planes")
	{
		constexpr auto minPlaneDistance = 20.f;
		constexpr auto planesCount = 100;
		constexpr auto pointSamplesCount = 100;
		constexpr auto allowSamplesOverPlane = true;

		auto planes = SPGMT::Debug::RandomPlaneSampling(planesCount, -400, 400);

		std::vector<SPGMT::Point3> points;

		const auto specialPoints = SPGMT::Debug::SamplePointsAlongPlaneIntersections(planes, 100);
		std::copy(specialPoints.begin(), specialPoints.end(), std::back_inserter(points));
		const auto specialPointsVertices = SPGMT::Debug::SampleTriplePlaneIntersectionPoints(planes, 100);
		std::copy(specialPointsVertices.begin(), specialPointsVertices.end(), std::back_inserter(points));
		CGAL::Random r{ };

		if (planes.size() > 0)
		{
			while (points.size() < 10000)
			{
				const auto i = r.get_int(0, planes.size());
				const auto& planePoints = SPGMT::Debug::RandomPointsPartitionedByPlane(
					pointSamplesCount, planes[i], minPlaneDistance - 1.f, allowSamplesOverPlane);
				std::copy(planePoints.mySamples.begin(), planePoints.mySamples.end(), std::back_inserter(points));
			}
		}

		using std::chrono::high_resolution_clock;
		using std::chrono::duration_cast;
		using std::chrono::duration;
		using std::chrono::milliseconds;

		auto t1 = high_resolution_clock::now();
		std::cout << "starting NOW!" << std::endl;
		// run function
		const auto result = SPGMT::BatchPointLocation(planes, points);

		auto t2 = high_resolution_clock::now();

		auto ms_int = duration_cast<milliseconds>(t2 - t1);

		std::cout << "Elapsed time: " << ms_int.count() << std::endl;

		REQUIRE(result.myRangeWrappers.size() == points.size());

		// Check ranges
		for (auto i = 0; i < result.myRangeWrappers.size(); ++i)
		{
			const auto& zoneRange = result.myRangeWrappers[i];
			const auto& sortedPlanesIndices =
				result.mySortedPlanesCache.at(zoneRange.myCacheIndex).at(zoneRange.myIndex);

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

/*

TEST_CASE("Benchmark ComputeLowerEnvelope-BruteForce", "[Benchmark-PCLE]")
{
	constexpr auto benchmarksCount = 6;
	constexpr auto inputCount = 10;
	for (auto i = 0; i < benchmarksCount; ++i)
	{
		std::string title = "PCLE-BruteForce (" + std::to_string(160 + i * inputCount) + ")";
		const auto& planes = SPGMT::Debug::RandomPlaneSampling(160 + i * inputCount, -300, 300);
		BENCHMARK(title.c_str()) {
			return SPGMT::ComputeLowerEnvelope<SPGMT::ExecutionPolicy::PAR_UNSEQ>(planes);
		};
	}
}

TEST_CASE("Benchmark ComputeLowerEnvelope-BruteForce", "[Benchmark-CLE]")
{
	constexpr auto benchmarksCount = 15;
	constexpr auto inputCount = 10;
	for (auto i = 0; i <= benchmarksCount; ++i)
	{
		std::string title = "CLE-BruteForce (" + std::to_string(i * inputCount) + ")";
		const auto& planes = SPGMT::Debug::RandomPlaneSampling(i * inputCount, -300, 300);
		BENCHMARK(title.c_str()) {
			return SPGMT::ComputeLowerEnvelope<SPGMT::ExecutionPolicy::SEQ>(planes);
		};
	}
}
*/
TEST_CASE("Benchmark BatchPointLocation", "[Benchmark-BPL]")
{
	constexpr auto minPlaneDistance = 20.f;
	constexpr auto pointSamplesCount = 200;
	constexpr auto allowSamplesOverPlane = true;
	constexpr auto pointsCount = 10000;
	constexpr auto benchmarksCount = 13;
	constexpr auto inputCount = 10;

	// Planes
	const auto& allPlanes = SPGMT::Debug::RandomPlaneSampling(benchmarksCount * inputCount);

	// Points
	std::vector<SPGMT::Point3> allPoints;
	const auto specialPoints = SPGMT::Debug::SamplePointsAlongPlaneIntersections(allPlanes, 1000);
	std::copy(specialPoints.begin(), specialPoints.end(), std::back_inserter(allPoints));
	const auto specialPointsVertices = SPGMT::Debug::SampleTriplePlaneIntersectionPoints(allPlanes, 1000);
	std::copy(specialPointsVertices.begin(), specialPointsVertices.end(), std::back_inserter(allPoints));

	if (allPlanes.size() > 0)
	{
		CGAL::Random r{ 0 };
		while (allPoints.size() < pointsCount)
		{
			const auto rand = r.get_int(0, allPlanes.size());
			const auto& planePoints = SPGMT::Debug::RandomPointsPartitionedByPlane(
				pointSamplesCount, allPlanes[rand], minPlaneDistance - 1.f, allowSamplesOverPlane);
			std::copy(planePoints.mySamples.begin(), planePoints.mySamples.end(), std::back_inserter(allPoints));
		}
	}

	for (auto cores = 1; cores < 15; ++cores)
	{
		// Limit the number of threads to two for all oneTBB parallel interfaces
		oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, cores);

		for (auto i = 0; i <= benchmarksCount; ++i)
		{
			const auto planesCount = i * inputCount;
			std::vector<SPGMT::Plane> planes(allPlanes.begin(), allPlanes.begin() + planesCount);
			const std::string title = "BPL (" + std::to_string(cores) + " cores, " + std::to_string(planesCount) + " planes, " + std::to_string(pointsCount) + " points)";

			BENCHMARK(title.c_str()) {
				return SPGMT::BatchPointLocation(planes, allPoints);
			};
		}
	}
}


TEST_CASE("ComputeLowerEnvelope with some parallel planes", "[ComputeLowerEnvelope]")
{
	//std::cout << tbb::TBB_malloc_replacement_log(nullptr) << std::endl;

	// Limit the number of threads to two for all oneTBB parallel interfaces
	oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, THREADS_COUNT);

/*
	SECTION("100 random, non-vertical parallel planes")
	{
		constexpr auto planesCount = 100;
		const auto& planes = SPGMT::Debug::RandomParallelPlanesSampling(planesCount);

		const auto result = SPGMT::ComputeLowerEnvelope(planes);

		REQUIRE(std::holds_alternative<size_t>(result));
		REQUIRE(SPGMT::Debug::IsLowerEnvelopeCorrect(result, planes));

		//SPGMT::Visualization::VisualizeLowerEnvelope(result);
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

		REQUIRE(std::holds_alternative<std::vector<Edge<Point3>>>(result));
		REQUIRE(SPGMT::Debug::IsLowerEnvelopeCorrect(result, planes));

		//SPGMT::Visualization::VisualizeLowerEnvelope(result);
	}
*/
	
	SECTION("20 random dual planes")
	{
		using namespace SPGMT;
		constexpr auto planesCount = 25;

		auto planes = Debug::RandomPlaneSamplingTest(planesCount);
		Utils::FlipPlaneNormalsIfFacingDownwards(planes);
		//const auto& planes = SPGMT::Debug::RandomPlaneSampling(planesCount, -300, 300);
		auto result = SPGMT::ComputeLowerEnvelope(planes);

		//Debug::PrintLowerEnvelope(result);

		REQUIRE(Debug::IsLowerEnvelopeCorrect(result, planes));

		TriangulateLowerEnvelope(result);

		constexpr auto showTriangleEdges = true;
		SPGMT::Visualization::VisualizeLowerEnvelope(result, showTriangleEdges);
	}
}

TEST_CASE("Testing duality map properties Point->Plane, Plane->Point", "[DualityMap]")
{
	SECTION("A point above a plane should become a plane below a point")
	{
		using namespace SPGMT;

		constexpr auto itemsCount = 1000;
		constexpr auto distance = 10.f;

		const auto& startingPlanes = Debug::RandomPlaneSamplingTest(itemsCount);
		std::vector<Point3> startingPoints;
		startingPoints.reserve(itemsCount);

		for (const auto& plane : startingPlanes)
		{
			startingPoints.emplace_back(plane.point() + plane.orthogonal_vector() * distance);
			REQUIRE(plane.has_on_positive_side(startingPoints.back()));
		}

		// Apply duality transform
		const auto& dualizedPlanes = Utils::DualMapping(startingPlanes);
		auto dualizedPoints = Utils::DualMapping(startingPoints);

		// Check reverse order property
		for (int i = 0; i < itemsCount; ++i)
		{
			// Reversed check because each plane normal is pointing downwards (-z)
			REQUIRE(!dualizedPoints[i].has_on_negative_side(dualizedPlanes[i]));
		}
	}
}


TEST_CASE("ComputeSmartLowerEnvelope with some random planes", "[ComputeSmartLowerEnvelope]")
{
	SECTION("20 random dual planes")
	{
		using namespace SPGMT;
		constexpr auto planesCount = 25;

		const auto& planes = Debug::RandomPlaneSamplingTest(planesCount);
		REQUIRE(Utils::AreItemsUnique(planes));

		auto result = SPGMT::ComputeLowerEnvelopeSmart(planes);

		REQUIRE(Debug::IsLowerEnvelopeCorrect(result, planes));

		//SPGMT::TriangulateLowerEnvelope(result);

		//SPGMT::Visualization::VisualizeLowerEnvelope(result);
	}
}