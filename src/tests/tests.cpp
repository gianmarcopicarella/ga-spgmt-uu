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

// Including it here otherwise the compiler moans
#include <hpx/hpx_init.hpp>

/**/
TEST_CASE("BatchPointLocation with no planes returns an empty list of ranges", "[BatchPointLocation]")
{
	SECTION("100 Random Points")
	{
		constexpr auto sampleCount = 100;
		std::vector<SPGMT::Plane> planes;
		std::vector<SPGMT::Point3> points = SPGMT::Debug::Uniform3DCubeSampling(100.f, sampleCount);

		const auto result = SPGMT::ParallelBatchPointLocation(planes, points);

		REQUIRE((result.myRangeWrappers.empty() && result.mySortedPlanesCache.empty()));
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
		const auto result = SPGMT::ParallelBatchPointLocation(planes, points);

		REQUIRE(result.myRangeWrappers.size() == (planeSamplesCount * pointSamplesCount));
		REQUIRE((result.mySortedPlanesCache.size() == 1 && result.mySortedPlanesCache[0].at(0).size() == planeSamplesCount));

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
					const int planeIdx = result.mySortedPlanesCache[0].at(0)[k];
					const auto requirement = planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}

				for (int k = firstPlaneAboveIdx; k < rangeEnd; ++k)
				{
					const int planeIdx = result.mySortedPlanesCache[0].at(0)[k];
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
			const auto result = SPGMT::ParallelBatchPointLocation({ plane }, points.mySamples);
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
		const auto result = SPGMT::ParallelBatchPointLocation(planes, points);

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
}

TEST_CASE("BatchPointLocation with some random planes", "[BatchPointLocation]")
{
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


	SECTION("random planes")
	{
		constexpr auto minPlaneDistance = 20.f;
		constexpr auto planesCount = 15;
		constexpr auto pointSamplesCount = 100;
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
		const auto result = SPGMT::ParallelBatchPointLocation(planes, points);

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
}

TEST_CASE("Benchmark ComputeLowerEnvelope-BruteForce", "[Benchmark-PCLE]")
{
	constexpr auto benchmarksCount = 15;
	constexpr auto inputCount = 10;
	for (auto i = 0; i < benchmarksCount; ++i)
	{
		std::string title = "PCLE-BruteForce (" + std::to_string(i * inputCount) + ")";
		const auto& planes = SPGMT::Debug::RandomPlaneSampling(i * inputCount, -300, 300);
		BENCHMARK(title.c_str()) {
			return SPGMT::ParallelComputeLowerEnvelope(planes);
		};
	}
}

TEST_CASE("Benchmark ComputeLowerEnvelope-BruteForce", "[Benchmark-CLE]")
{
	constexpr auto benchmarksCount = 15;
	constexpr auto inputCount = 10;
	for (auto i = 0; i < benchmarksCount; ++i)
	{
		std::string title = "CLE-BruteForce (" + std::to_string(i * inputCount) + ")";
		const auto& planes = SPGMT::Debug::RandomPlaneSampling(i * inputCount, -300, 300);
		BENCHMARK(title.c_str()) {
			return SPGMT::ComputeLowerEnvelope(planes);
		};
	}
}

TEST_CASE("Benchmark BatchPointLocation", "[Benchmark-PBPL]")
{
	constexpr auto minPlaneDistance = 20.f;
	constexpr auto pointSamplesCount = 1000;
	constexpr auto allowSamplesOverPlane = true;

	constexpr auto benchmarksCount = 5;
	constexpr auto inputCount = 10;

	for (auto i = 0; i < benchmarksCount; ++i)
	{
		std::string title = "PBPL (" + std::to_string(i * inputCount) + ")";

		// Compute input data
		auto planes = SPGMT::Debug::RandomPlaneSampling(i * inputCount);
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
		// End input data computation

		BENCHMARK(title.c_str()) {
			return SPGMT::ParallelBatchPointLocation(planes, points);
		};
	}
}

TEST_CASE("Benchmark BatchPointLocation", "[Benchmark-BPL]")
{
	constexpr auto minPlaneDistance = 20.f;
	constexpr auto pointSamplesCount = 1000;
	constexpr auto allowSamplesOverPlane = true;

	constexpr auto benchmarksCount = 5;
	constexpr auto inputCount = 10;

	for (auto i = 0; i < benchmarksCount; ++i)
	{
		std::string title = "BPL (" + std::to_string(i * inputCount) + ")";

		// Compute input data
		auto planes = SPGMT::Debug::RandomPlaneSampling(i * inputCount);
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
		// End input data computation

		BENCHMARK(title.c_str()) {
			return SPGMT::BatchPointLocation(planes, points);
		};
	}
}


TEST_CASE("ComputeLowerEnvelope with some parallel planes", "[ComputeLowerEnvelope]")
{

	SECTION("100 random, non-vertical parallel planes")
	{
		constexpr auto planesCount = 100;
		const auto& planes = SPGMT::Debug::RandomParallelPlanesSampling(planesCount);

		const auto result = SPGMT::ParallelComputeLowerEnvelope(planes);

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

		const auto result = SPGMT::ParallelComputeLowerEnvelope(planes);

		REQUIRE(std::holds_alternative<std::vector<Edge<Point3>>>(result));
		REQUIRE(SPGMT::Debug::IsLowerEnvelopeCorrect(result, planes));

		//SPGMT::Visualization::VisualizeLowerEnvelope(result);
	}

	SECTION("20 random dual planes")
	{
		using namespace SPGMT;
		constexpr auto planesCount = 25;
		const auto& planes = SPGMT::Debug::RandomPlaneSampling(planesCount, -300, 300);
		auto result = SPGMT::ParallelComputeLowerEnvelope(planes);

		REQUIRE(Debug::IsLowerEnvelopeCorrect(result, planes));

		TriangulateLowerEnvelope(result);

		SPGMT::Visualization::VisualizeLowerEnvelope(result);
	}
}

bool locIsPlaneFacingUpwards(const SPGMT::Plane& aPlane)
{
	using namespace SPGMT;
	static const Vec3 up{ 0,0,1 };
	return CGAL::sign(CGAL::scalar_product(up, aPlane.orthogonal_vector())) != CGAL::Sign::NEGATIVE;
}

TEST_CASE("Testing duality map properties Point->Plane, Plane->Point", "[DualityMap]")
{
	SECTION("A point above a plane should become a plane below a point")
	{
		using namespace SPGMT;

		constexpr auto planesCount = 100;

		const auto& planes = SPGMT::Debug::RandomPlaneSampling(planesCount);
		std::vector<Point3> points;

		for (const auto& plane : planes)
		{
			points.emplace_back(plane.point() + plane.orthogonal_vector() * 10.f);
			REQUIRE(plane.has_on_positive_side(points.back()));
		}

		// Apply duality
		const auto& dualPoints = Utils::DualMapping(points);
		const auto& dualPlanes = Utils::DualMapping(planes);

		// Check property
		for (int i = 0; i < planesCount; ++i)
		{
			std::cout << "Plane: " << dualPoints[i] << ", Point: " << dualPlanes[i] << std::endl;

			if (locIsPlaneFacingUpwards(dualPoints[i]))
			{
				REQUIRE(!dualPoints[i].has_on_positive_side(dualPlanes[i]));
			}
			else
			{
				REQUIRE(!dualPoints[i].has_on_negative_side(dualPlanes[i]));
			}


		}

	}
}

TEST_CASE("ComputeSmartLowerEnvelope with some random planes", "[ComputeSmartLowerEnvelope]")
{
	SECTION("20 random dual planes")
	{
		using namespace SPGMT;
		/*constexpr auto planesCount = 20;
		const auto& planes = SPGMT::Debug::RandomPlaneSampling(planesCount);
		auto result = SPGMT::ComputeLowerEnvelopeSmart(planes);

		REQUIRE(Debug::IsLowerEnvelopeCorrect(result, planes));
		*/
		//SPGMT::TriangulateLowerEnvelope(result);

		//SPGMT::Visualization::VisualizeLowerEnvelope(result);
	}
}


int RunTests(int argc, char* argv[])
{
	Catch::ConfigData settings;
	settings.benchmarkSamples = 3;

	std::vector<std::string> testOrTag;

	//testOrTag.push_back("[ComputeLowerEnvelope]");
	testOrTag.push_back("[BatchPointLocation]");//[Benchmark-BPL], [Benchmark-PBPL]");//"[BatchPointLocation]");//"[Benchmark-PCLE], [Benchmark-CLE]");//[BatchPointLocation]");//"[ComputeLowerEnvelope]");//[Benchmark - PCLE], [Benchmark - CLE]"); //,[Benchmark-CLE]

	settings.testsOrTags = testOrTag;

	Catch::Session session;

	session.useConfigData(settings);
	session.run();
	//Catch::Session().run(argc, argv);
	return hpx::finalize();
}

int main(int argc, char* argv[])
{
	return hpx::init(RunTests, argc, argv);
}