#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "../common/BatchPointLocation.h"
#include "../common/DebugUtils.h"

/*
Currently this function must behave correctly for the following cases:
1) No Plane intersections found -> Parallel planes (vertical/horizontal cases too) DONE
2) No Plane intersections found -> Only one plane. DONE
*/

/*
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
		constexpr auto planeSamplesCount = 1000;
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
				for (int k = 0; k < firstPlaneAboveIdx - 1; ++k)
				{
					const int planeIdx = result.mySortedPlanesIndices.at(0)[k];
					const auto requirement = planes[planeIdx].has_on_positive_side(points[i]);
					REQUIRE(requirement);
				}

				for (int k = firstPlaneAboveIdx; k < rangeEnd; ++k)
				{
					const int planeIdx = result.mySortedPlanesIndices.at(0)[k];
					const auto requirement = planes[planeIdx].has_on_negative_side(points[i]);
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
        constexpr auto planesCount = 68;
        constexpr auto pointSamplesCount = 0;
        constexpr auto specialPointSamplesCount = 10000;
        constexpr auto allowSamplesOverPlane = false;

        auto planes = SPGMT::Debug::RandomParallelPlanesSampling(planesCount, minPlaneDistance);
        {
            SPGMT::Plane horizontalPlane{ SPGMT::Point3{0,0,0}, SPGMT::Dir3{0,1,0} };
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
                for (int k = 0; k < firstPlaneAboveIdx - 1; ++k)
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
*/
TEST_CASE("BatchPointLocation with some random planes", "[BatchPointLocation]")
{
    SECTION("random planes")
    {
        constexpr auto minPlaneDistance = 3.f;
        constexpr auto planesCount = 25;
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

        // run function
        const auto result = SPGMT::BatchPointLocation(planes, points);

        REQUIRE(result.myRangeWrappers.size() == points.size());

        // Check ranges
        for (auto i = 0; i < result.myRangeWrappers.size(); ++i)
        {
            const auto& zoneRange = result.myRangeWrappers[i];

            // TEMPORARY SKIP
            if(zoneRange.myRefIndex == -1)
            {
                continue;
            }

            const auto& sortedPlanesIndices =
                    result.mySortedPlanesIndices.at(zoneRange.myRefIndex);

            // This is the index related to the sorted list of planes in unpackedResult, NOT the list "planes"
            const int firstPlaneAboveIdx = zoneRange.myRange.first;
            const int rangeEnd = zoneRange.myRange.second;

            if (firstPlaneAboveIdx != planes.size())
            {
                // Check that all planes before are below the point and all planes in the range are above the point
                for (int k = 0; k < firstPlaneAboveIdx - 1; ++k)
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