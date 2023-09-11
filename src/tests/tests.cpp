#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "../common/BatchPointLocation.h"
#include "../common/DebugUtils.h"

/*
Currently this function must behave correctly for the following cases:
1) No Plane intersections found -> Parallel planes (vertical/horizontal cases too) DONE
2) No Plane intersections found -> Only one plane. DONE
*/


TEST_CASE("BatchPointLocation with no planes returns an empty list of ranges", "[BatchPointLocation]")
{
    SPGMT::LocationResult result;
    std::vector<SPGMT::Plane> planes;

    SECTION("100 Random Points") 
    {
        constexpr auto sampleCount = 100;
        std::vector<SPGMT::Point3> points = SPGMT::Debug::Uniform3DCubeSampling(100.f, sampleCount);
        
        SPGMT::BatchPointLocation(planes, points, result);

        REQUIRE(result.myRanges.size() == 0);
        REQUIRE(result.myType == SPGMT::ResultType::NONE);
    }
}

TEST_CASE("BatchPointLocation with one plane returns a list with pair <0, -1> if point is above plane, else <-1, -1>", "[BatchPointLocation]")
{
    // Parallel Random planes
    SECTION("Multiple Parallel Planes and Random Points around them")
    {
        constexpr auto planeSamplesCount = 100;
        constexpr auto pointSamplesCount = 10;
        constexpr auto minPlaneDistance = 20.f;

        auto& planes = SPGMT::Debug::RandomParallelPlanesSampling(planeSamplesCount, minPlaneDistance);
        std::vector<SPGMT::Point3> points;
        
        for(auto i = 0; i < planes.size(); ++i)
        {
            const auto& planePoints = SPGMT::Debug::RandomPointsPartitionedByPlane(pointSamplesCount, planes[i], minPlaneDistance - 1.f);
            std::copy(planePoints.mySamples.begin(), planePoints.mySamples.end(), std::back_inserter(points));
        }

        // run function
        SPGMT::LocationResult result;
        SPGMT::BatchPointLocation(planes, points, result);

        REQUIRE(result.myType == SPGMT::ResultType::PARALLEL_PLANES);
        REQUIRE(result.myRanges.size() == (planeSamplesCount * pointSamplesCount));
        REQUIRE(result.mySortedPlanesIndices.size() == planeSamplesCount);

        // Sort planes
        std::vector<SPGMT::Plane> sortedPlanes;
        for (auto i = 0; i < planes.size(); ++i)
        {
            const auto planeIdx = result.mySortedPlanesIndices[i];
            sortedPlanes.push_back(planes[planeIdx]);
        }

        // check ranges
        for (auto i = 0; i < result.myRanges.size(); ++i)
        {
            const int startPlaneIdx = result.myRanges[i].first;
            const int endPlaneIdx = result.myRanges[i].second;
            if (startPlaneIdx != -1)
            {
                const auto requirement = (sortedPlanes[startPlaneIdx].has_on_negative_side(points[i]) ||
                    sortedPlanes[startPlaneIdx].has_on(points[i])) && endPlaneIdx > startPlaneIdx;

                REQUIRE(requirement);
            }
            else
            {
                const auto requirement = sortedPlanes.back().has_on_positive_side(points[i]) && endPlaneIdx == -1;
                REQUIRE(requirement);
            }
        }
    }
     
    
    // Random planes
    SECTION("Single Random Plane and Random Points around it")
    {
        constexpr auto planeSamplesCount = 100;
        constexpr auto pointSamplesCount = 1000;

        for (const auto& plane : SPGMT::Debug::RandomPlaneSampling(planeSamplesCount))
        {
            // Generate 1000 points partitioned by the plane
            const auto& points = SPGMT::Debug::RandomPointsPartitionedByPlane(pointSamplesCount, plane);

            // run function
            SPGMT::LocationResult result;
            SPGMT::BatchPointLocation({ plane }, points.mySamples, result);

            REQUIRE(result.myType == SPGMT::ResultType::SINGLE_PLANE);
            REQUIRE(result.myRanges.size() == pointSamplesCount);

            // compare results
            auto planeAboveFoundCount { 0 };

            for (auto i = 0; i < result.myRanges.size(); ++i)
            {
                const auto range = result.myRanges[i];
                if (range.first == 0 && range.second == -1)
                {
                    planeAboveFoundCount += 1;
                }
            }

            REQUIRE(planeAboveFoundCount == (points.myNegativeCount + points.myOverSurfaceCount));
        }
    }
}