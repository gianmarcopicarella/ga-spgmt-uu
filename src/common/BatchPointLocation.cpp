#include "BatchPointLocation.h"

#include <vector>

#include <CGAL/Arrangement_2.h>

namespace SPGMT
{
	namespace
    {
        struct PlanePlaneVisitor
        {
            typedef Line3 payload_type;
            typedef void result_type;
            void operator()(const Line3& aLine)
            {
                myOutResult.push_back(aLine);
            }
            void operator()(const Plane& aPlane)
            {
                CGAL_precondition(false);
            }

            std::vector<payload_type>& myOutResult;
        };

        struct LineParallelLineVisitor
        {
            struct Data
            {
                FT myDistance;
                int myIndex;
            };
            typedef Data payload_type;
            typedef void result_type;
            void operator()(const Point2& aPoint)
            {
                static const FT zero {0};
                const Vec2 vecToPoint{ myLine.point(), aPoint };
                const auto distanceSign = CGAL::scalar_product(vecToPoint, myLine.to_vector()) > zero ? 1.f : -1.f;
                myOutResult.push_back({distanceSign * vecToPoint.squared_length(), myPlaneIndex});
            }
            void operator()(const Line2& /**/)
            {
                CGAL_precondition(false);
            }

            const Line2& myLine;
            const int myPlaneIndex;
            std::vector<payload_type>& myOutResult;
        };

        struct LinePlaneVisitor
        {
            struct Data
            {
                FT myDistance;
                int myIndex;
            };
            typedef Data payload_type;
            typedef void result_type;
            void operator()(const Point3& aPoint)
            {
                static const FT zero {0};
                const Vec3 vecToPoint{ myLine.point(), aPoint };
                const auto distanceSign = CGAL::scalar_product(vecToPoint, myLine.to_vector()) > zero ? 1.f : -1.f;
                myOutResult.push_back({distanceSign * vecToPoint.squared_length(), myPlaneIndex});
            }
            void operator()(const Line3& /**/)
            {
                CGAL_precondition(false);
            }

            const Line3& myLine;
            const int myPlaneIndex;
            std::vector<payload_type>& myOutResult;
        };

        struct Line2Line2Visitor
        {
            typedef Point2 payload_type;
            typedef void result_type;
            void operator()(const Point2& aPoint)
            {
                myOutResult.push_back(aPoint);
            }
            void operator()(const Line2& /**/)
            {
                CGAL_precondition(false);
            }

            std::vector<payload_type>& myOutResult;
        };

        template<typename T>
        bool locAreItemsUnique(const std::vector<T>& someItems)
        {
            auto areItemsUnique{ true };
            for (int i = 0; i < someItems.size() && areItemsUnique; ++i)
            {
                for (int k = i + 1; k < someItems.size(); ++k)
                {
                    if (someItems[i] == someItems[k])
                    {
                        areItemsUnique = false;
                        break;
                    }
                }
            }
            return areItemsUnique;
        }
        bool locArePlanesNonVertical(const std::vector<Plane>& somePlanes)
        {
            static const Vec3 verticalAxis{ 0,1,0 };
            static const FT zero {0};
            auto arePlanesNonVertical{ true };
            for (int i = 0; i < somePlanes.size() && arePlanesNonVertical; ++i)
            {
                arePlanesNonVertical =
                        CGAL::scalar_product(somePlanes[i].orthogonal_vector(), verticalAxis) != zero;
            }
            return arePlanesNonVertical;
        }
        bool locArePlanesUniformlyOriented(const std::vector<Plane>& somePlanes)
        {
            static const Vec3 up {0, 1, 0};
            static const FT zero {0};
            auto isPointingDown {false}, isPointingUp {false};
            for(auto& plane : somePlanes)
            {
                const auto dot = CGAL::scalar_product(up, plane.orthogonal_vector());
                // Cannot handle vertical planes
                CGAL_precondition(dot != zero);
                isPointingUp |= dot > zero;
                isPointingDown |= dot < zero;

                if(isPointingUp && isPointingDown)
                {
                    return false;
                }
            }
            return true;
        }
        template <typename Strategy, typename T>
        std::vector<typename Strategy::payload_type> locFindIntersections(const std::vector<T>& someItems, const bool aKeepOnlyUnique = false)
        {
            std::vector<typename Strategy::payload_type> result;
            for(int i = 0; i < someItems.size(); ++i)
            {
                for(int k = i + 1; k < someItems.size(); ++k)
                {
                    const auto intersection = CGAL::intersection(someItems[i], someItems[k]);
                    if (const auto* data = intersection.get_ptr())
                    {
                        Strategy visitor{ result };
                        data->apply_visitor(visitor);
                    }
                }
            }
            if(aKeepOnlyUnique)
            {
                result.erase(std::unique(result.begin(), result.end(), CGAL::Equal_to<typename Strategy::payload_type,
                        typename Strategy::payload_type>()), result.end());
            }
            return result;
        }
        template <typename Strategy, typename T, typename K>
        std::vector<typename Strategy::payload_type> locFindIntersections(const std::vector<T>& someItems, const K& anotherItem)
        {
            std::vector<typename Strategy::payload_type> result;
            for(int i = 0; i < someItems.size(); ++i)
            {
                const auto intersection = CGAL::intersection(someItems[i], anotherItem);
                if (const auto* data = intersection.get_ptr())
                {
                    Strategy visitor{ anotherItem, i, result };
                    data->apply_visitor(visitor);
                }
            }
            return result;
        }
        template<typename T>
        void locSortByDistanceMember(std::vector<T>& someOutData)
        {
            const auto sort = [](auto& aFirst, auto& aSecond) {return aFirst.myDistance < aSecond.myDistance;};
            std::sort(someOutData.begin(), someOutData.end(), sort);
        }
        template<typename T>
        std::vector<int> locExtractIndicesMember(const std::vector<T>& someData)
        {
            const auto transform = [](const auto& aFirst) {return aFirst.myIndex;};
            std::vector<int> indices;
            indices.reserve(someData.size());
            std::transform(someData.begin(), someData.end(), std::back_inserter(indices), transform);
            return indices;
        }
        template <typename Strategy, typename T, typename K>
        std::vector<int> locGetSortedItemsIndicesFromIntersections(const std::vector<T>& someItems, const K& anotherItem)
        {
            auto planesLineIntersections = locFindIntersections<Strategy>(someItems, anotherItem);
            locSortByDistanceMember(planesLineIntersections);
            return locExtractIndicesMember(planesLineIntersections);
        }

        Line2 Create2DLine(const Point2& aFirst, const Point2& aSecond)
        {
            CGAL_precondition(aFirst != aSecond);

            if (aFirst < aSecond)
            {
                return Line2{ aFirst, aSecond };
            }
            else
            {
                return Line2{ aSecond, aFirst };
            }
        }

        std::vector<Line2> locProjectOnXZPlane(const std::vector<Line3>& someLines)
        {
            static const FT distance {1000.f};
            std::vector<Line2> result;
            for(const auto& line : someLines)
            {
                const auto& translation = line.to_vector() * distance;
                const auto firstPoint{ line.point() };
                const auto secondPoint{ firstPoint + translation };
                const auto line2d = Create2DLine(
                        Point2{ firstPoint.x(), firstPoint.z() },
                        Point2{ secondPoint.x(), secondPoint.z() });
                result.push_back(line2d);
            }
            return result;
        }

        std::vector<Line2> locProjectOnXZPlaneUnique(const std::vector<Line3>& someLines)
        {
            auto result = locProjectOnXZPlane(someLines);
            result.erase(std::unique(result.begin(), result.end(),
                                     CGAL::Equal_to<Line2, Line2>()), result.end());
            return result;
        }

        template<typename T, typename K>
        int BinarySearch(const std::vector<T> & someItems, const std::vector<int>& someSortedIndices, const K& aPoint, int aMin, int aMax)
        {
            const auto mid = (aMin + aMax) / 2;
            const auto& item = someItems[someSortedIndices[mid]];

            if (item.has_on_negative_side(aPoint))
            {
                if (mid == 0)
                {
                    return 0;
                }
                else if (someItems[someSortedIndices[mid - 1]].has_on_positive_side(aPoint))
                {
                    return mid;
                }
                else
                {
                    return BinarySearch(someItems, someSortedIndices, aPoint, aMin, mid - 1);
                }
            }
            else if (item.has_on_positive_side(aPoint))
            {
                if (mid == someItems.size() - 1)
                {
                    return someItems.size();
                }
                else if (someItems[someSortedIndices[mid + 1]].has_on_negative_side(aPoint))
                {
                    return mid + 1;
                }
                else
                {
                    return BinarySearch(someItems, someSortedIndices, aPoint, mid + 1, aMax);
                }
            }
            else
            {
                if(mid == 0 || (mid > 0 && !someItems[someSortedIndices[mid - 1]].has_on(aPoint)))
                {
                    return mid;
                }
                else
                {
                    return BinarySearch(someItems, someSortedIndices, aPoint, 0, mid - 1);
                }
            }
        }
        template<typename T, typename K>
        int BinarySearch(const std::vector<T> & someItems, const K& aPoint, int aMin, int aMax)
        {
            const auto mid = (aMin + aMax) / 2;
            const auto& item = someItems[mid];

            if (item.has_on_negative_side(aPoint))
            {
                if (mid == 0)
                {
                    return 0;
                }
                else if (someItems[mid - 1].has_on_positive_side(aPoint))
                {
                    return mid;
                }
                else
                {
                    return BinarySearch(someItems, aPoint, aMin, mid - 1);
                }
            }
            else if (item.has_on_positive_side(aPoint))
            {
                if (mid == someItems.size() - 1)
                {
                    return someItems.size();
                }
                else if (someItems[mid + 1].has_on_negative_side(aPoint))
                {
                    return mid + 1;
                }
                else
                {
                    return BinarySearch(someItems, aPoint, mid + 1, aMax);
                }
            }
            else
            {
                if(mid == 0 || (mid > 0 && !someItems[mid - 1].has_on(aPoint)))
                {
                    return mid;
                }
                else
                {
                    return BinarySearch(someItems, aPoint, 0, mid - 1);
                }
            }
        }

        int FindSlabIndex(const std::vector<FT> & someItems, const Point2& aPoint, int aMin, int aMax)
        {
            const auto mid = (aMin + aMax) / 2;
            const auto& x = someItems[mid];

            if (aPoint.x() < x)
            {
                if (mid == 0)
                {
                    return -1;
                }
                else if (aPoint.x() > someItems[mid - 1])
                {
                    return mid - 1;
                }
                else
                {
                    return FindSlabIndex(someItems, aPoint, aMin, mid - 1);
                }
            }
            else if (aPoint.x() > x)
            {
                if (mid == (someItems.size() - 1) ||
                    aPoint.x() < someItems[mid + 1])
                {
                    return mid;
                }
                else
                {
                    return FindSlabIndex(someItems, aPoint, mid + 1, aMax);
                }
            }
            else
            {
                return mid;
            }
        }

        struct LineInfo
        {
            FT myY;
            FT mySlope;
            int myIndex;
        };

        struct SlabPartition
        {
            std::vector<FT> myUniqueXValues; // N values with N = # unique x-values
            std::unordered_map<int, int> myVerticalLineIndex;
            std::vector<std::vector<int>> mySortedLines; // N+1 values ( considering -inf partition)
        };

        struct LineData
        {
            std::vector<Line2> myLines;
            std::vector<std::vector<FT>> mySortedVertices;
            std::vector<Point2> myUniqueVertices;

            std::function<FT(const Point2&)> GetGetter(const Line2& aLine)
            {
                static const std::array<std::function<FT(const Point2&)>, 2> getters = {
                        [](const Point2& aPoint) { return aPoint.x(); },
                        [](const Point2& aPoint) { return aPoint.y(); }
                };
                const int getterIndex = static_cast<int>(aLine.is_vertical());
                return getters[getterIndex];
            }

            int FindVertexIndex(const Point2& aPoint, int aLineIdx)
            {
                CGAL_precondition(myLines[aLineIdx].has_on(aPoint));
                const auto& values = mySortedVertices[aLineIdx];
                const auto& getter = GetGetter(myLines[aLineIdx]);
                int left = 0, right = values.size();
                while(left < right)
                {
                    const int mid = (left + right) / 2;
                    if(values[mid] == getter(aPoint))
                    {
                        return mid;
                    }
                    else if(getter(aPoint) < values[mid])
                    {
                        right = mid - 1;
                    }
                    else
                    {
                        left = mid + 1;
                    }
                }
                return -1;
            }

            int FindSegmentIndex(const Point2& aPoint, int aLineIdx)
            {
                CGAL_precondition(myLines[aLineIdx].has_on(aPoint));
                const auto& values = mySortedVertices[aLineIdx];
                const auto& getter = GetGetter(myLines[aLineIdx]);

                if(getter(aPoint) < values[0])
                {
                    return 0;
                }

                int left = 0, right = values.size();
                while(left < right)
                {
                    const int mid = (left + right) / 2;
                    if(getter(aPoint) < values[mid])
                    {
                        if(mid == 0) return 0;
                        else if(getter(aPoint) > values[mid - 1])
                        {
                            return mid;
                        }
                        else
                        {
                            right = mid - 1;
                        }
                    }
                    else if(getter(aPoint) > values[mid])
                    {
                        if(mid == values.size() - 1) return mid;
                        else if(getter(aPoint) < values[mid + 1])
                        {
                            return mid + 1;
                        }
                        else
                        {
                            left = mid + 1;
                        }
                    }
                    else
                    {
                        CGAL_precondition(false);
                    }
                }
                return -1;
            }
        };

        LineData locComputeLinesWrapper(const std::vector<Line3>& someLines)
        {
            LineData data;

            data.myLines = locProjectOnXZPlane(someLines);
            data.myLines.erase( std::unique(data.myLines.begin(), data.myLines.end(),
                                CGAL::Equal_to<Line2, Line2>()), data.myLines.end());

            for(int i = 0; i < data.myLines.size(); ++i)
            {
                const auto& line = data.myLines[i];
                const auto getter = data.GetGetter(line);

                std::vector<FT> values;

                for(int k = 0; k < data.myLines.size(); ++k)
                {
                    if(i != k)
                    {
                        const auto result = CGAL::intersection(data.myLines[i], data.myLines[k]);
                        if (result)
                        {
                            if (const Point2* p = boost::get<Point2>(&*result))
                            {
                                values.push_back(getter(*p));
                                data.myUniqueVertices.push_back(*p);
                            }
                            else
                            {
                                CGAL_precondition(false);
                            }
                        }
                    }
                }

                std::sort(values.begin(), values.end(), CGAL::Less<FT, FT>());
                data.mySortedVertices.push_back(values);
            }

            data.myUniqueVertices.erase( std::unique(data.myUniqueVertices.begin(), data.myUniqueVertices.end(),
                                            CGAL::Less<Point2, Point2>()), data.myUniqueVertices.end());

            return data;
        }

        SlabPartition locComputeSlabs(const std::vector<Point2> & someVertices, const std::vector<Line2>& someUniqueLines)
        {
            const auto slopeTransform = [](const auto& aLine) {
                std::optional<FT> result{ std::nullopt };
                if (!aLine.is_vertical())
                {
                    const auto& dir = aLine.direction();
                    result = dir.dy() / dir.dx();
                }
                return result;
            };
            const auto sort = [](const auto& aFirst, const auto& aSecond) {
                return aFirst.myY < aSecond.myY ||
                       (aFirst.myY == aSecond.myY && aFirst.mySlope < aSecond.mySlope);
            };
            const auto sortInverseSlope = [](const auto& aFirst, const auto& aSecond) {
                return aFirst.myY < aSecond.myY ||
                       (aFirst.myY == aSecond.myY && aFirst.mySlope > aSecond.mySlope);
            };

            // compute slope
            std::vector<std::optional<FT>> slopes;

            // Compute slopes for each line
            std::transform(someUniqueLines.begin(), someUniqueLines.end(), std::back_inserter(slopes), slopeTransform);

            std::vector<size_t> sortedVerticesIndices;
            sortedVerticesIndices.resize(someVertices.size());
            std::iota(sortedVerticesIndices.begin(), sortedVerticesIndices.end(), 0);

            // Sorted along x-axis and sorted from bottom to top along y-axis when the x value is the same
            std::sort(sortedVerticesIndices.begin(), sortedVerticesIndices.end(),
            [&someVertices](auto a, auto b) {
                return someVertices[a] < someVertices[b];
            });

            SlabPartition partition;

            for(int i = 0; i < someVertices.size(); ++i)
            {
                const auto& currentVertex = someVertices[sortedVerticesIndices[i]];

                if(!partition.myUniqueXValues.empty() && partition.myUniqueXValues.back() == currentVertex.x())
                {
                    //slabs.myInfo.back().mySortedVertices.push_back(sortedVerticesIndices[i]);
                    continue;
                }
                else
                {
                    // Get slab index
                    const int slabIdx = partition.myUniqueXValues.size();
                    // Add x value
                    partition.myUniqueXValues.push_back(currentVertex.x());

                    // Sort lines
                    std::vector<LineInfo> sortedLines;

                    for(int k = 0; k < someUniqueLines.size(); ++k)
                    {
                        const auto& currentLine = someUniqueLines[k];
                        if(currentLine.is_vertical())
                        {
                            partition.myVerticalLineIndex.insert(std::make_pair(slabIdx, k));
                        }
                        else
                        {
                            const auto& dir = currentLine.direction();
                            const auto slope = dir.dy() / dir.dx();
                            const auto Y = currentLine.y_at_x(currentVertex.x());
                            sortedLines.push_back({Y, slope, k});
                        }
                    }

                    std::vector<int> sortedIndices;
                    if(i == 0)
                    {
                        std::sort(sortedLines.begin(),sortedLines.end(), sortInverseSlope);
                        std::transform(sortedLines.begin(), sortedLines.end(), std::back_inserter(sortedIndices),
                                       [](auto& a) {return a.myIndex;});
                        partition.mySortedLines.push_back(sortedIndices);
                    }

                    sortedIndices.clear();
                    std::sort(sortedLines.begin(),sortedLines.end(), sort);
                    std::transform(sortedLines.begin(), sortedLines.end(), std::back_inserter(sortedIndices),
                                   [](auto& a) {return a.myIndex;});
                    partition.mySortedLines.push_back(sortedIndices);
                }
            }

            return partition;
        }

    }

    Result BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints)
    {
        if(somePlanes.empty() || somePoints.empty())
        {
            return Result{};
        }

        // Requirements check
        CGAL_precondition(locAreItemsUnique(somePlanes));
        CGAL_precondition(locArePlanesNonVertical(somePlanes));
        CGAL_precondition(locArePlanesUniformlyOriented(somePlanes));

        // Compute plane-plane intersections
        const auto lines3d = locFindIntersections<PlanePlaneVisitor>(somePlanes);

        // 1st edge case: single plane or parallel planes
        if(lines3d.empty())
        {
            const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
                                                                                       somePlanes[0].perpendicular_line(somePlanes[0].point()));
            Result result;
            result.mySortedPlanesIndices.insert(std::make_pair(0, sortedPlanesIndices));

            for(const auto& queryPoint : somePoints)
            {
                const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
                                                queryPoint, 0, sortedPlanesIndices.size());
                result.myRangeWrappers.push_back({Range{firstUpperPlaneIdx, somePlanes.size()}, 0});
            }

            return result;
        }

        // Project 3d lines and construct a list of enhanced lines


        // Compute 3d lines projection on xz-plane
        const auto lines2d = locProjectOnXZPlaneUnique(lines3d);
        constexpr auto uniqueVertices = true;
        const auto vertices = locFindIntersections<Line2Line2Visitor>(lines2d, uniqueVertices);

        // 2nd edge case: single or parallel lines
        if(vertices.empty())
        {
            const auto sortedLinesIndices = locGetSortedItemsIndicesFromIntersections<LineParallelLineVisitor>(lines2d,
                                                                                                               lines2d[0].perpendicular(lines2d[0].point()));
            Result result;
            for(const auto& queryPoint : somePoints)
            {
                const Point2 projectedPoint { queryPoint.x(), queryPoint.z() };
                const auto firstUpperLineIdx = BinarySearch<Line2, Point2>(lines2d, sortedLinesIndices,
                                                                           projectedPoint, 0, sortedLinesIndices.size());

                const auto isPointAlongLine = firstUpperLineIdx < sortedLinesIndices.size() &&
                        lines2d[sortedLinesIndices[firstUpperLineIdx]].has_on(projectedPoint);

                const int bucketIndex = isPointAlongLine ? firstUpperLineIdx : (firstUpperLineIdx + lines2d.size());

                if(result.mySortedPlanesIndices.find(bucketIndex) == result.mySortedPlanesIndices.end())
                {
                    const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
                                                                                                                 Line3{queryPoint, Vec3{0,1,0}});
                    result.mySortedPlanesIndices.insert(std::make_pair(bucketIndex, sortedPlanesIndices));
                }

                const auto sortedPlanesIndices = result.mySortedPlanesIndices.at(bucketIndex);
                const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
                                                                            queryPoint, 0, sortedPlanesIndices.size());
                result.myRangeWrappers.push_back({Range{firstUpperPlaneIdx, somePlanes.size()}, bucketIndex });
            }
            return result;
        }

        // Slab-based approach
        const auto partition = locComputeSlabs(vertices, lines2d);

        int cacheUsageCount = 0;

        Result result;
        for(const auto& queryPoint : somePoints)
        {
            const Point2 projectedPoint { queryPoint.x(), queryPoint.z() };
            int slabIdx = FindSlabIndex(partition.myUniqueXValues, projectedPoint, 0, partition.myUniqueXValues.size());

            // Found slab idx and it is not -inf
            if(slabIdx > -1)
            {
                // Along vertical axis
                if(partition.myUniqueXValues[slabIdx] == projectedPoint.x())
                {
                    // Use bruteforce
                    // Return range
                    result.myRangeWrappers.push_back({Range{}, -1});
                    continue;
                }

                slabIdx += 1;
            }
            else
            {
                slabIdx = 0;
            }

            const auto& sortedLines = partition.mySortedLines[slabIdx];
            const int lineIdx = BinarySearch(lines2d, sortedLines, projectedPoint, 0, sortedLines.size());

            if(lineIdx < sortedLines.size() && lines2d[sortedLines[lineIdx]].has_on(projectedPoint))
            {
                // Use bruteforce
                // Return range
                result.myRangeWrappers.push_back({Range{}, -1});
                continue;
            }

            const int globalFaceIdx = lineIdx + slabIdx * lines2d.size();

            if(result.mySortedPlanesIndices.find(globalFaceIdx) == result.mySortedPlanesIndices.end())
            {
                const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
                                                                                                             Line3{queryPoint, Vec3{0,1,0}});
                result.mySortedPlanesIndices.insert(std::make_pair(globalFaceIdx, sortedPlanesIndices));
            }
            else
            {
                // Temporary
                cacheUsageCount++;
            }

            const auto sortedPlanesIndices = result.mySortedPlanesIndices.at(globalFaceIdx);
            const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
                                                                        queryPoint, 0, sortedPlanesIndices.size());
            result.myRangeWrappers.push_back({Range{firstUpperPlaneIdx, somePlanes.size()}, globalFaceIdx });
        }

        std::cout << "Cache used: " << cacheUsageCount << " times";
        return result;
    }

}