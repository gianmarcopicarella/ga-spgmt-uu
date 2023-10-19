#include "BatchPointLocation.h"
#include "../common/Utils.h"

#include <vector>

#include <hpx/hpx.hpp>

namespace SPGMT
{

	namespace
	{
		template<ExecutionPolicy E, typename Action>
		auto BindExecutionPolicy(Action&& anOutAction)
		{
			switch (E)
			{
			case SPGMT::ExecutionPolicy::PARALLEL:
				return anOutAction(hpx::execution::par);
			default:
				return anOutAction(hpx::execution::seq);
			}
		}


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
				const FT zero{ 0 };
				const Vec2 vecToPoint{ myLine.point(), aPoint };
				const auto distanceSign = CGAL::scalar_product(vecToPoint, myLine.to_vector()) > zero ? 1.f : -1.f;
				myOutResult.push_back({ distanceSign * vecToPoint.squared_length(), myPlaneIndex });
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
				const FT zero{ 0 };
				const Vec3 vecToPoint{ myLine.point(), aPoint };
				const auto distanceSign = CGAL::scalar_product(vecToPoint, myLine.to_vector()) > zero ? 1.f : -1.f;
				myOutResult.push_back({ distanceSign * vecToPoint.squared_length(), myPlaneIndex });
			}
			void operator()(const Line3& /**/)
			{
				CGAL_precondition(false);
			}

			const Line3& myLine;
			const int myPlaneIndex;
			std::vector<payload_type>& myOutResult;
		};

		template<typename T>
		void locSortByDistanceMember(std::vector<T>& someOutData)
		{
			const auto sort = [](auto& aFirst, auto& aSecond) {return aFirst.myDistance < aSecond.myDistance; };
			std::sort(someOutData.begin(), someOutData.end(), sort);
		}
		template<typename T>
		std::vector<int> locExtractIndicesMember(const std::vector<T>& someData)
		{
			const auto transform = [](const auto& aFirst) {return aFirst.myIndex; };
			std::vector<int> indices;
			indices.reserve(someData.size());
			std::transform(someData.begin(), someData.end(), std::back_inserter(indices), transform);
			return indices;
		}
		template <typename Strategy, typename T, typename K>
		std::vector<int> locGetSortedItemsIndicesFromIntersections(const std::vector<T>& someItems, const K& anotherItem)
		{
			auto planesLineIntersections = Utils::FindIntersections<Strategy>(someItems, anotherItem);
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

		template<ExecutionPolicy E>
		std::vector<Line2> locProjectOnXYPlane(const std::vector<Line3>& someLines)
		{
			static const FT distance{ 1000.f };
			std::vector<Line2> result(someLines.size());

			BindExecutionPolicy<E>(std::bind(hpx::for_loop, std::placeholders::_1, 0, someLines.size(), [&](size_t aLineIdx) {
				const auto& translation = someLines[aLineIdx].to_vector() * distance;
				const auto firstPoint{ someLines[aLineIdx].point() };
				const auto secondPoint{ firstPoint + translation };
				result[aLineIdx] = Create2DLine(
					Point2{ firstPoint.x(), firstPoint.y() },
					Point2{ secondPoint.x(), secondPoint.y() });
				}));

			return result;
		}
		template<ExecutionPolicy E>
		std::vector<Line2> locProjectOnXYPlaneUnique(const std::vector<Line3>& someLines)
		{
			auto result = locProjectOnXYPlane<E>(someLines);
			auto uniqueIter =
				BindExecutionPolicy<E>(std::bind(hpx::unique, std::placeholders::_1,
					result.begin(), result.end(), CGAL::Equal_to<Line2, Line2>()));
			result.erase(uniqueIter, result.end());
			return result;
		}

		template<typename T, typename K>
		int BinarySearch(const std::vector<T>& someItems, const std::vector<int>& someSortedIndices, const K& aPoint, int aMin, int aMax)
		{

			const auto mid = aMin + (aMax - aMin) / 2;
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
				if (mid == 0 || (mid > 0 && !someItems[someSortedIndices[mid - 1]].has_on(aPoint)))
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
		int BinarySearch(const std::vector<T>& someItems, const K& aPoint, int aMin, int aMax)
		{
			const auto mid = aMin + (aMax - aMin) / 2;
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
				if (mid == 0 || (mid > 0 && !someItems[mid - 1].has_on(aPoint)))
				{
					return mid;
				}
				else
				{
					return BinarySearch(someItems, aPoint, 0, mid - 1);
				}
			}
		}

		int FindSlabIndex(const std::vector<FT>& someItems, const Point2& aPoint, int aMin, int aMax)
		{
			const auto mid = aMin + (aMax - aMin) / 2;
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
			std::unordered_map<int, int> myVerticalLineIndex;
			std::vector<std::vector<int>> mySortedLines; // N+1 values ( considering -inf partition)
		};

		struct LineData
		{
			std::vector<Line2> myLines;
			std::vector<std::vector<std::pair<FT, int>>> mySortedVertices;
			std::vector<FT> myUniqueAndSortedVerticesX;
			std::vector<int> mySegmentsCount;
			int myUniqueVerticesCount;

			const std::function<FT(const Point2&)>& GetGetter(const Line2& aLine) const
			{
				static const std::array<std::function<FT(const Point2&)>, 2> getters = {
						[](const Point2& aPoint) { return aPoint.x(); },
						[](const Point2& aPoint) { return aPoint.y(); }
				};
				const auto getterIndex = static_cast<int>(aLine.is_vertical());
				return getters[getterIndex];
			}
			struct LookupResult
			{
				int myIndex{ -1 };
				bool myIsVertex{ false };
			};
			LookupResult FindSegmentOrVertexIndex(const Point2& aPoint, int aLineIdx) const
			{
				CGAL_precondition(myLines[aLineIdx].has_on(aPoint));
				const auto& values = mySortedVertices[aLineIdx];
				const auto& getter = GetGetter(myLines[aLineIdx]);
				int left = 0, right = values.size();
				LookupResult result;
				while (left < right)
				{
					const int mid = left + (right - left) / 2;
					if (getter(aPoint) < values[mid].first)
					{
						if (mid == 0 || getter(aPoint) > values[mid - 1].first)
						{
							result.myIndex = mid;
							return result;
						}
						else
						{
							right = mid;
						}
					}
					else if (getter(aPoint) > values[mid].first)
					{
						if (mid == values.size() - 1 ||
							getter(aPoint) < values[mid + 1].first)
						{
							result.myIndex = mid + 1;
							return result;
						}
						else
						{
							left = mid + 1;
						}
					}
					else
					{
						result.myIndex = values[mid].second;
						result.myIsVertex = true;
						return result;
					}
				}

				return result;
			}
		};

		LineData locComputeLinesWrapper(const std::vector<Line3>& someLines)
		{
			struct PointIndexPair
			{
				Point2 myKey;
				int myIndex{ -1 };
			};

			LineData data;

			data.myLines = locProjectOnXYPlaneUnique<ExecutionPolicy::SEQUENTIAL>(someLines);
			data.mySortedVertices.resize(data.myLines.size());

			constexpr auto pointIndexPairCmp = [](auto& a, auto& b)
			{ return a.myKey < b.myKey; };

			data.myUniqueVerticesCount = 0;
			std::set<FT, CGAL::Less<FT, FT>> allUniqueVertices{ CGAL::Less<FT, FT>{} };

			for (int i = 0; i < data.myLines.size(); ++i)
			{
				std::set<PointIndexPair, decltype(pointIndexPairCmp)> uniqueVertices{ pointIndexPairCmp };

				for (int k = i + 1; k < data.myLines.size(); ++k)
				{
					const auto result = CGAL::intersection(data.myLines[i], data.myLines[k]);
					if (result)
					{
						if (const Point2* point = boost::get<Point2>(&*result))
						{
							auto vertexIter = uniqueVertices.find(PointIndexPair{ *point });
							if (vertexIter == uniqueVertices.end())
							{
								const int vertexIndex = data.myUniqueVerticesCount;
								vertexIter = uniqueVertices.insert(PointIndexPair{ *point, vertexIndex }).first;
								const auto& firstGetter = data.GetGetter(data.myLines[i]);
								data.mySortedVertices[i].push_back(std::make_pair(firstGetter(*point), vertexIndex));
								++data.myUniqueVerticesCount;
							}

							const int vertexIndex = vertexIter->myIndex;
							const auto& secondGetter = data.GetGetter(data.myLines[k]);
							data.mySortedVertices[k].push_back(std::make_pair(secondGetter(*point), vertexIndex));
						}
						else
						{
							CGAL_precondition(false);
						}
					}
				}

				for (const auto& uniqueVertex : uniqueVertices)
				{
					allUniqueVertices.insert(uniqueVertex.myKey.x());
				}
			}

			constexpr auto pairComparator = [](auto& aFirst, auto& aSecond) { return aFirst.first < aSecond.first; };

			data.mySegmentsCount.push_back(0);
			std::sort(data.mySortedVertices[0].begin(), data.mySortedVertices[0].end(), pairComparator);

			for (int i = 1; i < data.myLines.size(); ++i)
			{
				std::sort(data.mySortedVertices[i].begin(), data.mySortedVertices[i].end(), pairComparator);
				data.mySegmentsCount.push_back(data.mySortedVertices[i - 1].size() + 1 + data.mySegmentsCount.back());
			}

			data.mySegmentsCount.push_back(data.mySortedVertices.back().size() + 1 + data.mySegmentsCount.back());

			data.myUniqueAndSortedVerticesX.resize(allUniqueVertices.size());
			std::copy(allUniqueVertices.begin(), allUniqueVertices.end(), data.myUniqueAndSortedVerticesX.begin());

			return data;
		}

		namespace Common
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
		}

		namespace Parallel
		{
			SlabPartition locComputeSlabs(const LineData& aLineData)
			{
				using namespace Common;
				// compute slope
				std::vector<std::optional<FT>> slopes(aLineData.myLines.size());

				// Compute slopes for each line
				hpx::for_loop(hpx::execution::par_unseq, 0, aLineData.myLines.size(), [&](size_t aLineIdx) {
					slopes[aLineIdx] = slopeTransform(aLineData.myLines[aLineIdx]);
					});

				SlabPartition partition;
				partition.mySortedLines.resize(aLineData.myUniqueAndSortedVerticesX.size() + 1);

				std::mutex verticalLinesMutex;

				hpx::for_loop(hpx::execution::par_unseq, 0, aLineData.myUniqueAndSortedVerticesX.size(), [&](size_t aSlabIdx) {

					const auto& currentX = aLineData.myUniqueAndSortedVerticesX[aSlabIdx];
					std::vector<LineInfo> sortedLines(aLineData.myLines.size());
					std::atomic<size_t> atomicCounter{ 0 };

					hpx::for_loop(hpx::execution::par, 0, aLineData.myLines.size(), [&](size_t k) {
						const auto& currentLine = aLineData.myLines[k];
						if (currentLine.is_vertical())
						{
							std::lock_guard guard(verticalLinesMutex);
							CGAL_precondition(partition.myVerticalLineIndex.find(aSlabIdx) ==
								partition.myVerticalLineIndex.end());
							partition.myVerticalLineIndex.insert(std::make_pair((int)aSlabIdx, (int)k));
						}
						else
						{
							const auto& dir = currentLine.direction();
							const auto slope = dir.dy() / dir.dx();
							const auto Y = currentLine.y_at_x(currentX);
							sortedLines[atomicCounter.fetch_add(1)] = LineInfo{ Y, slope, (int)k };
						}
						});
					sortedLines.resize(atomicCounter.fetch_xor(atomicCounter.load()));

					if (aSlabIdx == 0)
					{
						hpx::sort(hpx::execution::par_unseq, sortedLines.begin(), sortedLines.end(), sortInverseSlope);

						partition.mySortedLines[aSlabIdx].resize(sortedLines.size());

						hpx::transform(hpx::execution::par_unseq, sortedLines.begin(), sortedLines.end(), 
							partition.mySortedLines[aSlabIdx].begin(), [](const auto& aLine) { return aLine.myIndex; });
					}

					hpx::sort(hpx::execution::par_unseq, sortedLines.begin(), sortedLines.end(), sort);
					
					partition.mySortedLines[aSlabIdx + 1].resize(sortedLines.size());

					hpx::transform(hpx::execution::par_unseq, sortedLines.begin(), sortedLines.end(),
						partition.mySortedLines[aSlabIdx + 1].begin(), [](const auto& aLine) { return aLine.myIndex; });
					});

				return partition;
			}
		}


		SlabPartition locComputeSlabs(const LineData& aLineData)
		{
			using namespace Common;
			// compute slope
			std::vector<std::optional<FT>> slopes(aLineData.myLines.size());

			// Compute slopes for each line
			hpx::for_loop(hpx::execution::seq, 0, aLineData.myLines.size(), [&](size_t aLineIdx) {
				slopes[aLineIdx] = slopeTransform(aLineData.myLines[aLineIdx]);
				});

			SlabPartition partition;

			for (int slabIdx = 0; slabIdx < aLineData.myUniqueAndSortedVerticesX.size(); ++slabIdx)
			{
				const auto& currentX = aLineData.myUniqueAndSortedVerticesX[slabIdx];
				std::vector<LineInfo> sortedLines;

				for (int k = 0; k < aLineData.myLines.size(); ++k)
				{
					const auto& currentLine = aLineData.myLines[k];
					if (currentLine.is_vertical())
					{
						CGAL_precondition(partition.myVerticalLineIndex.find(slabIdx) ==
							partition.myVerticalLineIndex.end());
						partition.myVerticalLineIndex.insert(std::make_pair(slabIdx, k));
					}
					else
					{
						const auto& dir = currentLine.direction();
						const auto slope = dir.dy() / dir.dx();
						const auto Y = currentLine.y_at_x(currentX);
						sortedLines.push_back({ Y, slope, k });
					}
				}

				std::vector<int> sortedIndices;
				if (slabIdx == 0)
				{
					std::sort(sortedLines.begin(), sortedLines.end(), sortInverseSlope);
					std::transform(sortedLines.begin(), sortedLines.end(), std::back_inserter(sortedIndices),
						[](auto& a) {return a.myIndex; });
					partition.mySortedLines.push_back(sortedIndices);
				}

				sortedIndices.clear();
				std::sort(sortedLines.begin(), sortedLines.end(), sort);
				std::transform(sortedLines.begin(), sortedLines.end(), std::back_inserter(sortedIndices),
					[](auto& a) {return a.myIndex; });
				partition.mySortedLines.push_back(sortedIndices);
			}

			return partition;
		}
	}

	// Skeleton
	BatchPointResult ParallelBatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints)
	{
		constexpr auto policy = ExecutionPolicy::PARALLEL;

		if (somePlanes.empty() || somePoints.empty())
		{
			return BatchPointResult{};
		}

		// Requirements check
		CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));
		CGAL_precondition(Utils::ArePlanesUniformlyOriented(somePlanes));

		// Find 3d lines
		const auto maxLinesCount = somePlanes.size() * (somePlanes.size() - 1) / 2;
		std::vector<Line3> lines3d(maxLinesCount);
		std::atomic<size_t> atomicCount{ 0 };

		hpx::for_loop(hpx::execution::par, 0, somePlanes.size(), [&](size_t aPlaneIdx) {

			hpx::for_loop(hpx::execution::par, aPlaneIdx + 1, somePlanes.size(), [&](size_t anoterPlaneIdx) {
				const auto intersection = CGAL::intersection(somePlanes[aPlaneIdx], somePlanes[anoterPlaneIdx]);
				if (const auto* data = intersection.get_ptr())
				{
					const Line3* line = boost::get<Line3>(&*intersection);
					CGAL_precondition(line != nullptr);
					lines3d[atomicCount.fetch_add(1)] = *line;
				}
				});

			});

		lines3d.resize(atomicCount.fetch_xor(atomicCount.load()));

		// 1st edge case: single plane or parallel planes
		if (lines3d.empty())
		{
			BatchPointResult result;
			std::vector<int> sortedPlanesIndices(somePlanes.size());

			// Compute sorted list of plane indices
			{
				std::vector<std::pair<FT, size_t>> localPlanes(somePlanes.size());
				const auto line = somePlanes[0].perpendicular_line(somePlanes[0].point());
				hpx::for_loop(hpx::execution::par_unseq, 0, somePlanes.size(), [&](size_t aPlaneIdx) {
					const auto intersection = CGAL::intersection(somePlanes[aPlaneIdx], line);
					const auto* data = intersection.get_ptr();
					const Point3* point = boost::get<Point3>(&*intersection);
					CGAL_precondition(point != nullptr);
					localPlanes[aPlaneIdx] = std::make_pair(point->z(), aPlaneIdx);
					});

				hpx::sort(hpx::execution::par_unseq, localPlanes.begin(), localPlanes.end(),
					[&](const auto& aFirst, const auto& aSecond) { return aFirst.first < aSecond.first; });
				hpx::for_loop(hpx::execution::par_unseq, 0, somePlanes.size(), [&](size_t aPlaneIdx) {
					sortedPlanesIndices[aPlaneIdx] = localPlanes[aPlaneIdx].second;
					});

			}

			result.mySortedPlanesIndices.insert(std::make_pair(0, sortedPlanesIndices));
			result.myRangeWrappers.resize(somePoints.size());

			hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {
				const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
				somePoints[aPointIdx], 0, sortedPlanesIndices.size());
			result.myRangeWrappers[aPointIdx] = { Range{firstUpperPlaneIdx, somePlanes.size()}, 0 };
				});

			return result;
		}

		// Compute lines data
		// TODO: Find a way to parallelize efficiently
		const auto linesData = locComputeLinesWrapper(lines3d);

		// 2nd edge case: single or parallel lines
		if (linesData.myUniqueVerticesCount == 0)
		{
			BatchPointResult result;
			std::vector<int> sortedLinesIndices(linesData.myLines.size());

			// Compute sorted list of lines indices
			{
				std::vector<std::pair<FT, size_t>> localLines(linesData.myLines.size());
				const auto line = linesData.myLines[0].perpendicular(linesData.myLines[0].point());
				hpx::for_loop(hpx::execution::par_unseq, 0, linesData.myLines.size(), [&](size_t aLineIdx) {
					const auto intersection = CGAL::intersection(linesData.myLines[aLineIdx], line);
					const auto* data = intersection.get_ptr();
					const Point2* point = boost::get<Point2>(&*intersection);
					CGAL_precondition(point != nullptr);

					const Vec2 vecToPoint{ line.point(), *point };
					const auto distanceSign = CGAL::sign(CGAL::scalar_product(vecToPoint, line.to_vector()))
						== CGAL::Sign::POSITIVE ? FT(1.f) : FT(-1.f);

					localLines[aLineIdx] = std::make_pair(distanceSign * vecToPoint.squared_length(), aLineIdx);
					});

				hpx::sort(hpx::execution::par_unseq, localLines.begin(), localLines.end(),
					[&](const auto& aFirst, const auto& aSecond) { return aFirst.first < aSecond.first; });
				hpx::for_loop(hpx::execution::par_unseq, 0, localLines.size(), [&](size_t aLineIdx) {
					sortedLinesIndices[aLineIdx] = localLines[aLineIdx].second;
					});
			}

			std::mutex planesIndicesMutex;

			result.myRangeWrappers.resize(somePoints.size());

			hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {

				const Point2 projectedPoint{ somePoints[aPointIdx].x(), somePoints[aPointIdx].y() };
				const auto firstUpperLineIdx = BinarySearch<Line2, Point2>(linesData.myLines, sortedLinesIndices,
					projectedPoint, 0, sortedLinesIndices.size());

				const auto isPointAlongLine = firstUpperLineIdx < sortedLinesIndices.size() &&
					linesData.myLines[sortedLinesIndices[firstUpperLineIdx]].has_on(projectedPoint);

				const int bucketIndex = isPointAlongLine ? firstUpperLineIdx : (firstUpperLineIdx + linesData.myLines.size());

				std::vector<int> sortedPlanesIndices;

				{
					std::lock_guard guard(planesIndicesMutex);
					if (result.mySortedPlanesIndices.find(bucketIndex) == result.mySortedPlanesIndices.end())
					{
						const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
							Line3{ somePoints[aPointIdx], Vec3{0,0,1} });
						result.mySortedPlanesIndices.insert(std::make_pair(bucketIndex, sortedPlanesIndices));
					}

					sortedPlanesIndices = result.mySortedPlanesIndices.at(bucketIndex);
				}

				const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
					somePoints[aPointIdx], 0, sortedPlanesIndices.size());
				result.myRangeWrappers[aPointIdx] = { Range{firstUpperPlaneIdx, somePlanes.size()}, bucketIndex };


				});

			return result;
		}

		// Slab-based approach
		const auto partition = Parallel::locComputeSlabs(linesData);

		BatchPointResult result;

		result.myRangeWrappers.resize(somePoints.size());

		std::mutex planesIndicesMutex;

		hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {
			const Point2 projectedPoint{ somePoints[aPointIdx].x(), somePoints[aPointIdx].y()};
			int slabIdx = FindSlabIndex(linesData.myUniqueAndSortedVerticesX, 
				projectedPoint, 0, linesData.myUniqueAndSortedVerticesX.size());
			int globalBucketIndex = -1;

			// Found slab idx and it is not -inf
			if (slabIdx > -1)
			{
				CGAL_precondition(linesData.myUniqueAndSortedVerticesX[slabIdx] <= projectedPoint.x());
				// Check if point is along vertical slab line
				if (linesData.myUniqueAndSortedVerticesX[slabIdx] == projectedPoint.x())
				{
					const auto vLineIter = partition.myVerticalLineIndex.find(slabIdx);
					if (vLineIter != partition.myVerticalLineIndex.end())
					{
						const int vLineIdx = vLineIter->second;
						const auto lineSearchResult = linesData.FindSegmentOrVertexIndex(projectedPoint, vLineIdx);
						if (!lineSearchResult.myIsVertex && lineSearchResult.myIndex > -1)
						{
							// Compute global index
							globalBucketIndex = linesData.myUniqueVerticesCount + // Vertices count
								linesData.mySegmentsCount[vLineIdx] + lineSearchResult.myIndex;
						}
					}
				}
				slabIdx += 1;
			}
			else
			{
				slabIdx = 0;
			}

			const auto& sortedLines = partition.mySortedLines[slabIdx];
			const int lineIdx = BinarySearch(linesData.myLines, sortedLines, projectedPoint, 0, sortedLines.size());

			if (lineIdx < sortedLines.size() && linesData.myLines[sortedLines[lineIdx]].has_on(projectedPoint))
			{
				const auto lineSearchResult = linesData.FindSegmentOrVertexIndex(projectedPoint, sortedLines[lineIdx]);
				CGAL_precondition(lineSearchResult.myIndex > -1);

				if (lineSearchResult.myIsVertex)
				{
					// Compute global index
					globalBucketIndex = lineSearchResult.myIndex;
					//CGAL_precondition(lineSearchResult.myIndex < linesData.myUniqueVerticesCount);
				}
				else
				{
					globalBucketIndex = linesData.myUniqueVerticesCount + // Vertices count
						linesData.mySegmentsCount[sortedLines[lineIdx]] + lineSearchResult.myIndex;
				}
			}

			if (globalBucketIndex == -1)
			{
				globalBucketIndex = linesData.myUniqueVerticesCount + // Vertices count
					linesData.mySegmentsCount.back() + // Segments count
					lineIdx + slabIdx * (linesData.myLines.size() + 1); // Faces count
			}

			std::vector<int> sortedPlanesIndices;

			{
				std::lock_guard guard(planesIndicesMutex);
				if (result.mySortedPlanesIndices.find(globalBucketIndex) == result.mySortedPlanesIndices.end())
				{
					const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
						Line3{ somePoints[aPointIdx], Vec3{0,0,1}});
					result.mySortedPlanesIndices.insert(std::make_pair(globalBucketIndex, sortedPlanesIndices));
				}

				sortedPlanesIndices = result.mySortedPlanesIndices.at(globalBucketIndex);
			}

			//CGAL_precondition(SameLists(sortedPlanesIndices, sortedPlanesIndices_cache));
			const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
				somePoints[aPointIdx], 0, sortedPlanesIndices.size());
			result.myRangeWrappers[aPointIdx] = { Range{firstUpperPlaneIdx, somePlanes.size()}, globalBucketIndex };

			});

		return result;
	}


	BatchPointResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints)
	{
		if (somePlanes.empty() || somePoints.empty())
		{
			return BatchPointResult{};
		}

		// Requirements check
		CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));
		CGAL_precondition(Utils::ArePlanesUniformlyOriented(somePlanes));

		// Compute plane-plane intersections
		const auto& lines3d = Utils::FindIntersections<Utils::PlanePlaneVisitor>(somePlanes);

		// 1st edge case: single plane or parallel planes
		if (lines3d.empty())
		{
			const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
				somePlanes[0].perpendicular_line(somePlanes[0].point()));
			BatchPointResult result;
			result.mySortedPlanesIndices.insert(std::make_pair(0, sortedPlanesIndices));

			for (const auto& queryPoint : somePoints)
			{
				const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
					queryPoint, 0, sortedPlanesIndices.size());
				result.myRangeWrappers.push_back({ Range{firstUpperPlaneIdx, somePlanes.size()}, 0 });
			}

			return result;
		}

		// Compute lines data
		const auto linesData = locComputeLinesWrapper(lines3d);

		// 2nd edge case: single or parallel lines
		if (linesData.myUniqueVerticesCount == 0)
		{
			const auto sortedLinesIndices = locGetSortedItemsIndicesFromIntersections<LineParallelLineVisitor>(linesData.myLines,
				linesData.myLines[0].perpendicular(linesData.myLines[0].point()));
			BatchPointResult result;
			for (const auto& queryPoint : somePoints)
			{
				const Point2 projectedPoint{ queryPoint.x(), queryPoint.y() };
				const auto firstUpperLineIdx = BinarySearch<Line2, Point2>(linesData.myLines, sortedLinesIndices,
					projectedPoint, 0, sortedLinesIndices.size());

				const auto isPointAlongLine = firstUpperLineIdx < sortedLinesIndices.size() &&
					linesData.myLines[sortedLinesIndices[firstUpperLineIdx]].has_on(projectedPoint);

				const int bucketIndex = isPointAlongLine ? firstUpperLineIdx : (firstUpperLineIdx + linesData.myLines.size());

				if (result.mySortedPlanesIndices.find(bucketIndex) == result.mySortedPlanesIndices.end())
				{
					const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
						Line3{ queryPoint, Vec3{0,0,1} });
					result.mySortedPlanesIndices.insert(std::make_pair(bucketIndex, sortedPlanesIndices));
				}

				const auto& sortedPlanesIndices = result.mySortedPlanesIndices.at(bucketIndex);
				const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
					queryPoint, 0, sortedPlanesIndices.size());
				result.myRangeWrappers.push_back({ Range{firstUpperPlaneIdx, somePlanes.size()}, bucketIndex });
			}
			return result;
		}

		// Slab-based approach
		const auto partition = locComputeSlabs(linesData);

		int cacheUsageCount = 0;

		BatchPointResult result;

		int counter = 0;

		for (const auto& queryPoint : somePoints)
		{
			const Point2 projectedPoint{ queryPoint.x(), queryPoint.y() };
			int slabIdx = FindSlabIndex(linesData.myUniqueAndSortedVerticesX, projectedPoint, 0, linesData.myUniqueAndSortedVerticesX.size());
			int globalBucketIndex = -1;

			// Found slab idx and it is not -inf
			if (slabIdx > -1)
			{
				CGAL_precondition(linesData.myUniqueAndSortedVerticesX[slabIdx] <= projectedPoint.x());
				// Check if point is along vertical slab line
				if (linesData.myUniqueAndSortedVerticesX[slabIdx] == projectedPoint.x())
				{
					const auto vLineIter = partition.myVerticalLineIndex.find(slabIdx);
					if (vLineIter != partition.myVerticalLineIndex.end())
					{
						const int vLineIdx = vLineIter->second;
						const auto lineSearchResult = linesData.FindSegmentOrVertexIndex(projectedPoint, vLineIdx);
						if (!lineSearchResult.myIsVertex && lineSearchResult.myIndex > -1)
						{
							// Compute global index
							globalBucketIndex = linesData.myUniqueVerticesCount + // Vertices count
								linesData.mySegmentsCount[vLineIdx] + lineSearchResult.myIndex;
						}
					}
				}
				slabIdx += 1;
			}
			else
			{
				slabIdx = 0;
			}

			const auto& sortedLines = partition.mySortedLines[slabIdx];
			const int lineIdx = BinarySearch(linesData.myLines, sortedLines, projectedPoint, 0, sortedLines.size());

			if (lineIdx < sortedLines.size() && linesData.myLines[sortedLines[lineIdx]].has_on(projectedPoint))
			{
				const auto lineSearchResult = linesData.FindSegmentOrVertexIndex(projectedPoint, sortedLines[lineIdx]);
				CGAL_precondition(lineSearchResult.myIndex > -1);

				if (lineSearchResult.myIsVertex)
				{
					// Compute global index
					globalBucketIndex = lineSearchResult.myIndex;
					//CGAL_precondition(lineSearchResult.myIndex < linesData.myUniqueVerticesCount);
				}
				else
				{
					globalBucketIndex = linesData.myUniqueVerticesCount + // Vertices count
						linesData.mySegmentsCount[sortedLines[lineIdx]] + lineSearchResult.myIndex;
				}
			}

			if (globalBucketIndex == -1)
			{
				globalBucketIndex = linesData.myUniqueVerticesCount + // Vertices count
					linesData.mySegmentsCount.back() + // Segments count
					lineIdx + slabIdx * (linesData.myLines.size() + 1); // Faces count
			}

			if (result.mySortedPlanesIndices.find(globalBucketIndex) == result.mySortedPlanesIndices.end())
			{
				const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor>(somePlanes,
					Line3{ queryPoint, Vec3{0,0,1} });
				result.mySortedPlanesIndices.insert(std::make_pair(globalBucketIndex, sortedPlanesIndices));
			}
			else
			{
				// Temporary
				cacheUsageCount++;
			}

			const auto& sortedPlanesIndices = result.mySortedPlanesIndices.at(globalBucketIndex);
			//CGAL_precondition(SameLists(sortedPlanesIndices, sortedPlanesIndices_cache));
			const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
				queryPoint, 0, sortedPlanesIndices.size());
			result.myRangeWrappers.push_back({ Range{firstUpperPlaneIdx, somePlanes.size()}, globalBucketIndex });

			// DEBUG ONLY
			counter++;
		}

		std::cout << "Cache used: " << cacheUsageCount << " times";
		return result;
	}

}