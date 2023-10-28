#include "BatchPointLocation.h"
#include "../common/Utils.h"

#include <vector>

namespace SPGMT
{
	namespace
	{
		template<ExecutionPolicy E = ExecutionPolicy::SEQ>
		struct LineParallelLineVisitor
		{
			constexpr static ExecutionPolicy GetExecutionPolicy()
			{
				return E;
			}

			struct Data
			{
				FT myDistance;
				size_t myIndex;
			};
			typedef Data payload_type;
			typedef void result_type;
			void operator()(const Point2& aPoint)
			{
				const FT zero{ 0 };
				const Vec2 vecToPoint{ myLine.point(), aPoint };
				const auto distanceSign = CGAL::scalar_product(vecToPoint, myLine.to_vector()) > zero ? 1.f : -1.f;

				if (E == ExecutionPolicy::SEQ)
				{
					myOutResult.emplace_back(Data{ distanceSign * vecToPoint.squared_length(), myPlaneIndex });
				}
				else
				{
					myOutResult[myPlaneIndex] = Data{ distanceSign * vecToPoint.squared_length(), myPlaneIndex };
				}
			}
			void operator()(const Line2& /**/)
			{
				CGAL_precondition(false);
			}

			const Line2& myLine;
			const size_t myPlaneIndex;
			std::vector<payload_type>& myOutResult;
		};

		template<ExecutionPolicy E = ExecutionPolicy::SEQ>
		struct LinePlaneVisitor
		{
			constexpr static ExecutionPolicy GetExecutionPolicy()
			{
				return E;
			}

			struct Data
			{
				FT myDistance;
				size_t myIndex;
			};
			typedef Data payload_type;
			typedef void result_type;
			void operator()(const Point3& aPoint)
			{
				const FT zero{ 0 };
				const Vec3 vecToPoint{ myLine.point(), aPoint };
				const auto distanceSign = CGAL::scalar_product(vecToPoint, myLine.to_vector()) > zero ? 1.f : -1.f;

				if (E == ExecutionPolicy::SEQ)
				{
					myOutResult.emplace_back(Data{ distanceSign * vecToPoint.squared_length(), myPlaneIndex });
				}
				else
				{
					myOutResult[myPlaneIndex] = Data{ distanceSign * vecToPoint.squared_length(), myPlaneIndex };
				}

			}
			void operator()(const Line3& /**/)
			{
				CGAL_precondition(false);
			}

			const Line3& myLine;
			const size_t myPlaneIndex;
			std::vector<payload_type>& myOutResult;
		};

		template<ExecutionPolicy E, typename T>
		void locSortByDistanceMember(std::vector<T>& someOutData)
		{
			const auto sort = [](const auto& aFirst, const auto& aSecond) {return aFirst.myDistance < aSecond.myDistance; };
			BindExecutionPolicy<E>(hpx::sort, someOutData.begin(), someOutData.end(), sort);
		}
		template<ExecutionPolicy E, typename T>
		std::vector<size_t> locExtractIndicesMember(const std::vector<T>& someData)
		{
			const auto transform = [](const auto& aFirst) {return aFirst.myIndex; };
			std::vector<size_t> indices(someData.size());
			BindExecutionPolicy<E>(hpx::transform, someData.begin(), someData.end(), indices.begin(), transform);
			return indices;
		}
		template <typename Strategy, typename T, typename K>
		std::vector<size_t> locGetSortedItemsIndicesFromIntersections(const std::vector<T>& someItems, const K& anotherItem)
		{
			auto planesLineIntersections = Utils::FindIntersections<Strategy>(someItems, anotherItem);
			locSortByDistanceMember<Strategy::GetExecutionPolicy()>(planesLineIntersections);
			return locExtractIndicesMember<Strategy::GetExecutionPolicy()>(planesLineIntersections);
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

			BindExecutionPolicy<E>(hpx::for_loop, 0, someLines.size(), [&](size_t aLineIdx) {
				const auto& translation = someLines[aLineIdx].to_vector() * distance;
				const auto firstPoint{ someLines[aLineIdx].point() };
				const auto secondPoint{ firstPoint + translation };
				result[aLineIdx] = Create2DLine(
					Point2{ firstPoint.x(), firstPoint.y() },
					Point2{ secondPoint.x(), secondPoint.y() });
				});

			return result;
		}
		template<ExecutionPolicy E>
		std::vector<Line2> locProjectOnXYPlaneUnique(const std::vector<Line3>& someLines)
		{
			auto result = locProjectOnXYPlane<E>(someLines);
			auto uniqueIter =
				BindExecutionPolicy<E>(hpx::unique, result.begin(), result.end(), CGAL::Equal_to<Line2, Line2>());
			result.erase(uniqueIter, result.end());
			return result;
		}

		template<typename T, typename K>
		int BinarySearch(const std::vector<T>& someItems, const std::vector<size_t>& someSortedIndices, const K& aPoint, int aMin, int aMax)
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
			size_t myIndex;
		};

		struct SlabPartition
		{
			std::unordered_map<int, int> myVerticalLineIndex;
			std::vector<std::vector<size_t>> mySortedLines; // N+1 values ( considering -inf partition)
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

		// new data
		struct BatchPointLocationData
		{
			std::vector<Point2> myVertices;
			std::vector<Line2> myLines;
			size_t myNonVerticalLinesCount;
			std::optional<BatchPointResult> myBatchPointResultOpt;
		};

		template<ExecutionPolicy E>
		std::vector<Line2> locSortLines(
			const FT& aSlabX,
			const std::vector<Line2>& someLines,
			const size_t aNonVerticalLinesCount,
			bool& anOutHasVerticalLine)
		{

			// The number of vertical lines is, on average, very rare and so it makes sense to run this code sequentially for the sake of efficiency
			anOutHasVerticalLine =
				hpx::find_if(hpx::execution::seq, someLines.begin() + aNonVerticalLinesCount, someLines.end(), [&](const auto& aVerticalLine) {
				return aVerticalLine.point().x() == aSlabX; }) != someLines.end();

				std::vector<Line2> sortedLines(someLines.begin(), someLines.begin() + aNonVerticalLinesCount);
				BindExecutionPolicy<E>(hpx::sort, sortedLines.begin(), sortedLines.end(),
					[&](const auto& a, const auto& b) {
						const auto& ya = a.y_at_x(aSlabX);
						const auto& yb = b.y_at_x(aSlabX);
						return ya < yb || (ya == yb && CGAL::compare_slope(a, b) == CGAL::Comparison_result::SMALLER);
					});

				return sortedLines;
		}

		void locComputeData(
			const std::vector<Plane>& somePlanes,
			const std::vector<Point3>& somePoints,
			BatchPointLocationData& anOutData)
		{
			const static Point3 zero{ 0,0,0 };
			constexpr auto policy = ExecutionPolicy::PAR_UNSEQ;
			// Compute all 3D lines and map to 2D
			const auto maxLinesCount = somePlanes.size() * somePlanes.size() / 2;

			std::atomic<size_t> counter{ 0 };
			std::vector<std::tuple<bool, Point2, Point2>> lines(maxLinesCount);

			hpx::for_loop(hpx::execution::par_unseq, 0, somePlanes.size(), [&](size_t aPlaneIdx) {
				hpx::for_loop(hpx::execution::par_unseq, aPlaneIdx + 1, somePlanes.size(), [&](size_t anotherPlaneIdx) {
					const auto intersection = CGAL::intersection(somePlanes[aPlaneIdx], somePlanes[anotherPlaneIdx]);
					if (intersection)
					{
						const Line3* line = boost::get<Line3>(&*intersection);
						CGAL_precondition(line != nullptr);
						const auto outIdx = counter.fetch_add(1, std::memory_order_relaxed);

						auto start = line->projection(zero);
						auto end = start + line->to_vector() * FT(1000.f);
						if (end < start) std::swap(start, end);

						lines[outIdx] = std::make_tuple(
							start.y() == end.y(),
							Point2{ start.x(), start.y() },
							Point2{ end.x(), end.y() });
					}
					});
				});

			lines.resize(counter.fetch_xor(counter.load()));
			hpx::sort(hpx::execution::par_unseq, lines.begin(), lines.end());

			const auto linesEndIter = hpx::unique(hpx::execution::par_unseq, lines.begin(), lines.end());
			lines.erase(linesEndIter, lines.end());

			anOutData.myNonVerticalLinesCount = hpx::count_if(hpx::execution::par_unseq, lines.begin(), lines.end(),
				[](const auto& a) {return !std::get<0>(a); });

			anOutData.myLines.resize(lines.size());
			hpx::transform(hpx::execution::par_unseq, lines.begin(), lines.end(), anOutData.myLines.begin(),
				[](const auto& a) { return Line2{ std::get<1>(a), std::get<2>(a) }; });

			// If all the lines are parallel then compute result and return early
			const auto areLinesParallel = Utils::AreItemsParallel<policy>(anOutData.myLines,
				[](const auto& aLine) { return aLine.direction(); });

			if (areLinesParallel)
			{
				BatchPointResult result;

				const auto line = anOutData.myLines[0].perpendicular(anOutData.myLines[0].point());
				const auto& sortedLinesIndices =
					locGetSortedItemsIndicesFromIntersections<LineParallelLineVisitor<policy>>(anOutData.myLines, line);

				std::mutex planesIndicesMutex;

				result.myRangeWrappers.resize(somePoints.size());

				hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {

					const Point2 projectedPoint{ somePoints[aPointIdx].x(), somePoints[aPointIdx].y() };
					const auto firstUpperLineIdx = BinarySearch<Line2, Point2>(anOutData.myLines, sortedLinesIndices,
						projectedPoint, 0, sortedLinesIndices.size());

					const auto isPointAlongLine = firstUpperLineIdx < sortedLinesIndices.size() &&
						anOutData.myLines[sortedLinesIndices[firstUpperLineIdx]].has_on(projectedPoint);

					const int bucketIndex = isPointAlongLine ? firstUpperLineIdx : (firstUpperLineIdx + anOutData.myLines.size());

					std::vector<size_t> sortedPlanesIndices;

					{
						std::lock_guard guard(planesIndicesMutex);
						if (result.mySortedPlanesIndices.find(bucketIndex) == result.mySortedPlanesIndices.end())
						{
							const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
								Line3{ somePoints[aPointIdx], Vec3{0,0,1} });
							result.mySortedPlanesIndices.insert(std::make_pair(bucketIndex, sortedPlanesIndices));
						}

						sortedPlanesIndices = result.mySortedPlanesIndices.at(bucketIndex);
					}

					const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
						somePoints[aPointIdx], 0, sortedPlanesIndices.size());
					result.myRangeWrappers[aPointIdx] = { Range{firstUpperPlaneIdx, somePlanes.size()}, bucketIndex };

					});

				anOutData.myBatchPointResultOpt = result;
				return;
			}

			// Compute all intersection points
			const auto maxVerticesCount = anOutData.myLines.size() * (anOutData.myLines.size() - 1) / 2;
			anOutData.myVertices.resize(maxVerticesCount);

			hpx::for_loop(hpx::execution::par_unseq, 0, anOutData.myLines.size(), [&](size_t aLineIdx) {
				hpx::for_loop(hpx::execution::par_unseq, aLineIdx + 1, anOutData.myLines.size(), [&](size_t anotherLineIdx) {
					const auto intersection = CGAL::intersection(anOutData.myLines[aLineIdx], anOutData.myLines[anotherLineIdx]);
					if (intersection)
					{
						const Point2* point = boost::get<Point2>(&*intersection);
						CGAL_precondition(point != nullptr);
						const auto outIdx = counter.fetch_add(1, std::memory_order_relaxed);
						anOutData.myVertices[outIdx] = *point;
					}
					});
				});

			anOutData.myVertices.resize(counter.load());
			hpx::sort(hpx::execution::par_unseq, anOutData.myVertices.begin(), anOutData.myVertices.end(), CGAL::Less<Point2, Point2>());

			// This line is the bottleneck (memory)
			const auto verticesEndIter = hpx::unique(hpx::execution::par_unseq, anOutData.myVertices.begin(),
				anOutData.myVertices.end(), CGAL::Equal_to<Point2, Point2>());

			anOutData.myVertices.erase(verticesEndIter, anOutData.myVertices.end());
		}

		LineData locComputeLinesWrapper(const std::vector<Line3>& someLines)
		{


			struct PointIndexPair
			{
				Point2 myKey;
				int myIndex{ -1 };
			};

			LineData data;

			data.myLines = locProjectOnXYPlaneUnique<ExecutionPolicy::SEQ>(someLines);
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

			for (size_t slabIdx = 0; slabIdx < aLineData.myUniqueAndSortedVerticesX.size(); ++slabIdx)
			{
				const auto& currentX = aLineData.myUniqueAndSortedVerticesX[slabIdx];
				std::vector<LineInfo> sortedLines;

				for (size_t k = 0; k < aLineData.myLines.size(); ++k)
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

				std::vector<size_t> sortedIndices;
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

	std::tuple<size_t, bool> locFindSegmentOrVertexIndex(
		const std::vector<Point2>::const_iterator aVerticesStart,
		const std::vector<Point2>::const_iterator aVerticesEnd,
		const Point2& aQuery)
	{
		const auto verticesCount = std::distance(aVerticesStart, aVerticesEnd);
		CGAL_precondition(verticesCount > 0);
		int left = 0, right = verticesCount;
		while (left < right)
		{
			const size_t mid = left + (right - left) / 2;
			if (aQuery.y() < (aVerticesStart + mid)->y())
			{
				if (mid == 0 || aQuery.y() < (aVerticesStart + mid - 1)->y())
				{
					return std::make_tuple(mid, false);
				}
				else
				{
					right = mid;
				}
			}
			else if (aQuery.y() > (aVerticesStart + mid)->y())
			{
				if (mid == verticesCount - 1 ||
					aQuery.y() < (aVerticesStart + mid + 1)->y())
				{
					return std::make_tuple(mid + 1, false);
				}
				else
				{
					left = mid + 1;
				}
			}
			else
			{
				return std::make_tuple(mid, true);
			}
		}

		return std::make_tuple(-1, false);
	}

	// Skeleton
	BatchPointResult ParallelBatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints)
	{
		using namespace Parallel;

		constexpr auto policy = ExecutionPolicy::PAR_UNSEQ;

		if (somePlanes.empty() || somePoints.empty())
		{
			return BatchPointResult{};
		}

		// Requirements check (TEMPORARY DISABLED)
		// CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		// CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));
		// CGAL_precondition(Utils::ArePlanesUniformlyOriented(somePlanes));

		// 1st edge case: single or parallel planes
		const auto arePlanesParallel = Utils::AreItemsParallel<policy>(somePlanes,
			[](const Plane& aPlane) { return aPlane.orthogonal_direction(); });

		if (arePlanesParallel)
		{
			BatchPointResult result;

			const auto line = somePlanes[0].perpendicular_line(somePlanes[0].point());
			const auto& sortedPlanesIndices =
				locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<policy>>(somePlanes, line);

			result.mySortedPlanesIndices.insert(std::make_pair(0, sortedPlanesIndices));
			result.myRangeWrappers.resize(somePoints.size());

			hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {
				const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
				somePoints[aPointIdx], 0, sortedPlanesIndices.size());
			result.myRangeWrappers[aPointIdx] = { Range{firstUpperPlaneIdx, somePlanes.size()}, 0 };
				});

			return result;
		}

		BatchPointLocationData newData;

		locComputeData(somePlanes, somePoints, newData);

		if (newData.myBatchPointResultOpt.has_value())
		{
			return newData.myBatchPointResultOpt.value();
		}

		BatchPointResult result;

		// New method

		// sort query points
#if 1
		std::atomic<size_t> cacheIndex{ 0 };

		std::vector<Point3> queryPoints(somePoints.size());
		std::copy(somePoints.begin(), somePoints.end(), queryPoints.begin());
		std::sort(queryPoints.begin(), queryPoints.end());

		// Add 2 fake points
		// 1) is equal to vertex[0] - offset
		// 2) queryPoint.last() + offset
		{
			const Vec2 offset{ 1000.f, 0.f };
			const auto fakeStart = newData.myVertices[0] - offset;
			const auto fakeEnd = Point2{ queryPoints.back().x(), queryPoints.back().y() } + offset;

			newData.myVertices.insert(newData.myVertices.begin(), fakeStart);
			newData.myVertices.emplace_back(fakeEnd);
		}

		// Pre-allocate cache
		result.mySortedPlanesCache.resize(newData.myVertices.size() - 1);

		size_t i = 0, k = 0;
		while (i < newData.myVertices.size() - 1 && k < queryPoints.size())
		{
			const auto& slabX = newData.myVertices[i].x();

			const auto startVerticesIter = newData.myVertices.begin() + i;
			while (slabX == newData.myVertices[++i].x());

			const auto nextX = newData.myVertices[i].x();

			if (queryPoints[k].x() < nextX)
			{
				const auto endVerticesIter = newData.myVertices.begin() + i;
				const auto verticesCount = std::distance(startVerticesIter, endVerticesIter);
				const auto claimedIdx = cacheIndex.fetch_add(1, std::memory_order_relaxed);
				bool hasVerticalLine{ false };
				const auto& sortedLines = locSortLines<policy>(slabX, newData.myLines, newData.myNonVerticalLinesCount, hasVerticalLine);
				auto& cache = result.mySortedPlanesCache[claimedIdx];

				while (k < queryPoints.size() && queryPoints[k].x() < nextX)
				{
					const Point2 projectedPoint{ queryPoints[k].x(), queryPoints[k].y() };
					auto cacheIdx = -1;

					if (projectedPoint.x() == slabX)
					{
						const auto test = locFindSegmentOrVertexIndex(startVerticesIter, endVerticesIter, projectedPoint);
						if (std::get<1>(test))
						{
							// Vertex
							cacheIdx = std::get<0>(test);
						}
						else if (hasVerticalLine)
						{
							// Vertical segment
							cacheIdx = cache.size() - (verticesCount + 1) + std::get<0>(test);
						}
					}
					else
					{
						const auto lineIdx = BinarySearch(sortedLines, projectedPoint, 0, sortedLines.size());
						if (lineIdx < sortedLines.size() && sortedLines[lineIdx].has_on(projectedPoint))
						{
							// Line
							cacheIdx = verticesCount + lineIdx;
						}
						else
						{
							// Area
							cacheIdx = verticesCount + sortedLines.size() + lineIdx;
						}
					}

					CGAL_postcondition(cacheIdx > -1);

					if (cache.find(cacheIdx) == cache.end())
					{
						const auto& sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
							Line3{ queryPoints[k], Vec3{0,0,1} });
						cache.insert(std::make_pair(cacheIdx, sortedPlanesIndices));
					}

					const auto& sortedPlanesIndices = cache.at(cacheIdx);
					const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
						queryPoints[k], 0, sortedPlanesIndices.size());

					using Res = BatchPointResult::NewRangeWrapper;
					result.myNewRangeWrappers.emplace_back(Res{ Range{firstUpperPlaneIdx, somePlanes.size()}, claimedIdx, static_cast<size_t>(cacheIdx) });

					++k;
				}
			}
		}

		// Resize cache
		result.mySortedPlanesCache.resize(cacheIndex.load());

		// remove 2 fake points
		{
			newData.myVertices.erase(newData.myVertices.begin(), newData.myVertices.begin() + 1);
			newData.myVertices.pop_back();
		}

#endif

#if 0
		// Slab-based approach
		const auto partition = Parallel::locComputeSlabs(linesData);

		result.myRangeWrappers.resize(somePoints.size());

		std::mutex planesIndicesMutex;

		for (size_t aPointIdx = 0; aPointIdx < somePoints.size(); ++aPointIdx)
		{

			//hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {

			if (aPointIdx == 237)//somePoints[aPointIdx] == Point3{ -114.345f, 25.3222f, 81.2164f })
			{
				std::cout << "same" << std::endl;
			}

			const Point2 projectedPoint{ somePoints[aPointIdx].x(), somePoints[aPointIdx].y() };
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

			std::vector<size_t> sortedPlanesIndices;

			{
				std::lock_guard guard(planesIndicesMutex);
				if (result.mySortedPlanesIndices.find(globalBucketIndex) == result.mySortedPlanesIndices.end())
				{
					const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
						Line3{ somePoints[aPointIdx], Vec3{0,0,1} });
					result.mySortedPlanesIndices.insert(std::make_pair(globalBucketIndex, sortedPlanesIndices));
				}

				sortedPlanesIndices = result.mySortedPlanesIndices.at(globalBucketIndex);
			}

			//CGAL_precondition(SameLists(sortedPlanesIndices, sortedPlanesIndices_cache));
			const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
				somePoints[aPointIdx], 0, sortedPlanesIndices.size());
			result.myRangeWrappers[aPointIdx] = { Range{firstUpperPlaneIdx, somePlanes.size()}, globalBucketIndex };

		}//);
#endif

		//std::cout << "DONE" << std::endl;

		return result;
	}


	BatchPointResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints)
	{
		if (somePlanes.empty() || somePoints.empty())
		{
			return BatchPointResult{};
		}

		// Requirements check (TEMPORARY DISABLED)
		//CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		//CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));
		//CGAL_precondition(Utils::ArePlanesUniformlyOriented(somePlanes));

		// Compute plane-plane intersections
		const auto& lines3d = Utils::FindIntersections<Utils::PlanePlaneVisitor>(somePlanes);

		// 1st edge case: single plane or parallel planes
		if (lines3d.empty())
		{
			const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
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
			const auto sortedLinesIndices = locGetSortedItemsIndicesFromIntersections<LineParallelLineVisitor<ExecutionPolicy::SEQ>>(linesData.myLines,
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
					const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
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
				const auto sortedPlanesIndices = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
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