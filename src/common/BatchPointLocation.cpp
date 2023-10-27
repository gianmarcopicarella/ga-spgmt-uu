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
			std::vector<std::tuple<Point2, size_t>> myVertices;
			std::vector<Line2> myLines;
			std::optional<BatchPointResult> myBatchPointResultOpt;
		};

		void locComputeData(
			const std::vector<Plane>& somePlanes,
			const std::vector<Point3>& somePoints,
			BatchPointLocationData& anOutData,
			LineData& anOutOldData)
		{
			constexpr auto policy = ExecutionPolicy::PAR_UNSEQ;
			// Compute all 3D lines and map to 2D
			using LineWrapper = std::tuple<STATUS, Line2>;
			const auto maxLinesCount = somePlanes.size() * somePlanes.size();

			std::vector<LineWrapper> lines(maxLinesCount);

			hpx::for_loop(hpx::execution::par_unseq, 0, somePlanes.size(), [&](size_t aPlaneIdx) {
				hpx::for_loop(hpx::execution::par_unseq, aPlaneIdx + 1, somePlanes.size(), [&](size_t anotherPlaneIdx) {
					const auto intersection = CGAL::intersection(somePlanes[aPlaneIdx], somePlanes[anotherPlaneIdx]);
					if (intersection)
					{
						const Line3* line = boost::get<Line3>(&*intersection);
						CGAL_precondition(line != nullptr);
						const auto outIdx = aPlaneIdx * somePlanes.size() + anotherPlaneIdx;

						auto start = line->point();
						auto end = start + line->to_vector() * FT(1000.f);
						if (end < start) std::swap(start, end);

						lines[outIdx] = std::make_tuple(STATUS::INIT, Line2{
							Point2{ start.x(), start.y() },
							Point2{ end.x(), end.y() } });
					}
					});
				});

			hpx::sort(hpx::execution::par_unseq, lines.begin(), lines.end(),
				[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
			Utils::TrimToLastValidItem(lines);
			const auto linesEndIter = hpx::unique(hpx::execution::par_unseq, lines.begin(), lines.end(),
				[](const auto& a, const auto& b) { return std::get<1>(a) == std::get<1>(b); });
			lines.erase(linesEndIter, lines.end());

			// If all the lines are parallel then compute result and return early
			const auto areLinesParallel = Utils::AreItemsParallel<policy>(lines,
				[](const auto& aTuple) { return std::get<1>(aTuple).direction(); });

			if (areLinesParallel)
			{
				BatchPointResult result;

				// Compute sorted list of lines indices
				std::vector<Line2> plainLines(lines.size());
				hpx::transform(hpx::execution::par_unseq, lines.begin(), lines.end(), plainLines.begin(),
					[](const auto& a) { return std::get<1>(a); });


				const auto line = plainLines[0].perpendicular(plainLines[0].point());
				const auto& sortedLinesIndices =
					locGetSortedItemsIndicesFromIntersections<LineParallelLineVisitor<policy>>(plainLines, line);

				std::mutex planesIndicesMutex;

				result.myRangeWrappers.resize(somePoints.size());

				hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {

					const Point2 projectedPoint{ somePoints[aPointIdx].x(), somePoints[aPointIdx].y() };
					const auto firstUpperLineIdx = BinarySearch<Line2, Point2>(plainLines, sortedLinesIndices,
						projectedPoint, 0, sortedLinesIndices.size());

					const auto isPointAlongLine = firstUpperLineIdx < sortedLinesIndices.size() &&
						plainLines[sortedLinesIndices[firstUpperLineIdx]].has_on(projectedPoint);

					const int bucketIndex = isPointAlongLine ? firstUpperLineIdx : (firstUpperLineIdx + plainLines.size());

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
			const auto maxVerticesCount = lines.size() * lines.size();
			using VertexWrapper = std::tuple<STATUS, size_t, Point2, size_t>;
			std::vector<VertexWrapper> vertices(maxVerticesCount);

			hpx::for_loop(hpx::execution::par_unseq, 0, lines.size(), [&](size_t aLineIdx) {
				hpx::for_loop(hpx::execution::par_unseq, aLineIdx + 1, lines.size(), [&](size_t anotherLineIdx) {
					const auto intersection = CGAL::intersection(std::get<1>(lines[aLineIdx]), std::get<1>(lines[anotherLineIdx]));
					if (intersection)
					{
						const Point2* point = boost::get<Point2>(&*intersection);
						CGAL_precondition(point != nullptr);
						const auto firstOutIdx = aLineIdx * lines.size() + anotherLineIdx;
						const auto secondOutIdx = anotherLineIdx * lines.size() + aLineIdx;

						CGAL_precondition(std::get<0>(vertices[firstOutIdx]) == STATUS::NONE);
						CGAL_precondition(std::get<0>(vertices[secondOutIdx]) == STATUS::NONE);

						vertices[firstOutIdx] = std::make_tuple(STATUS::INIT, aLineIdx, *point, 0);
						vertices[secondOutIdx] = std::make_tuple(STATUS::INIT, anotherLineIdx, *point, 0);
					}
					});
				});

			hpx::sort(hpx::execution::par_unseq, vertices.begin(), vertices.end(),
				[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
			Utils::TrimToLastValidItem(vertices);
			const auto verticesEndIter = hpx::unique(hpx::execution::par_unseq, vertices.begin(), vertices.end(),
				[](const auto& a, const auto& b) { return std::get<2>(a) == std::get<2>(b); });
			vertices.erase(verticesEndIter, vertices.end());

			// STORE IN NEW DATA
			anOutData.myVertices.resize(vertices.size());
			std::transform(vertices.begin(), vertices.end(), anOutData.myVertices.begin(),
				[](const auto& a) { return std::make_tuple(std::get<2>(a), std::get<1>(a)); });
			std::sort(anOutData.myVertices.begin(), anOutData.myVertices.end());

			anOutData.myLines.resize(lines.size());
			std::transform(lines.begin(), lines.end(), anOutData.myLines.begin(),
				[](const auto& a) { return std::get<1>(a); });

			// TEMPORARY (JUST FOR COMPATIBILITY)
			// Build back to old data format

			constexpr auto pairComparator = [](auto& aFirst, auto& aSecond) { return aFirst.first < aSecond.first; };

			anOutOldData.myLines.resize(lines.size());
			std::transform(lines.begin(), lines.end(), anOutOldData.myLines.begin(), [](const auto& a) { return std::get<1>(a); });

			anOutOldData.mySortedVertices.resize(lines.size());
			anOutOldData.myUniqueVerticesCount = vertices.size();

			// Giving unique id to each vertex
			hpx::sort(hpx::execution::par_unseq, vertices.begin(), vertices.end(),
				[](const auto& a, const auto& b) { return std::get<2>(a) < std::get<2>(b); });

			// Sequentially assign an id to every vertex
			std::get<3>(vertices[0]) = 0;
			for (size_t i = 1, id = 0; i < vertices.size(); ++i)
			{
				if (std::get<2>(vertices[i - 1]) < std::get<2>(vertices[i]))
				{
					++id;
				}
				std::get<3>(vertices[i]) = id;
			}

			for (const auto& t : vertices)
			{
				const auto& getter = anOutOldData.GetGetter(std::get<1>(lines[std::get<1>(t)]));
				anOutOldData.mySortedVertices[std::get<1>(t)].push_back(std::make_pair(getter(std::get<2>(t)), std::get<3>(t)));

				anOutOldData.myUniqueAndSortedVerticesX.push_back(std::get<2>(t).x());
			}

			std::sort(anOutOldData.myUniqueAndSortedVerticesX.begin(), anOutOldData.myUniqueAndSortedVerticesX.end(), CGAL::Less<FT, FT>());

			anOutOldData.mySegmentsCount.push_back(0);

			for (int i = 1; i < anOutOldData.mySortedVertices.size(); ++i)
			{
				std::sort(anOutOldData.mySortedVertices[i].begin(), anOutOldData.mySortedVertices[i].end(), pairComparator);
				anOutOldData.mySegmentsCount.push_back(anOutOldData.mySortedVertices[i - 1].size() + 1 + anOutOldData.mySegmentsCount.back());
			}
			anOutOldData.mySegmentsCount.push_back(anOutOldData.mySortedVertices.back().size() + 1 + anOutOldData.mySegmentsCount.back());

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
			void locSortLines(const FT& aSlabX, std::vector<Line2>& someOutLines, std::optional<size_t>& anOutVerticalLineIdxOpt)
			{
				using LineWrapper = std::tuple<STATUS, FT, FT, size_t>;
				std::vector<LineWrapper> lines(someOutLines.size());

				hpx::for_loop(hpx::execution::par, 0, someOutLines.size(), [&](size_t aLineIdx) {
					if (someOutLines[aLineIdx].is_vertical() &&
						someOutLines[aLineIdx].point().x() == aSlabX)
					{
						anOutVerticalLineIdxOpt = aLineIdx;
					}
					else
					{
						const auto Y = someOutLines[aLineIdx].y_at_x(aSlabX);
						const auto& dir = someOutLines[aLineIdx].direction();
						const auto slope = dir.dy() / dir.dx();
						lines[aLineIdx] = std::make_tuple(STATUS::INIT, Y, slope, aLineIdx);
					}
					});

				hpx::sort(hpx::execution::par, lines.begin(), lines.end());
				Utils::TrimToLastValidItem(lines);

				std::vector<Line2> sortedResult(lines.size());
				hpx::transform(hpx::execution::par, lines.begin(), lines.end(), sortedResult.begin(),
					[&](const auto& a) { return someOutLines[std::get<3>(a)]; });
				someOutLines = sortedResult;
			}


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
							partition.myVerticalLineIndex.insert(std::make_pair(aSlabIdx, k));
						}
						else
						{
							const auto& dir = currentLine.direction();
							const auto slope = dir.dy() / dir.dx();
							const auto Y = currentLine.y_at_x(currentX);
							sortedLines[atomicCounter.fetch_add(1)] = LineInfo{ Y, slope, k };
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

			void locComputeUniqueLines(
				const std::vector<Plane>& somePlanes,
				std::vector<Line3>& someOutLines)
			{
				using LineWrapper = std::tuple<STATUS, Line3>;
				const auto maxLinesCount = somePlanes.size() * somePlanes.size();
				std::vector<LineWrapper> lines(maxLinesCount);

				hpx::for_loop(hpx::execution::par_unseq, 0, somePlanes.size(), [&](size_t aPlaneIdx) {
					hpx::for_loop(hpx::execution::par_unseq, aPlaneIdx + 1, somePlanes.size(), [&](size_t anotherPlaneIdx) {
						const auto intersection = CGAL::intersection(somePlanes[aPlaneIdx], somePlanes[anotherPlaneIdx]);
						if (intersection)
						{
							const Line3* line = boost::get<Line3>(&*intersection);
							CGAL_precondition(line != nullptr);
							const auto outIdx = aPlaneIdx * somePlanes.size() + anotherPlaneIdx;
							lines[outIdx] = std::make_tuple(STATUS::INIT, *line);
						}
						});
					});

				hpx::sort(hpx::execution::par_unseq, lines.begin(), lines.end(),
					[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
				Utils::TrimToLastValidItem(lines);
				someOutLines.resize(lines.size());
				hpx::transform(hpx::execution::par_unseq, lines.begin(), lines.end(), someOutLines.begin(),
					[](const auto& anItem) { return std::get<1>(anItem); });
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

	void locComputeSlabData(
		const FT aCurrentX,
		const FT aNextX,
		const std::vector<Point2>& someVertices,
		const std::vector<Plane>& somePlanes,
		const std::vector<Line2>& someSortedLines,
		const bool anIsVerticalLineOnSlab,
		std::vector<std::vector<size_t>>& someOutData)
	{
		const auto cacheItemsCount = someVertices.size() + 2 * someSortedLines.size() + 1 +
			static_cast<int>(anIsVerticalLineOnSlab) * (someVertices.size() + 1);
		someOutData.resize(cacheItemsCount);

		const auto midX = (aCurrentX + aNextX) * FT(0.5f);
		size_t outIdx = 0;

		// Vertices
		for (size_t i = 0; i < someVertices.size(); ++i)
		{
			const Point2& point = someVertices[i];
			someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
				Line3{ Point3 {point.x(), point.y(), 0}, Vec3{0,0,1} });
		}

		// Lines
		for (size_t i = 0; i < someSortedLines.size(); ++i)
		{
			const Point3 sampleAlongLine{ midX, someSortedLines[i].y_at_x(midX), 0 };
			someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
				Line3{ sampleAlongLine, Vec3{0,0,1} });
		}

		// First area (Bottom-most one)
		{
			const Point3 sampleAlongLine{ midX, someSortedLines[0].y_at_x(aCurrentX) - FT(1000.f), 0 };
			someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
				Line3{ sampleAlongLine, Vec3{0,0,1} });
		}

		for (size_t i = 1; i < someSortedLines.size(); ++i)
		{
			const auto areaMidY = (someSortedLines[i - 1].y_at_x(midX) + someSortedLines[i].y_at_x(midX)) * FT(0.5f);
			const Point3 sampleAlongLine{ midX, areaMidY, 0 };
			someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
				Line3{ sampleAlongLine, Vec3{0,0,1} });
		}

		// Last area (Top-most one)
		{
			const Point3 sampleAlongLine{ midX, someSortedLines.back().y_at_x(aCurrentX) + FT(1000.f), 0 };
			someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
				Line3{ sampleAlongLine, Vec3{0,0,1} });
		}

		if (anIsVerticalLineOnSlab)
		{
			// First segment (Bottom-most one)
			{
				const Point3 point{ aCurrentX, someVertices[0].y() - FT(1000.f), 0.f };
				someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
					Line3{ Point3 {point.x(), point.y(), 0}, Vec3{0,0,1} });
			}

			for (size_t i = 1; i < someVertices.size(); ++i)
			{
				const auto midY = (someVertices[i].y() + someVertices[i - 1].y()) * FT(0.5f);
				const Point3 point{ aCurrentX, midY, 0.f };
				someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
					Line3{ Point3 {point.x(), point.y(), 0}, Vec3{0,0,1} });
			}

			// Last segment (Bottom-most one)
			{
				const Point3 point{ aCurrentX, someVertices.back().y() + FT(1000.f), 0.f };
				someOutData[outIdx++] = locGetSortedItemsIndicesFromIntersections<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
					Line3{ Point3 {point.x(), point.y(), 0}, Vec3{0,0,1} });
			}
		}

		CGAL_postcondition(outIdx == cacheItemsCount);
	}

	std::tuple<size_t, bool> locFindSegmentOrVertexIndex(const std::vector<Point2>& someVertices, const Point2& aQuery)
	{
		CGAL_precondition(someVertices.size() > 0);
		int left = 0, right = someVertices.size();
		while (left < right)
		{
			const int mid = left + (right - left) / 2;
			if (aQuery.y() < someVertices[mid].y())
			{
				if (mid == 0 || aQuery.y() < someVertices[mid - 1].y())
				{
					return std::make_tuple(mid, false);
				}
				else
				{
					right = mid;
				}
			}
			else if (aQuery.y() > someVertices[mid].y())
			{
				if (mid == someVertices.size() - 1 ||
					aQuery.y() < someVertices[mid + 1].y())
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

		constexpr auto policy = ExecutionPolicy::PAR;

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

		LineData linesData;

		locComputeData(somePlanes, somePoints, newData, linesData);

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
			const auto fakeStart = std::make_tuple(std::get<0>(newData.myVertices[0]) - offset, 0);
			const auto fakeEnd = std::make_tuple(Point2{ queryPoints.back().x(), queryPoints.back().y() } + offset, 0);

			newData.myVertices.insert(newData.myVertices.begin(), fakeStart);
			newData.myVertices.emplace_back(fakeEnd);
		}

		// Pre-allocate cache
		result.mySortedPlanesCache.resize(newData.myVertices.size() - 1);

		std::vector<Point2> vertices;
		size_t i = 0, k = 0;
		while (i < newData.myVertices.size() - 1 && k < queryPoints.size())
		{
			vertices.clear();
			const auto slabX = std::get<0>(newData.myVertices[i]).x();
			
			{
				vertices.emplace_back(std::get<0>(newData.myVertices[i]));
				while (++i < newData.myVertices.size() && slabX == std::get<0>(newData.myVertices[i]).x())
				{
					if (vertices.back().y() != std::get<0>(newData.myVertices[i]).y())
					{
						vertices.emplace_back(std::get<0>(newData.myVertices[i]));
					}
				}
			}

			const auto nextX = std::get<0>(newData.myVertices[i]).x();

			const auto isThereAtLeastOneQuery = queryPoints[k].x() < nextX;
			if (isThereAtLeastOneQuery)
			{
				const auto claimedIdx = cacheIndex.fetch_add(1, std::memory_order_relaxed);
				std::vector<Line2> sortedLines(newData.myLines);
				std::optional<size_t> verticalLineIdx;
				locSortLines(slabX, sortedLines, verticalLineIdx);
				locComputeSlabData(slabX, nextX, vertices, somePlanes, sortedLines, verticalLineIdx.has_value(), result.mySortedPlanesCache[claimedIdx]);

				const auto& cache = result.mySortedPlanesCache[claimedIdx];

				//std::cout << "sLAB x: " << slabX << std::endl;

				while (k < queryPoints.size() && queryPoints[k].x() < nextX)
				{
					const Point2 projectedPoint{ queryPoints[k].x(), queryPoints[k].y() };

					//std::cout << "qpoint: " << queryPoints[k] << std::endl;

					auto cacheIdx = -1;

					if (projectedPoint.x() == slabX)
					{
						const auto test = locFindSegmentOrVertexIndex(vertices, projectedPoint);
						if (std::get<1>(test))
						{
							// Vertex
							cacheIdx = std::get<0>(test);
						}
						else if (verticalLineIdx.has_value())
						{
							// Vertical segment
							cacheIdx = cache.size() - (vertices.size() + 1) + std::get<0>(test);
						}
					}
					else
					{
						const auto lineIdx = BinarySearch(sortedLines, projectedPoint, 0, sortedLines.size());
						if (lineIdx < sortedLines.size() && sortedLines[lineIdx].has_on(projectedPoint))
						{
							// Line
							cacheIdx = vertices.size() + lineIdx;
						}
						else
						{
							// Area
							cacheIdx = vertices.size() + sortedLines.size() + lineIdx;
						}
					}

					CGAL_postcondition((cacheIdx > -1 && cacheIdx < cache.size()));

					const auto& sortedPlanesIndices = cache[cacheIdx];
					const auto firstUpperPlaneIdx = BinarySearch<Plane, Point3>(somePlanes, sortedPlanesIndices,
						queryPoints[k], 0, sortedPlanesIndices.size());

					using Res = BatchPointResult::NewRangeWrapper;
					result.myNewRangeWrappers.emplace_back(Res{ Range{firstUpperPlaneIdx, somePlanes.size()}, claimedIdx, static_cast<size_t>(cacheIdx) });

					++k;
				}
			}
		}

		// remove 2 fake points
		{
			newData.myVertices.erase(newData.myVertices.begin(), newData.myVertices.begin()+1);
			newData.myVertices.pop_back();
		}

#endif

		// Slab-based approach
		const auto partition = Parallel::locComputeSlabs(linesData);

		result.myRangeWrappers.resize(somePoints.size());

		std::mutex planesIndicesMutex;

		for(size_t aPointIdx = 0; aPointIdx < somePoints.size(); ++aPointIdx)
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