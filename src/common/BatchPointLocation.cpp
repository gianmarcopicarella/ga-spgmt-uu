#include "BatchPointLocation.h"
#include "../common/Utils.h"

#include <vector>
#include <functional>

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

		template <typename Strategy, typename T, typename K>
		std::vector<size_t> locGetSortedItemsByIntersection(const std::vector<T>& someItems, const K& anotherItem)
		{
			const auto sort = [](const auto& aFirst, const auto& aSecond) {return aFirst.myDistance < aSecond.myDistance; };
			const auto transform = [](const auto& aFirst) {return aFirst.myIndex; };
			constexpr auto E = Strategy::GetExecutionPolicy();

			auto items = Utils::FindIntersections<Strategy>(someItems, anotherItem);
			BindExecutionPolicy<E>(hpx::sort, items.begin(), items.end(), sort);

			std::vector<size_t> indices(items.size());
			BindExecutionPolicy<E>(hpx::transform, items.begin(), items.end(), indices.begin(), transform);

			return indices;
		}

		bool locHasVerticalLine(
			const FT& aSlabX,
			const std::vector<Line2>& someLines,
			const size_t aNonVerticalLinesCount)
		{
			// The number of vertical lines is, on average, very rare and so it makes sense to run this code sequentially for the sake of efficiency
			return hpx::find_if(hpx::execution::seq, someLines.begin() + aNonVerticalLinesCount, someLines.end(), [&](const auto& aVerticalLine) {
				return aVerticalLine.point().x() == aSlabX; }) != someLines.end();
		}

		template<ExecutionPolicy E>
		std::vector<Line2> locGetSortedLines(
			const FT& aSlabX,
			const std::vector<Line2>& someLines,
			const size_t aNonVerticalLinesCount)
		{
			std::vector<Line2> sortedLines(someLines.begin(), someLines.begin() + aNonVerticalLinesCount);
			BindExecutionPolicy<E>(hpx::sort, sortedLines.begin(), sortedLines.end(),
				[&](const auto& a, const auto& b) {
					const auto& ya = a.y_at_x(aSlabX);
					const auto& yb = b.y_at_x(aSlabX);
					return ya < yb || (ya == yb && CGAL::compare_slope(a, b) == CGAL::Comparison_result::SMALLER);
				});

			return sortedLines;
		}

		struct BatchPointLocationData
		{
			std::vector<std::tuple<FT, std::vector<FT>>> mySlabs;
			std::vector<Line2> myLines;
			size_t myNonVerticalLinesCount;
			std::optional<BatchPointResult> myBatchPointResultOpt;
		};

		template<ExecutionPolicy E, typename T, typename F>
		void locFillRequiredCacheEntries(
			const std::vector<Plane>& somePlanes,
			const std::vector<T>& somePoints,
			const size_t aStartIdx,
			const size_t anEndIdx,
			const F& anAreaIndexGetter,
			std::vector<size_t>& someOutAreaIndices,
			std::unordered_map<size_t, std::vector<size_t>>& anOutCache)
		{
			const auto verticesCount = anEndIdx - aStartIdx;
			someOutAreaIndices.resize(verticesCount);
			BindExecutionPolicy<E>(hpx::for_loop, 0, verticesCount, [&](size_t aPointIdx) {
				someOutAreaIndices[aPointIdx] = anAreaIndexGetter(aStartIdx + aPointIdx);
				});

			std::vector<std::tuple<size_t, size_t>> uniqueAreaIndices(someOutAreaIndices.size());
			BindExecutionPolicy<E>(hpx::for_loop, 0, someOutAreaIndices.size(),
				[&](size_t aPointIdx) { uniqueAreaIndices[aPointIdx] = std::make_tuple(someOutAreaIndices[aPointIdx], aStartIdx + aPointIdx); });

			BindExecutionPolicy<E>(hpx::sort, uniqueAreaIndices.begin(), uniqueAreaIndices.end(),
				[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
			const auto itemsEndIter = BindExecutionPolicy<E>(hpx::unique, uniqueAreaIndices.begin(), uniqueAreaIndices.end(),
				[](const auto& a, const auto& b) { return std::get<0>(a) == std::get<0>(b); });
			uniqueAreaIndices.erase(itemsEndIter, uniqueAreaIndices.end());

			// Fill cache (KEEP THIS LOOP SEQUENTIAL)
			hpx::for_loop(hpx::execution::seq, uniqueAreaIndices.begin(), uniqueAreaIndices.end(),
				[&](const auto& anAreaTuple) { anOutCache.insert(std::make_pair(std::get<0>(*anAreaTuple), std::vector<size_t>{})); });

			BindExecutionPolicy<E>(hpx::for_loop, uniqueAreaIndices.begin(), uniqueAreaIndices.end(),
				[&](const auto& anAreaTuple) {
					if constexpr (std::is_same<T, Point3>::value)
					{
						anOutCache.find(std::get<0>(*anAreaTuple))->second =
							locGetSortedItemsByIntersection<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
								Line3{ somePoints[std::get<1>(*anAreaTuple)], Vec3{0,0,1} });
					}
					else
					{
						anOutCache.find(std::get<0>(*anAreaTuple))->second =
							locGetSortedItemsByIntersection<LinePlaneVisitor<ExecutionPolicy::SEQ>>(somePlanes,
								Line3{ std::get<0>(somePoints[std::get<1>(*anAreaTuple)]), Vec3{0,0,1} });
					}
				});
		}

		template<typename K, typename F>
		size_t locGeometricBinarySearch(
			const K& aPoint,
			const size_t anItemsCount,
			F&& anItemGetter)
		{
			CGAL_precondition(anItemsCount > 0);
			int left = 0, right = anItemsCount;
			while (left < right)
			{
				const auto mid = left + (right - left) / 2;
				const auto& item = anItemGetter(mid);
				if (item.has_on_negative_side(aPoint))
				{
					if (mid == 0 || anItemGetter(mid - 1).has_on_positive_side(aPoint))
					{
						return mid;
					}
					else
					{
						right = mid;
					}
				}
				else if (item.has_on_positive_side(aPoint))
				{
					if (mid == anItemsCount - 1 ||
						anItemGetter(mid + 1).has_on_negative_side(aPoint))
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
					if (mid == 0 || (mid > 0 && !anItemGetter(mid - 1).has_on(aPoint)))
					{
						return mid;
					}
					else
					{
						right = mid;
					}
				}
			}

			// MUST NEVER REACH THIS AREA
			CGAL_precondition(false);
			return SIZE_MAX;
		}

		enum class SlabSearchType
		{
			INTERNAL, EXTERNAL
		};

		template<SlabSearchType T, typename F>
		std::tuple<size_t, bool> locSlabBinarySearch(
			const FT& aValue,
			const size_t anItemsCount,
			F&& anItemGetter)
		{
			CGAL_precondition(anItemsCount > 0);
			int left = 0, right = anItemsCount;
			while (left < right)
			{
				const size_t mid = left + (right - left) / 2;
				const auto& item = anItemGetter(mid);
				if (aValue < item)
				{
					if (mid == 0 || aValue > anItemGetter(mid - 1))
					{
						if constexpr (T == SlabSearchType::INTERNAL)
						{
							return std::make_tuple(mid, false);
						}
						else
						{
							return std::make_tuple(mid - 1, false);
						}
					}
					else
					{
						right = mid;
					}
				}
				else if (aValue > item)
				{
					if (mid == anItemsCount - 1 ||
						aValue < anItemGetter(mid + 1))
					{
						if constexpr (T == SlabSearchType::INTERNAL)
						{
							return std::make_tuple(mid + 1, false);
						}
						else
						{
							return std::make_tuple(mid, false);
						}
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

			// It should never reach this area
			CGAL_postcondition(false);
			return std::make_tuple(SIZE_MAX, false);
		}

		template<ExecutionPolicy E>
		void locComputeData(
			const std::vector<Plane>& somePlanes,
			const std::vector<Point3>& somePoints,
			BatchPointLocationData& anOutData)
		{
			const static Point3 zero{ 0,0,0 };

			std::vector<std::tuple<bool, Point2, Point2>> lines;

			const auto PlanePlaneSeq = [&](const size_t aFirstPlaneIdx, const size_t aSecondPlaneIdx)
			{
				const auto intersection = CGAL::intersection(somePlanes[aFirstPlaneIdx], somePlanes[aSecondPlaneIdx]);
				if (intersection)
				{
					const Line3* line = boost::get<Line3>(&*intersection);
					CGAL_precondition(line != nullptr);

					auto start = line->projection(zero);
					auto end = start + line->to_vector() * FT(1000.f);
					if (end < start) std::swap(start, end);

					lines.emplace_back(std::make_tuple(
						start.y() == end.y(), Point2{ start.x(), start.y() }, Point2{ end.x(), end.y() }));
				}
			};
			const auto PlanePlanePar = [&](const size_t aFirstPlaneIdx, const size_t aSecondPlaneIdx, std::atomic<size_t>& anOutCounter)
			{
				const auto intersection = CGAL::intersection(somePlanes[aFirstPlaneIdx], somePlanes[aSecondPlaneIdx]);
				if (intersection)
				{
					const Line3* line = boost::get<Line3>(&*intersection);
					CGAL_precondition(line != nullptr);

					auto start = line->projection(zero);
					auto end = start + line->to_vector() * FT(1000.f);
					if (end < start) std::swap(start, end);

					const auto outIdx = anOutCounter.fetch_add(1, std::memory_order_relaxed);

					lines[outIdx] = std::make_tuple(
						start.y() == end.y(),
						Point2{ start.x(), start.y() },
						Point2{ end.x(), end.y() });
				}
			};

			if constexpr (E == ExecutionPolicy::SEQ)
			{
				hpx::for_loop(hpx::execution::seq, 0, somePlanes.size(), [&](size_t aPlaneIdx) {
					hpx::for_loop(hpx::execution::seq, aPlaneIdx + 1, somePlanes.size(), [&](size_t anotherPlaneIdx) {
						PlanePlaneSeq(aPlaneIdx, anotherPlaneIdx);
						});
					});
			}
			else
			{
				std::atomic<size_t> counter{ 0 };
				const auto maxLinesCount = somePlanes.size() * (somePlanes.size() - 1) / 2;
				lines.resize(maxLinesCount);

				BindExecutionPolicy<E>(hpx::for_loop, 0, somePlanes.size(), [&](size_t aPlaneIdx) {
					BindExecutionPolicy<E>(hpx::for_loop, aPlaneIdx + 1, somePlanes.size(), [&](size_t anotherPlaneIdx) {
						PlanePlanePar(aPlaneIdx, anotherPlaneIdx, counter);
						});
					});

				lines.resize(counter.fetch_xor(counter.load()));
			}

			BindExecutionPolicy<E>(hpx::sort, lines.begin(), lines.end());
			const auto linesEndIter = BindExecutionPolicy<E>(hpx::unique, lines.begin(), lines.end());
			lines.erase(linesEndIter, lines.end());

			anOutData.myNonVerticalLinesCount = BindExecutionPolicy<E>(hpx::count_if, lines.begin(), lines.end(),
				[](const auto& a) {return !std::get<0>(a); });

			anOutData.myLines.resize(lines.size());
			BindExecutionPolicy<E>(hpx::transform, lines.begin(), lines.end(), anOutData.myLines.begin(),
				[](const auto& a) { return Line2{ std::get<1>(a), std::get<2>(a) }; });

			// If all the lines are parallel then compute result and return early
			const auto areLinesParallel = Utils::AreItemsParallel<E>(anOutData.myLines,
				[](const auto& aLine) { return aLine.direction(); });

			if (areLinesParallel)
			{
				BatchPointResult result;

				const auto line = anOutData.myLines[0].perpendicular(anOutData.myLines[0].point());
				const auto& sortedLinesIndices =
					locGetSortedItemsByIntersection<LineParallelLineVisitor<E>>(anOutData.myLines, line);

				const auto GetAreaIdx = [&](size_t aPointIdx)
				{
					const Point2 projectedPoint{ somePoints[aPointIdx].x(), somePoints[aPointIdx].y() };
					const size_t firstUpperLineIdx = locGeometricBinarySearch(projectedPoint, sortedLinesIndices.size(), [&](size_t i) {
						return anOutData.myLines[sortedLinesIndices[i]];
						});
					const auto isPointAlongLine = firstUpperLineIdx < sortedLinesIndices.size() &&
						anOutData.myLines[sortedLinesIndices[firstUpperLineIdx]].has_on(projectedPoint);

					if (isPointAlongLine)
					{
						return firstUpperLineIdx;
					}
					else
					{
						return firstUpperLineIdx + anOutData.myLines.size();
					}
				};

				std::vector<size_t> areaIndices;
				result.mySortedPlanesCache.resize(1);

				locFillRequiredCacheEntries<E>(somePlanes, somePoints, 0, somePoints.size(), GetAreaIdx, areaIndices, result.mySortedPlanesCache[0]);

				result.myRangeWrappers.resize(somePoints.size());
				BindExecutionPolicy<E>(hpx::for_loop, 0, somePoints.size(), [&](size_t aPointIdx) {
					const auto cacheIdx = areaIndices[aPointIdx];
					const auto& sortedPlanesIndices = result.mySortedPlanesCache[0].at(cacheIdx);
					const auto firstUpperPlaneIdx = locGeometricBinarySearch(somePoints[aPointIdx], sortedPlanesIndices.size(),
						[&](size_t i) { return somePlanes[sortedPlanesIndices[i]]; });
					using Res = BatchPointResult::RangeWrapper;
					result.myRangeWrappers[aPointIdx] = Res{ Range{firstUpperPlaneIdx, somePlanes.size()}, 0, cacheIdx };
					});

				anOutData.myBatchPointResultOpt = result;
				return;
			}

			std::vector<Point2> vertices;

			const auto LineLineSeq = [&](size_t aFirstLineIdx, size_t aSecondLineIdx)
			{
				const auto intersection = CGAL::intersection(anOutData.myLines[aFirstLineIdx], anOutData.myLines[aSecondLineIdx]);
				if (intersection)
				{
					const Point2* point = boost::get<Point2>(&*intersection);
					CGAL_precondition(point != nullptr);
					vertices.emplace_back(*point);
				}
			};
			const auto LineLinePar = [&](size_t aFirstLineIdx, size_t aSecondLineIdx, std::atomic<size_t>& anOutCounter)
			{
				const auto intersection = CGAL::intersection(anOutData.myLines[aFirstLineIdx], anOutData.myLines[aSecondLineIdx]);
				if (intersection)
				{
					const Point2* point = boost::get<Point2>(&*intersection);
					CGAL_precondition(point != nullptr);
					const auto outIdx = anOutCounter.fetch_add(1, std::memory_order_relaxed);
					vertices[outIdx] = *point;
				}
			};

			if constexpr (E == ExecutionPolicy::SEQ)
			{
				hpx::for_loop(hpx::execution::seq, 0, anOutData.myLines.size(), [&](size_t aLineIdx) {
					hpx::for_loop(hpx::execution::seq, aLineIdx + 1, anOutData.myLines.size(), [&](size_t anotherLineIdx) {
						LineLineSeq(aLineIdx, anotherLineIdx);
						});
					});
			}
			else
			{
				const auto maxVerticesCount = anOutData.myLines.size() * (anOutData.myLines.size() - 1) / 2;
				vertices.resize(maxVerticesCount);
				std::atomic<size_t> counter{ 0 };
				BindExecutionPolicy<E>(hpx::for_loop, 0, anOutData.myLines.size(), [&](size_t aLineIdx) {
					BindExecutionPolicy<E>(hpx::for_loop, aLineIdx + 1, anOutData.myLines.size(), [&](size_t anotherLineIdx) {
						LineLinePar(aLineIdx, anotherLineIdx, counter);
						});
					});
				vertices.resize(counter.load());
			}

			BindExecutionPolicy<E>(hpx::sort, vertices.begin(), vertices.end(), CGAL::Less<Point2, Point2>());

			if constexpr (E == ExecutionPolicy::SEQ)
			{
				std::vector<FT> uniqueYValues(1, vertices[0].y());
				FT lastSlabX = vertices[0].x();

				hpx::for_loop(hpx::execution::seq, 1, vertices.size(), [&](size_t aVertexIdx) {
					if (vertices[aVertexIdx].x() == lastSlabX)
					{
						if (vertices[aVertexIdx].y() != uniqueYValues.back())
						{
							uniqueYValues.emplace_back(vertices[aVertexIdx].y());
						}
					}
					else
					{
						anOutData.mySlabs.emplace_back(std::make_tuple(lastSlabX, uniqueYValues));

						lastSlabX = vertices[aVertexIdx].x();
						uniqueYValues.clear();
						uniqueYValues.emplace_back(vertices[aVertexIdx].y());
					}
					});

				anOutData.mySlabs.emplace_back(std::make_tuple(lastSlabX, uniqueYValues));
			}
			else
			{
				std::atomic<size_t> counter{ 0 };
				anOutData.mySlabs.resize(vertices.size());

				BindExecutionPolicy<E>(hpx::for_loop, 0, vertices.size(), [&](size_t aVertexIdx) {
					const auto& slabX = vertices[aVertexIdx].x();
					if (aVertexIdx == 0 ||
						slabX > vertices[aVertexIdx - 1].x())
					{
						std::vector<FT> uniqueYValues(1, vertices[aVertexIdx].y());

						for (auto i = aVertexIdx + 1; i < vertices.size() && slabX == vertices[i].x(); ++i)
						{
							if (vertices[i].y() != uniqueYValues.back())
							{
								uniqueYValues.emplace_back(vertices[i].y());
							}
						}

						anOutData.mySlabs[counter.fetch_add(1, std::memory_order_relaxed)] = std::make_tuple(slabX, uniqueYValues);
					}
					});

				anOutData.mySlabs.resize(counter.load());

				BindExecutionPolicy<E>(hpx::sort, anOutData.mySlabs.begin(), anOutData.mySlabs.end(),
					[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
			}

			const auto minMaxItems =
				BindExecutionPolicy<E>(hpx::minmax_element, somePoints.begin(), somePoints.end(), CGAL::Less<Point3, Point3>());

			// Maybe add 2 fake boundaries
			static const FT offset = 1000.f;

			if (minMaxItems.min->x() < std::get<0>(anOutData.mySlabs[0]))
			{
				std::vector<FT> uniqueY(1, minMaxItems.min->y());
				const auto& startEntry = std::make_tuple(minMaxItems.min->x() - offset, uniqueY);
				anOutData.mySlabs.insert(anOutData.mySlabs.begin(), startEntry);
			}

			if (minMaxItems.max->x() >= std::get<0>(anOutData.mySlabs.back()))
			{
				std::vector<FT> uniqueY(1, minMaxItems.max->y());
				const auto& endEntry = std::make_tuple(minMaxItems.max->x() + offset, uniqueY);
				anOutData.mySlabs.emplace_back(endEntry);
			}
		}
	}

	template<ExecutionPolicy E>
	BatchPointResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints)
	{
		if (somePlanes.empty() || somePoints.empty())
		{
			return BatchPointResult{};
		}

		// Requirements check (TEMPORARY DISABLED)
		// CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		// CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));
		// CGAL_precondition(Utils::ArePlanesUniformlyOriented(somePlanes));

		// 1st edge case: single or parallel planes
		const auto arePlanesParallel = Utils::AreItemsParallel<E>(somePlanes,
			[](const Plane& aPlane) { return aPlane.orthogonal_direction(); });

		if (arePlanesParallel)
		{
			BatchPointResult result;

			const auto line = somePlanes[0].perpendicular_line(somePlanes[0].point());
			const auto& sortedPlanesIndices =
				locGetSortedItemsByIntersection<LinePlaneVisitor<E>>(somePlanes, line);

			result.mySortedPlanesCache.emplace_back(std::unordered_map<size_t, std::vector<size_t>>{});
			result.mySortedPlanesCache[0].insert(std::make_pair(0, sortedPlanesIndices));
			result.myRangeWrappers.resize(somePoints.size());

			hpx::for_loop(hpx::execution::par_unseq, 0, somePoints.size(), [&](size_t aPointIdx) {
				const auto firstUpperPlaneIdx = locGeometricBinarySearch(somePoints[aPointIdx], sortedPlanesIndices.size(),
				[&](size_t i) { return somePlanes[sortedPlanesIndices[i]]; });
			result.myRangeWrappers[aPointIdx] = { Range{firstUpperPlaneIdx, somePlanes.size()}, 0, 0 };
				});

			return result;
		}

		BatchPointLocationData data;

		locComputeData<E>(somePlanes, somePoints, data);

		if (data.myBatchPointResultOpt.has_value())
		{
			return data.myBatchPointResultOpt.value();
		}

		BatchPointResult result;

		std::vector<std::tuple<Point3, size_t, size_t>> queryPoints(somePoints.size());

		BindExecutionPolicy<E>(hpx::for_loop, 0, somePoints.size(),
			[&](size_t aPointIdx) { queryPoints[aPointIdx] = std::make_tuple(somePoints[aPointIdx], aPointIdx, 0); });

		BindExecutionPolicy<E>(hpx::sort, queryPoints.begin(), queryPoints.end());

		if constexpr (E == ExecutionPolicy::SEQ)
		{
			const Point2 projectedPoint{ std::get<0>(queryPoints[0]).x(), std::get<0>(queryPoints[0]).y() };
			const auto startSlabIdx = std::get<0>(locSlabBinarySearch<SlabSearchType::EXTERNAL>(projectedPoint.x(), data.mySlabs.size(),
				[&](size_t i) { return std::get<0>(data.mySlabs[i]); }));
			for (size_t i = startSlabIdx, k = 0; i < data.mySlabs.size() - 1 && k < queryPoints.size(); ++i)
			{
				while (k < queryPoints.size() && std::get<0>(queryPoints[k]).x() < std::get<0>(data.mySlabs[i + 1]))
				{
					std::get<2>(queryPoints[k++]) = i;
				}
			}
		}
		else
		{
			std::atomic<size_t> counter{ 0 };

			const auto runningThreadsCount = hpx::resource::get_num_threads();
			const auto chunkSize = static_cast<size_t>(std::ceil(queryPoints.size() / (float)runningThreadsCount));

			// Additional logic on chunk size in order to take the most appropriate decision here

			BindExecutionPolicy<E>(hpx::for_loop_strided, 0, queryPoints.size(), chunkSize,
				[&](size_t aChunkStartIdx) {
					const Point2 projectedPoint{ std::get<0>(queryPoints[aChunkStartIdx]).x(), std::get<0>(queryPoints[aChunkStartIdx]).y() };
					const auto startSlabIdx = std::get<0>(locSlabBinarySearch<SlabSearchType::EXTERNAL>(projectedPoint.x(), data.mySlabs.size(),
						[&](size_t i) { return std::get<0>(data.mySlabs[i]); }));
					const auto chunkEnd = std::min<size_t>(aChunkStartIdx + chunkSize, queryPoints.size());

					for (auto i = startSlabIdx, k = aChunkStartIdx; i < data.mySlabs.size() - 1 && k < chunkEnd; ++i)
					{
						counter.fetch_add(1, std::memory_order_relaxed);
						while (k < chunkEnd && std::get<0>(queryPoints[k]).x() < std::get<0>(data.mySlabs[i + 1]))
						{
							std::get<2>(queryPoints[k++]) = i;
						}
					}
				});

			result.mySortedPlanesCache.resize(counter.load());
		}

		BindExecutionPolicy<E>(hpx::sort, queryPoints.begin(), queryPoints.end(),
			[](const auto& a, const auto& b) { return std::get<2>(a) < std::get<2>(b); });

		result.myRangeWrappers.resize(queryPoints.size());

		const auto queryPointsProcessing = [&](size_t aPointIdx, size_t claimedIdx)
		{
			const size_t rangeStartIdx = aPointIdx;
			while (++aPointIdx < queryPoints.size() && std::get<2>(queryPoints[aPointIdx]) == std::get<2>(queryPoints[rangeStartIdx]));
			const size_t rangeEndIdx = aPointIdx;

			const auto& slab = data.mySlabs[std::get<2>(queryPoints[rangeStartIdx])];

			bool hasVerticalLine = locHasVerticalLine(std::get<0>(slab), data.myLines, data.myNonVerticalLinesCount);
			const auto& sortedLines = locGetSortedLines<ExecutionPolicy::SEQ>(std::get<0>(slab), data.myLines, data.myNonVerticalLinesCount);

			auto& cache = result.mySortedPlanesCache[claimedIdx];

			const auto GetAreaIdx = [&](const size_t aPointIdx)
			{
				const Point2 projectedPoint{ std::get<0>(queryPoints[aPointIdx]).x(), std::get<0>(queryPoints[aPointIdx]).y() };
				const auto verticesCount = std::get<1>(slab).size();
				auto cacheIdx = -1;

				if (projectedPoint.x() == std::get<0>(slab))
				{
					const auto test = locSlabBinarySearch<SlabSearchType::INTERNAL>(projectedPoint.y(), std::get<1>(slab).size(),
						[&](size_t i) { return std::get<1>(slab)[i]; });

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
					const auto lineIdx = locGeometricBinarySearch(projectedPoint, sortedLines.size(), [&](size_t i) { return sortedLines[i]; });
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
				return cacheIdx;
			};

			std::vector<size_t> areaIndices;
			locFillRequiredCacheEntries<ExecutionPolicy::SEQ>(somePlanes, queryPoints, rangeStartIdx, rangeEndIdx, GetAreaIdx, areaIndices, cache);

			hpx::for_loop(hpx::execution::seq, rangeStartIdx, rangeEndIdx, [&](size_t aPointIdx) {
				const auto cacheIdx = areaIndices[aPointIdx - rangeStartIdx];
				const auto& sortedPlanesIndices = cache.at(cacheIdx);
				const auto firstUpperPlaneIdx = locGeometricBinarySearch(std::get<0>(queryPoints[aPointIdx]), sortedPlanesIndices.size(),
					[&](size_t i) { return somePlanes[sortedPlanesIndices[i]]; });
				using Res = BatchPointResult::RangeWrapper;
				result.myRangeWrappers[std::get<1>(queryPoints[aPointIdx])] =
					Res{ Range{firstUpperPlaneIdx, somePlanes.size()}, claimedIdx, cacheIdx };
				});
			return rangeEndIdx;
		};

		if constexpr (E == ExecutionPolicy::SEQ)
		{
			size_t i = 0;
			while (i < queryPoints.size())
			{
				const auto claimedIdx = result.mySortedPlanesCache.size();
				result.mySortedPlanesCache.emplace_back(std::unordered_map<size_t, std::vector<size_t>>{});
				i = queryPointsProcessing(i, claimedIdx);
			}
		}
		else
		{
			std::atomic<size_t> counter{ 0 };

			BindExecutionPolicy<E>(hpx::for_loop, 0, queryPoints.size(),
				[&](size_t aPointIdx) {
					if (aPointIdx == 0 ||
						std::get<2>(queryPoints[aPointIdx]) != std::get<2>(queryPoints[aPointIdx - 1]))
					{
						queryPointsProcessing(aPointIdx, counter.fetch_add(1, std::memory_order_relaxed));
					}
				});

			result.mySortedPlanesCache.resize(counter.load());
		}

		return result;
	}

	// ------------------------------------------------
	// Explicit instantiation for each execution policy
	template BatchPointResult BatchPointLocation<ExecutionPolicy::SEQ>(const std::vector<Plane>&, const std::vector<Point3>&);
	template BatchPointResult BatchPointLocation<ExecutionPolicy::PAR>(const std::vector<Plane>&, const std::vector<Point3>&);
	template BatchPointResult BatchPointLocation<ExecutionPolicy::PAR_UNSEQ>(const std::vector<Plane>&, const std::vector<Point3>&);
}