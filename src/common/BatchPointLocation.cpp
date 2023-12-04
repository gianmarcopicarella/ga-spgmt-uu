#include "BatchPointLocation.h"
#include "../common/Utils.h"

#include <vector>
#include <functional>

// just for debugging
#include <iostream>
#include <chrono>

#include <CGAL/Compact_container.h>
#include <CGAL/Concurrent_compact_container.h>

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/blocked_range2d.h>

#include <oneapi/tbb/parallel_sort.h>

#include <oneapi/tbb/task_scheduler_observer.h>

#include <functional>
#include <oneapi/tbb/task_arena.h>

#include <oneapi/tbb/parallel_reduce.h>

//#include <CGAL/hilbert_sort.h>

namespace SPGMT
{
	namespace
	{
		template<typename T, typename K>
		void locGetSortedItemsByLine3(const std::vector<T>& someItems, const K& aLine, std::vector<size_t>& someOutItemsIndices)
		{
			struct P { size_t myIdx; FT myValue; };
			std::vector<P> planes;
			planes.reserve(someItems.size());
			for (size_t i = 0; i < someItems.size(); ++i)
			{
				const auto& intersection = CGAL::intersection(someItems[i], aLine);
				if (intersection)
				{
					const Point3* point = boost::get<Point3>(&*intersection);
					CGAL_precondition(point != nullptr);
					const Vec3 vec{ aLine.point(), *point };
					if (CGAL::sign(CGAL::scalar_product(vec, aLine.to_vector())) == CGAL::Sign::POSITIVE)
					{
						planes.emplace_back(P{ i, vec.squared_length() });
					}
					else
					{
						planes.emplace_back(P{ i, -vec.squared_length() });
					}
				}
			}

			std::sort(planes.begin(), planes.end(), [](const auto& a, const auto& b) { return a.myValue < b.myValue; });
			someOutItemsIndices.resize(planes.size());
			std::transform(planes.begin(), planes.end(), someOutItemsIndices.begin(), [](const auto& a) { return a.myIdx; });
		}

		template<typename T, typename K>
		void locGetSortedItemsByLine2(const std::vector<T>& someItems, const K& aLine, std::vector<size_t>& someOutItemsIndices)
		{
			struct P { size_t myIdx; FT myValue; };
			std::vector<P> planes;
			planes.reserve(someItems.size());
			for (size_t i = 0; i < someItems.size(); ++i)
			{
				const auto& intersection = CGAL::intersection(someItems[i], aLine);
				if (intersection)
				{
					const Point2* point = boost::get<Point2>(&*intersection);
					CGAL_precondition(point != nullptr);
					const Vec2 vec{ aLine.point(), *point };
					if (CGAL::sign(CGAL::scalar_product(vec, aLine.to_vector())) == CGAL::Sign::POSITIVE)
					{
						planes.emplace_back(P{ i, vec.squared_length() });
					}
					else
					{
						planes.emplace_back(P{ i, -vec.squared_length() });
					}
				}
			}

			std::sort(planes.begin(), planes.end(), [](const auto& a, const auto& b) { return a.myValue < b.myValue; });
			someOutItemsIndices.resize(planes.size());
			std::transform(planes.begin(), planes.end(), someOutItemsIndices.begin(), [](const auto& a) { return a.myIdx; });
		}

		struct BatchPointLocationData
		{
			std::vector<Point2> myVertices;
			//std::vector</*std::tuple<FT, std::vector<FT>>*/Slab> mySlabs;
			std::vector<Line2> myLines;
			size_t myNonVerticalLinesCount;
			std::optional<BatchPointResult> myBatchPointResultOpt;
		};

		template<typename K, typename F>
		int locGeometricBinarySearch(
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
			return -1;
		}

		int locSlabBinarySearch(
			const Point2& aPoint,
			const std::vector<Point2>& someVertices)
		{
			int left = 0, right = someVertices.size();
			int mid;
			while (left < right)
			{
				mid = left + (right - left) / 2;
				const auto& item = someVertices[mid];
				if (aPoint.x() < item.x())
				{
					if (mid == 0 || aPoint.x() > someVertices[mid - 1].x())
					{
						mid = mid - 1;
						break;
					}
					else
					{
						right = mid;
					}
				}
				else if (aPoint.x() > item.x())
				{
					if (mid == someVertices.size() - 1 ||
						aPoint.x() < someVertices[mid + 1].x())
					{
						break;
					}
					else
					{
						left = mid + 1;
					}
				}
				else
				{
					break;
				}
			}

			// Go to first slab vertex
			while (--mid > 0 && someVertices[mid].x() == someVertices[mid + 1].x());
			return mid + 1;

			// It should never reach this area
			CGAL_postcondition(false);
			return -1;
		}

		std::tuple<size_t, bool> locTestBinarySearch(
			const Point2& aPoint,
			const std::vector<Point2>& someVertices)
		{
			CGAL_precondition(someVertices.size() > 0);
			int left = 0, right = someVertices.size();
			while (left < right)
			{
				const size_t mid = left + (right - left) / 2;
				const auto& item = someVertices[mid];
				if (aPoint.y() < item.y())
				{
					if (mid == 0 || aPoint.y() > someVertices[mid - 1].y())
					{
						return std::make_tuple(mid, false);
					}
					else
					{
						right = mid;
					}
				}
				else if (aPoint.y() > item.y())
				{
					if (mid == someVertices.size() - 1 ||
						aPoint.y() < someVertices[mid + 1].y())
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

			// It should never reach this area
			CGAL_postcondition(false);
			return std::make_tuple(SIZE_MAX, false);
		}

		void locComputeData(
			const std::vector<Plane>& somePlanes,
			const std::vector<Point3>& somePoints,
			BatchPointLocationData& anOutData)
		{
			// Start lines computation
			{
				using LW = std::tuple<bool, Point2, Point2>;
				tbb::enumerable_thread_specific<std::vector<LW>> data;
				tbb::enumerable_thread_specific<size_t> count{ 0 };

				{
					static tbb::affinity_partitioner ap;
					tbb::parallel_for(tbb::blocked_range2d<int, int>(0, somePlanes.size(), 0, somePlanes.size()),
						[&](tbb::blocked_range2d<int, int>& r) {
							static const Point3 zero{ 0,0,0 };
							if (r.rows().end() > r.cols().begin())
							{
								for (int i = r.rows().begin(); i < r.rows().end(); ++i)
								{
									for (int k = r.cols().begin(); k < r.cols().end(); ++k)
									{
										if (i <= k) continue;
										const auto intersection = CGAL::intersection(somePlanes[k], somePlanes[i]);
										if (intersection)
										{
											const Line3* line = boost::get<Line3>(&*intersection);
											CGAL_precondition(line != nullptr);

											auto start = line->projection(zero);
											auto end = start + line->to_vector() * FT(1000.f);
											if (end < start) std::swap(start, end);

											data.local().emplace_back(
												std::make_tuple(start.y() == end.y(), Point2{ start.x(), start.y() }, Point2{ end.x(), end.y() }));
											++count.local();
										}
									}
								}
							}
						});
				}

				std::vector<LW> lines;
				lines.resize(count.combine(std::plus<size_t>()));
				size_t copyStartIdx = 0;
				for (auto it = data.begin(); it != data.end(); copyStartIdx += it->size(), ++it)
				{
					std::copy(it->begin(), it->end(), lines.begin() + copyStartIdx);
				}

				{
					static tbb::affinity_partitioner ap;
					tbb::parallel_sort(lines.begin(), lines.end());
				}

				// Sequential copy unique
				const auto& newEnd = std::unique(lines.begin(), lines.end());
				anOutData.myLines.reserve(std::distance(lines.begin(), newEnd));
				for (auto it = lines.begin(); it != newEnd; ++it)
				{
					anOutData.myLines.emplace_back(Line2{ std::get<1>(*it), std::get<2>(*it) });
				}

				anOutData.myNonVerticalLinesCount = std::count_if(lines.begin(), newEnd, [](const auto& aLine) {
					return !std::get<0>(aLine);
					});
			}
			// End lines computation

			// TOADD: Edge case parallel lines
			{
				const auto& refDir = anOutData.myLines[0].direction();
				const auto areParallel = std::all_of(anOutData.myLines.begin(), anOutData.myLines.end(), [&](const auto& aLine) {
					return refDir == aLine.direction() || refDir == -aLine.direction();
					});
				if (areParallel)
				{
					BatchPointResult result;

					const auto& line = anOutData.myLines[0].perpendicular(anOutData.myLines[0].point());

					std::vector<size_t> sortedLineindices;
					locGetSortedItemsByLine2(anOutData.myLines, line, sortedLineindices);

					// Start from here <--

					struct S { Point3 myQueryPoint; size_t mySlabIdx; };
					std::vector<S> slabs;
					slabs.resize(somePoints.size());

					{
						static tbb::affinity_partitioner ap;
						tbb::parallel_for(tbb::blocked_range<int>(0, somePoints.size()),
							[&](tbb::blocked_range<int>& r) {
								std::vector<S> local;
								for (auto i = r.begin(); i < r.end(); ++i)
								{
									const Point2 projectedPoint{ somePoints[i].x(), somePoints[i].y() };
									const auto lineIdx = locGeometricBinarySearch(projectedPoint, sortedLineindices.size(),
										[&](size_t k) { return anOutData.myLines[sortedLineindices[k]]; });

									const auto isPointAlongLine = lineIdx < sortedLineindices.size() &&
										anOutData.myLines[sortedLineindices[lineIdx]].has_on(projectedPoint);
									const auto finalLineIdx = isPointAlongLine ? lineIdx : (lineIdx + anOutData.myLines.size());
									local.emplace_back(S{ somePoints[i], finalLineIdx });
								}
								std::copy(local.begin(), local.end(), slabs.begin() + r.begin());
							});
					}

					// Sequential Slab filtering
					{
						std::vector<S> cp(slabs);
						std::sort(cp.begin(), cp.end(), [](const auto& a, const auto& b) {
							return a.mySlabIdx < b.mySlabIdx;
							});
						const auto newEnd = std::unique(cp.begin(), cp.end(), [](const auto& a, const auto& b) {
							return a.mySlabIdx == b.mySlabIdx;
							});
						cp.erase(newEnd, cp.end());

						result.mySortedPlanesCache.insert(std::make_pair(0, BatchPointResult::SlabCache{}));
						for (const auto slab : cp)
						{
							result.mySortedPlanesCache.at(0).insert(std::make_pair(slab.mySlabIdx, std::vector<size_t>{}));
						}

						static tbb::affinity_partitioner ap;
						tbb::parallel_for(tbb::blocked_range<int>(0, cp.size()),
							[&](tbb::blocked_range<int>& r) {
								for (auto i = r.begin(); i < r.end(); ++i)
								{
									locGetSortedItemsByLine3(somePlanes, Line3{ cp[i].myQueryPoint, Vec3{0,0,1} },
										result.mySortedPlanesCache.at(0).find(cp[i].mySlabIdx)->second);
								}
							});
					}


					result.myRangeWrappers.resize(somePoints.size());

					{
						static tbb::affinity_partitioner ap;
						tbb::parallel_for(tbb::blocked_range<int>(0, somePoints.size()),
							[&](tbb::blocked_range<int>& r) {
								for (auto i = r.begin(); i < r.end(); ++i)
								{
									const auto& sortedPlanesIndices = result.mySortedPlanesCache[0].at(slabs[i].mySlabIdx);
									const auto firstUpperPlaneIdx = locGeometricBinarySearch(somePoints[i], sortedPlanesIndices.size(),
										[&](size_t i) { return somePlanes[sortedPlanesIndices[i]]; });
									using R = BatchPointResult::RangeWrapper;
									result.myRangeWrappers[i] = R{ Range{firstUpperPlaneIdx, somePlanes.size()}, 0, slabs[i].mySlabIdx };
								}
							});
					}

					anOutData.myBatchPointResultOpt = result;
					return;
				}
			}

			// Start vertices computation
			{
				anOutData.myVertices.resize(anOutData.myLines.size() * (anOutData.myLines.size() - 1));

				class t : public tbb::task_scheduler_observer
				{

				public:
					tbb::enumerable_thread_specific<std::vector<Point2>> data;

					t(std::vector<Point2>& abc, size_t& finalCount)
						: tbb::task_scheduler_observer( /*local=*/ /*true*/), out(abc), fc(finalCount) {
						observe(true); // activate the observer
					}
					~t() {
						observe(false); // deactivate the observer

						for (auto& t : data) {
							fc += t.size();
						}
					}
					/*override*/ void on_scheduler_entry(bool worker) {
						CGAL_precondition(data.local().empty());
					}
					/*override*/ void on_scheduler_exit(bool worker) {
						const int start = c.fetch_add(data.local().size());
						std::copy(data.local().begin(), data.local().end(), out.begin() + start);
					}
					size_t& fc;
					std::atomic<size_t> c{ 0 };
					std::vector<Point2>& out;

				};

				size_t finalSize{ 0 };

				{
					t tso(anOutData.myVertices, finalSize);
					static tbb::affinity_partitioner ap;
					tbb::parallel_for(tbb::blocked_range2d<int, int>(0, anOutData.myLines.size(), 0, anOutData.myLines.size()),
						[&](tbb::blocked_range2d<int, int>& r) {
							if (r.rows().end() > r.cols().begin())
							{
								for (int i = r.rows().begin(); i < r.rows().end(); ++i)
								{
									for (int k = r.cols().begin(); k < r.cols().end(); ++k)
									{
										if (i <= k) continue;
										const auto intersection = CGAL::intersection(anOutData.myLines[k], anOutData.myLines[i]);
										if (intersection)
										{
											const Point2* point = boost::get<Point2>(&*intersection);
											CGAL_precondition(point != nullptr);
											tso.data.local().emplace_back(*point);
										}
									}
								}
							}

						});

					const int start = tso.c.fetch_add(tso.data.local().size());
					std::copy(tso.data.local().begin(), tso.data.local().end(), anOutData.myVertices.begin() + start);
				}

				anOutData.myVertices.resize(finalSize);
				const auto& extremes = std::minmax_element(somePoints.begin(), somePoints.end());
				anOutData.myVertices.emplace_back(Point2{ extremes.first->x(), extremes.first->y() });
				anOutData.myVertices.emplace_back(Point2{ extremes.second->x(), extremes.second->y() });

				{
					static tbb::affinity_partitioner ap;
					tbb::parallel_sort(anOutData.myVertices.begin(), anOutData.myVertices.end(), CGAL::Less<Point2, Point2>());
				}
			}
			// Done vertices computation
		}
	}

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
		{
			const auto& refDir = somePlanes[0].orthogonal_direction();
			const auto arePlanesParallel = std::all_of(somePlanes.begin(), somePlanes.end(), [&](const auto& aPlane) {
				return refDir == aPlane.orthogonal_direction() || refDir == -aPlane.orthogonal_direction();
				});
			if (arePlanesParallel)
			{
				BatchPointResult result;
				const auto line = somePlanes[0].perpendicular_line(somePlanes[0].point());

				result.mySortedPlanesCache.insert(std::make_pair(0, BatchPointResult::SlabCache{}));
				auto& sortedPlanesIndices = result.mySortedPlanesCache.at(0).insert(std::make_pair(0, std::vector<size_t> {})).first->second;
				locGetSortedItemsByLine3(somePlanes, line, sortedPlanesIndices);
				result.myRangeWrappers.resize(somePoints.size());

				tbb::parallel_for(tbb::blocked_range<int>(0, somePoints.size()), [&](tbb::blocked_range<int>& r) {
					using RW = BatchPointResult::RangeWrapper;
					std::vector<RW> local;
					for (auto i = r.begin(); i < r.end(); ++i)
					{
						const auto firstUpperPlaneIdx = locGeometricBinarySearch(somePoints[i], sortedPlanesIndices.size(),
							[&](size_t k) { return somePlanes[sortedPlanesIndices[k]]; });
						local.emplace_back(RW{ Range{firstUpperPlaneIdx, somePlanes.size()}, 0, 0 });
					}
					std::copy(local.begin(), local.end(), result.myRangeWrappers.begin() + r.begin());
					});
				return result;
			}
		}

		BatchPointLocationData data;

		locComputeData(somePlanes, somePoints, data);

		if (data.myBatchPointResultOpt.has_value())
		{
			return data.myBatchPointResultOpt.value();
		}

		//std::cout << "START" << std::endl;

		BatchPointResult result;
		result.myRangeWrappers.resize(somePoints.size());

		//std::cout << "START slabs lookup" << std::endl;

		using P = std::tuple<Point3, int, int>;
		std::vector<P> queryPoints;
		queryPoints.resize(somePoints.size());

		// Find slab index for each query point
		{
			tbb::enumerable_thread_specific<std::vector<P>> local;
			static tbb::affinity_partitioner ap;
			tbb::parallel_for(tbb::blocked_range<int>(0, somePoints.size()), [&](tbb::blocked_range<int>& r) {
				for (auto i = r.begin(); i < r.end(); ++i)
				{
					const Point2 projectedPoint{ somePoints[i].x(), somePoints[i].y() };
					local.local().emplace_back(std::make_tuple(somePoints[i], locSlabBinarySearch(projectedPoint, data.myVertices), i));
				}
				});

			// Sequential copy (Potential optimization)
			size_t startIdx = 0;
			for (auto it = local.begin(); it != local.end(); startIdx += it->size(), ++it)
			{
				std::copy(it->begin(), it->end(), queryPoints.begin() + startIdx);
			}
		}
		// End slab lookup
		//std::cout << "DONE slabs lookup" << std::endl;

		//std::cout << "START query point sort" << std::endl;
		{
			static tbb::affinity_partitioner ap;
			tbb::parallel_sort(queryPoints.begin(), queryPoints.end(),
				[](const auto& a, const auto& b) {return std::get<1>(a) < std::get<1>(b); });
		}
		//std::cout << "DONE query point sort" << std::endl;

		// Find result
		// 1. Find lines slopes
		std::vector<FT> slopes;
		slopes.reserve(data.myNonVerticalLinesCount);
		for (auto i = 0; i < data.myNonVerticalLinesCount; ++i)
		{
			const auto& dir = data.myLines[i].direction();
			slopes.emplace_back(dir.dy() / dir.dx());
		}

		// 2. Find ranges sequentially
		std::vector<size_t> slabsRanges;

		{
			//std::cout << "START fill slabs ranges and cache" << std::endl;
			using C = BatchPointResult::SlabCache;
			std::vector <std::pair<size_t, C>> cacheEntries;
			size_t i = 0;
			slabsRanges.emplace_back(i);
			do
			{
				cacheEntries.emplace_back(std::make_pair(std::get<1>(queryPoints[i]), C{}));
				while (++i < queryPoints.size() && std::get<1>(queryPoints[i]) == std::get<1>(queryPoints[i - 1]));
				slabsRanges.emplace_back(i);

			} while (i < queryPoints.size());
			//slabsRanges.pop_back();
			result.mySortedPlanesCache.insert(cacheEntries.begin(), cacheEntries.end());
			//std::cout << "END fill slabs ranges and cache" << std::endl;
		}

		// 3. Compute results
		{
			tbb::enumerable_thread_specific<std::vector<Line2>> lines{ data.myLines.begin(), data.myLines.begin() + data.myNonVerticalLinesCount };
			static tbb::affinity_partitioner ap;
			tbb::parallel_for(tbb::blocked_range<int>(0, slabsRanges.size() - 1), [&](tbb::blocked_range<int>& r) {

				for (auto i = r.begin(); i < r.end(); ++i)
				{
					// Common slab work
					const auto pointStart = slabsRanges[i];
					const auto pointEnd = slabsRanges[i + 1];
					const auto& slabIdx = std::get<1>(queryPoints[pointStart]);
					const auto slabStartVertex = data.myVertices[slabIdx];
					const auto slabX = slabStartVertex.x();

					// Count slab vertices
					std::vector<Point2> uniqueY;
					{
						auto z = slabIdx;
						while (z < data.myVertices.size() && data.myVertices[z].x() == slabX)
						{
							uniqueY.emplace_back(data.myVertices[z]);
							while (++z < data.myVertices.size() && data.myVertices[z].y() == uniqueY.back().y());
						}
					}

					// Check for vertical line
					const auto hasVerticalLine = std::find_if(data.myLines.begin() + data.myNonVerticalLinesCount, data.myLines.end(),
						[&](const auto& aVerticalLine) { return aVerticalLine.point().x() == slabX; }) != data.myLines.end();

					auto& cache = result.mySortedPlanesCache.at(slabIdx);

					// Compute sorted lines
					std::vector<FT> yIntercept;
					using L = std::tuple<Line2, size_t>;
					std::vector<L> sortedLines;
					sortedLines.reserve(lines.local().size());

					yIntercept.reserve(data.myNonVerticalLinesCount);
					sortedLines.reserve(data.myNonVerticalLinesCount);
					for (size_t k = 0; k < data.myNonVerticalLinesCount; ++k)
					{
						yIntercept.emplace_back(lines.local()[k].y_at_x(slabX));
						sortedLines.emplace_back(L{ lines.local()[k], k });
					}

					std::sort(sortedLines.begin(), sortedLines.end(), [&](const auto& a, const auto& b) {
						const auto m = std::get<1>(a);
						const auto k = std::get<1>(b);
						return yIntercept[m] < yIntercept[k] || (yIntercept[m] == yIntercept[k] && slopes[m] < slopes[k]);
						});

					// Fill required cache entries
					using A = std::tuple<size_t, size_t>;
					std::vector<A> areaIndices;
					areaIndices.reserve(pointEnd - pointStart);

					{
						// Define area index getter
						const auto GetAreaIdx = [&](const size_t aPointIdx)
						{
							const Point2 projectedPoint{ std::get<0>(queryPoints[aPointIdx]).x(), std::get<0>(queryPoints[aPointIdx]).y() };
							auto cacheIdx = -1;
							const auto verticesCount = uniqueY.size();

							if (projectedPoint.x() == slabX)
							{
								const auto test = locTestBinarySearch(projectedPoint, uniqueY);

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
								const auto lineIdx = locGeometricBinarySearch(projectedPoint, sortedLines.size(), [&](size_t i) { return std::get<0>(sortedLines[i]); });
								if (lineIdx < sortedLines.size() && std::get<0>(sortedLines[lineIdx]).has_on(projectedPoint))
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

						for (auto k = pointStart; k < pointEnd; ++k)
						{
							areaIndices.emplace_back(std::make_pair(GetAreaIdx(k), k));
						}

						std::vector<A> uniqueAreaIndices(areaIndices);
						std::sort(uniqueAreaIndices.begin(), uniqueAreaIndices.end());

						const auto newEnd = std::unique(uniqueAreaIndices.begin(), uniqueAreaIndices.end(),
							[](const auto& a, const auto& b) { return std::get<0>(a) == std::get<0>(b); });
						uniqueAreaIndices.erase(newEnd, uniqueAreaIndices.end());

						for (const auto& a : uniqueAreaIndices)
						{
							CGAL_precondition(cache.find(std::get<0>(a)) == cache.end());
							auto& sortedPlanesIndices = cache.insert(std::make_pair(std::get<0>(a), std::vector<size_t>{})).first->second;
							const auto& queryPoint = std::get<0>(queryPoints[std::get<1>(a)]);
							locGetSortedItemsByLine3(somePlanes, Line3{ queryPoint, Vec3{0,0,1} }, sortedPlanesIndices);
						}
					}

					// Fill result
					for (auto k = pointStart; k < pointEnd; ++k)
					{
						const auto& queryPoint = queryPoints[k];
						const auto& slabIdx = (size_t)std::get<1>(queryPoint);
						const auto cacheIdx = std::get<0>(areaIndices[k - pointStart]);
						const auto& sortedPlanesIndices = cache.at(cacheIdx);
						const auto firstUpperPlaneIdx = locGeometricBinarySearch(std::get<0>(queryPoint), sortedPlanesIndices.size(),
							[&](size_t i) { return somePlanes[sortedPlanesIndices[i]]; });
						using Res = BatchPointResult::RangeWrapper;
						result.myRangeWrappers[std::get<2>(queryPoint)] =
							Res{ Range{firstUpperPlaneIdx, somePlanes.size()}, slabIdx, cacheIdx };
					}

				}
				});
		}

		//std::cout << "Done" << std::endl;

		return result;
	}
}