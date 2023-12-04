#include "BruteForce.h"
#include "Utils.h"
#include "DebugUtils.h"

#include <CGAL/Polygon_2_algorithms.h>

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

namespace SPGMT
{

	namespace
	{
		bool locIsVertexInLowerEnvelope(
			const std::vector<Plane>& somePlanes,
			const Point3& aPoint)
		{
			CGAL_precondition(somePlanes.size() > 0);
			return std::all_of(somePlanes.begin(), somePlanes.end(), [&](const auto& aPlane) {
				return !aPlane.has_on_positive_side(aPoint); });
		}

		std::tuple<size_t, size_t, FT, FT> locGetMinMaxPlaneAtPoint(
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& somePlaneIndices,
			const Point3& aPoint)
		{
			CGAL_precondition(somePlanes.size() > 0);

			static const Vec3 up{ 0,0,1 };
			const Line3 upLine{ aPoint, up };

			const auto computePlaneHeight = [&](const size_t aPlaneIdx)
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes[aPlaneIdx]);
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				return point->z();
			};

			size_t minIdx = 0, maxIdx = 0;
			FT minZ, maxZ;

			{
				const auto planeZ = computePlaneHeight(somePlaneIndices[0]);
				minZ = planeZ;
				maxZ = planeZ;
			}

			for (auto i = 1; i < somePlaneIndices.size(); ++i)
			{
				const auto& planeZ = computePlaneHeight(somePlaneIndices[i]);
				if (planeZ < minZ)
				{
					minZ = planeZ;
					minIdx = somePlaneIndices[i];
				}
				if (planeZ > maxZ)
				{
					maxZ = planeZ;
					maxIdx = somePlaneIndices[i];
				}
			}

			return std::make_tuple(minIdx, maxIdx, minZ, maxZ);
		}

		// Visual Help. https://www.falstad.com/dotproduct/ | https://www.geogebra.org/3d?lang=en
		std::pair<size_t, size_t> locGetMinMaxSteepPlaneIndices(
			const Line3& aLine,
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& somePlanesIndices)
		{
			static const Vec3 up{ 0,0,1 };

			CGAL_precondition(somePlanes.size() > 0);
			CGAL_precondition(somePlanesIndices.size() > 0);

			const auto d = aLine.direction();
			const Vec3 leftDir{ -CGAL::cross_product(Vec3{d.dx(), d.dy(), d.dz()}, up) };
			const auto sampleLoc = aLine.point() + leftDir * 1000.f;

			const auto minMaxPlanes = locGetMinMaxPlaneAtPoint(somePlanes, somePlanesIndices, sampleLoc);
			return std::make_pair(std::get<0>(minMaxPlanes), std::get<1>(minMaxPlanes));
		}

		void locSortEdgesCCW(std::vector<Edge<Point3>>& someOutEdges)
		{
			static const Vec3 ref{ 1, 0, 0 };
			static tbb::affinity_partitioner ap;
			tbb::parallel_sort(someOutEdges.begin(), someOutEdges.end(), [&](const auto& aFirstEdge, const auto& aSecondEdge) {
				if (aFirstEdge.myStart < aSecondEdge.myStart)
				{
					return true;
				}
				else if (aFirstEdge.myStart == aSecondEdge.myStart)
				{
					const auto firstAngle = CGAL::approximate_angle(ref, Vec3{ aFirstEdge.myStart, aFirstEdge.myEnd });
					const auto secondAngle = CGAL::approximate_angle(ref, Vec3{ aFirstEdge.myStart, aSecondEdge.myEnd });
					return firstAngle < secondAngle;
				}
				return false;
				});
		}

		struct BruteForceData
		{
			std::vector<std::tuple<size_t, size_t>> myMinMaxPlanes;
			std::vector<std::tuple<size_t, Point3>> myVertices;
			std::vector<std::tuple<size_t, size_t>> myRanges;
			std::optional<LowerEnvelope3d> myLowerEnvelopeOpt;
		};

		void locComputeData(
			const std::vector<Plane>& somePlanes,
			BruteForceData& anOutData)
		{
			static const Point3 zero{ 0,0,0 };
			static const Vec3 up{ 0,0,1 };
			constexpr auto distance = 1000.f;

			std::vector<Line3> uniqueLines;

			// Unique lines and MinMax planes

			{
				using LW = std::tuple<Point3, Point3, FT, size_t>;
				std::vector<LW> lines;

				tbb::enumerable_thread_specific<std::vector<LW>> data;
				tbb::enumerable_thread_specific<size_t> count{ 0 };
				static tbb::affinity_partitioner ap;

				tbb::parallel_for(tbb::blocked_range2d<int, int>(0, somePlanes.size(), 0, somePlanes.size()),
					[&](tbb::blocked_range2d<int, int>& r) {
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

										// Sort line segment
										auto start = line->projection(zero);
										auto end = start + line->to_vector() * FT(1000.f);
										if (end < start) std::swap(start, end);

										// Compute left planes Z
										const auto& d = Vec3{ start, end }.direction();
										const Vec3 leftDir{ -CGAL::cross_product(Vec3{d.dx(), d.dy(), d.dz()}, up) };
										const auto sampleLoc = start + leftDir * FT(1000.f);
										const Line3 vLine{ sampleLoc, up };

										const auto firstPlaneIntersection = CGAL::intersection(vLine, somePlanes[i]);
										const auto secondPlaneIntersection = CGAL::intersection(vLine, somePlanes[k]);

										CGAL_precondition(firstPlaneIntersection && secondPlaneIntersection);

										const Point3* firstPoint = boost::get<Point3>(&*firstPlaneIntersection);
										CGAL_precondition(firstPoint != nullptr);

										const Point3* secondPoint = boost::get<Point3>(&*secondPlaneIntersection);
										CGAL_precondition(secondPoint != nullptr);

										data.local().emplace_back(std::make_tuple(start, end, firstPoint->z(), i));
										data.local().emplace_back(std::make_tuple(start, end, secondPoint->z(), k));

										count.local() += 2;
									}
								}
							}
						}
					});

				// Copy data back sequentially (possible optimization)
				lines.resize(count.combine(std::plus<size_t>()));
				auto itemsCount = 0;
				for (const auto& chunk : data)
				{
					std::copy(chunk.begin(), chunk.end(), lines.begin() + itemsCount);
					itemsCount += chunk.size();
				}
				lines.resize(itemsCount);

				// Sort data
				{
					static tbb::affinity_partitioner ap;
					tbb::parallel_sort(lines.begin(), lines.end());
				}

				// Pick unique lines and their min-max planes
				size_t i = 0;
				while (i < lines.size())
				{
					std::tuple<size_t, size_t> minMax;
					uniqueLines.emplace_back(Line3{ std::get<0>(lines[i]), std::get<1>(lines[i]) });

					std::get<0>(minMax) = std::get<3>(lines[i]);

					const auto k = i;
					while (++i < lines.size() &&
						std::get<0>(lines[k]) == std::get<0>(lines[i]) &&
						std::get<0>(lines[k]) == std::get<0>(lines[i]));

					std::get<1>(minMax) = std::get<3>(lines[i - 1]);

					anOutData.myMinMaxPlanes.emplace_back(minMax);
				}
			}
			// END

			// If all the lines are parallel then return lower envelope
			const auto areAllLinesParallel = Utils::AreItemsParallel(uniqueLines,
				[](const auto& aLine) { return aLine.direction(); });
			if (areAllLinesParallel)
			{
				static tbb::affinity_partitioner ap;
				std::atomic<int> index{ -1 };
				tbb::parallel_for(tbb::blocked_range<int>(0, uniqueLines.size()),
					[&](tbb::blocked_range<int>& r) {
						for (auto i = r.begin(); i < r.end(); ++i)
						{
							if (locIsVertexInLowerEnvelope(somePlanes, uniqueLines[i].point()))
							{
								index.store(i);
								tbb::task::current_context()->cancel_group_execution();
							}
						}
					});

				CGAL_postcondition(index.load() > -1);

				const auto& minMaxPlanes = anOutData.myMinMaxPlanes[index.load()];
				const auto& start = uniqueLines[index.load()].projection(zero);
				const auto& end = start + uniqueLines[index.load()].to_vector() * distance;

				std::vector<Edge<Point3>> edges(2);
				edges[0].myStart = start < end ? start : end;
				edges[0].myEnd = start < end ? end : start;
				edges[0].myLowestLeftPlane = std::get<0>(minMaxPlanes);
				edges[0].myType = EdgeType::LINE;
				edges[1].myStart = edges[0].myEnd;
				edges[1].myEnd = edges[0].myStart;
				edges[1].myType = EdgeType::LINE;
				edges[1].myLowestLeftPlane = std::get<1>(minMaxPlanes);

				anOutData.myLowerEnvelopeOpt = edges;
				return;
			}


			// Compute vertices
			{
				using VW = std::tuple<size_t, Point3>;

				tbb::enumerable_thread_specific<std::vector<VW>> data;
				tbb::enumerable_thread_specific<size_t> count{ 0 };
				static tbb::affinity_partitioner ap;

				tbb::parallel_for(tbb::blocked_range2d<int, int>(0, uniqueLines.size(), 0, uniqueLines.size()),
					[&](tbb::blocked_range2d<int, int>& r) {
						if (r.rows().end() > r.cols().begin())
						{
							for (int i = r.rows().begin(); i < r.rows().end(); ++i)
							{
								for (int k = r.cols().begin(); k < r.cols().end(); ++k)
								{
									if (i <= k) continue;
									const auto intersection = CGAL::intersection(uniqueLines[k], uniqueLines[i]);
									if (intersection)
									{
										const Point3* point = boost::get<Point3>(&*intersection);
										CGAL_precondition(point != nullptr);

										data.local().emplace_back(std::make_pair(i, *point));
										data.local().emplace_back(std::make_pair(k, *point));

										count.local() += 2;
									}
								}
							}
						}
					});

				// Copy data back sequentially (possible optimization)
				anOutData.myVertices.resize(count.combine(std::plus<size_t>()));
				auto itemsCount = 0;
				for (const auto& chunk : data)
				{
					std::copy(chunk.begin(), chunk.end(), anOutData.myVertices.begin() + itemsCount);
					itemsCount += chunk.size();
				}
				anOutData.myVertices.resize(itemsCount);

				// Sort data
				{
					static tbb::affinity_partitioner ap;
					tbb::parallel_sort(anOutData.myVertices.begin(), anOutData.myVertices.end());
				}

				// Remove duplicate vertices
				anOutData.myVertices.erase(
					std::unique(anOutData.myVertices.begin(), anOutData.myVertices.end()),
					anOutData.myVertices.end());
			}
			// END

			// Add vertices at infinity and compute ranges
			{
				size_t i = 0;
				const auto prevVerticesSize = anOutData.myVertices.size();
				while (i < prevVerticesSize)
				{
					const auto lineIdx = std::get<0>(anOutData.myVertices[i]);
					const auto& minPoint = std::get<1>(anOutData.myVertices[i]);

					while (++i < prevVerticesSize &&
						std::get<0>(anOutData.myVertices[i]) == std::get<0>(anOutData.myVertices[i - 1]));

					const auto& maxPoint = std::get<1>(anOutData.myVertices[i - 1]);

					const auto& minusInfPoint = minPoint - uniqueLines[lineIdx].to_vector() * FT(distance);
					const auto& plusInfPoint = maxPoint + uniqueLines[lineIdx].to_vector() * FT(distance);

					anOutData.myVertices.emplace_back(std::make_tuple(lineIdx, minusInfPoint));
					anOutData.myVertices.emplace_back(std::make_tuple(lineIdx, plusInfPoint));
				}

				// Sort data
				{
					static tbb::affinity_partitioner ap;
					tbb::parallel_sort(anOutData.myVertices.begin(), anOutData.myVertices.end());
				}

				// Compute ranges
				i = 0;
				while (i < anOutData.myVertices.size())
				{
					std::tuple<size_t, size_t> range;
					std::get<0>(range) = i;

					while (++i < anOutData.myVertices.size() &&
						std::get<0>(anOutData.myVertices[i]) == std::get<0>(anOutData.myVertices[i - 1]));

					std::get<1>(range) = i;
					anOutData.myRanges.emplace_back(range);
				}
			}
			// END
		}


	}

	namespace Observers
	{
		/*
		template<typename T>
		class TaskObserver : public tbb::task_scheduler_observer
		{

		public:
			tbb::enumerable_thread_specific<std::vector<T>> myData;

			t(std::vector<Point2>& abc, size_t& finalCount)
				: tbb::task_scheduler_observer(), out(abc), fc(finalCount) {
				observe(true); // activate the observer
			}
			~t() {
				observe(false); // deactivate the observer

				for (auto& t : myData) {
					fc += t.size();
				}
			}
			void on_scheduler_entry(bool worker) {
				CGAL_precondition(myData.local().empty());
			}
			void on_scheduler_exit(bool worker) {
				const int start = c.fetch_add(myData.local().size());
				std::copy(myData.local().begin()myData.local().end(), out.begin() + start);
			}
			size_t& fc;
			std::atomic<size_t> c{ 0 };
			std::vector<Point2>& out;

		};
	*/
	}


	bool locIsEdgeInLowerEnvelope(
		const std::vector<Plane>& somePlanes,
		const BruteForceData& aData,
		const size_t aRangeIdx,
		const size_t aCurrIdx,
		Edge<Point3>& aFirstOutEdge,
		Edge<Point3>& aSecondOutEdge,
		int& anOutSegmentStart)
	{
		const auto lineIdx = std::get<0>(aData.myVertices[aCurrIdx]);
		const auto& midpoint = CGAL::midpoint(
			std::get<1>(aData.myVertices[aCurrIdx]),
			std::get<1>(aData.myVertices[aCurrIdx + 1]));

		const auto isSegmentGood = locIsVertexInLowerEnvelope(somePlanes, midpoint);

		if (isSegmentGood && anOutSegmentStart == -1)
		{
			anOutSegmentStart = aCurrIdx;
		}

		const auto rangeStart = std::get<0>(aData.myRanges[aRangeIdx]);
		const auto rangeEnd = std::get<1>(aData.myRanges[aRangeIdx]);

		const auto isStartIdxCase = anOutSegmentStart != -1 && !isSegmentGood;
		const auto isEndIdxCase = anOutSegmentStart != -1 && isSegmentGood && aCurrIdx == (rangeEnd - 2);

		if (isStartIdxCase || isEndIdxCase)
		{
			const int segmentEndIdx = isStartIdxCase ? aCurrIdx : (aCurrIdx + 1);
			const auto caseIdx =
				static_cast<size_t>(anOutSegmentStart == rangeStart) |
				static_cast<size_t>(segmentEndIdx == rangeEnd - 1) << 1;

			switch (caseIdx)
			{
			case 0:
				aFirstOutEdge.myType = EdgeType::SEGMENT;
				aSecondOutEdge.myType = EdgeType::SEGMENT;
				break;
			case 1:
				aFirstOutEdge.myType = EdgeType::HALF_EDGE_EF;
				aSecondOutEdge.myType = EdgeType::HALF_EDGE_SF;
				break;
			case 2:
				aFirstOutEdge.myType = EdgeType::HALF_EDGE_SF;
				aSecondOutEdge.myType = EdgeType::HALF_EDGE_EF;
				break;
			default:
				aFirstOutEdge.myType = EdgeType::LINE;
				aSecondOutEdge.myType = EdgeType::LINE;
				break;
			}

			const auto& startVertex = std::get<1>(aData.myVertices[anOutSegmentStart]);
			const auto& endVertex = std::get<1>(aData.myVertices[segmentEndIdx]);

			aFirstOutEdge.myStart = startVertex;
			aFirstOutEdge.myEnd = endVertex;
			aFirstOutEdge.myLowestLeftPlane = std::get<0>(aData.myMinMaxPlanes[lineIdx]);

			aSecondOutEdge.myStart = endVertex;
			aSecondOutEdge.myEnd = startVertex;
			aSecondOutEdge.myLowestLeftPlane = std::get<1>(aData.myMinMaxPlanes[lineIdx]);
			anOutSegmentStart = -1;
			return true;
		}
		return false;
	}



	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes)
	{
		if (somePlanes.empty())
		{
			return LowerEnvelope3d{};
		}

		constexpr auto distance = 1000.f;

		// Requirements check (TEMPORARY DISABLED)
		//CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		//CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));

		const auto arePlanesParallel = Utils::AreItemsParallel(somePlanes,
			[](const Plane& aPlane) { return aPlane.orthogonal_direction(); });

		// Edge case: all planes are parallel
		if (arePlanesParallel)
		{
			static const Point3 origin{ 0,0,0 };
			std::vector<size_t> indices(somePlanes.size());
			std::iota(indices.begin(), indices.end(), 0);
			const auto minMaxPlanes = locGetMinMaxPlaneAtPoint(somePlanes, indices, origin);
			CGAL_precondition(std::get<0>(minMaxPlanes) != -1);
			return std::get<0>(minMaxPlanes);
		}

		BruteForceData data;
		locComputeData(somePlanes, data);

		if (data.myLowerEnvelopeOpt.has_value())
		{
			return data.myLowerEnvelopeOpt.value();
		}

		std::vector<Edge<Point3>> result;

		{
			tbb::enumerable_thread_specific<std::vector<Edge<Point3>>> edges;
			tbb::enumerable_thread_specific<size_t> count{ 0 };
			static tbb::affinity_partitioner ap;

			tbb::parallel_for(tbb::blocked_range<int>(0, data.myRanges.size()),
				[&](tbb::blocked_range<int>& r) {
					for (int i = r.begin(); i < r.end(); ++i)
					{
						int segmentStartIdx = -1;
						Edge<Point3> edge, oppEdge;

						for (auto k = std::get<0>(data.myRanges[i]); k < std::get<1>(data.myRanges[i]) - 1; ++k)
						{
							if (locIsEdgeInLowerEnvelope(somePlanes, data, i, k, edge, oppEdge, segmentStartIdx))
							{
								edges.local().emplace_back(edge);
								edges.local().emplace_back(oppEdge);
								count.local() += 2;
							}
						}
					}
				});

			// Copy data back sequentially (possible optimization)
			result.resize(count.combine(std::plus<size_t>()));
			auto itemsCount = 0;
			for (const auto& chunk : edges)
			{
				std::copy(chunk.begin(), chunk.end(), result.begin() + itemsCount);
				itemsCount += chunk.size();
			}
			result.resize(itemsCount);
		}

		locSortEdgesCCW(result);

		return result;
	}

	void TriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope)
	{
		CGAL_precondition(std::holds_alternative<std::vector<Edge<Point3>>>(anOutLowerEnvelope));
		auto& edges = std::get<std::vector<Edge<Point3>>>(anOutLowerEnvelope);
		CGAL_precondition(edges.size() > 2);

		tbb::parallel_sort(edges.begin(), edges.end(), [](const auto& aFirstEdge, const auto& aSecondEdge)
			{
				return aFirstEdge.myLowestLeftPlane < aSecondEdge.myLowestLeftPlane ||
					(aFirstEdge.myLowestLeftPlane == aSecondEdge.myLowestLeftPlane && aFirstEdge.myType < aSecondEdge.myType);
			});

		for (int i = 1, lastStart = 0, prevStartIdx = -1, edgesOldSize = edges.size(); i <= edgesOldSize; ++i)
		{
			if (i == edgesOldSize ||
				edges[lastStart].myLowestLeftPlane != edges[i].myLowestLeftPlane)
			{
				for (auto k = lastStart + 1; k < i; ++k)
				{
					if (prevStartIdx > -1 && (k == prevStartIdx || edges[k].myEnd == edges[prevStartIdx].myStart))
					{
						continue;
					}

					CGAL_precondition(edges[lastStart].myStart != edges[k].myEnd);

					edges.emplace_back(Edge<Point3>{ edges[lastStart].myStart, edges[k].myEnd, EdgeType::SEGMENT_TRIANGLE });
					edges.emplace_back(Edge<Point3>{ edges[k].myEnd, edges[lastStart].myStart, EdgeType::SEGMENT_TRIANGLE });
				}

				lastStart = i;
				prevStartIdx = -1;
			}
			else if (edges[i].myEnd == edges[lastStart].myStart)
			{
				CGAL_precondition(prevStartIdx == -1);
				prevStartIdx = i;
			}
		}

		locSortEdgesCCW(edges);
	}

	size_t CountVerticesInLowerEnvelope(const LowerEnvelope3d& aLowerEnvelope)
	{
		if (std::holds_alternative<std::monostate>(aLowerEnvelope))
		{
			return 0;
		}
		else if (std::holds_alternative<size_t>(aLowerEnvelope))
		{
			return 1;
		}
		else
		{
			// Sequential is already very fast, no parallelization needed
			using EdgesList = std::vector<Edge<Point3>>;
			CGAL_precondition(std::holds_alternative<EdgesList>(aLowerEnvelope));
			const auto& edges = std::get<EdgesList>(aLowerEnvelope);
			size_t verticesCount = 0, i = 0;
			while (i < edges.size() && ++verticesCount)
			{
				while (++i < edges.size() && edges[i].myStart == edges[i - 1].myStart);
			}
			return verticesCount;
		}
	}

	std::vector<size_t> GetLowerEnvelopePlanesIndices(const LowerEnvelope3d& aLowerEnvelope)
	{
		return std::vector<size_t>{};

		if (std::holds_alternative<std::monostate>(aLowerEnvelope))
		{
			return std::vector<size_t>{};
		}
		else if (std::holds_alternative<size_t>(aLowerEnvelope))
		{
			std::vector<size_t> indices;
			indices.emplace_back(std::get<size_t>(aLowerEnvelope));
			return indices;
		}
		else
		{
			CGAL_precondition(std::holds_alternative<std::vector<Edge<Point3>>>(aLowerEnvelope));
			const auto& edges = std::get<std::vector<Edge<Point3>>>(aLowerEnvelope);
			constexpr auto edgeMask = (int)EdgeType::LINE | (int)EdgeType::HALF_EDGE_SF | (int)EdgeType::HALF_EDGE_EF | (int)EdgeType::SEGMENT;
			std::vector<size_t> planesIndices;
			planesIndices.reserve(edges.size());

			for (auto i = 0; i < edges.size(); ++i)
			{
				if ((int)edges[i].myType & edgeMask)
				{
					CGAL_precondition(edges[i].myLowestLeftPlane > -1);
					planesIndices.emplace_back(edges[i].myLowestLeftPlane);
				}
			}

			planesIndices.reserve(planesIndices.size());

			static tbb::affinity_partitioner ap;
			tbb::parallel_sort(planesIndices.begin(), planesIndices.end());

			planesIndices.erase(std::unique(planesIndices.begin(), planesIndices.end()), planesIndices.end());

			return planesIndices;
		}
	}
}