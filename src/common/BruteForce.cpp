#include "BruteForce.h"
#include "Utils.h"
#include "DebugUtils.h"

#include <CGAL/Polygon_2_algorithms.h>

#include <atomic>

namespace SPGMT
{
	namespace
	{
		template<ExecutionPolicy E>
		bool locIsVertexInLowerEnvelope(
			const std::vector<Plane>& somePlanes,
			const Point3& aPoint)
		{
			CGAL_precondition(somePlanes.size() > 0);
			return BindExecutionPolicy<E>(hpx::all_of, somePlanes.begin(), somePlanes.end(), [&](const auto& aPlane) {
				return !aPlane.has_on_positive_side(aPoint); });
		}

		template<ExecutionPolicy E>
		std::tuple<size_t, size_t, FT, FT> locGetMinMaxPlaneAtPoint(
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& somePlaneIndices,
			const Point3& aPoint)
		{
			CGAL_precondition(somePlanes.size() > 0);

			static const Vec3 up{ 0,0,1 };
			const Line3 upLine{ aPoint, up };

			std::vector<std::pair<FT, size_t>> pairs(somePlaneIndices.size());
			std::atomic<size_t> pairsCount{ 0 };

			const auto computePlaneHeight = [&](const auto aPlaneIdx)
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes[*aPlaneIdx]);
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				pairs[pairsCount.fetch_add(1, std::memory_order_relaxed)] = std::make_pair(point->z(), *aPlaneIdx);
			};

			constexpr auto pairsReductionMin = [](const auto& aFirst, const auto& aSecond)
			{
				return aFirst.first < aSecond.first ? aFirst : aSecond;
			};

			constexpr auto pairsReductionMax = [](const auto& aFirst, const auto& aSecond)
			{
				return aFirst.first < aSecond.first ? aSecond : aFirst;
			};

			BindExecutionPolicy<E>(hpx::for_loop, somePlaneIndices.begin(), somePlaneIndices.end(), computePlaneHeight);
			const auto minPlane = BindExecutionPolicy<E>(hpx::reduce, std::next(pairs.begin()), pairs.end(), pairs[0], pairsReductionMin);
			const auto maxPlane = BindExecutionPolicy<E>(hpx::reduce, std::next(pairs.begin()), pairs.end(), pairs[0], pairsReductionMax);
			return std::make_tuple(minPlane.second, maxPlane.second, minPlane.first, maxPlane.first);
		}

		// Visual Help. https://www.falstad.com/dotproduct/ | https://www.geogebra.org/3d?lang=en
		template<ExecutionPolicy E>
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

			const auto minMaxPlanes = locGetMinMaxPlaneAtPoint<E>(somePlanes, somePlanesIndices, sampleLoc);
			return std::make_pair(std::get<0>(minMaxPlanes), std::get<1>(minMaxPlanes));
		}

		template<ExecutionPolicy E>
		void locSortEdgesCCW(std::vector<Edge<Point3>>& someOutEdges)
		{
			static const Vec3 ref{ 1, 0, 0 };
			BindExecutionPolicy<E>(hpx::sort, someOutEdges.begin(), someOutEdges.end(),
				[&](const auto& aFirst, const auto& aSecond) {
					if (aFirst.myStart < aSecond.myStart)
					{
						return true;
					}
					else if (aFirst.myStart == aSecond.myStart)
					{
						const auto firstAngle = CGAL::approximate_angle(ref, Vec3{ aFirst.myStart, aFirst.myEnd });
						const auto secondAngle = CGAL::approximate_angle(ref, Vec3{ aFirst.myStart, aSecond.myEnd });
						return firstAngle < secondAngle;
					}
					return false;
				});
		}

		void locInsertInfinityPoints(
			const std::tuple<size_t, Point3>& aMinPoint,
			const std::tuple<size_t, Point3>& aMaxPoint,
			const std::vector<Line3>& someUniqueLines,
			const size_t aFiniteVerticesCount,
			std::vector<std::tuple<size_t, Point3>>& someOutVertices)
		{
			constexpr auto distance = 1000.f;

			const auto minLineIdx = std::get<0>(aMinPoint);
			const auto maxLineIdx = std::get<0>(aMaxPoint);

			const auto& minusInfPoint = std::get<1>(aMinPoint) - someUniqueLines[minLineIdx].to_vector() * distance;
			const auto& plusInfPoint = std::get<1>(aMaxPoint) + someUniqueLines[maxLineIdx].to_vector() * distance;

			CGAL_precondition(someUniqueLines[minLineIdx].has_on(minusInfPoint));
			CGAL_precondition(someUniqueLines[maxLineIdx].has_on(plusInfPoint));

			const auto minusInfOutIdx = aFiniteVerticesCount + 2 * minLineIdx;
			const auto plusInfOutIdx = aFiniteVerticesCount + 2 * maxLineIdx + 1;

			someOutVertices[minusInfOutIdx] = std::make_tuple(minLineIdx, minusInfPoint);
			someOutVertices[plusInfOutIdx] = std::make_tuple(maxLineIdx, plusInfPoint);
		}

		struct BruteForceData
		{
			std::vector<std::tuple<size_t, size_t>> myMinMaxPlanes;
			std::vector<std::tuple<size_t, Point3>> myVertices;
			std::vector<size_t> myRanges;
			std::optional<LowerEnvelope3d> myLowerEnvelopeOpt;
		};

		template<ExecutionPolicy E>
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

			const auto isSegmentGood = locIsVertexInLowerEnvelope<E>(somePlanes, midpoint);

			if (isSegmentGood && anOutSegmentStart == -1)
			{
				anOutSegmentStart = aCurrIdx;
			}

			const auto isStartIdxCase = anOutSegmentStart != -1 && !isSegmentGood;
			const auto isEndIdxCase = anOutSegmentStart != -1 && isSegmentGood && aCurrIdx == (aData.myRanges[aRangeIdx + 1] - 2);

			if (isStartIdxCase || isEndIdxCase)
			{
				const int segmentEndIdx = isStartIdxCase ? aCurrIdx : (aCurrIdx + 1);
				const auto caseIdx =
					static_cast<size_t>(anOutSegmentStart == aData.myRanges[aRangeIdx]) |
					static_cast<size_t>(segmentEndIdx == aData.myRanges[aRangeIdx + 1] - 1) << 1;

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

		template<ExecutionPolicy E = ExecutionPolicy::SEQ>
		void locComputeData(const std::vector<Plane>& somePlanes, BruteForceData& anOutData)
		{
			static const Point3 zero{ 0,0,0 };
			constexpr auto policy = ExecutionPolicy::SEQ;
			constexpr auto distance = 1000.f;
			using LineWrapper = std::tuple<Point3, Point3, size_t>;
			std::vector<LineWrapper> lines;

			const auto planePlaneIntersection = [&](size_t i, size_t k, std::atomic<size_t>* anOutCounterPtr)
			{
				const auto intersection =
					CGAL::intersection(somePlanes[i], somePlanes[k]);
				if (intersection)
				{
					const Line3* line = boost::get<Line3>(&*intersection);
					CGAL_precondition(line != nullptr);

					auto start = line->projection(zero);
					auto end = start + line->to_vector() * FT(distance);
					if (start > end) std::swap(start, end);
					CGAL_precondition(start < end);

					if constexpr (E == ExecutionPolicy::SEQ)
					{
						lines.emplace_back(std::make_tuple(start, end, i));
						lines.emplace_back(std::make_tuple(start, end, k));
					}
					else
					{
						CGAL_precondition(anOutCounterPtr != nullptr);
						const auto outIdx = anOutCounterPtr->fetch_add(2, std::memory_order_relaxed);
						lines[outIdx] = std::make_tuple(start, end, i);
						lines[outIdx + 1] = std::make_tuple(start, end, k);
					}
				}
			};

			// Find intersection lines
			if constexpr (E == ExecutionPolicy::SEQ)
			{
				for (size_t i = 0; i < somePlanes.size(); ++i)
				{
					for (size_t k = i + 1; k < somePlanes.size(); ++k)
					{
						planePlaneIntersection(i, k, nullptr);
					}
				}
			}
			else
			{
				lines.resize(somePlanes.size() * (somePlanes.size() - 1));
				std::atomic<size_t> counter{ 0 };
				BindExecutionPolicy<E>(hpx::for_loop, 0, somePlanes.size(), [&](size_t i) {
					BindExecutionPolicy<E>(hpx::for_loop, i + 1, somePlanes.size(), [&](size_t k) {
						planePlaneIntersection(i, k, &counter);
						});
					});
				lines.resize(counter.load());
			}

			BindExecutionPolicy<E>(hpx::sort, lines.begin(), lines.end());

			lines.erase(BindExecutionPolicy<E>(hpx::unique, lines.begin(), lines.end()), lines.end());

			// If all the lines are parallel then return lower envelope
			const auto areAllLinesParallel = Utils::AreItemsParallel<E>(lines,
				[](const auto& aLine) { return Dir3{ std::get<1>(aLine) - std::get<0>(aLine) }; });
			if (areAllLinesParallel)
			{
				auto lineIter = BindExecutionPolicy<E>(hpx::find_if, lines.begin(), lines.end(),
					[&](const auto& aLine) { return locIsVertexInLowerEnvelope<E>(somePlanes, std::get<1>(aLine)); });

				std::vector<size_t> planesIndices;

				for (auto it = lineIter;
					it != lines.end() &&
					std::get<0>(*lineIter) == std::get<0>(*it) &&
					std::get<1>(*lineIter) == std::get<1>(*it);
					++it)
				{
					planesIndices.emplace_back(std::get<2>(*it));
				}

				if constexpr (E != ExecutionPolicy::SEQ)
				{
					using RI = std::vector<LineWrapper>::reverse_iterator;
					for (RI it = std::make_reverse_iterator(std::prev(lineIter));
						it != lines.rend() &&
						std::get<0>(*lineIter) == std::get<0>(*it) &&
						std::get<1>(*lineIter) == std::get<1>(*it);
						++it)
					{
						planesIndices.emplace_back(std::get<2>(*it));
					}
				}

				const Line3 lowerEnvelopeLine{ std::get<0>(*lineIter), std::get<1>(*lineIter) };
				const auto minMaxPlanes = locGetMinMaxSteepPlaneIndices<E>(lowerEnvelopeLine, somePlanes, planesIndices);
				const auto& start = lowerEnvelopeLine.projection(zero);
				const auto& end = start + lowerEnvelopeLine.to_vector() * distance;

				std::vector<Edge<Point3>> edges(2);
				edges[0].myStart = start < end ? start : end;
				edges[0].myEnd = start < end ? end : start;
				edges[0].myLowestLeftPlane = minMaxPlanes.first;
				edges[0].myType = EdgeType::LINE;
				edges[1].myStart = edges[0].myEnd;
				edges[1].myEnd = edges[0].myStart;
				edges[1].myType = EdgeType::LINE;
				edges[1].myLowestLeftPlane = minMaxPlanes.second;

				anOutData.myLowerEnvelopeOpt = edges;
				return;
			}

			// Compute min-max planes and unique lines
			std::vector<Line3> uniqueLines;

			const auto uniqueLinesSelection = [&](size_t i, std::atomic<size_t>* aCounterPtr)
			{
				std::vector<size_t> indices;

				size_t k = i;
				while (k < lines.size() &&
					std::get<0>(lines[k]) == std::get<0>(lines[i]) &&
					std::get<1>(lines[k]) == std::get<1>(lines[i]))
				{
					indices.emplace_back(std::get<2>(lines[k++]));
				}

				constexpr auto innerPolicy = ExecutionPolicy::SEQ;

				if constexpr (E == ExecutionPolicy::SEQ)
				{
					uniqueLines.emplace_back(Line3{ std::get<0>(lines[i]), std::get<1>(lines[i]) });
					anOutData.myMinMaxPlanes.emplace_back(
						locGetMinMaxSteepPlaneIndices<innerPolicy>(uniqueLines.back(), somePlanes, indices));
				}
				else
				{
					CGAL_precondition(aCounterPtr != nullptr);
					const auto outIdx = aCounterPtr->fetch_add(1, std::memory_order_relaxed);
					uniqueLines[outIdx] = Line3{ std::get<0>(lines[i]), std::get<1>(lines[i]) };
					anOutData.myMinMaxPlanes[outIdx] =
						locGetMinMaxSteepPlaneIndices<innerPolicy>(uniqueLines[outIdx], somePlanes, indices);
				}

				return k;
			};

			if constexpr (E == ExecutionPolicy::SEQ)
			{
				size_t i = 0;
				while ((i = uniqueLinesSelection(i, nullptr)) < lines.size());
			}
			else
			{
				uniqueLines.resize(lines.size() / 2);
				anOutData.myMinMaxPlanes.resize(lines.size() / 2);
				std::atomic<size_t> counter{ 0 };

				BindExecutionPolicy<E>(hpx::for_loop, uniqueLinesSelection(0, &counter), lines.size(), [&](size_t i) {
					if (std::get<0>(lines[i - 1]) != std::get<0>(lines[i]) ||
						std::get<1>(lines[i - 1]) != std::get<1>(lines[i]))
					{
						uniqueLinesSelection(i, &counter);
					}
					});

				anOutData.myMinMaxPlanes.resize(counter.load());
				uniqueLines.resize(counter.load());
			}

			const auto lineLineIntersection = [&](size_t i, size_t k, std::atomic<size_t>* aCounterPtr)
			{
				const auto intersection =
					CGAL::intersection(uniqueLines[i], uniqueLines[k]);
				if (intersection)
				{
					const Point3* point = boost::get<Point3>(&*intersection);
					CGAL_precondition(point != nullptr);

					if constexpr (E == ExecutionPolicy::SEQ)
					{
						anOutData.myVertices.emplace_back(std::make_tuple(i, *point));
						anOutData.myVertices.emplace_back(std::make_tuple(k, *point));
					}
					else
					{
						CGAL_precondition(aCounterPtr != nullptr);
						const auto outIdx = aCounterPtr->fetch_add(2, std::memory_order_relaxed);
						anOutData.myVertices[outIdx] = std::make_tuple(i, *point);
						anOutData.myVertices[outIdx + 1] = std::make_tuple(k, *point);
					}
				}
			};

			// Compute vertices
			if constexpr (E == ExecutionPolicy::SEQ)
			{
				for (size_t i = 0; i < uniqueLines.size(); ++i)
				{
					for (size_t k = i + 1; k < uniqueLines.size(); ++k)
					{
						lineLineIntersection(i, k, nullptr);
					}
				}
			}
			else
			{
				anOutData.myVertices.resize(uniqueLines.size() * (uniqueLines.size() - 1));
				std::atomic<size_t> counter{ 0 };

				BindExecutionPolicy<E>(hpx::for_loop, 0, uniqueLines.size(), [&](size_t i) {
					BindExecutionPolicy<E>(hpx::for_loop, i + 1, uniqueLines.size(), [&](size_t k) {
						lineLineIntersection(i, k, &counter);
						});
					});

				anOutData.myVertices.resize(counter.load());
			}

			// Sort and remove duplicate vertices
			BindExecutionPolicy<E>(hpx::sort, anOutData.myVertices.begin(), anOutData.myVertices.end());
			anOutData.myVertices.erase(
				BindExecutionPolicy<E>(hpx::unique, anOutData.myVertices.begin(), anOutData.myVertices.end()), anOutData.myVertices.end());

			const auto finiteVerticesCount = anOutData.myVertices.size();
			anOutData.myVertices.resize(finiteVerticesCount + 2 * uniqueLines.size());

			anOutData.myRanges.resize(uniqueLines.size() + 1);
			anOutData.myRanges.front() = 0;
			anOutData.myRanges.back() = anOutData.myVertices.size();

			// Compute points at infinity for first and last vertex in vector
			locInsertInfinityPoints(
				anOutData.myVertices[0],
				anOutData.myVertices[finiteVerticesCount - 1],
				uniqueLines, finiteVerticesCount, anOutData.myVertices);

			BindExecutionPolicy<E>(hpx::for_loop, 1, finiteVerticesCount - 1, [&](size_t i) {
				if (std::get<0>(anOutData.myVertices[i]) <
					std::get<0>(anOutData.myVertices[i + 1]))
				{
					locInsertInfinityPoints(
						anOutData.myVertices[i + 1],
						anOutData.myVertices[i],
						uniqueLines, finiteVerticesCount, anOutData.myVertices);

					const auto nextLineIdx = std::get<0>(anOutData.myVertices[i + 1]);
					anOutData.myRanges[nextLineIdx] = (i + 1) + 2 * nextLineIdx;
				}
				});

			// Sort vertices
			BindExecutionPolicy<E>(hpx::sort, anOutData.myVertices.begin(), anOutData.myVertices.end());
			BindExecutionPolicy<E>(hpx::sort, anOutData.myRanges.begin(), anOutData.myRanges.end());
		}
	}

	template<ExecutionPolicy E>
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

		const auto arePlanesParallel = Utils::AreItemsParallel<E>(somePlanes,
			[](const Plane& aPlane) { return aPlane.orthogonal_direction(); });

		// Edge case: all planes are parallel
		if (arePlanesParallel)
		{
			static const Point3 origin{ 0,0,0 };
			std::vector<size_t> indices(somePlanes.size());
			std::iota(indices.begin(), indices.end(), 0);
			const auto minMaxPlanes = locGetMinMaxPlaneAtPoint<E>(somePlanes, indices, origin);
			CGAL_precondition(std::get<0>(minMaxPlanes) != -1);
			return std::get<0>(minMaxPlanes);
		}

		BruteForceData data;
		locComputeData<E>(somePlanes, data);

		if (data.myLowerEnvelopeOpt.has_value())
		{
			return data.myLowerEnvelopeOpt.value();
		}

		std::vector<Edge<Point3>> edges;

		const auto buildEdges = [&](size_t aRangeIdx, std::atomic<size_t>* aCounterPtr)
		{
			int segmentStartIdx = -1;
			Edge<Point3> edge, oppEdge;
			for (size_t i = data.myRanges[aRangeIdx]; i < data.myRanges[aRangeIdx + 1] - 1; ++i)
			{
				if (locIsEdgeInLowerEnvelope<ExecutionPolicy::SEQ>(somePlanes, data, aRangeIdx, i, edge, oppEdge, segmentStartIdx))
				{
					if constexpr (E == ExecutionPolicy::SEQ)
					{
						edges.emplace_back(edge);
						edges.emplace_back(oppEdge);
					}
					else
					{
						CGAL_precondition(aCounterPtr != nullptr);
						const auto outIdx = aCounterPtr->fetch_add(2, std::memory_order_relaxed);
						edges[outIdx] = edge;
						edges[outIdx + 1] = oppEdge;
					}
				}
			}
		};

		if constexpr (E == ExecutionPolicy::SEQ)
		{
			for (size_t rangeIdx = 0; rangeIdx < data.myRanges.size() - 1; ++rangeIdx)
			{
				buildEdges(rangeIdx, nullptr);
			}
		}
		else
		{
			edges.resize(2 * data.myVertices.size());
			std::atomic<size_t> counter{ 0 };

			BindExecutionPolicy<E>(hpx::for_loop, 0, data.myRanges.size() - 1, [&](const size_t aRangeIdx)
				{
					buildEdges(aRangeIdx, &counter);
				});
			edges.resize(counter.load());
		}

		// Sort edges
		locSortEdgesCCW<E>(edges);

		return edges;
	}

	template<ExecutionPolicy E>
	void TriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope)
	{
		CGAL_precondition(std::holds_alternative<std::vector<Edge<Point3>>>(anOutLowerEnvelope));
		auto& edges = std::get<std::vector<Edge<Point3>>>(anOutLowerEnvelope);
		CGAL_precondition(edges.size() > 2);

		BindExecutionPolicy<E>(hpx::sort, edges.begin(), edges.end(), [](const auto& aFirstEdge, const auto& aSecondEdge)
			{
				return aFirstEdge.myLowestLeftPlane < aSecondEdge.myLowestLeftPlane ||
					(aFirstEdge.myLowestLeftPlane == aSecondEdge.myLowestLeftPlane && aFirstEdge.myType < aSecondEdge.myType);
			});

		if constexpr (E == ExecutionPolicy::SEQ)
		{
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
		}
		else
		{
			const auto edgesOldSize = edges.size();
			// Dummy upper-bound
			edges.resize(edgesOldSize * 4);
			std::atomic<size_t> atomicCounter{ edgesOldSize };

			BindExecutionPolicy<E>(hpx::for_loop, 0, edgesOldSize, [&](size_t i) {
				if (i == 0 ||
					edges[i].myLowestLeftPlane != edges[i - 1].myLowestLeftPlane)
				{
					const auto startIdx = i;
					auto prevStartIdx = -1;
					while (++i < edgesOldSize &&
						edges[i].myLowestLeftPlane == edges[startIdx].myLowestLeftPlane)
					{
						if (edges[i].myEnd == edges[startIdx].myStart)
						{
							CGAL_precondition(prevStartIdx == -1);
							prevStartIdx = i;
						}
					}

					const auto newEdgesCount = i - startIdx - 1 - 2 * static_cast<int>(prevStartIdx > -1);
					auto outIdx = atomicCounter.fetch_add(newEdgesCount, std::memory_order_relaxed);
					i = startIdx;

					while (++i < edgesOldSize && edges[i].myLowestLeftPlane == edges[startIdx].myLowestLeftPlane)
					{
						if (prevStartIdx > -1 && (i == prevStartIdx || edges[i].myEnd == edges[prevStartIdx].myStart))
						{
							continue;
						}
						edges[outIdx++] = Edge<Point3>{ edges[startIdx].myStart, edges[i].myEnd, EdgeType::SEGMENT_TRIANGLE };
						edges[outIdx++] = Edge<Point3>{ edges[i].myEnd, edges[startIdx].myStart, EdgeType::SEGMENT_TRIANGLE };
					}
				}
				});
			edges.resize(atomicCounter.load());
		}

		locSortEdgesCCW<E>(edges);
	}

	template<ExecutionPolicy E>
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
			using EdgesList = std::vector<Edge<Point3>>;
			CGAL_precondition(std::holds_alternative<EdgesList>(aLowerEnvelope));
			const auto& edges = std::get<EdgesList>(aLowerEnvelope);

			// Dummy upper bound
			std::vector<Point3> vertices(edges.size() / 2);

			BindExecutionPolicy<E>(hpx::transform, edges.begin(), edges.end(), vertices.begin(),
				[](const auto& anEdge) { return anEdge.myStart; });

			BindExecutionPolicy<E>(hpx::sort, vertices.begin(), vertices.end(), CGAL::Less<Point3, Point3>());

			const auto uniqueEndIter =
				BindExecutionPolicy<E>(hpx::unique, vertices.begin(), vertices.end(), CGAL::Equal_to<Point3, Point3>());

			return static_cast<size_t>(std::distance(vertices.begin(), uniqueEndIter));
		}
	}

	template<ExecutionPolicy E>
	std::vector<size_t> GetLowerEnvelopePlanesIndices(const LowerEnvelope3d& aLowerEnvelope)
	{
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
			std::vector<size_t> indices(edges.size());
			std::atomic<size_t> counter{ 0 };
			BindExecutionPolicy<E>(hpx::for_loop, 0, edges.size(), [&](size_t i) {
				if (edges[i].myType == EdgeType::LINE ||
					edges[i].myType == EdgeType::HALF_EDGE_SF ||
					edges[i].myType == EdgeType::HALF_EDGE_EF ||
					edges[i].myType == EdgeType::SEGMENT)
				{
					CGAL_precondition(edges[i].myLowestLeftPlane > -1);
					indices[counter.fetch_add(1, std::memory_order_relaxed)] = edges[i].myLowestLeftPlane;
				}
				});
			indices.resize(counter.load());
			BindExecutionPolicy<E>(hpx::sort, indices.begin(), indices.end());
			indices.erase(BindExecutionPolicy<E>(hpx::unique, indices.begin(), indices.end()), indices.end());
			return indices;
		}
	}


	// ------------------------------------------------
	// Explicit instantiation for each execution policy
	template LowerEnvelope3d ComputeLowerEnvelope<ExecutionPolicy::SEQ>(const std::vector<Plane>& somePlanes);
	template LowerEnvelope3d ComputeLowerEnvelope<ExecutionPolicy::PAR>(const std::vector<Plane>& somePlanes);
	template LowerEnvelope3d ComputeLowerEnvelope<ExecutionPolicy::PAR_UNSEQ>(const std::vector<Plane>& somePlanes);

	// Explicit instantiation for each execution policy
	template void TriangulateLowerEnvelope<ExecutionPolicy::SEQ>(LowerEnvelope3d& anOutLowerEnvelope);
	template void TriangulateLowerEnvelope<ExecutionPolicy::PAR>(LowerEnvelope3d& anOutLowerEnvelope);
	template void TriangulateLowerEnvelope<ExecutionPolicy::PAR_UNSEQ>(LowerEnvelope3d& anOutLowerEnvelope);

	// Explicit instantiation for each execution policy
	template size_t CountVerticesInLowerEnvelope<ExecutionPolicy::SEQ>(const LowerEnvelope3d& aLowerEnvelope);
	template size_t CountVerticesInLowerEnvelope<ExecutionPolicy::PAR>(const LowerEnvelope3d& aLowerEnvelope);
	template size_t CountVerticesInLowerEnvelope<ExecutionPolicy::PAR_UNSEQ>(const LowerEnvelope3d& aLowerEnvelope);

	// Explicit instantiation for each execution policy
	template std::vector<size_t> GetLowerEnvelopePlanesIndices<ExecutionPolicy::SEQ>(const LowerEnvelope3d& aLowerEnvelope);
	template std::vector<size_t> GetLowerEnvelopePlanesIndices<ExecutionPolicy::PAR>(const LowerEnvelope3d& aLowerEnvelope);
	template std::vector<size_t> GetLowerEnvelopePlanesIndices<ExecutionPolicy::PAR_UNSEQ>(const LowerEnvelope3d& aLowerEnvelope);
}