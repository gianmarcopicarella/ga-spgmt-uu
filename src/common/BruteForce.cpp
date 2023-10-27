#include "BruteForce.h"
#include "Utils.h"
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

		template<ExecutionPolicy E>
		void locTriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope)
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
				auto boundaryStartIt = edges.begin();
				std::vector<Edge<Point3>> temporaryEdges;

				for (size_t i = 1, lastStart = 0, edgesOldSize = edges.size(); i < edgesOldSize; ++i)
				{
					if (edges[lastStart].myLowestLeftPlane != edges[i].myLowestLeftPlane ||
						i == edgesOldSize - 1)
					{
						for (size_t k = lastStart + 1; k < i; ++k)
						{
							edges.emplace_back(Edge<Point3>{ edges[lastStart].myStart, edges[i].myEnd, EdgeType::SEGMENT_TRIANGLE });
							edges.emplace_back(Edge<Point3>{ edges[i].myEnd, edges[lastStart].myStart, EdgeType::SEGMENT_TRIANGLE });
						}
						lastStart = i;
					}
				}
			}
			else
			{
				// Dummy upper bound
				std::vector<size_t> ranges(edges.size() / 2);
				std::atomic<size_t> atomicCounter{ 1 };
				ranges.front() = 0;
				BindExecutionPolicy<E>(hpx::for_loop, 1, edges.size(), [&](size_t anEdgeIdx) {
					if (edges[anEdgeIdx].myLowestLeftPlane != edges[anEdgeIdx - 1].myLowestLeftPlane)
					{
						ranges[atomicCounter.fetch_add(1, std::memory_order_relaxed)] = anEdgeIdx;
					}
					});

				ranges.resize(atomicCounter.fetch_xor(atomicCounter.load()) + 1);
				ranges.back() = edges.size();

				const auto prevEdgesCount = edges.size();
				atomicCounter.store(prevEdgesCount);

				// Dummy upper-bound
				edges.resize(prevEdgesCount * 3);

				BindExecutionPolicy<E>(hpx::for_loop, 0, ranges.size() - 1, [&](size_t aRangeIdx) {
					BindExecutionPolicy<E>(hpx::for_loop, ranges[aRangeIdx] + 1, ranges[aRangeIdx + 1], [&](size_t anEdgeIdx) {
						const auto& startEdge = edges[ranges[aRangeIdx]];
						const auto storeIdx = atomicCounter.fetch_add(2, std::memory_order_relaxed);

						edges[storeIdx].myType = EdgeType::SEGMENT_TRIANGLE;
						edges[storeIdx].myStart = startEdge.myStart;
						edges[storeIdx].myEnd = edges[anEdgeIdx].myEnd;

						edges[storeIdx + 1].myType = EdgeType::SEGMENT_TRIANGLE;
						edges[storeIdx + 1].myStart = edges[anEdgeIdx].myEnd;
						edges[storeIdx + 1].myEnd = startEdge.myStart;
						});
					});

				edges.resize(atomicCounter.fetch_xor(atomicCounter.load()) - prevEdgesCount);
			}

			locSortEdgesCCW<E>(edges);
		}

		void locInsertInfinityPoints(
			const std::tuple<size_t, Point3>& aMinPoint,
			const std::tuple<size_t, Point3>& aMaxPoint,
			const Line3& aFirstLine,
			const Line3& aSecondLine,
			std::vector<std::tuple<size_t, Point3>>::iterator anOutIter)
		{
			constexpr auto distance = 1000.f;
			const auto& minusInfPoint = std::get<1>(aMinPoint) - aFirstLine.to_vector() * distance;
			const auto& plusInfPoint = std::get<1>(aMaxPoint) + aSecondLine.to_vector() * distance;

			CGAL_precondition(aFirstLine.has_on(minusInfPoint));
			CGAL_precondition(aSecondLine.has_on(plusInfPoint));

			*anOutIter = std::make_tuple(std::get<0>(aMinPoint), minusInfPoint);
			*(anOutIter + 1) = std::make_tuple(std::get<0>(aMaxPoint), plusInfPoint);
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
	}

	namespace Parallel
	{
		void locComputeData(const std::vector<Plane>& somePlanes, BruteForceData& anOutData)
		{
			static const Point3 zero{ 0,0,0 };
			constexpr auto policy = ExecutionPolicy::PAR_UNSEQ;
			const auto maxLinesCount = somePlanes.size() * (somePlanes.size());
			using LineWrapper = std::tuple<STATUS, Point3, Point3, size_t>;

			std::vector<LineWrapper> lines(maxLinesCount);

			// Find intersection lines
			hpx::for_loop(hpx::execution::par_unseq, 0, somePlanes.size(), [&](size_t aPlaneIdx)
				{
					//auto startOutIdx = outStart[aPlaneIdx];
					hpx::for_loop(hpx::execution::par_unseq, aPlaneIdx + 1, somePlanes.size(), [&](size_t anotherPlaneIdx)
						{
							const auto intersection =
								CGAL::intersection(somePlanes[aPlaneIdx], somePlanes[anotherPlaneIdx]);
							if (intersection)
							{
								const Line3* line = boost::get<Line3>(&*intersection);
								CGAL_precondition(line != nullptr);

								auto start = line->projection(zero);
								auto end = start + line->to_vector() * FT(1000.f);
								if (start > end) std::swap(start, end);
								CGAL_precondition(start < end);

								const auto outFirstIdx = aPlaneIdx * somePlanes.size() + anotherPlaneIdx;
								const auto outSecondIdx = anotherPlaneIdx * somePlanes.size() + aPlaneIdx;

								CGAL_precondition(std::get<0>(lines[outFirstIdx]) == STATUS::NONE);
								CGAL_precondition(std::get<0>(lines[outSecondIdx]) == STATUS::NONE);

								lines[outFirstIdx] = std::make_tuple(STATUS::INIT, start, end, aPlaneIdx);
								lines[outSecondIdx] = std::make_tuple(STATUS::INIT, start, end, anotherPlaneIdx);
							}
						});
				});

			// Sort and keep only unique lines
			hpx::sort(hpx::execution::par_unseq, lines.begin(), lines.end());

			Utils::TrimToLastValidItem(lines);

			const auto uniqueEntriesCount =
				static_cast<size_t>(std::distance(
					lines.begin(),
					hpx::unique(hpx::execution::par_unseq, lines.begin(), lines.end()))) - 1;

			// If all the lines are parallel then return lower envelope
			const auto areAllLinesParallel = Utils::AreItemsParallel<policy>(lines,
				[](const auto& aLine) { return Dir3{ std::get<2>(aLine) - std::get<1>(aLine) }; });
			if (areAllLinesParallel)
			{
				static const Point3 zero{ 0,0,0 };

				const auto virtualLinesEnd = lines.begin() + uniqueEntriesCount;

				auto lineIter =
					hpx::find_if(hpx::execution::par_unseq, lines.begin(), virtualLinesEnd,
						[&](const auto& aLine) { return locIsVertexInLowerEnvelope<policy>(somePlanes, std::get<1>(aLine)); });

				while (lineIter != lines.begin() && std::get<3>(*lineIter) == std::get<3>(*(lineIter - 1))) { --lineIter; }

				CGAL_precondition(lineIter != virtualLinesEnd);
				const auto& lowerEnvelopeLine = Line3{ std::get<1>(*lineIter), std::get<2>(*lineIter) };

				std::vector<size_t> planesIndices;
				for (auto iter = lineIter;
					iter != virtualLinesEnd &&
					std::get<3>(*iter) == std::get<3>(*lineIter); ++iter)
				{
					planesIndices.emplace_back(std::get<3>(*iter));
				}

				const auto minMaxPlanes = locGetMinMaxSteepPlaneIndices<policy>(lowerEnvelopeLine, somePlanes, planesIndices);

				constexpr auto distance = 1000.f;
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

			// Collect min-max planes and unique lines
			using UniqueLineWrapper = std::tuple<STATUS, Line3, std::tuple<size_t, size_t>>;
			std::vector<UniqueLineWrapper> uniqueLines(uniqueEntriesCount);

			hpx::for_loop(hpx::execution::par_unseq, 0, uniqueEntriesCount, [&](size_t aLineIdx) {
				if (aLineIdx == 0 ||
					std::get<1>(lines[aLineIdx]) != std::get<1>(lines[aLineIdx - 1]) ||
					std::get<2>(lines[aLineIdx]) != std::get<2>(lines[aLineIdx - 1]))
				{
					std::vector<size_t> planesIndices;
					planesIndices.emplace_back(std::get<3>(lines[aLineIdx]));

					for (size_t i = aLineIdx + 1;
						i < lines.size() &&
						std::get<1>(lines[aLineIdx]) == std::get<1>(lines[i]) &&
						std::get<2>(lines[aLineIdx]) == std::get<2>(lines[i]); ++i)
					{
						planesIndices.emplace_back(std::get<3>(lines[i]));
					}

					const auto& outLine = Line3{ std::get<1>(lines[aLineIdx]), std::get<2>(lines[aLineIdx]) };
					const auto& minMaxPlanes = locGetMinMaxSteepPlaneIndices<ExecutionPolicy::SEQ>(outLine, somePlanes, planesIndices);

					CGAL_precondition(std::get<0>(uniqueLines[aLineIdx]) == STATUS::NONE);

					uniqueLines[aLineIdx] = std::make_tuple(STATUS::INIT, outLine, minMaxPlanes);
				}
				});

			// Sort and keep only unique lines
			hpx::sort(hpx::execution::par_unseq, uniqueLines.begin(), uniqueLines.end(),
				[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });

			Utils::TrimToLastValidItem(uniqueLines);

			// Compute vertices
			const auto maxVerticesCount = uniqueLines.size() * (uniqueLines.size());

			using VertexWrapper = std::tuple<STATUS, size_t, Point3>;
			std::vector<VertexWrapper> vertices(maxVerticesCount);

			hpx::for_loop(hpx::execution::par_unseq, 0, uniqueLines.size(), [&](size_t aLineIdx) {
				hpx::for_loop(hpx::execution::par_unseq, aLineIdx + 1, uniqueLines.size(), [&](size_t anotherLineIdx) {
					const auto intersection = CGAL::intersection(std::get<1>(uniqueLines[aLineIdx]), std::get<1>(uniqueLines[anotherLineIdx]));
					if (intersection)
					{
						const Point3* point = boost::get<Point3>(&*intersection);
						CGAL_precondition(point != nullptr);
						//const auto outIndex = atomicCount.fetch_add(2, std::memory_order_relaxed);

						const auto outFirstIdx = aLineIdx * uniqueLines.size() + anotherLineIdx;
						const auto outSecondIdx = anotherLineIdx * uniqueLines.size() + aLineIdx;

						vertices[outFirstIdx] = std::make_tuple(STATUS::INIT, aLineIdx, *point);
						vertices[outSecondIdx] = std::make_tuple(STATUS::INIT, anotherLineIdx, *point);
					}
					});
				});

			// Sort vertices NEW
			hpx::sort(hpx::execution::par_unseq, vertices.begin(), vertices.end());

			Utils::TrimToLastValidItem(vertices);

			const auto endUniqueIter =
				hpx::unique(hpx::execution::par_unseq, vertices.begin(), vertices.end());
			vertices.erase(endUniqueIter, vertices.end());


			anOutData.myVertices.resize(vertices.size() + 2 * uniqueLines.size());
			hpx::transform(hpx::execution::par_unseq, vertices.begin(), vertices.end(), anOutData.myVertices.begin(),
				[](const auto& a) { return std::make_tuple(std::get<1>(a), std::get<2>(a)); });

			//std::cout << "Done sorting" << std::endl;

			anOutData.myRanges.resize(uniqueLines.size() + 1);
			anOutData.myRanges.front() = 0;
			anOutData.myRanges.back() = anOutData.myVertices.size();

			const auto finiteVerticesCount = vertices.size();

			// Compute points at infinity
			{
				constexpr auto distance = 1000.f;
				const auto& minusInfPoint = std::get<1>(anOutData.myVertices[0]) - std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[0])]).to_vector() * distance;
				const auto& plusInfPoint = std::get<1>(anOutData.myVertices[vertices.size() - 1]) + std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[vertices.size() - 1])]).to_vector() * distance;

				CGAL_precondition(std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[0])]).has_on(minusInfPoint));
				CGAL_precondition(std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[vertices.size() - 1])]).has_on(plusInfPoint));

				const auto minusInfOutIdx = finiteVerticesCount + 2 * std::get<0>(anOutData.myVertices[0]);
				const auto plusInfOutIdx = finiteVerticesCount + 2 * std::get<0>(anOutData.myVertices[vertices.size() - 1]) + 1;

				anOutData.myVertices[minusInfOutIdx] = std::make_tuple(std::get<0>(anOutData.myVertices[0]), minusInfPoint);
				anOutData.myVertices[plusInfOutIdx] = std::make_tuple(std::get<0>(anOutData.myVertices[vertices.size() - 1]), plusInfPoint);
			}

			hpx::for_loop(hpx::execution::par_unseq, 1, finiteVerticesCount - 1, [&](size_t aVertexIdx) {
				if (std::get<0>(anOutData.myVertices[aVertexIdx]) <
					std::get<0>(anOutData.myVertices[aVertexIdx + 1]))
				{
					constexpr auto distance = 1000.f;
					const auto& minusInfPoint = std::get<1>(anOutData.myVertices[aVertexIdx + 1]) - std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[aVertexIdx + 1])]).to_vector() * distance;
					const auto& plusInfPoint = std::get<1>(anOutData.myVertices[aVertexIdx]) + std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[aVertexIdx])]).to_vector() * distance;

					CGAL_precondition(std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[aVertexIdx + 1])]).has_on(minusInfPoint));
					CGAL_precondition(std::get<1>(uniqueLines[std::get<0>(anOutData.myVertices[aVertexIdx])]).has_on(plusInfPoint));

					const auto minusInfOutIdx = finiteVerticesCount + 2 * std::get<0>(anOutData.myVertices[aVertexIdx + 1]);
					const auto plusInfOutIdx = finiteVerticesCount + 2 * std::get<0>(anOutData.myVertices[aVertexIdx]) + 1;

					anOutData.myVertices[minusInfOutIdx] = std::make_tuple(std::get<0>(anOutData.myVertices[aVertexIdx + 1]), minusInfPoint);
					anOutData.myVertices[plusInfOutIdx] = std::make_tuple(std::get<0>(anOutData.myVertices[aVertexIdx]), plusInfPoint);

					anOutData.myRanges[std::get<0>(anOutData.myVertices[aVertexIdx + 1])] = aVertexIdx;
				}
				});

			//std::cout << "1" << std::endl;

			anOutData.myMinMaxPlanes.resize(uniqueLines.size());
			hpx::transform(hpx::execution::par_unseq, uniqueLines.begin(), uniqueLines.end(), anOutData.myMinMaxPlanes.begin(),
				[](const auto& a) { return std::get<2>(a); });

			hpx::sort(hpx::execution::par_unseq, anOutData.myVertices.begin(), anOutData.myVertices.end());
			hpx::sort(hpx::execution::par_unseq, anOutData.myRanges.begin(), anOutData.myRanges.end());
		}
	}

	namespace Sequential
	{
		void locComputeData(const std::vector<Plane>& somePlanes, BruteForceData& anOutData)
		{
			static const Point3 zero{ 0,0,0 };
			constexpr auto policy = ExecutionPolicy::SEQ;
			constexpr auto distance = 1000.f;
			using LineWrapper = std::tuple<Point3, Point3, size_t>;
			std::vector<LineWrapper> lines;


			// Find intersection lines
			for (size_t i = 0; i < somePlanes.size(); ++i)
			{
				for (size_t k = i + 1; k < somePlanes.size(); ++k)
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

						lines.emplace_back(std::make_tuple(start, end, i));
						lines.emplace_back(std::make_tuple(start, end, k));
					}
				}
			}

			std::sort(lines.begin(), lines.end());
			const auto uniqueEntriesCount =
				static_cast<size_t>(std::distance(lines.begin(), std::unique(lines.begin(), lines.end())));

			// If all the lines are parallel then return lower envelope
			const auto areAllLinesParallel = Utils::AreItemsParallel<ExecutionPolicy::SEQ>(lines,
				[](const auto& aLine) { return Dir3{ std::get<1>(aLine) - std::get<0>(aLine) }; });
			if (areAllLinesParallel)
			{
				const auto virtualLinesEnd = lines.begin() + uniqueEntriesCount;

				auto lineIter =
					std::find_if(lines.begin(), virtualLinesEnd,
						[&](const auto& aLine) { return locIsVertexInLowerEnvelope<policy>(somePlanes, std::get<0>(aLine)); });

				CGAL_precondition(lineIter != virtualLinesEnd);
				const auto& lowerEnvelopeLine = Line3{ std::get<0>(*lineIter), std::get<1>(*lineIter) };

				std::vector<size_t> planesIndices;
				for (auto iter = lineIter;
					iter != virtualLinesEnd &&
					std::get<2>(*iter) == std::get<2>(*lineIter); ++iter)
				{
					planesIndices.emplace_back(std::get<2>(*iter));
				}

				const auto minMaxPlanes = locGetMinMaxSteepPlaneIndices<policy>(lowerEnvelopeLine, somePlanes, planesIndices);

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
			uniqueLines.emplace_back(Line3{ std::get<0>(lines[0]), std::get<1>(lines[1]) });

			{
				std::vector<size_t> uniquePlanes;
				uniquePlanes.emplace_back(std::get<2>(lines[2]));

				for (size_t i = 1, lastIdx = 0; i < uniqueEntriesCount; ++i)
				{
					if (std::get<0>(lines[lastIdx]) != std::get<0>(lines[i]) ||
						std::get<1>(lines[lastIdx]) != std::get<1>(lines[i]))
					{
						anOutData.myMinMaxPlanes.emplace_back(
							locGetMinMaxSteepPlaneIndices<policy>(uniqueLines.back(), somePlanes, uniquePlanes)
						);
						uniqueLines.emplace_back(Line3{ std::get<0>(lines[i]), std::get<1>(lines[i]) });
						uniquePlanes.clear();
						lastIdx = i;
					}

					uniquePlanes.emplace_back(std::get<2>(lines[i]));
				}

				anOutData.myMinMaxPlanes.emplace_back(
					locGetMinMaxSteepPlaneIndices<policy>(uniqueLines.back(), somePlanes, uniquePlanes)
				);
			}

			// Compute vertices
			for (size_t i = 0; i < uniqueLines.size(); ++i)
			{
				for (size_t k = i + 1; k < uniqueLines.size(); ++k)
				{
					const auto intersection =
						CGAL::intersection(uniqueLines[i], uniqueLines[k]);
					if (intersection)
					{
						const Point3* point = boost::get<Point3>(&*intersection);
						CGAL_precondition(point != nullptr);
						anOutData.myVertices.emplace_back(std::make_tuple(i, *point));
						anOutData.myVertices.emplace_back(std::make_tuple(k, *point));
					}
				}
			}

			// Sort vertices
			std::sort(anOutData.myVertices.begin(), anOutData.myVertices.begin());

			// Compute points at infinity
			{
				anOutData.myVertices.emplace_back(); anOutData.myVertices.emplace_back();
				locInsertInfinityPoints(
					anOutData.myVertices[0],
					*(anOutData.myVertices.end() - 3),
					uniqueLines[std::get<0>(anOutData.myVertices[0])],
					uniqueLines[std::get<0>(*(anOutData.myVertices.end() - 3))],
					anOutData.myVertices.end() - 2);
			}

			const auto prevSize = anOutData.myVertices.size() - 3;
			for (size_t i = 1; i < prevSize; ++i)
			{
				if (std::get<0>(anOutData.myVertices[i]) <
					std::get<0>(anOutData.myVertices[i + 1]))
				{
					anOutData.myVertices.emplace_back(); anOutData.myVertices.emplace_back();
					locInsertInfinityPoints(
						anOutData.myVertices[i + 1],
						anOutData.myVertices[i],
						uniqueLines[std::get<0>(anOutData.myVertices[i + 1])],
						uniqueLines[std::get<0>(anOutData.myVertices[i])],
						anOutData.myVertices.end() - 2);
				}
			}

			std::sort(anOutData.myVertices.begin(), anOutData.myVertices.end());

			// Remove duplicates and store line ranges
			auto first = anOutData.myVertices.begin();
			const auto last = anOutData.myVertices.end();
			auto result = first;
			auto lastLineIdx = std::get<0>(*first);
			auto vertexIdx = 0;
			anOutData.myRanges.emplace_back(0);

			while (++first != last)
			{
				if (!(*result == *first))
				{
					++vertexIdx;
					if (++result != first)
					{
						*result = std::move(*first);
					}

					const auto currLineIdx = std::get<0>(*result);
					if (lastLineIdx != currLineIdx)
					{
						lastLineIdx = currLineIdx;
						anOutData.myRanges.emplace_back(vertexIdx);
					}
				}
			}

			anOutData.myVertices.erase(++result, anOutData.myVertices.end());
			anOutData.myRanges.emplace_back(anOutData.myVertices.size());
		}
	}

	LowerEnvelope3d ParallelComputeLowerEnvelope(const std::vector<Plane>& somePlanes)
	{
		using namespace Parallel;
		constexpr auto policy = ExecutionPolicy::PAR_UNSEQ;
		constexpr auto distance = 1000.f;

		if (somePlanes.empty())
		{
			return LowerEnvelope3d{};
		}

		// Requirements check (TEMPORARY DISABLED)
		//CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		//CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));

		const auto arePlanesParallel = Utils::AreItemsParallel<policy>(somePlanes,
			[](const Plane& aPlane) { return aPlane.orthogonal_direction(); });

		// Edge case: all planes are parallel
		if (arePlanesParallel)
		{
			static const Point3 origin{ 0,0,0 };
			std::vector<size_t> indices(somePlanes.size());
			std::iota(indices.begin(), indices.end(), 0);
			const auto minMaxPlanes = locGetMinMaxPlaneAtPoint<policy>(somePlanes, indices, origin);
			CGAL_precondition(std::get<0>(minMaxPlanes) != -1);
			return std::get<0>(minMaxPlanes);
		}

		BruteForceData data;
		locComputeData(somePlanes, data);

		if (data.myLowerEnvelopeOpt.has_value())
		{
			return data.myLowerEnvelopeOpt.value();
		}

		// Construct result
		using EdgeWrapper = std::tuple<STATUS, Edge<Point3>>;
		std::vector<EdgeWrapper> edgesWrap(2 * data.myVertices.size());

		hpx::for_loop(hpx::execution::par_unseq, 0, data.myRanges.size() - 1, [&](const size_t aRangeIdx)
			{
				int segmentStartIdx = -1;
				Edge<Point3> edge, oppEdge;
				for (size_t i = data.myRanges[aRangeIdx]; i < data.myRanges[aRangeIdx + 1] - 1; ++i)
				{
					if (locIsEdgeInLowerEnvelope<ExecutionPolicy::SEQ>(somePlanes, data, aRangeIdx, i, edge, oppEdge, segmentStartIdx))
					{
						edgesWrap[2 * i] = std::make_tuple(STATUS::INIT, edge);
						edgesWrap[2 * i + 1] = std::make_tuple(STATUS::INIT, oppEdge);
					}
				}
			});

		// Sort vertices NEW
		hpx::sort(hpx::execution::par_unseq, edgesWrap.begin(), edgesWrap.end(),
			[](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });

		Utils::TrimToLastValidItem(edgesWrap);

		std::vector<Edge<Point3>> edges(edgesWrap.size());
		hpx::transform(hpx::execution::par_unseq, edgesWrap.begin(), edgesWrap.end(), edges.begin(),
			[](const auto& a) { return std::get<1>(a); });

		// Sort edges
		locSortEdgesCCW<policy>(edges);

		return edges;
	}

	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes)
	{
		using namespace Sequential;
		constexpr auto policy = ExecutionPolicy::SEQ;
		constexpr auto distance = 1000.f;

		if (somePlanes.empty())
		{
			return LowerEnvelope3d{};
		}

		// Requirements check (TEMPORARY)
		//CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		//CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));

		const auto arePlanesParallel = Utils::AreItemsParallel<policy>(somePlanes,
			[](const Plane& aPlane) { return aPlane.orthogonal_direction(); });

		// Edge case: all planes are parallel
		if (arePlanesParallel)
		{
			static const Point3 origin{ 0,0,0 };
			std::vector<size_t> indices(somePlanes.size());
			std::iota(indices.begin(), indices.end(), 0);
			const auto minMaxPlanes = locGetMinMaxPlaneAtPoint<policy>(somePlanes, indices, origin);
			CGAL_precondition(std::get<0>(minMaxPlanes) != -1);
			return std::get<0>(minMaxPlanes);
		}

		BruteForceData data;
		locComputeData(somePlanes, data);

		// Edge case: no triple of planes intersect
		if (data.myLowerEnvelopeOpt.has_value())
		{
			return data.myLowerEnvelopeOpt.value();
		}

		// Construct result
		std::vector<Edge<Point3>> edges;

		for (size_t rangeIdx = 0; rangeIdx < data.myRanges.size() - 1; ++rangeIdx)
		{
			int segmentStartIdx = -1;
			Edge<Point3> edge, oppEdge;
			for (size_t i = data.myRanges[rangeIdx]; i < data.myRanges[rangeIdx + 1] - 1; ++i)
			{
				if (locIsEdgeInLowerEnvelope<policy>(somePlanes, data, rangeIdx, i, edge, oppEdge, segmentStartIdx))
				{
					edges.emplace_back(edge);
					edges.emplace_back(oppEdge);
				}
			}
		}

		// Sort edges
		locSortEdgesCCW<policy>(edges);

		return edges;
	}

	void TriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope)
	{
		constexpr auto policy = ExecutionPolicy::SEQ;
		locTriangulateLowerEnvelope<policy>(anOutLowerEnvelope);
	}

	void ParallelTriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope)
	{
		constexpr auto policy = ExecutionPolicy::PAR_UNSEQ;
		locTriangulateLowerEnvelope<policy>(anOutLowerEnvelope);
	}

	size_t CountUniqueVertices(const LowerEnvelope3d& aLowerEnvelope)
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

			hpx::transform(hpx::execution::seq, edges.begin(), edges.end(), vertices.begin(),
				[](const auto& anEdge) { return anEdge.myStart; });
			hpx::sort(hpx::execution::seq, vertices.begin(), vertices.end(), CGAL::Less<Point3, Point3>());
			const auto uniqueEndIter =
				hpx::unique(hpx::execution::seq, vertices.begin(), vertices.end(), CGAL::Equal_to<Point3, Point3>());

			return static_cast<size_t>(std::distance(vertices.begin(), uniqueEndIter));
		}
	}
}