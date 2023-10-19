#include "BruteForce.h"
#include "Utils.h"
#include <CGAL/Polygon_2_algorithms.h>

#include <hpx/hpx.hpp>

#include <atomic>
namespace SPGMT
{
	namespace Common
	{
		struct LineWrapper
		{
			Point3 myPoint;
			Dir3 myDir;
			size_t myPlanes[2];
		};

		template<typename LData>
		Vec3 locGetUniformLineVector(const LData& aLineData)
		{
			constexpr auto distance = 1000.f;
			const auto& a = aLineData.myLine.point();
			const auto& b = a + aLineData.myLine.to_vector() * distance;
			if (a < b)
			{
				return Vec3{ a, b };
			}
			else
			{
				return Vec3{ b, a };
			}
		}

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

		template<ExecutionPolicy E>
		std::pair<FT, size_t> locGetLowestPlaneAtPoint(const std::vector<Plane>& somePlanes, const Point3& aPoint)
		{
			CGAL_precondition(somePlanes.size() > 0);

			static const Vec3 up{ 0,0,1 };
			const Line3 upLine{ aPoint, up };

			using ValueIndexPair = std::pair<FT, size_t>;
			std::vector<ValueIndexPair> pairs(somePlanes.size());

			const auto computePlaneHeight = [&](const size_t& aPlaneIdx)
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes[aPlaneIdx]);
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				pairs[aPlaneIdx] = std::make_pair(point->z(), aPlaneIdx);
			};

			constexpr auto pairsReduction = [](const auto& aFirst, const auto& aSecond)
			{
				return aFirst.first < aSecond.first ? aFirst : aSecond;
			};

			BindExecutionPolicy<E>(
				std::bind(hpx::for_loop, std::placeholders::_1, 0, somePlanes.size(), computePlaneHeight));

			const auto lowestPlane = BindExecutionPolicy<E>(
				std::bind(hpx::reduce, std::placeholders::_1, std::next(pairs.begin()), pairs.end(), pairs[0], pairsReduction));

			return lowestPlane;
		}

		template<ExecutionPolicy E>
		std::pair<FT, size_t> locGetLowestPlaneAtPointFromIndices(
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& somePlaneIndices,
			const Point3& aPoint)
		{
			CGAL_precondition(somePlanes.size() > 0);

			std::vector<Plane> selectedPlanes;
			BindExecutionPolicy<E>(std::bind(hpx::transform, somePlaneIndices.begin(), somePlaneIndices.end(),
				std::back_inserter(selectedPlanes), [&](size_t aPlaneIdx) {
					return somePlanes[aPlaneIdx];
				}));

			return locGetLowestPlaneAtPoint<E>(selectedPlanes, aPoint);
		}

		template<ExecutionPolicy E>
		bool locIsVertexInLowerEnvelope(const std::vector<Plane>& somePlanes, const Point3& aPoint)
		{
			const auto& lowestPlane = locGetLowestPlaneAtPoint<E>(somePlanes, aPoint);
			const auto comparisonResult = CGAL::compare(aPoint.z(), lowestPlane.first);
			return comparisonResult <= CGAL::Comparison_result::EQUAL;
		}

		// Visual Help. https://www.falstad.com/dotproduct/ | https://www.geogebra.org/3d?lang=en
		template<ExecutionPolicy E>
		size_t locMinSteepPlaneIndexThroughSegment(const Point3& aStart, const Point3& anEnd, const std::vector<Plane>& somePlanes, const std::vector<size_t>& somePlanesIndices)
		{
			static const Vec3 up{ 0,0,1 };

			CGAL_precondition(somePlanes.size() > 0);
			CGAL_precondition(somePlanesIndices.size() > 0);

			const Vec3 segmentDir{ aStart, anEnd };
			const Vec3 leftDir{ -CGAL::cross_product(segmentDir, up) };
			const auto sampleLoc = CGAL::midpoint(aStart, anEnd) + leftDir * 1000.f;

			CGAL_precondition((somePlanes.size() > 0 && somePlanesIndices.size() > 0));

			return locGetLowestPlaneAtPointFromIndices<E>(somePlanes, somePlanesIndices, sampleLoc).second;
		}

		void locSortNeighboursCCW(
			const std::vector<Edge<Point3>>::iterator aStartIter,
			const std::vector<Edge<Point3>>::iterator anEndIter)
		{
			static const Vec3 ref{ 1, 0, 0 };
			std::sort(aStartIter, anEndIter,
				[&](auto& aFirst, auto& aSecond) {
					CGAL_precondition(aFirst.myStart == aSecond.myStart);
					const auto firstAngle = CGAL::approximate_angle(ref, Vec3{ aFirst.myStart, aFirst.myEnd });
					const auto secondAngle = CGAL::approximate_angle(ref, Vec3{ aFirst.myStart, aSecond.myEnd });
					return firstAngle < secondAngle;
				});
		}

		template<ExecutionPolicy E>
		void locSortLowerEnvelopeCCW(std::vector<Edge<Point3>>& someEdges)
		{
			const auto sortEdges = [](const auto& aFirstEdge, const auto& aSecondEdge)
			{
				return aFirstEdge.myStart < aSecondEdge.myStart;
			};

			BindExecutionPolicy<E>(
				std::bind(hpx::sort, std::placeholders::_1, someEdges.begin(), someEdges.end(), sortEdges));

			auto iterStart = someEdges.begin();

			for (auto iter = someEdges.begin(); iter != someEdges.end(); ++iter)
			{
				if (iterStart->myStart != iter->myStart)
				{
					locSortNeighboursCCW(iterStart, iter);
					iterStart = iter;
				}
			}

			locSortNeighboursCCW(iterStart, someEdges.end());
		}
	}

	namespace Parallel
	{
		struct LineData
		{
			Line3 myLine;
			std::vector<size_t> myPlanesIndices;
			std::vector<Point3> myVertices;
		};
		struct LinesAndVerticesData
		{
			std::vector<LineData> myUniqueLines;
			bool myDoLinesIntersect{ false };
		};

		using namespace Common;
		void locComputeLines(const std::vector<Plane>& somePlanes, LinesAndVerticesData& anOutResult)
		{
			static const Point3 zero{ 0,0,0 };
			const auto maximumLinesCount = somePlanes.size() * (somePlanes.size() - 1) / 2;
			std::vector<LineWrapper> lines(maximumLinesCount);
			std::atomic<size_t> atomicCount{ 0 };

			// Find intersection lines
			hpx::for_loop(hpx::execution::par, 0, somePlanes.size(), [&](size_t aPlaneIdx)
				{
					hpx::for_loop(hpx::execution::par, aPlaneIdx + 1, somePlanes.size(), [&](size_t anotherPlaneIdx)
						{
							const auto intersection =
								CGAL::intersection(somePlanes[aPlaneIdx], somePlanes[anotherPlaneIdx]);
							if (intersection)
							{
								const Line3* line = boost::get<Line3>(&*intersection);
								CGAL_precondition(line != nullptr);
								auto& lineData = lines[atomicCount.fetch_add(1)];
								lineData.myPoint = line->projection(zero);
								lineData.myDir = line->direction();
								lineData.myPlanes[0] = aPlaneIdx;
								lineData.myPlanes[1] = anotherPlaneIdx;
							}
						});
				});
			lines.resize(atomicCount.fetch_xor(atomicCount.load()));

			// Sort intersected lines
			hpx::sort(hpx::execution::par, lines.begin(), lines.end(), [&](const auto& aFirst, const auto& aSecond) {
				return aFirst.myPoint < aSecond.myPoint ||
					(aFirst.myPoint == aSecond.myPoint && aFirst.myDir == aSecond.myDir);
				});

			// Extract ranges
			std::vector<size_t> ranges(lines.size());
			hpx::for_loop(hpx::execution::par, 0, lines.size(), [&](size_t aLineIdx) {
				if (aLineIdx == 0 ||
					Line3{ lines[aLineIdx - 1].myPoint, lines[aLineIdx - 1].myDir } !=
					Line3{ lines[aLineIdx].myPoint, lines[aLineIdx].myDir })
				{
					ranges[atomicCount.fetch_add(1)] = aLineIdx;
				}
				});
			ranges.resize(atomicCount.fetch_xor(atomicCount.load()) + 1);
			ranges.back() = lines.size();

			// Sort ranges
			hpx::sort(hpx::execution::par, ranges.begin(), ranges.end());

			// Process ranges
			anOutResult.myUniqueLines.resize(ranges.size() - 1);
			hpx::for_loop(hpx::execution::par, 0, ranges.size() - 1, [&](size_t aRangeIdx) {
				const auto startLineIdx = ranges[aRangeIdx];
				auto& outLine = anOutResult.myUniqueLines[atomicCount.fetch_add(1)];
				outLine.myLine = Line3{ lines[startLineIdx].myPoint, lines[startLineIdx].myDir };
				std::unordered_set<int> uniquePlanesIndices;
				for (size_t i = startLineIdx;
					i < ranges[aRangeIdx + 1]; ++i)
				{
					uniquePlanesIndices.insert(lines[i].myPlanes[0]);
					uniquePlanesIndices.insert(lines[i].myPlanes[1]);
				}
				outLine.myPlanesIndices.insert(
					outLine.myPlanesIndices.end(), uniquePlanesIndices.begin(), uniquePlanesIndices.end());
				});
		}
		void locComputeVertices(LinesAndVerticesData& anOutResult)
		{
			std::atomic<size_t> atomicCount{ 0 };

			hpx::for_loop(hpx::execution::par, 0, anOutResult.myUniqueLines.size(), [&](size_t aLineIdx) {
				auto& vertices = anOutResult.myUniqueLines[aLineIdx].myVertices;
				const auto lineInnerLoop = [&](const size_t anotherLineIdx)
				{
					const auto intersection = CGAL::intersection(
						anOutResult.myUniqueLines[aLineIdx].myLine, anOutResult.myUniqueLines[anotherLineIdx].myLine);
					if (intersection)
					{
						const Point3* point = boost::get<Point3>(&*intersection);
						CGAL_precondition(point != nullptr);
						vertices.emplace_back(*point);
						atomicCount.fetch_add(1);
					}
				};

				hpx::for_loop(hpx::execution::seq, 0, aLineIdx, lineInnerLoop);
				hpx::for_loop(hpx::execution::seq, aLineIdx + 1, anOutResult.myUniqueLines.size(), lineInnerLoop);
				hpx::sort(hpx::execution::seq, vertices.begin(), vertices.end(), CGAL::Less<Point3, Point3>());

				const auto newEndIter = hpx::unique(hpx::execution::seq, vertices.begin(), vertices.end());
				vertices.erase(newEndIter, vertices.end());
				});

			anOutResult.myDoLinesIntersect = atomicCount.load() > 0;
		}
	}

	namespace Sequential
	{
		struct LineData
		{
			Line3 myLine;
			std::vector<size_t> myPlanesIndices;
			std::vector<size_t> myVerticesIndices;
		};
		struct LinesAndVerticesData
		{
			std::vector<LineData> myUniqueLines;
			std::vector<Point3> myUniqueVertices;
		};

		std::vector<Edge<Point3>> locTriangulateConvexFace(
			const std::vector<Edge<Point3>>::iterator aStartIter,
			const std::vector<Edge<Point3>>::iterator anEndIter)
		{
			CGAL_precondition(std::distance(aStartIter, anEndIter) > 1);

			std::vector<Edge<Point3>> edges;

			for (auto it = aStartIter + 1; it != anEndIter; ++it)
			{
				if (aStartIter->myStart == it->myEnd)
				{
#ifndef NDEBUG
					const auto isBoundedFace = aStartIter->myType == EdgeType::SEGMENT;
					CGAL_postcondition(isBoundedFace);
#endif			
					continue;
				}

				Edge<Point3> edge, oppositeEdge;

				edge.myType = EdgeType::SEGMENT_TRIANGLE;
				edge.myStart = aStartIter->myStart;
				edge.myEnd = it->myEnd;

				oppositeEdge.myType = EdgeType::SEGMENT_TRIANGLE;
				oppositeEdge.myStart = it->myEnd;
				oppositeEdge.myEnd = aStartIter->myStart;

				edges.emplace_back(edge);
				edges.emplace_back(oppositeEdge);
			}
			return edges;
		}

		using namespace Common;
		void locComputeLines(const std::vector<Plane>& somePlanes, LinesAndVerticesData& anOutResult)
		{
			static const Point3 zero{ 0,0,0 };
			std::vector<LineWrapper> lines;
			for (size_t i = 0; i < somePlanes.size(); ++i)
			{
				for (size_t k = i + 1; k < somePlanes.size(); ++k)
				{
					const auto intersection = CGAL::intersection(somePlanes[i], somePlanes[k]);
					if (intersection)
					{
						const Line3* line = boost::get<Line3>(&*intersection);
						CGAL_precondition(line != nullptr);
						lines.emplace_back(LineWrapper{ line->projection(zero), line->direction(), i, k });
					}
				}
			}

			hpx::sort(hpx::execution::seq, lines.begin(), lines.end(), [&](const auto& aFirst, const auto& aSecond) {
				return aFirst.myPoint < aSecond.myPoint ||
					(aFirst.myPoint == aSecond.myPoint && aFirst.myDir == aSecond.myDir);
				});

			if (lines.size() > 0)
			{
				std::unordered_set<size_t> uniquePlanes;
				Line3 uniqueLine{ lines.front().myPoint, lines.front().myDir };
				anOutResult.myUniqueLines.emplace_back(LineData{ uniqueLine });

				for (const auto& line : lines)
				{
					const Line3 maybeUniqueLine{ line.myPoint, line.myDir };
					if (anOutResult.myUniqueLines.back().myLine != maybeUniqueLine)
					{
						{
							auto& planesIndices = anOutResult.myUniqueLines.back().myPlanesIndices;
							planesIndices.insert(planesIndices.end(), uniquePlanes.begin(), uniquePlanes.end());
							uniquePlanes.clear();
						}
						anOutResult.myUniqueLines.emplace_back(LineData{ maybeUniqueLine });
					}

					uniquePlanes.insert(line.myPlanes[0]);
					uniquePlanes.insert(line.myPlanes[1]);
				}

				auto& planesIndices = anOutResult.myUniqueLines.back().myPlanesIndices;
				planesIndices.insert(planesIndices.end(), uniquePlanes.begin(), uniquePlanes.end());
			}
		}
		void locComputeVertices(LinesAndVerticesData& anOutResult)
		{
			constexpr auto pointIndexPairCmp = [](auto& a, auto& b) { return a.myKey < b.myKey; };
			struct PointIndexPair { Point3 myKey; size_t myIndex; };
			std::set<PointIndexPair, decltype(pointIndexPairCmp)> uniqueVertices{ pointIndexPairCmp };

			for (int i = 0; i < anOutResult.myUniqueLines.size(); ++i)
			{
				for (int k = i + 1; k < anOutResult.myUniqueLines.size(); ++k)
				{
					const auto intersection =
						CGAL::intersection(anOutResult.myUniqueLines[i].myLine, anOutResult.myUniqueLines[k].myLine);
					if (intersection)
					{
						const Point3* point = boost::get<Point3>(&*intersection);
						CGAL_precondition(point != nullptr);
						auto vertexIter = uniqueVertices.find(PointIndexPair{ *point });
						if (vertexIter == uniqueVertices.end())
						{
							vertexIter = uniqueVertices.insert(PointIndexPair{ *point, anOutResult.myUniqueVertices.size() }).first;
							anOutResult.myUniqueVertices.emplace_back(*point);
							anOutResult.myUniqueLines[i].myVerticesIndices.emplace_back(vertexIter->myIndex);
						}
						auto& verticesIndicesK = anOutResult.myUniqueLines[k].myVerticesIndices;
						if (std::find(verticesIndicesK.begin(), verticesIndicesK.end(), vertexIter->myIndex) == verticesIndicesK.end())
						{
							verticesIndicesK.emplace_back(vertexIter->myIndex);
						}
					}
				}
			}
		}
	}

	LowerEnvelope3d ParallelComputeLowerEnvelope(const std::vector<Plane>& somePlanes)
	{
		using namespace Parallel;
		constexpr auto policy = ExecutionPolicy::PARALLEL;
		constexpr auto distance = 1000.f;

		if (somePlanes.empty())
		{
			return LowerEnvelope3d{};
		}

		// Requirements check
		CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));

		LinesAndVerticesData linesVerticesData;
		locComputeLines(somePlanes, linesVerticesData);

		// Edge case: all planes are parallel
		if (linesVerticesData.myUniqueLines.empty())
		{
			const Point3 origin{ 0,0,0 };
			const auto lowestPlaneIndex = locGetLowestPlaneAtPoint<policy>(somePlanes, origin).second;
			CGAL_precondition(lowestPlaneIndex != -1);
			return lowestPlaneIndex;
		}

		locComputeVertices(linesVerticesData);

		// Edge case: no triple of planes intersect
		if (!linesVerticesData.myDoLinesIntersect)
		{
			static const Point3 zero{ 0,0,0 };
			constexpr auto innerPolicy = ExecutionPolicy::SEQUENTIAL;
			auto lineIter =
				hpx::find_if(hpx::execution::par, linesVerticesData.myUniqueLines.begin(), linesVerticesData.myUniqueLines.end(),
					[&](const auto& aLine) { return locIsVertexInLowerEnvelope<innerPolicy>(somePlanes, aLine.myLine.point()); });

			CGAL_precondition(lineIter != linesVerticesData.myUniqueLines.end());
			std::vector<Edge<Point3>> edges;

			Edge<Point3> firstEdge;
			firstEdge.myStart = lineIter->myLine.projection(zero);
			firstEdge.myEnd = firstEdge.myStart + lineIter->myLine.to_vector() * distance;
			firstEdge.myType = EdgeType::LINE;
			firstEdge.myLowestLeftPlane =
				locMinSteepPlaneIndexThroughSegment<policy>(firstEdge.myStart, firstEdge.myEnd, somePlanes, lineIter->myPlanesIndices);

			Edge<Point3> secondEdge;
			secondEdge.myStart = firstEdge.myEnd;
			secondEdge.myEnd = firstEdge.myStart;
			secondEdge.myType = EdgeType::LINE;
			secondEdge.myLowestLeftPlane =
				locMinSteepPlaneIndexThroughSegment<policy>(secondEdge.myStart, secondEdge.myEnd, somePlanes, lineIter->myPlanesIndices);

			edges.emplace_back(firstEdge);
			edges.emplace_back(secondEdge);

			return edges;
		}

		// Construct result
		constexpr auto innerPolicy = ExecutionPolicy::PARALLEL;
		std::vector<Edge<Point3>> edges;

		std::mutex edgesMutex;

		hpx::for_loop(hpx::execution::par, 0, linesVerticesData.myUniqueLines.size(), [&](const size_t aLineIdx)
			{
				auto& currentLine = linesVerticesData.myUniqueLines[aLineIdx];
				const auto& uniformLineVector = locGetUniformLineVector(currentLine);
				const auto& positive = currentLine.myVertices.back() + uniformLineVector * distance;
				const auto& negative = currentLine.myVertices.front() - uniformLineVector * distance;
				CGAL_precondition(currentLine.myLine.has_on(positive));
				CGAL_precondition(currentLine.myLine.has_on(negative));
				currentLine.myVertices.emplace_back(positive);
				currentLine.myVertices.insert(currentLine.myVertices.begin(), negative);

				// Build result
				for (int k = 0, segmentStartIdx = -1; k < currentLine.myVertices.size() - 1; ++k)
				{
					const auto& midpoint = CGAL::midpoint(currentLine.myVertices[k], currentLine.myVertices[k + 1]);
					CGAL_precondition(currentLine.myLine.has_on(midpoint));
					const auto isSegmentGood = locIsVertexInLowerEnvelope<innerPolicy>(somePlanes, midpoint);

					if (isSegmentGood && segmentStartIdx == -1)
					{
						segmentStartIdx = k;
					}

					const auto isStartIdxCase = segmentStartIdx != -1 && !isSegmentGood;
					const auto isEndIdxCase = segmentStartIdx != -1 && isSegmentGood && k == (currentLine.myVertices.size() - 2);

					if (isStartIdxCase || isEndIdxCase)
					{
						const int segmentEndIdx = isStartIdxCase ? k : (k + 1);
						const auto& startVertex = currentLine.myVertices[segmentStartIdx];
						const auto& endVertex = currentLine.myVertices[segmentEndIdx];

						Edge<Point3> edge, oppositeEdge;

						const auto caseIdx =
							static_cast<size_t>(segmentStartIdx == 0) |
							static_cast<size_t>(segmentEndIdx == currentLine.myVertices.size() - 1) << 1;

						switch (caseIdx)
						{
						case 0:
							edge.myType = EdgeType::SEGMENT;
							oppositeEdge.myType = EdgeType::SEGMENT;
							break;
						case 1:
							edge.myType = EdgeType::HALF_EDGE_EF;
							oppositeEdge.myType = EdgeType::HALF_EDGE_SF;
							break;
						case 2:
							edge.myType = EdgeType::HALF_EDGE_SF;
							oppositeEdge.myType = EdgeType::HALF_EDGE_EF;
							break;
						default:
							edge.myType = EdgeType::LINE;
							oppositeEdge.myType = EdgeType::LINE;
							break;
						}

						edge.myStart = startVertex;
						edge.myEnd = endVertex;
						edge.myLowestLeftPlane =
							locMinSteepPlaneIndexThroughSegment<innerPolicy>(edge.myStart, edge.myEnd, somePlanes, currentLine.myPlanesIndices);

						oppositeEdge.myStart = endVertex;
						oppositeEdge.myEnd = startVertex;
						oppositeEdge.myLowestLeftPlane =
							locMinSteepPlaneIndexThroughSegment<innerPolicy>(oppositeEdge.myStart, oppositeEdge.myEnd, somePlanes, currentLine.myPlanesIndices);

						{
							std::lock_guard guard(edgesMutex);
							edges.emplace_back(edge);
							edges.emplace_back(oppositeEdge);
						}

						segmentStartIdx = -1;
					}
				}
			});

		// Sort edges
		locSortLowerEnvelopeCCW<policy>(edges);

		return edges;
	}

	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes)
	{
		using namespace Sequential;
		constexpr auto policy = ExecutionPolicy::SEQUENTIAL;
		constexpr auto distance = 1000.f;

		if (somePlanes.empty())
		{
			return LowerEnvelope3d{};
		}

		// Requirements check
		CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));

		LinesAndVerticesData linesVerticesData;
		locComputeLines(somePlanes, linesVerticesData);

		// Edge case: all planes are parallel
		if (linesVerticesData.myUniqueLines.empty())
		{
			const Point3 origin{ 0,0,0 };
			const auto lowestPlaneIndex = locGetLowestPlaneAtPoint<policy>(somePlanes, origin).second;
			CGAL_precondition(lowestPlaneIndex != -1);
			return lowestPlaneIndex;
		}

		locComputeVertices(linesVerticesData);

		// Edge case: no triple of planes intersect
		if (linesVerticesData.myUniqueVertices.empty())
		{
			static const Point3 zero{ 0,0,0 };
			std::vector<Edge<Point3>> edges;

			for (size_t i = 0; i < linesVerticesData.myUniqueLines.size(); ++i)
			{
				const auto& currentLineData = linesVerticesData.myUniqueLines[i];
				const auto& currentLinePoint = currentLineData.myLine.point();
				if (locIsVertexInLowerEnvelope<policy>(somePlanes, currentLinePoint))
				{
					Edge<Point3> firstEdge;
					firstEdge.myStart = currentLineData.myLine.projection(zero);
					firstEdge.myEnd = firstEdge.myStart + currentLineData.myLine.to_vector() * distance;
					firstEdge.myType = EdgeType::LINE;
					firstEdge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment<policy>(firstEdge.myStart, firstEdge.myEnd, somePlanes, currentLineData.myPlanesIndices);

					Edge<Point3> secondEdge;
					secondEdge.myStart = firstEdge.myEnd;
					secondEdge.myEnd = firstEdge.myStart;
					secondEdge.myType = EdgeType::LINE;
					secondEdge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment<policy>(secondEdge.myStart, secondEdge.myEnd, somePlanes, currentLineData.myPlanesIndices);

					edges.emplace_back(firstEdge);
					edges.emplace_back(secondEdge);
				}
			}

			return edges;
		}

		// Construct result
		std::vector<Edge<Point3>> edges;

		for (size_t i = 0; i < linesVerticesData.myUniqueLines.size(); ++i)
		{
			auto& currentLine = linesVerticesData.myUniqueLines[i];
			const auto& uniformLineVector = locGetUniformLineVector(currentLine);

			hpx::sort(hpx::execution::seq, currentLine.myVerticesIndices.begin(), currentLine.myVerticesIndices.end(),
				[&](size_t aFirstIdx, size_t aSecondIdx) {
					return linesVerticesData.myUniqueVertices[aFirstIdx] < linesVerticesData.myUniqueVertices[aSecondIdx]; });

			const auto& positive =
				linesVerticesData.myUniqueVertices[currentLine.myVerticesIndices.back()] + uniformLineVector * distance;
			const auto& negative =
				linesVerticesData.myUniqueVertices[currentLine.myVerticesIndices.front()] - uniformLineVector * distance;

			CGAL_precondition(currentLine.myLine.has_on(positive));
			CGAL_precondition(currentLine.myLine.has_on(negative));

			linesVerticesData.myUniqueVertices.emplace_back(positive);
			linesVerticesData.myUniqueVertices.emplace_back(negative);

			currentLine.myVerticesIndices.emplace_back(linesVerticesData.myUniqueVertices.size() - 2);
			currentLine.myVerticesIndices.insert(currentLine.myVerticesIndices.begin(), linesVerticesData.myUniqueVertices.size() - 1);

			// Segments selection
			for (int k = 0, segmentStartIdx = -1; k < currentLine.myVerticesIndices.size() - 1; ++k)
			{
				const auto startIdx = currentLine.myVerticesIndices[k];
				const auto endIdx = currentLine.myVerticesIndices[k + 1];
				const auto& startVertex = linesVerticesData.myUniqueVertices[startIdx];
				const auto& endVertex = linesVerticesData.myUniqueVertices[endIdx];
				const auto& midpoint = CGAL::midpoint(startVertex, endVertex);
				CGAL_precondition(currentLine.myLine.has_on(midpoint));
				const auto isSegmentGood = locIsVertexInLowerEnvelope<policy>(somePlanes, midpoint);
				if (isSegmentGood && segmentStartIdx == -1)
				{
					segmentStartIdx = startIdx;
				}
				const auto isStartIdxCase = segmentStartIdx != -1 && !isSegmentGood;
				const auto isEndIdxCase = segmentStartIdx != -1 && isSegmentGood && k == (currentLine.myVerticesIndices.size() - 2);

				if (isStartIdxCase || isEndIdxCase)
				{
					const auto segmentEndIdx = isStartIdxCase ? startIdx : endIdx;
					const auto& startVertex = linesVerticesData.myUniqueVertices[segmentStartIdx];
					const auto& endVertex = linesVerticesData.myUniqueVertices[segmentEndIdx];

					Edge<Point3> edge, oppositeEdge;

					const auto caseIdx =
						static_cast<size_t>(segmentStartIdx == currentLine.myVerticesIndices.front()) |
						static_cast<size_t>(segmentEndIdx == currentLine.myVerticesIndices.back()) << 1;

					switch (caseIdx)
					{
					case 0:
						edge.myType = EdgeType::SEGMENT;
						oppositeEdge.myType = EdgeType::SEGMENT;
						break;
					case 1:
						edge.myType = EdgeType::HALF_EDGE_EF;
						oppositeEdge.myType = EdgeType::HALF_EDGE_SF;
						break;
					case 2:
						edge.myType = EdgeType::HALF_EDGE_SF;
						oppositeEdge.myType = EdgeType::HALF_EDGE_EF;
						break;
					default:
						edge.myType = EdgeType::LINE;
						oppositeEdge.myType = EdgeType::LINE;
						break;
					}

					edge.myStart = startVertex;
					edge.myEnd = endVertex;
					edge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment<policy>(edge.myStart, edge.myEnd, somePlanes, currentLine.myPlanesIndices);

					oppositeEdge.myStart = endVertex;
					oppositeEdge.myEnd = startVertex;
					oppositeEdge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment<policy>(oppositeEdge.myStart, oppositeEdge.myEnd, somePlanes, currentLine.myPlanesIndices);

					edges.emplace_back(edge);
					edges.emplace_back(oppositeEdge);

					segmentStartIdx = -1;
				}
			}
		}

		// Sort edges
		locSortLowerEnvelopeCCW<policy>(edges);

		return edges;
	}

	void TriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope)
	{
		using namespace Sequential;
		constexpr auto policy = ExecutionPolicy::SEQUENTIAL;
		CGAL_precondition(std::holds_alternative<std::vector<Edge<Point3>>>(anOutLowerEnvelope));
		auto& edges = std::get<std::vector<Edge<Point3>>>(anOutLowerEnvelope);
		CGAL_precondition(edges.size() > 2);

		hpx::sort(hpx::execution::seq, edges.begin(), edges.end(), [](const auto& aFirstEdge, const auto& aSecondEdge)
			{
				return aFirstEdge.myLowestLeftPlane < aSecondEdge.myLowestLeftPlane ||
					(aFirstEdge.myLowestLeftPlane == aSecondEdge.myLowestLeftPlane && aFirstEdge.myType < aSecondEdge.myType);
			});

		auto boundaryStartIt = edges.begin();

		std::vector<Edge<Point3>> temporaryEdges;

		for (auto it = edges.begin(); it != edges.end(); ++it)
		{
			if (boundaryStartIt->myLowestLeftPlane != it->myLowestLeftPlane)
			{
				// Triangles building function
				const auto& faceDiagonals = locTriangulateConvexFace(boundaryStartIt, it);
				temporaryEdges.insert(temporaryEdges.end(), faceDiagonals.begin(), faceDiagonals.end());
				boundaryStartIt = it;
			}
		}

		{
			const auto& faceDiagonals = locTriangulateConvexFace(boundaryStartIt, edges.end());
			temporaryEdges.insert(temporaryEdges.end(), faceDiagonals.begin(), faceDiagonals.end());
		}

		edges.insert(edges.end(), temporaryEdges.begin(), temporaryEdges.end());
		locSortLowerEnvelopeCCW<policy>(edges);
	}

	void ParallelTriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope)
	{
		using namespace Common;
		constexpr auto policy = ExecutionPolicy::PARALLEL;
		CGAL_precondition(std::holds_alternative<std::vector<Edge<Point3>>>(anOutLowerEnvelope));
		auto& edges = std::get<std::vector<Edge<Point3>>>(anOutLowerEnvelope);
		CGAL_precondition(edges.size() > 2);

		hpx::sort(hpx::execution::par, edges.begin(), edges.end(), [](const auto& aFirstEdge, const auto& aSecondEdge)
			{
				return aFirstEdge.myLowestLeftPlane < aSecondEdge.myLowestLeftPlane ||
					(aFirstEdge.myLowestLeftPlane == aSecondEdge.myLowestLeftPlane && aFirstEdge.myType < aSecondEdge.myType);
			});

		// Dummy upper bound
		std::vector<size_t> ranges(edges.size() / 2);
		std::atomic<size_t> atomicCounter{ 0 };

		hpx::for_loop(hpx::execution::par_unseq, 0, edges.size(), [&](size_t anEdgeIdx) {
			if (anEdgeIdx == 0 ||
				edges[anEdgeIdx].myLowestLeftPlane != edges[anEdgeIdx - 1].myLowestLeftPlane)
			{
				ranges[atomicCounter.fetch_add(1)] = anEdgeIdx;
			}
			});

		ranges.resize(atomicCounter.fetch_xor(atomicCounter.load()) + 1);
		ranges.back() = edges.size();

		const auto prevEdgesCount = edges.size();
		atomicCounter.store(prevEdgesCount);

		// Dummy upper-bound
		edges.resize(prevEdgesCount * 3);

		hpx::for_loop(hpx::execution::par_unseq, 0, ranges.size() - 1, [&](size_t aRangeIdx) {
			hpx::for_loop(hpx::execution::par_unseq, ranges[aRangeIdx] + 1, ranges[aRangeIdx + 1], [&](size_t anEdgeIdx) {
				const auto& startEdge = edges[ranges[aRangeIdx]];
				Edge<Point3> edge, oppositeEdge;

				edge.myType = EdgeType::SEGMENT_TRIANGLE;
				edge.myStart = startEdge.myStart;
				edge.myEnd = edges[anEdgeIdx].myEnd;

				oppositeEdge.myType = EdgeType::SEGMENT_TRIANGLE;
				oppositeEdge.myStart = edges[anEdgeIdx].myEnd;
				oppositeEdge.myEnd = startEdge.myStart;

				const auto storeIdx = atomicCounter.fetch_add(2);
				edges[storeIdx] = edge;
				edges[storeIdx + 1] = oppositeEdge;
				});
			});

		edges.resize(atomicCounter.fetch_xor(atomicCounter.load()) - prevEdgesCount);
		locSortLowerEnvelopeCCW<policy>(edges);
	}
}