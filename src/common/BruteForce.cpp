#include "BruteForce.h"
#include "Utils.h"
#include <CGAL/Polygon_2_algorithms.h>

namespace SPGMT
{
	namespace
	{
		int locGetLowestPlaneAtOrigin(const std::vector<Plane>& somePlanes, const Point3& aPoint = Point3{ 0,0,0 })
		{
			CGAL_precondition(somePlanes.size() > 0);
			const Line3 upLine{ aPoint, Vec3{0,0,1} };
			// Initialize minimum z with the first plane
			FT minZ;
			int lowestPlaneIdx = 0;
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes.front());
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				minZ = point->z();
			}
			// Check all the remaining planes' z values
			for (int i = 1; i < somePlanes.size(); ++i)
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes[i]);
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				if (point->z() < minZ)
				{
					minZ = point->z();
					lowestPlaneIdx = i;
				}
			}
			return lowestPlaneIdx;
		}
		bool locIsVertexInLowerEnvelope(const std::vector<Plane>& somePlanes, const Point3& aPoint)
		{
			const Line3 upLine{ aPoint, Vec3{0,0,1} };
			FT minZ = aPoint.z();
			for (int i = 0; i < somePlanes.size(); ++i)
			{
				const auto inter = CGAL::intersection(somePlanes[i], upLine);
				CGAL_precondition(inter.has_value());
				const Point3* point = boost::get<Point3>(&*inter);
				CGAL_precondition(point != nullptr);
				if (point->z() < minZ)
				{
					return false;
				}
			}
			return true;
		}

		struct LineData
		{
			Line3 myLine;
			std::vector<int> myPlanesIndices;
			std::vector<int> mySortedVerticesIndices;
		};

		struct VertexData
		{
			Point3 myPoint;
			bool myIsAtInfinity{ false };
		};

		struct LinesAndVerticesData
		{
			std::vector<LineData> myUniqueLines;
			std::vector<VertexData> myUniqueVertices;
		};

		Point3 locProject(const Point3& aPoint)
		{
			return Point3{ aPoint.x(), aPoint.y(), 0 };
		}

		LinesAndVerticesData locComputeLinesAndVertices(const std::vector<Plane>& somePlanes)
		{
			static const Point3 zero{ 0,0,0 };

			struct LineWrapper
			{
				Point3 myPoint;
				Dir3 myDir;
				int myFirstPlane{ -1 }, mySecondPlane{ -1 };
			};

			std::cout << "START line computation" << std::endl;

			const auto maximumLinesCount = somePlanes.size() * (somePlanes.size() - 1) / 2;
			std::vector<LineWrapper> lines;
			lines.reserve(maximumLinesCount);

			for (int i = 0; i < somePlanes.size(); ++i)
			{
				for (int k = i + 1; k < somePlanes.size(); ++k)
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

			lines.shrink_to_fit();
			const auto linesSort = [&](auto& aFirstWrapper, auto& aSecondWrapper)
			{
				return aFirstWrapper.myPoint < aSecondWrapper.myPoint ||
					(aFirstWrapper.myPoint == aSecondWrapper.myPoint && aFirstWrapper.myDir == aSecondWrapper.myDir);
			};
			std::sort(lines.begin(), lines.end(), linesSort);

			// Store unique lines
			LinesAndVerticesData result;

			if (lines.size() > 0)
			{
				Line3 uniqueLine{ lines.front().myPoint, lines.front().myDir };
				result.myUniqueLines.emplace_back(LineData{ std::move(uniqueLine) });
				result.myUniqueLines.back().myPlanesIndices.push_back(lines.front().myFirstPlane);
				result.myUniqueLines.back().myPlanesIndices.push_back(lines.front().mySecondPlane);
			}

			for (const auto& line : lines)
			{
				Line3 maybeUniqueLine{ line.myPoint, line.myDir };
				if (result.myUniqueLines.back().myLine != maybeUniqueLine)
				{
					result.myUniqueLines.emplace_back(LineData{ std::move(maybeUniqueLine) });
				}
				result.myUniqueLines.back().myPlanesIndices.push_back(line.myFirstPlane);
				result.myUniqueLines.back().myPlanesIndices.push_back(line.mySecondPlane);
			}

			// unique lines stored

			std::cout << "DONE line computation" << std::endl;
			// -----------------------------------------------------


			std::cout << "START vertices computation" << std::endl;
			// store unique vertices
			struct VertexWrapper
			{
				Point3 myPoint;
				int myFirstLine{ -1 }, mySecondLine{ -1 };
			};

			const auto maximumVerticesCount = result.myUniqueLines.size() * (result.myUniqueLines.size() - 1) / 2;
			std::vector<VertexWrapper> vertices;
			vertices.reserve(maximumVerticesCount);

			for (int i = 0; i < result.myUniqueLines.size(); ++i)
			{
				for (int k = i + 1; k < result.myUniqueLines.size(); ++k)
				{
					const auto intersection = CGAL::intersection(result.myUniqueLines[i].myLine, result.myUniqueLines[k].myLine);
					if (intersection)
					{
						const Point3* point = boost::get<Point3>(&*intersection);
						CGAL_precondition(point != nullptr);
						vertices.emplace_back(VertexWrapper{ *point, i, k });
					}
				}
			}

			vertices.shrink_to_fit();
			const auto verticesSort = [&](auto& aFirstWrapper, auto& aSecondWrapper)
			{
				return aFirstWrapper.myPoint < aSecondWrapper.myPoint;
			};
			std::sort(vertices.begin(), vertices.end(), verticesSort);

			if (vertices.size() > 0)
			{
				result.myUniqueVertices.emplace_back(VertexData{ vertices.front().myPoint });
				constexpr auto vertexIdx = 0;
				result.myUniqueLines[vertices.front().myFirstLine].mySortedVerticesIndices.push_back(vertexIdx);
				result.myUniqueLines[vertices.front().mySecondLine].mySortedVerticesIndices.push_back(vertexIdx);
			}

			for (auto& vertex : vertices)
			{
				if (result.myUniqueVertices.back().myPoint != vertex.myPoint)
				{
					result.myUniqueVertices.emplace_back(VertexData{ vertex.myPoint });
				}

				const auto vertexIdx = result.myUniqueVertices.size() - 1;
				result.myUniqueLines[vertex.myFirstLine].mySortedVerticesIndices.push_back(vertexIdx);
				result.myUniqueLines[vertex.mySecondLine].mySortedVerticesIndices.push_back(vertexIdx);
			}
			std::cout << "DONE vertices computation" << std::endl;
			//End vertices computation
			// 
			// Keep only unique
			std::cout << "START keep unique planes and vertices indices" << std::endl;

			for (auto& line : result.myUniqueLines)
			{
				std::sort(line.myPlanesIndices.begin(), line.myPlanesIndices.end());
				std::sort(line.mySortedVerticesIndices.begin(), line.mySortedVerticesIndices.end());
				line.myPlanesIndices.erase(std::unique(line.myPlanesIndices.begin(), line.myPlanesIndices.end()), line.myPlanesIndices.end());
				line.mySortedVerticesIndices.erase(std::unique(line.mySortedVerticesIndices.begin(), line.mySortedVerticesIndices.end()), line.mySortedVerticesIndices.end());
			}

			std::cout << "DONE keep unique planes and vertices indices" << std::endl;
			return result;
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

		// Visual Help. https://www.falstad.com/dotproduct/ | https://www.geogebra.org/3d?lang=en
		int locMinSteepPlaneIndexThroughSegment(const Point3& aStart, const Point3& anEnd, const std::vector<Plane>& somePlanes, const std::vector<int>& somePlanesIndices)
		{
			const Vec3 upDir{ 0,0,1 };

			CGAL_precondition(somePlanes.size() > 0);
			CGAL_precondition(somePlanesIndices.size() > 0);

			const Vec3 segmentDir{ aStart, anEnd };
			const Vec3 leftDir{ -CGAL::cross_product(segmentDir, upDir) };
			const auto sampleLoc = CGAL::midpoint(aStart, anEnd) + leftDir * 1000.f;

			CGAL_precondition((somePlanes.size() > 0 && somePlanesIndices.size() > 0));
			const Line3 upLine{ sampleLoc, upDir };
			// Initialize minimum z with the first plane
			FT minZ;
			int lowestPlaneIdx = somePlanesIndices[0];
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes[somePlanesIndices.front()]);
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				minZ = point->z();
			}
			// Check all the remaining planes' z values
			for (int i = 1; i < somePlanesIndices.size(); ++i)
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes[somePlanesIndices[i]]);
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				if (point->z() < minZ)
				{
					minZ = point->z();
					lowestPlaneIdx = somePlanesIndices[i];
				}
			}
			return lowestPlaneIdx;
		}
		Vec3 locGetUniformLineVector(const LineData& aLineData)
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
		
		std::vector<Point3> locTriangulateConvexFace(
			const std::vector<Edge<Point3>>::iterator aStartIter,
			const std::vector<Edge<Point3>>::iterator anEndIter)
		{
			std::vector<Point3> vertices;

			CGAL_precondition(std::distance(aStartIter, anEndIter) > 1);

			// The triangulation is guaranteed to be CCW because the first and second vertices always bound the face to their left
			// So every triangulation will pick the vertices in CCW order by design
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

				vertices.emplace_back(aStartIter->myStart);
				vertices.emplace_back(it->myStart);
				vertices.emplace_back(it->myEnd);
#ifndef NDEBUG
				{
					std::vector<Point2> vertices2d;
					vertices2d.emplace_back(Point2{ aStartIter->myStart.x(), aStartIter->myStart.y() });
					vertices2d.emplace_back(Point2{ it->myStart.x(), it->myStart.y() });
					vertices2d.emplace_back(Point2{ it->myEnd.x(), it->myEnd.y() });
					CGAL_postcondition(CGAL::orientation_2(vertices2d.begin(), vertices2d.end()) == CGAL::Orientation::COUNTERCLOCKWISE);
				}
#endif
			}
#ifndef NDEBUG
			const auto isUnboundedFace = aStartIter->myType != EdgeType::SEGMENT;
			const auto edgesCount = std::distance(aStartIter, anEndIter) + (isUnboundedFace ? 1 : 0);
			CGAL_postcondition((vertices.size() / 3) == edgesCount - 2);
#endif
			return vertices;
		}
	}

	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes)
	{
		if (somePlanes.empty())
		{
			return LowerEnvelope3d{};
		}

		// Requirements check
		CGAL_precondition(Utils::AreItemsUnique(somePlanes));
		CGAL_precondition(Utils::ArePlanesNonVertical(somePlanes));
		//CGAL_precondition(locArePlanesUniformlyOrientedUp(somePlanes));

		auto& linesVerticesData = locComputeLinesAndVertices(somePlanes);

		// Edge case: all planes are parallel
		if (linesVerticesData.myUniqueLines.empty())
		{
			const auto lowestPlaneIndex = locGetLowestPlaneAtOrigin(somePlanes);
			CGAL_precondition(lowestPlaneIndex != -1);
			return lowestPlaneIndex;
		}

		std::vector<Edge<Point3>> edges;

		// Edge case: no triple of planes intersect
		if (linesVerticesData.myUniqueVertices.empty())
		{
			for (int i = 0; i < linesVerticesData.myUniqueLines.size(); ++i)
			{
				const auto& currentLineData = linesVerticesData.myUniqueLines[i];
				const auto& currentLinePoint = currentLineData.myLine.point();
				if (locIsVertexInLowerEnvelope(somePlanes, currentLinePoint))
				{
					constexpr auto distance = 1000.f;
					static const Point3 zero{ 0,0,0 };

					Edge<Point3> firstEdge;
					firstEdge.myStart = currentLineData.myLine.projection(zero);
					firstEdge.myEnd = firstEdge.myStart + currentLineData.myLine.to_vector() * distance;
					firstEdge.myType = EdgeType::LINE;
					firstEdge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment(firstEdge.myStart, firstEdge.myEnd, somePlanes, currentLineData.myPlanesIndices);

					Edge<Point3> secondEdge;
					secondEdge.myStart = firstEdge.myEnd;
					secondEdge.myEnd = firstEdge.myStart;
					secondEdge.myType = EdgeType::LINE;
					secondEdge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment(secondEdge.myStart, secondEdge.myEnd, somePlanes, currentLineData.myPlanesIndices);

					edges.emplace_back(firstEdge);
					edges.emplace_back(secondEdge);

					break;
				}
			}
			return edges;
		}

		// Base case
		for (int i = 0; i < linesVerticesData.myUniqueLines.size(); ++i)
		{
			// First sort lexicographically vertices along the line
			auto& currentLine = linesVerticesData.myUniqueLines[i];
			auto& vertexIndices = currentLine.mySortedVerticesIndices;
			CGAL_precondition(!vertexIndices.empty());
			std::sort(vertexIndices.begin(), vertexIndices.end(), [&linesVerticesData](auto& a, auto& b) {
				return linesVerticesData.myUniqueVertices[a].myPoint < linesVerticesData.myUniqueVertices[b].myPoint;
				});

			// Add two custom vertices at each edge of the line representing +/- infinity
			constexpr auto offset = 1000.f;
			const auto& uniformLineVector = locGetUniformLineVector(currentLine);
			const auto& positive = linesVerticesData.myUniqueVertices[vertexIndices.back()].myPoint + uniformLineVector * offset;
			const auto& negative = linesVerticesData.myUniqueVertices[vertexIndices.front()].myPoint - uniformLineVector * offset;

			CGAL_precondition(currentLine.myLine.has_on(positive));
			CGAL_precondition(currentLine.myLine.has_on(negative));

			constexpr auto isAtInfinity = true;
			linesVerticesData.myUniqueVertices.push_back(VertexData{ positive, isAtInfinity });
			vertexIndices.push_back(linesVerticesData.myUniqueVertices.size() - 1);

			linesVerticesData.myUniqueVertices.push_back(VertexData{ negative, isAtInfinity });
			vertexIndices.insert(vertexIndices.begin(), linesVerticesData.myUniqueVertices.size() - 1);
		}

		// Construct result
		for (int i = 0; i < linesVerticesData.myUniqueLines.size(); ++i)
		{
			auto& currentLine = linesVerticesData.myUniqueLines[i];
			int segmentStartIdx = -1;
			for (int k = 0; k < currentLine.mySortedVerticesIndices.size() - 1; ++k)
			{
				const auto startIdx = currentLine.mySortedVerticesIndices[k];
				const auto endIdx = currentLine.mySortedVerticesIndices[k + 1];
				const auto& startVertex = linesVerticesData.myUniqueVertices[startIdx];
				const auto& endVertex = linesVerticesData.myUniqueVertices[endIdx];
				const auto& midpoint = CGAL::midpoint(startVertex.myPoint, endVertex.myPoint);
				CGAL_precondition(currentLine.myLine.has_on(midpoint));
				const auto isSegmentGood = locIsVertexInLowerEnvelope(somePlanes, midpoint);

				if (isSegmentGood && segmentStartIdx == -1)
				{
					segmentStartIdx = startIdx;
				}

				const auto isStartIdxCase = segmentStartIdx != -1 && !isSegmentGood;
				const auto isEndIdxCase = segmentStartIdx != -1 && isSegmentGood && k == (currentLine.mySortedVerticesIndices.size() - 2);

				if (isStartIdxCase || isEndIdxCase)
				{
					const int segmentEndIdx = isStartIdxCase ? startIdx : endIdx;
					const auto& startVertex = linesVerticesData.myUniqueVertices[segmentStartIdx];
					const auto& endVertex = linesVerticesData.myUniqueVertices[segmentEndIdx];

					Edge<Point3> edge, oppositeEdge;

					if (startVertex.myIsAtInfinity && endVertex.myIsAtInfinity)
					{
						edge.myType = EdgeType::LINE;
						oppositeEdge.myType = EdgeType::LINE;
					}
					else if (!startVertex.myIsAtInfinity && endVertex.myIsAtInfinity)
					{
						edge.myType = EdgeType::HALF_EDGE_SF;
						oppositeEdge.myType = EdgeType::HALF_EDGE_EF;
					}
					else if (startVertex.myIsAtInfinity && !endVertex.myIsAtInfinity)
					{
						edge.myType = EdgeType::HALF_EDGE_EF;
						oppositeEdge.myType = EdgeType::HALF_EDGE_SF;
					}
					else
					{
						edge.myType = EdgeType::SEGMENT;
						oppositeEdge.myType = EdgeType::SEGMENT;
					}

					edge.myStart = startVertex.myPoint;
					edge.myEnd = endVertex.myPoint;
					edge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment(edge.myStart, edge.myEnd, somePlanes, currentLine.myPlanesIndices);

					oppositeEdge.myStart = endVertex.myPoint;
					oppositeEdge.myEnd = startVertex.myPoint;
					oppositeEdge.myLowestLeftPlane =
						locMinSteepPlaneIndexThroughSegment(oppositeEdge.myStart, oppositeEdge.myEnd, somePlanes, currentLine.myPlanesIndices);

					edges.emplace_back(edge);
					edges.emplace_back(oppositeEdge);

					segmentStartIdx = -1;
				}
			}
		}

		const auto sortEdges = [](auto& aFirstEdge, auto& aSecondEdge)
		{
			return aFirstEdge.myStart < aSecondEdge.myStart;
		};

		std::sort(edges.begin(), edges.end(), sortEdges);

		auto iterStart = edges.begin();

		for (auto iter = edges.begin(); iter != edges.end(); ++iter)
		{
			if (iterStart->myStart != iter->myStart)
			{
				locSortNeighboursCCW(iterStart, iter);
				iterStart = iter;
			}
		}

		locSortNeighboursCCW(iterStart, edges.end());

		return edges;
	}

	Triangles3d TriangulateLowerEnvelope(const LowerEnvelope3d& aLowerEnvelope)
	{
		CGAL_precondition(std::holds_alternative<std::vector<Edge<Point3>>>(aLowerEnvelope));

		if (!std::holds_alternative<std::vector<Edge<Point3>>>(aLowerEnvelope))
		{
			return Triangles3d{};
		}

		auto edges = std::get<std::vector<Edge<Point3>>>(aLowerEnvelope);
		CGAL_precondition(edges.size() > 2);

		if (edges.size() < 3)
		{
			return Triangles3d{};
		}

		const auto sortEdges = [](auto& aFirstEdge, auto& aSecondEdge)
		{
			return aFirstEdge.myLowestLeftPlane < aSecondEdge.myLowestLeftPlane ||
				(aFirstEdge.myLowestLeftPlane == aSecondEdge.myLowestLeftPlane && aFirstEdge.myType < aSecondEdge.myType);
		};

		std::sort(edges.begin(), edges.end(), sortEdges);
		std::vector<Point3> result;
		auto boundaryStartIt = edges.begin();
		
		for (auto it = edges.begin(); it != edges.end(); ++it)
		{
			if (boundaryStartIt->myLowestLeftPlane != it->myLowestLeftPlane)
			{
				// Triangles building function
				const auto& faceTriangles = locTriangulateConvexFace(boundaryStartIt, it);
				result.insert(result.end(), faceTriangles.begin(), faceTriangles.end());
				boundaryStartIt = it;
			}
		}

		{
			const auto& faceTriangles = locTriangulateConvexFace(boundaryStartIt, edges.end());
			result.insert(result.end(), faceTriangles.begin(), faceTriangles.end());
		}

		return result;
	}
}