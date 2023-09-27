#include "BruteForce.h"
#include "Utils.h"

namespace SPGMT
{
	namespace
	{
		template<typename T>
		bool locAreItemsUnique(const std::vector<T>& someItems)
		{
			auto areItemsUnique{ true };
			for (int i = 0; i < someItems.size() && areItemsUnique; ++i)
			{
				for (int k = i + 1; k < someItems.size(); ++k)
				{
					if (someItems[i] == someItems[k])
					{
						areItemsUnique = false;
						break;
					}
				}
			}
			return areItemsUnique;
		}
		bool locArePlanesNonVertical(const std::vector<Plane>& somePlanes)
		{
			const Vec3 up{ 0,0,1 };
			const FT zero{ 0 };
			auto arePlanesNonVertical{ true };
			for (int i = 0; i < somePlanes.size() && arePlanesNonVertical; ++i)
			{
				arePlanesNonVertical =
					CGAL::scalar_product(somePlanes[i].orthogonal_vector(), up) != zero;
			}
			return arePlanesNonVertical;
		}
		bool locArePlanesUniformlyOrientedUp(const std::vector<Plane>& somePlanes)
		{
			const Vec3 up{ 0, 0, 1 };
			const FT zero{ 0 };
			auto isPointingUp{ false };
			for (auto& plane : somePlanes)
			{
				const auto dot = CGAL::scalar_product(up, plane.orthogonal_vector());
				// Cannot handle vertical planes
				CGAL_precondition(dot != zero);
				isPointingUp |= dot > zero;
			}
			return isPointingUp;
		}
		int locGetLowestPlaneAtOrigin(const std::vector<Plane>& somePlanes)
		{
			CGAL_precondition(somePlanes.size() > 0);
			const Line3 upLine{ Point3{0,0,0}, Vec3{0,0,1} };
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
			const Line3 upLine { aPoint, Vec3{0,0,1} };
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

		int locFindLineIndex(const std::vector<LineData>& someLines, const Line3& aLine)
		{
			const auto lineIter = std::find_if(someLines.begin(), someLines.end(), [&aLine](const auto& aLineData) {
				return aLine == aLineData.myLine;
				});

			if (lineIter == someLines.end())
			{
				return -1;
			}
			else
			{
				return std::distance(someLines.begin(), lineIter);
			}
		}
		LinesAndVerticesData locComputeLinesAndVertices(const std::vector<Plane>& somePlanes)
		{
			LinesAndVerticesData result;
			// Intersection lines computation
			for (int i = 0; i < somePlanes.size(); ++i)
			{
				for (int k = i + 1; k < somePlanes.size(); ++k)
				{
					const auto intersection = CGAL::intersection(somePlanes[i], somePlanes[k]);
					if (intersection)
					{
						const Line3* line = boost::get<Line3>(&*intersection);
						CGAL_precondition(line != nullptr);
						auto lineIdx = locFindLineIndex(result.myUniqueLines, *line);
						if (lineIdx == -1)
						{
							LineData data{ *line };
							data.myPlanesIndices.push_back(i);
							result.myUniqueLines.push_back(data);
							lineIdx = result.myUniqueLines.size() - 1;
						}
						result.myUniqueLines[lineIdx].myPlanesIndices.push_back(k);
					}
				}
			}
			// End intersection lines computation
			// Vertices computation
			constexpr auto pointIndexPairCmp = [](auto& a, auto& b) { return a.myKey < b.myKey; };
			struct PointIndexPair { Point3 myKey; size_t myIndex; };
			std::set<PointIndexPair, decltype(pointIndexPairCmp)> uniqueVertices{ pointIndexPairCmp };
			for (int i = 0; i < result.myUniqueLines.size(); ++i)
			{
				for (int k = i + 1; k < result.myUniqueLines.size(); ++k)
				{
					const auto intersection = CGAL::intersection(result.myUniqueLines[i].myLine, result.myUniqueLines[k].myLine);
					if (intersection)
					{
						const Point3* point = boost::get<Point3>(&*intersection);
						CGAL_precondition(point != nullptr);
						auto vertexIter = uniqueVertices.find(PointIndexPair{ *point });
						if (vertexIter == uniqueVertices.end())
						{
							vertexIter = uniqueVertices.insert(PointIndexPair{ *point, result.myUniqueVertices.size() }).first;
							result.myUniqueVertices.push_back(VertexData{ *point });
							result.myUniqueLines[i].mySortedVerticesIndices.push_back(vertexIter->myIndex);
						}
						auto& verticesIndicesK = result.myUniqueLines[k].mySortedVerticesIndices;
						if (std::find(verticesIndicesK.begin(), verticesIndicesK.end(), vertexIter->myIndex) == verticesIndicesK.end())
						{
							verticesIndicesK.push_back(vertexIter->myIndex);
						}
					}
				}
			}
			//End vertices computation
			return result;
		}

		template<typename T>
		T locNormalize(const T& aVec)
		{
			auto const slen = aVec.squared_length();
			auto const d = CGAL::sqrt(slen);
			return aVec / d;
		}

		void locSortClockWise(Vertex& anOutVertex, const std::vector<Vertex>& someVertices)
		{
			struct InternalWrapper
			{
				FT myAngle{ 0 };
				int myVertexIdx{ -1 }, myPlaneIdx{ -1 };
			};

			std::vector<InternalWrapper> wrappers;

			for (int i = 0; i < anOutVertex.mySortedNeighboursIndices.size(); ++i)
			{
				const auto lowestPlaneIdx = anOutVertex.myLowestLeftPlanes[i];
				const auto targetVertexIdx = anOutVertex.mySortedNeighboursIndices[i];

				const Vec3 segment{ anOutVertex.myPoint, someVertices[targetVertexIdx].myPoint };
				const auto angleDegrees = CGAL::approximate_angle(Vec3{ 0,0,1 }, segment);

				wrappers.push_back(InternalWrapper{ angleDegrees, targetVertexIdx, lowestPlaneIdx });
			}

			anOutVertex.mySortedNeighboursIndices.clear();
			anOutVertex.myLowestLeftPlanes.clear();

			std::sort(wrappers.begin(), wrappers.end(), [](auto& aFirst, auto& aSecond) {
				return aFirst.myAngle > aSecond.myAngle;
				});

			std::for_each(wrappers.begin(), wrappers.end(), [&anOutVertex](auto& aWrapper) {
				anOutVertex.mySortedNeighboursIndices.push_back(aWrapper.myVertexIdx);
				anOutVertex.myLowestLeftPlanes.push_back(aWrapper.myPlaneIdx);
				});
		}

		// Visual Help. https://www.falstad.com/dotproduct/
		std::pair<int, int> locMinMaxSteepPlaneIndexThroughSegment(const Point3& aStart, const Point3& anEnd, const std::vector<Plane>& somePlanes, const std::vector<int>& somePlanesIndices)
		{
			const Vec3 upDir{ 0,0,1 };

			CGAL_precondition(somePlanes.size() > 0);
			CGAL_precondition(somePlanesIndices.size() > 0);

			Vec3 segmentDir{ aStart, anEnd };
			segmentDir = locNormalize(segmentDir);

			const Vec3 leftDir{ CGAL::cross_product(segmentDir, upDir) };

			FT maxSteep = CGAL::scalar_product(somePlanes[somePlanesIndices[0]].orthogonal_vector(), leftDir);
			FT minSteep = maxSteep;
			auto maxPlaneIdx = somePlanesIndices[0];
			auto minPlaneIdx = maxPlaneIdx;

			for (int i = 1; i < somePlanesIndices.size(); ++i)
			{
				const auto currentPlaneIdx = somePlanesIndices[i];
				const auto dot = CGAL::scalar_product(somePlanes[currentPlaneIdx].orthogonal_vector(), leftDir);

				if (dot > maxSteep)
				{
					maxSteep = dot;
					maxPlaneIdx = currentPlaneIdx;
				}
				else if (dot < minSteep)
				{
					minSteep = dot;
					minPlaneIdx = currentPlaneIdx;
				}
			}

			return std::make_pair(minPlaneIdx, maxPlaneIdx);
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
	}

	std::vector<Vertex> ComputeLowerEnvelope(const std::vector<Plane>& somePlanes)
	{
		std::vector<Vertex> result;

		if (somePlanes.empty())
		{
			return result;
		}

		// Requirements check
		CGAL_precondition(locAreItemsUnique(somePlanes));
		CGAL_precondition(locArePlanesNonVertical(somePlanes));
		CGAL_precondition(locArePlanesUniformlyOrientedUp(somePlanes));

		auto& linesVerticesData = locComputeLinesAndVertices(somePlanes);

		// Edge case: all planes are parallel
		if (linesVerticesData.myUniqueLines.empty())
		{
			const auto lowestPlaneIndex = locGetLowestPlaneAtOrigin(somePlanes);
			CGAL_precondition(lowestPlaneIndex != -1);
			Vertex infinity{ VertexType::INFINITE };
			infinity.myLowestLeftPlanes.push_back(lowestPlaneIndex);
			result.push_back(infinity);
			return result;
		}

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
					Vertex start, end;

					start.myPoint = currentLineData.myLine.point();
					start.mySortedNeighboursIndices.push_back(1);
					start.myType = VertexType::INFINITE;

					end.myPoint = start.myPoint + currentLineData.myLine.to_vector() * distance;
					end.mySortedNeighboursIndices.push_back(0);
					end.myType = VertexType::INFINITE;

					const auto minMaxPlaneIndices =
						locMinMaxSteepPlaneIndexThroughSegment(start.myPoint, end.myPoint, somePlanes, currentLineData.myPlanesIndices);

					start.myLowestLeftPlanes.push_back(minMaxPlaneIndices.second);
					end.myLowestLeftPlanes.push_back(minMaxPlaneIndices.first);

					result.push_back(start);
					result.push_back(end);
					break;
				}
			}
			return result;
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
		constexpr auto pointIndexPairCmp = [](auto& a, auto& b) { return a.myKey < b.myKey; };
		struct PointIndexPair { Point3 myKey; size_t myIndex; };
		std::set<PointIndexPair, decltype(pointIndexPairCmp)> uniqueVertices{ pointIndexPairCmp };

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

					auto startIter = uniqueVertices.find(PointIndexPair{ startVertex.myPoint });
					auto endIter = uniqueVertices.find(PointIndexPair{ endVertex.myPoint });

					if (startIter == uniqueVertices.end())
					{
						startIter = uniqueVertices.insert(PointIndexPair{ startVertex.myPoint, uniqueVertices.size() }).first;
						Vertex start;
						start.myPoint = startIter->myKey;
						start.myType = startVertex.myIsAtInfinity ? VertexType::INFINITE : VertexType::FINITE;
						result.push_back(start);
					}

					if (endIter == uniqueVertices.end())
					{
						endIter = uniqueVertices.insert(PointIndexPair{ endVertex.myPoint, uniqueVertices.size() }).first;
						Vertex end;
						end.myPoint = endIter->myKey;
						end.myType = endVertex.myIsAtInfinity ? VertexType::INFINITE : VertexType::FINITE;
						result.push_back(end);
					}

					auto& start = result[startIter->myIndex];
					auto& end = result[endIter->myIndex];

					const auto minMaxPlaneIndices = locMinMaxSteepPlaneIndexThroughSegment(start.myPoint, end.myPoint,
						somePlanes, currentLine.myPlanesIndices);

					start.mySortedNeighboursIndices.push_back(endIter->myIndex);
					start.myLowestLeftPlanes.push_back(minMaxPlaneIndices.second);

					end.mySortedNeighboursIndices.push_back(startIter->myIndex);
					end.myLowestLeftPlanes.push_back(minMaxPlaneIndices.first);

					segmentStartIdx = -1;
				}
			}
		}

		// Sort each edge counter-clock wise
		for (int i = 0; i < result.size(); ++i)
		{
			locSortClockWise(result[i], result);

			// Set vertex to infinity if it has only one edge
			if (result[i].mySortedNeighboursIndices.size() == 1)
			{
				result[i].myType = VertexType::INFINITE;
			}
		}

		return result;
	}
}