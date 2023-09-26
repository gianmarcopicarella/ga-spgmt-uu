#include "BruteForce.h"
#include "Utils.h"

#include <CGAL/bounding_box.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/centroid.h>

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
			const Vec3 verticalAxis{ 0,1,0 };
			const FT zero{ 0 };
			auto arePlanesNonVertical{ true };
			for (int i = 0; i < somePlanes.size() && arePlanesNonVertical; ++i)
			{
				arePlanesNonVertical =
					CGAL::scalar_product(somePlanes[i].orthogonal_vector(), verticalAxis) != zero;
			}
			return arePlanesNonVertical;
		}
		bool locArePlanesUniformlyOrientedUp(const std::vector<Plane>& somePlanes)
		{
			const Vec3 up{ 0, 1, 0 };
			const FT zero{ 0 };
			auto isPointingDown{ false }, isPointingUp{ false };
			for (auto& plane : somePlanes)
			{
				const auto dot = CGAL::scalar_product(up, plane.orthogonal_vector());
				// Cannot handle vertical planes
				CGAL_precondition(dot != zero);
				isPointingUp |= dot > zero;
				isPointingDown |= dot < zero;

				if (isPointingUp && isPointingDown)
				{
					return false;
				}
			}

			if (!isPointingUp)
			{
				return false;
			}

			return true;
		}

		int locGetLowestPlaneAtOrigin(const std::vector<Plane>& somePlanes)
		{
			if (somePlanes.empty())
			{
				return -1;
			}
			const Line3 upLine{ Point3{0,0,0}, Vec3{0,0,1} };
			std::vector<std::pair<FT, int>> planeDepths;
			for (int i = 0; i < somePlanes.size(); ++i)
			{
				const auto intersection = CGAL::intersection(upLine, somePlanes[i]);
				CGAL_precondition(intersection.has_value());
				const Point3* point = boost::get<Point3>(&*intersection);
				CGAL_precondition(point != nullptr);
				planeDepths.push_back(std::make_pair(point->z(), i));
			}
			return CGAL::min_max_element(planeDepths.begin(), planeDepths.end()).first->second;
		}
		bool locIsVertexInLowerEnvelope(const std::vector<Plane>& somePlanes, const Point3& aPoint)
		{
			auto isValid{ true };

			std::vector<FT> zs;
			zs.push_back(aPoint.z());
			std::cout << "--Y: " << zs.back() << std::endl;
			for (int i = 0; i < somePlanes.size() && isValid; ++i)
			{
				const auto inter = CGAL::intersection(somePlanes[i], Line3{ aPoint, Vec3{0,0,1} });
				
				CGAL_precondition(inter.has_value());
				const Point3* point = boost::get<Point3>(&*inter);
				CGAL_precondition(point != nullptr);
				//planeDepths.push_back(std::make_pair(point->y(), i));
				zs.push_back(point->z());
				/*if (somePlanes[i].has_on_positive_side(aPoint))
				{
					isValid = false;
				}*/
				std::cout << "Y: " << zs.back() << std::endl;
			}

			isValid = *CGAL::min_max_element(zs.begin(), zs.end()).first == aPoint.z();

			return isValid;
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
				const auto angleDegrees = CGAL::approximate_angle(Vec3{ 1,0,0 }, segment);

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
		int locLowestPlaneIndexThroughSegment(const Point3& aStart, const Point3& anEnd, const std::vector<Plane>& somePlanes, const std::vector<int>& somePlanesIndices)
		{
			const Vec3 upDir{ 0,1,0 };

			CGAL_precondition(somePlanes.size() > 0);
			CGAL_precondition(somePlanesIndices.size() > 0);

			Vec3 segmentDir{ aStart, anEnd };
			segmentDir = locNormalize(segmentDir);

			const Vec3 leftDir{ CGAL::cross_product(segmentDir, upDir) };

			FT maxSteepness = CGAL::scalar_product(somePlanes[somePlanesIndices[0]].orthogonal_vector(), leftDir);
			int planeIdx = somePlanesIndices[0];

			for (int i = 1; i < somePlanesIndices.size(); ++i)
			{
				const auto currentPlaneIdx = somePlanesIndices[i];
				const auto dot = CGAL::scalar_product(somePlanes[currentPlaneIdx].orthogonal_vector(), leftDir);

				if (dot > maxSteepness)
				{
					maxSteepness = dot;
					planeIdx = currentPlaneIdx;
				}
			}

			return planeIdx;
		}

		CGAL::Bbox_3 ComputeBBox(const std::vector<Point3>& somePoints, const FT aPadding = 0.f)
		{
			std::vector<Point3> copy{ somePoints };
			const auto minMaxPair = CGAL::min_max_element(somePoints.begin(), somePoints.end());
			copy.push_back(Point3{
				minMaxPair.first->x() - aPadding,
				minMaxPair.first->y() - aPadding,
				minMaxPair.first->z() - aPadding });
			copy.push_back(Point3{
				minMaxPair.second->x() + aPadding,
				minMaxPair.second->y() + aPadding,
				minMaxPair.second->z() + aPadding });
			return CGAL::bbox_3(copy.begin(), copy.end());
		}
		Sphere3 ComputeBoundingSphere(const std::vector<Point3>& somePoints, const FT aPadding = 0.f)
		{
			CGAL_precondition(somePoints.size() > 0);
			const auto& sphereOrigin = CGAL::centroid(somePoints.begin(), somePoints.end());
			std::vector<FT> distances;
			distances.reserve(somePoints.size());
			const auto transform = [&sphereOrigin](const auto& aPoint) { return CGAL::squared_distance(aPoint, sphereOrigin); };
			std::transform(somePoints.begin(), somePoints.end(), std::back_inserter(distances), transform);
			const auto& minMaxDistances = CGAL::min_max_element(distances.begin(), distances.end());
			const auto sphereRadius = *minMaxDistances.second + aPadding;
			return Sphere3(sphereOrigin, sphereRadius);
		}

		Vec3 locGetLineDirection(const LineData& aLineData)
		{
			constexpr auto distance = 1000.f;
			const auto& a = aLineData.myLine.point();
			const auto& b = a + aLineData.myLine.to_vector() * distance;
			if (a < b)
			{
				return locNormalize(Vec3{ a, b });
			}
			else
			{
				return locNormalize(Vec3{ b, a });
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

					start.myLowestLeftPlanes.push_back(locLowestPlaneIndexThroughSegment(start.myPoint, end.myPoint, somePlanes, currentLineData.myPlanesIndices));
					end.myLowestLeftPlanes.push_back(locLowestPlaneIndexThroughSegment(end.myPoint, start.myPoint, somePlanes, currentLineData.myPlanesIndices));

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
			constexpr auto distance = 1000.f;
			const auto& uniformLineDirection = locGetLineDirection(currentLine);
			const auto& positive = linesVerticesData.myUniqueVertices[vertexIndices.back()].myPoint + uniformLineDirection * distance;
			const auto& negative = linesVerticesData.myUniqueVertices[vertexIndices.front()].myPoint - uniformLineDirection * distance;

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

				std::cout << "Segment: " << startVertex.myPoint << "  -  " << endVertex.myPoint << " | midpoint: " << midpoint << std::endl;

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

					if (!endVertex.myIsAtInfinity && !startVertex.myIsAtInfinity)
					{
						std::cout << "aaa" << std::endl;
					}

					

					auto& start = result[startIter->myIndex];
					auto& end = result[endIter->myIndex];

					std::cout << start.myPoint << "  --  " << end.myPoint << std::endl;

					start.mySortedNeighboursIndices.push_back(endIter->myIndex);
					start.myLowestLeftPlanes.push_back(locLowestPlaneIndexThroughSegment(start.myPoint, end.myPoint, 
						somePlanes, currentLine.myPlanesIndices));

					end.mySortedNeighboursIndices.push_back(startIter->myIndex);
					end.myLowestLeftPlanes.push_back(locLowestPlaneIndexThroughSegment(end.myPoint, start.myPoint, 
						somePlanes, currentLine.myPlanesIndices));

					segmentStartIdx = -1;
				}
			}
		}

		// Sort each edge counter-clock wise
		for (int i = 0; i < result.size(); ++i)
		{
			//locSortClockWise(result[i], result);

			// Set vertex to infinity if it has only one edge
			if (result[i].mySortedNeighboursIndices.size() == 1)
			{
				result[i].myType = VertexType::INFINITE;
			}
		}

		return result;
	}
}