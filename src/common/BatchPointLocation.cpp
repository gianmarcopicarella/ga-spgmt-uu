#include "BatchPointLocation.h"

#include <vector>

#include <CGAL/Arrangement_2.h>

namespace SPGMT
{
	struct LineInfo
	{
		int myLineIndex;
		FT myY;
		FT mySlope;
	};

	struct Slab
	{
		FT myValue;
		std::vector<LineInfo> myLinesInfo;
		int myVerticalLineIdx{ -1 };
		std::unordered_map<int, std::vector<int>> mySortedPlanesIndices_BS;
		std::vector<int> mySortedPlanesIndicesOnVertex_BS;
		std::vector<int> myIntersectionsIndices;
	};

	namespace Utils
	{
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
	}

	namespace
	{
		struct PlaneIntersectionStrategy
		{
			typedef struct PlaneIntersectionData
			{
				std::unordered_set<int> myCreationPlanesIndices;
				Line3 myLine;

				Line2 GetLine2() const
				{
					// DO NOT USE a static member to avoid the creation of the 2d line everytime.
					// In that case the static member would be initialized once for all instances of this struct
					constexpr auto distance = 1000.f;
					const auto& translation = myLine.to_vector() * distance;
					const auto firstPoint{ myLine.point() };
					const auto secondPoint{ firstPoint + translation };

					const auto line2d = Utils::Create2DLine(
						Point2{ firstPoint.x(), firstPoint.z() },
						Point2{ secondPoint.x(), secondPoint.z() });

					return line2d;
				}

			} Data;

			typedef void result_type;
			void operator()(const Line3& aLine)
			{
				Data data;
				data.myCreationPlanesIndices.insert(myFirstPlaneIdx);
				data.myCreationPlanesIndices.insert(mySecondPlaneIdx);
				data.myLine = aLine;

				myOutData.push_back(data);
			}

			void operator()(const Plane& /*aPlane*/)
			{
				constexpr auto nonEqualPlanes = false;
				CGAL_precondition(nonEqualPlanes);
			}

			int myFirstPlaneIdx;
			int mySecondPlaneIdx;
			std::vector<Data>& myOutData;
		};

		template<typename T, typename K, typename L>
		struct ParallelItemsStrategy
		{
			using Data = std::pair<FT, int>;
			typedef void result_type;
			void operator()(const K& aPoint)
			{
				const L vecToPoint{ myLine.point(), aPoint };
				const auto distanceSign = CGAL::scalar_product(vecToPoint, myLine.to_vector()) > 0.f ? 1.f : -1.f;
				myDistancePlanesIdxPairs.push_back(std::make_pair(distanceSign * vecToPoint.squared_length(), myPlaneIndex));
			}

			void operator()(const T& /*aLine*/)
			{
				constexpr auto isPlaneOrthogonal = false;
				CGAL_precondition(isPlaneOrthogonal);
			}

			const int myPlaneIndex;
			const T& myLine;
			std::vector<std::pair<FT, int>>& myDistancePlanesIdxPairs;
		};

		struct LineIntersectionStrategy
		{
			using Data = Point2;
			typedef void result_type;
			void operator()(const Point2& aPoint)
			{
				myIntersections.push_back(aPoint);
			}

			void operator()(const Line2& /*aLine*/)
			{
				constexpr auto nonEqualLines = false;
				CGAL_precondition(nonEqualLines);
			}

			int myFirstLineIdx;
			int mySecondLineIdx;
			std::vector<Point2>& myIntersections;
		};

		template<typename Strategy, typename InType>
		std::vector<typename Strategy::Data> HandleIntersections(const std::vector<InType>& someItems)
		{
			std::vector<typename Strategy::Data> result;
			for (int i = 0; i < someItems.size(); ++i)
			{
				for (int k = i + 1; k < someItems.size(); ++k)
				{
					const auto& intersection = CGAL::intersection(someItems[i], someItems[k]);
					if (const auto* data = intersection.get_ptr())
					{
						Strategy visitor{ i, k, result };
						data->apply_visitor(visitor);
					}
				}
			}
			return result;
		}

		template<typename Strategy, typename InFirstType, typename InSecondType>
		std::vector<typename Strategy::Data> HandleIntersectionsWith(const std::vector<InFirstType>& someItems, const InSecondType& anotherItem)
		{
			std::vector<typename Strategy::Data> result;
			for (int i = 0; i < someItems.size(); ++i)
			{
				const auto& intersection = CGAL::intersection(someItems[i], anotherItem);
				if (const auto* data = intersection.get_ptr())
				{
					Strategy visitor{ i, anotherItem, result };
					data->apply_visitor(visitor);
				}
			}
			return result;
		}

		Point2 Project3DPoint(const Point3& aPoint)
		{
			return Point2{ aPoint.x(), aPoint.z() };
		}

	}

	std::vector<Slab> ComputeSortedSlabs(const std::vector<Line2>& someLines, const std::vector<Point2>& someSortedPoints)
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

		std::vector<Slab> slabs;
		std::vector<std::optional<FT>> slopes;

		// Compute slopes for each line
		std::transform(someLines.begin(), someLines.end(), std::back_inserter(slopes), slopeTransform);

		// Insert slab at -inf (always at index 0)
		slabs.push_back(Slab{});

		// Compute remaining slabs
		for (auto i = 0; i < someSortedPoints.size(); ++i)
		{
			if (slabs.size() > 1 && slabs.back().myValue == someSortedPoints[i].x())
			{
				slabs.back().myIntersectionsIndices.push_back(i);
			}
			else
			{
				Slab currentSlab{ someSortedPoints[i].x() };
				currentSlab.myIntersectionsIndices.push_back(i);
				// This inner loop could be optimized by sorting the lines for the first slab and 
				// then flip the intersecting lines indices at each new slab
				for (int k = 0; k < someLines.size(); ++k)
				{
					const auto& currentLine = someLines[k];
					if (currentLine.is_vertical())
					{
						// A vertical line must not have a value
						CGAL_precondition(!slopes[k].has_value());
						currentSlab.myVerticalLineIdx = k;
					}
					else
					{
						// A non-vertical line must have a value
						CGAL_precondition(slopes[k].has_value());
						currentSlab.myLinesInfo.push_back(LineInfo{ k, currentLine.y_at_x(currentSlab.myValue), *slopes[k] });
					}
				}

				std::sort(currentSlab.myLinesInfo.begin(), currentSlab.myLinesInfo.end(), sort);
				slabs.push_back(currentSlab);
			}
		}

		// -Infinity slab custom sort case
		// Copy the lines order from the first slab but in case of same Y-intercept the slope sign is flipped while sorting
		std::copy(slabs[1].myLinesInfo.begin(), slabs[1].myLinesInfo.end(), std::back_inserter(slabs[0].myLinesInfo));
		std::sort(slabs[0].myLinesInfo.begin(), slabs[0].myLinesInfo.end(), sortInverseSlope);

		return slabs;
	}

	struct SlabSearchResult
	{
		int myIndex{ -1 };
		bool myIsOverLine{ false };
	};

	int BinarySearchSlab2(
		const std::vector<Slab>& someItems,
		const Point2& aPoint,
		const int aLow,
		const int anHigh)
	{
		const int middle = (aLow + anHigh) / 2;
		const auto& currentSlab = someItems[middle];

		if (currentSlab.myValue > aPoint.x())
		{
			if (middle == 0)
			{
				return 0;
			}
			else if (someItems[middle - 1].myValue < aPoint.x())
			{
				return middle;
			}
			else
			{
				return BinarySearchSlab2(someItems, aPoint, aLow, middle - 1);
			}
		}
		else if (currentSlab.myValue < aPoint.x())
		{
			if (middle == someItems.size() - 1)
			{
				return someItems.size() - 1;
			}
			else if (someItems[middle + 1].myValue > aPoint.x())
			{
				return middle + 1;
			}
			else
			{
				return BinarySearchSlab2(someItems, aPoint, middle + 1, anHigh);
			}
		}
		else
		{
			return middle + 1;
		}
	}

	
	int BinarySearchSlab(
		const std::vector<Slab>& someItems,
		const Point2& aPoint,
		const int aLow,
		const int anHigh)
	{
		const int middle = (aLow + anHigh) / 2;
		const auto& currentSlab = someItems[middle];

		if (currentSlab.myValue > aPoint.x())
		{
			if (middle == 0)
			{
				return 0;
			}
			else if (someItems[middle - 1].myValue < aPoint.x())
			{
				return middle;
			}
			else
			{
				return BinarySearchSlab(someItems, aPoint, aLow, middle - 1);
			}
		}
		else if (currentSlab.myValue < aPoint.x())
		{
			if (middle == someItems.size() - 1)
			{
				return someItems.size();
			}
			else if (someItems[middle + 1].myValue > aPoint.x())
			{
				return middle + 1;
			}
			else
			{
				return BinarySearchSlab(someItems, aPoint, middle + 1, anHigh);
			}
		}
		else
		{
			return middle + 1;
		}
	}

	SlabSearchResult BinarySearchSlab(
		const std::vector<Line2>& someItems,
		const std::vector<int>& someSortedItemIndices,
		const Point2& aPoint,
		const int aLow,
		const int anHigh)
	{
		const int middle = (aLow + anHigh) / 2;
		const auto& currentLine = someItems[someSortedItemIndices[middle]];

		if (currentLine.has_on_negative_side(aPoint))
		{
			if (middle == 0)
			{
				return { 0, false };
			}
			else if (someItems[someSortedItemIndices[middle - 1]].has_on_positive_side(aPoint))
			{
				return { middle, false };
			}
			else
			{
				return BinarySearchSlab(someItems, someSortedItemIndices, aPoint, aLow, middle - 1);
			}
		}
		else if (currentLine.has_on_positive_side(aPoint))
		{
			if (middle == someItems.size() - 1)
			{
				return { static_cast<int>(someItems.size()), false };
			}
			else if (someItems[someSortedItemIndices[middle + 1]].has_on_negative_side(aPoint))
			{
				return { middle + 1, false };
			}
			else
			{
				return BinarySearchSlab(someItems, someSortedItemIndices, aPoint, middle + 1, anHigh);
			}
		}
		else
		{
			// Lies over a slab boundary, that is, a line zone
			return { middle + 1, true };
		}
	}

	std::pair<int, int> BinarySearchItems(
		const std::vector<Plane>& someItems,
		const std::vector<int>& someSortedItemIndices,
		const Point3& aPoint,
		const int aLow,
		const int anHigh)
	{
		const int middle = (aLow + anHigh) / 2;
		const auto& currentPlane = someItems[someSortedItemIndices[middle]];

		if (currentPlane.has_on_negative_side(aPoint))
		{
			if (middle == 0)
			{
				return std::make_pair(0, static_cast<int>(someItems.size()));
			}
			else if (someItems[someSortedItemIndices[middle - 1]].has_on_positive_side(aPoint))
			{
				return std::make_pair(middle, static_cast<int>(someItems.size()));
			}
			else
			{
				return BinarySearchItems(someItems, someSortedItemIndices, aPoint, aLow, middle - 1);
			}
		}
		else if (currentPlane.has_on_positive_side(aPoint))
		{
			if (middle == someItems.size() - 1)
			{
				return std::make_pair(-1, -1);
			}
			else if (someItems[someSortedItemIndices[middle + 1]].has_on_negative_side(aPoint))
			{
				return std::make_pair(middle + 1, static_cast<int>(someItems.size()));
			}
			else
			{
				return BinarySearchItems(someItems, someSortedItemIndices, aPoint, middle + 1, anHigh);
			}
		}
		else if (middle < someItems.size() - 1)
		{
			return std::make_pair(middle + 1, static_cast<int>(someItems.size()));
		}
		else
		{
			return std::make_pair(-1, -1);
		}
	}

	bool ArePlanesNonVertical(const std::vector<Plane>& somePlanes)
	{
		static const Vec3 verticalAxis{ 0,1,0 };
		auto arePlanesNonVertical{ true };
		for (int i = 0; i < somePlanes.size() && arePlanesNonVertical; ++i)
		{
			arePlanesNonVertical =
				CGAL::scalar_product(somePlanes[i].orthogonal_vector(), verticalAxis) != 0.f;
		}
		return arePlanesNonVertical;
	}

	/*
	Assumptions:
	0) No duplicate planes
	1) No vertical planes -> should simplify the implementation
	2) A plane is above p if p is not contained in the plane and plane(p.x, p.z).y > p.y
	3) Every plane should always point upwards
	*/
	LocationResult BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints)
	{
		const auto sort = [](const auto& aFirst, const auto& aSecond) { return aFirst.first < aSecond.first; };
		const auto transform = [](const auto& anItem) { return anItem.second; };

		LocationResult result;

		if (somePlanes.size() == 0 || somePoints.size() == 0)
		{
			return result;
		}

		// No duplicate planes
		CGAL_precondition(Debug::AreItemsUnique(somePlanes));

		// No vertical planes
		CGAL_precondition(ArePlanesNonVertical(somePlanes));

		// 1) Edge case: single plane
		if (somePlanes.size() == 1)
		{
			OnePlaneResult resultPayload;
			resultPayload.myIsPointCovered.reserve(somePoints.size());

			for (const auto& queryPoint : somePoints)
			{
				const auto isPointCovered = somePlanes.back().has_on_negative_side(queryPoint);
				resultPayload.myIsPointCovered.push_back(isPointCovered);
			}

			result = resultPayload;
			return result;
		}

		// Compute plane-plane intersection lines
		const auto& linesData3D = HandleIntersections<PlaneIntersectionStrategy>(somePlanes);

		// 2) Edge case: parallel planes
		if (linesData3D.size() == 0)
		{
			ParallelPlanesResult resultPayload;

			auto distanceIndexPlanes =
				HandleIntersectionsWith<ParallelItemsStrategy<Line3, Point3, Vec3>>(somePlanes,
					somePlanes[0].perpendicular_line(somePlanes[0].point()));
			std::sort(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), sort);
			std::transform(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), std::back_inserter(resultPayload.mySortedPlanesIndices), transform);

			for (const auto& queryPoint : somePoints)
			{
				const auto range = BinarySearchItems(somePlanes, resultPayload.mySortedPlanesIndices, queryPoint, 0, somePlanes.size());
				resultPayload.myRanges.push_back(range);
			}

			result = resultPayload;
			return result;
		}

		// Get 2d unique lines
		std::vector<Line2> lines2D;
		std::transform(linesData3D.begin(), linesData3D.end(), std::back_inserter(lines2D),
			std::bind(&PlaneIntersectionStrategy::PlaneIntersectionData::GetLine2, std::placeholders::_1));
		lines2D.erase(std::unique(lines2D.begin(), lines2D.end()), lines2D.end());

		// Compute intersection points
		auto intersectionPoints = HandleIntersections<LineIntersectionStrategy>(lines2D);

		// 3) Edge case: parallel 2d lines
		if (intersectionPoints.empty())
		{
			BaseResult resultPayload;

			// Sort 2D lines
			std::vector<int> linesIndices_BS;
			auto distanceIndexLines =
				HandleIntersectionsWith<ParallelItemsStrategy<Line2, Point2, Vec2>>(lines2D,
					lines2D.front().perpendicular(lines2D.front().point()));
			std::sort(distanceIndexLines.begin(), distanceIndexLines.end(), sort);
			std::transform(distanceIndexLines.begin(), distanceIndexLines.end(), std::back_inserter(linesIndices_BS), transform);

			// Required to consider, among all vertically intersecting planes, only one
			// +1 because the number of faces is equal to lineCount + 1
			std::vector<std::vector<int>> lineZones_BS{ lines2D.size() + 1, std::vector<int>{} };

			for (const auto& queryPoint : somePoints)
			{
				const auto& projectedPoint = Project3DPoint(queryPoint);
				const auto zoneResult = BinarySearchSlab(lines2D, linesIndices_BS, projectedPoint, 0, linesIndices_BS.size());

				// Sort planes indices for this zone
				if (resultPayload.mySortedPlanesPerFaceZone.find(zoneResult.myIndex) ==
					resultPayload.mySortedPlanesPerFaceZone.end())
				{
					std::vector<int> sortedIndices;
					auto distanceIndexPlanes =
						HandleIntersectionsWith<ParallelItemsStrategy<Line3, Point3, Vec3>>(somePlanes, Line3{ queryPoint, Vec3{0,1,0} });
					std::sort(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), sort);
					std::transform(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), std::back_inserter(sortedIndices), transform);

					// Store planes in the output
					resultPayload.mySortedPlanesPerFaceZone.insert(std::make_pair(zoneResult.myIndex, sortedIndices));

					// No line to the left of zone 0 will ever ask for the list of z-sorted planes!
					if (zoneResult.myIsOverLine)
					{
						// From the sorted list of indices, remove all the planes with same z by 
						// keeping only the first one encountered in the list
						auto& sortedIndices_BS = lineZones_BS[zoneResult.myIndex];
						const auto equal = [](const auto& aFirst, const auto& aSecond) { return aFirst.first == aSecond.first; };
						distanceIndexPlanes.erase(std::unique(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), equal),
							distanceIndexPlanes.end());
						std::transform(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), std::back_inserter(sortedIndices_BS), transform);
					}
				}

				// Fetch z-sorted planes indices for current slab and face/line
				const auto& sortedIndices = zoneResult.myIsOverLine ? lineZones_BS[zoneResult.myIndex] :
					resultPayload.mySortedPlanesPerFaceZone.at(zoneResult.myIndex);

				// Search for the first plane above the query point
				const auto planesRange = BinarySearchItems(somePlanes, sortedIndices, queryPoint, 0, sortedIndices.size());
				resultPayload.myZoneRangesPairs.push_back(BaseResult::ZoneRange{ planesRange, zoneResult.myIndex, !zoneResult.myIsOverLine });
			}

			result = resultPayload;
			return result;
		}

		// Keep only unique intersection points
		intersectionPoints.erase(std::unique(intersectionPoints.begin(), intersectionPoints.end()), intersectionPoints.end());

		// Sort unique intersection points along x-axis
		std::sort(intersectionPoints.begin(), intersectionPoints.end());

		// Compute slabs
		auto slabs = ComputeSortedSlabs(lines2D, intersectionPoints);

		BaseResult resultPayload;

		//int c = -1;

		for (const auto& queryPoint : somePoints)
		{
			/*if (++c == 7)
			{
				std::cout << queryPoint << std::endl;
			}*/

			const auto& projectedQueryPoint = Project3DPoint(queryPoint);
			const auto slabIndex = BinarySearchSlab2(slabs, projectedQueryPoint, 0, slabs.size());
			auto& currentSlab = slabs[slabIndex];

			const bool pointOnVerticalSlabLine = currentSlab.myVerticalLineIdx > -1 &&
				lines2D[currentSlab.myVerticalLineIdx].has_on(projectedQueryPoint);

			const bool isQueryPointOnVertex =
				std::find_if(currentSlab.myIntersectionsIndices.begin(), currentSlab.myIntersectionsIndices.end(), 
				[&projectedQueryPoint, &intersectionPoints](const auto& aPointIdx) {
					return projectedQueryPoint == intersectionPoints[aPointIdx];
				}) != currentSlab.myIntersectionsIndices.end();

			// 1) Query point is a vertex
			if (isQueryPointOnVertex)
			{
				CGAL_precondition(false);
			}

			// Create a list of sorted lines indices based on the current slab sorting and find the top enclosing line of query point
			std::vector<int> sortedLinesIndices;
			const auto linesTransform = [](const auto& aLineInfo) {return aLineInfo.myLineIndex; };
			std::transform(currentSlab.myLinesInfo.begin(), currentSlab.myLinesInfo.end(),
				std::back_inserter(sortedLinesIndices), linesTransform);
			
			const auto topLineInfo = BinarySearchSlab(lines2D, sortedLinesIndices, projectedQueryPoint, 0, sortedLinesIndices.size());
			const int globalVerticalFaceIdx = slabIndex * (lines2D.size() + 1) + topLineInfo.myIndex;
			const int sortedPlanesIdx_BS = pointOnVerticalSlabLine ? sortedLinesIndices.size() : topLineInfo.myIndex;


			if (resultPayload.mySortedPlanesPerFaceZone.find(globalVerticalFaceIdx) ==
				resultPayload.mySortedPlanesPerFaceZone.end())
			{
				std::vector<int> sortedIndices;
				auto distanceIndexPlanes =
					HandleIntersectionsWith<ParallelItemsStrategy<Line3, Point3, Vec3>>(somePlanes, Line3{ queryPoint, Vec3{0,1,0} });
				std::sort(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), sort);
				std::transform(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), std::back_inserter(sortedIndices), transform);

				// Store planes in the output
				resultPayload.mySortedPlanesPerFaceZone.insert(std::make_pair(globalVerticalFaceIdx, sortedIndices));

				if (topLineInfo.myIsOverLine || pointOnVerticalSlabLine)
				{
					// From the sorted list of indices, remove all the planes with same z by 
					// keeping only the first one encountered
					std::vector<int> sortedIndices_BS;
					const auto equal = [](const auto& aFirst, const auto& aSecond) { return aFirst.first == aSecond.first; };
					distanceIndexPlanes.erase(std::unique(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), equal),
						distanceIndexPlanes.end());
					std::transform(distanceIndexPlanes.begin(), distanceIndexPlanes.end(), std::back_inserter(sortedIndices_BS), transform);
					currentSlab.mySortedPlanesIndices_BS.insert(std::make_pair(sortedPlanesIdx_BS, sortedIndices_BS));
				}
			}

			CGAL_precondition(!pointOnVerticalSlabLine);
			CGAL_precondition(!topLineInfo.myIsOverLine);

			// Fetch z-sorted planes indices for current slab and face/line
			const auto& sortedIndices = topLineInfo.myIsOverLine ? 
				currentSlab.mySortedPlanesIndices_BS[sortedPlanesIdx_BS] :
				resultPayload.mySortedPlanesPerFaceZone.at(globalVerticalFaceIdx);

			// Search for the first plane above the query point
			const auto planesRange = BinarySearchItems(somePlanes, sortedIndices, queryPoint, 0, sortedIndices.size());
			resultPayload.myZoneRangesPairs.push_back(BaseResult::ZoneRange{ planesRange, globalVerticalFaceIdx, !topLineInfo.myIsOverLine });
		}

		result = resultPayload;
		return result;
	}
}