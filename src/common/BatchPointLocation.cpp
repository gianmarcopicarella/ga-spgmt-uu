#include "BatchPointLocation.h"

#include <vector>

namespace SPGMT
{
	struct LineInfo
	{
		int myLineIndex;
		FT myYIntercept;
		FT mySlope;
	};

	struct Slab
	{
		FT myValue;
		std::vector<LineInfo> myLinesInfo;
		bool myStartsAtMinusInfinity{ false };
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
		struct PlaneIntersectionVisitor
		{
			typedef void result_type;
			void operator()(const Line3& aLine)
			{
				myIntersections.push_back(aLine);
			}

			void operator()(const Plane& /*aPlane*/)
			{
				constexpr auto nonEqualPlanes = false;
				CGAL_precondition(nonEqualPlanes);
			}

			std::vector<Line3>& myIntersections;
		};

		struct LineIntersectionVisitor
		{
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

			std::vector<Point2>& myIntersections;
		};

		struct Line3DIntersectionVisitor
		{
			typedef void result_type;
			void operator()(const Point3& aPoint)
			{
				const auto& lineDir = myLine.direction();
				const Vec3 lineVec{ lineDir.dx(), lineDir.dy(), lineDir.dz() };

				const Vec3 vec{ myLine.point(), aPoint };
				const Dir3 dir{ vec };

				if (CGAL::scalar_product(Vec3{ dir.dx(), dir.dy(), dir.dz() }, lineVec) < 0.f)
				{
					myIntersections.push_back(std::make_pair(-vec.squared_length(), myPlaneIndex));
				}
				else
				{
					myIntersections.push_back(std::make_pair(vec.squared_length(), myPlaneIndex));
				}
			}

			void operator()(const Line3& /*aLine*/)
			{
				constexpr auto isPlaneOrthogonal = false;
				CGAL_precondition(isPlaneOrthogonal);
			}

			const int myPlaneIndex;
			const Line3& myLine;
			std::vector<std::pair<FT, int>>& myIntersections;
		};

		std::vector<Line3> FindPlaneIntersectionLines(const std::vector<Plane>& somePlanes)
		{
			std::vector<Line3> intersections;

			for (int i = 0; i < somePlanes.size(); ++i)
			{
				for (int k = i + 1; k < somePlanes.size(); ++k)
				{
					const auto& firstPlane = somePlanes[i];
					const auto& secondPlane = somePlanes[k];
					const auto& result = CGAL::intersection(firstPlane, secondPlane);

					if (const auto* data = result.get_ptr())
					{
						PlaneIntersectionVisitor visitor{ intersections };
						data->apply_visitor(visitor);
					}
				}
			}

			return intersections;
		}

		Point2 ProjectToXZ(const Point3& aPoint)
		{
			return Point2{ aPoint.x(), aPoint.z() };
		}

		std::vector<Line2> ProjectToXZ(const std::vector<Line3>& someLines)
		{
			std::vector<Line2> projections;

			for (const auto& line : someLines)
			{
				const auto& lineDir = line.direction();
				const Vec3 lineVec{ lineDir.dx(), lineDir.dy(), lineDir.dz() };
				const Vec3 scaledLineVec = lineVec * 1000.f;
				const Point3 secondPoint{ scaledLineVec.x(), scaledLineVec.y(), scaledLineVec.z() };

				projections.push_back(Line2{ ProjectToXZ(line.point()), ProjectToXZ(secondPoint) });
			}

			return projections;
		}

		std::vector<Point2> FindLineIntersections(const std::vector<Line2>& someLines)
		{
			std::vector<Point2> intersections;

			for (int i = 0; i < someLines.size(); ++i)
			{
				for (int k = i + 1; k < someLines.size(); ++k)
				{
					const auto& firstLine = someLines[i];
					const auto& secondLine = someLines[k];

					const auto& result = CGAL::intersection(firstLine, secondLine);

					if (const auto* data = result.get_ptr())
					{
						LineIntersectionVisitor visitor{ intersections };
						data->apply_visitor(visitor);
					}
				}
			}
			return intersections;
		}

	}

	std::vector<Slab> ComputeSortedSlabs(const std::vector<Line2>& someLines, const std::vector<Point2>& someSortedPoints)
	{
		std::vector<Slab> slabs;

		// Compute slopes for each line
		std::vector<FT> slopes;
		std::transform(someLines.begin(), someLines.end(), std::back_inserter(slopes), [](const auto& aLine) {
			const auto& dir = aLine.direction();
			return dir.dy() / dir.dx();
			});

		// Insert -Infinity slab
		slabs.push_back(Slab{});
		slabs[0].myStartsAtMinusInfinity = true;

		// Compute remaining slabs
		for (auto i = 0; i < someSortedPoints.size(); ++i)
		{
			Slab currentSlab{ someSortedPoints[i].x() };

			for (int k = 0; k < someLines.size(); ++k)
			{
				const auto& currentLine = someLines[k];
				currentSlab.myLinesInfo.push_back(LineInfo{ k, currentLine.y_at_x(currentSlab.myValue), slopes[k] });
			}

			std::sort(currentSlab.myLinesInfo.begin(), currentSlab.myLinesInfo.end(), [](const auto& a, const auto& b) {
				return a.myYIntercept < b.myYIntercept || (a.myYIntercept == b.myYIntercept && a.mySlope < b.mySlope);
				});

			slabs.push_back(currentSlab);
		}

		// -Infinity slab custom sort case
		// Copy the lines order from the first slab but in case of same Y-intercept the slope sign is flipped while sorting
		std::copy(slabs[1].myLinesInfo.begin(), slabs[1].myLinesInfo.end(), std::back_inserter(slabs[0].myLinesInfo));
		std::sort(slabs[0].myLinesInfo.begin(), slabs[0].myLinesInfo.end(), [](const auto& a, const auto& b) {
			return a.myYIntercept < b.myYIntercept || (a.myYIntercept == b.myYIntercept && a.mySlope > b.mySlope);
			});

		return slabs;
	}

	std::pair<int, int> BinarySearchPlanes(const std::vector<Plane>& someSortedPlanes, const Point3& aPoint, const int aLow, const int anHigh, const bool anAreVerticalPlanes)
	{
		const int middle = (aLow + anHigh) / 2;
		const auto& currentPlane = someSortedPlanes[middle];
		
		if (currentPlane.has_on_negative_side(aPoint))
		{
			if (middle == 0)
			{
				if (anAreVerticalPlanes)
				{
					return std::make_pair(-1, -1);
				}
				else
				{
					return std::make_pair(0, someSortedPlanes.size());
				}
			}
			else if (someSortedPlanes[middle - 1].has_on_positive_side(aPoint))
			{
				if (anAreVerticalPlanes)
				{
					return std::make_pair(-1, -1);
				}
				else
				{
					return std::make_pair(middle, someSortedPlanes.size());
				}
			}
			else
			{
				return BinarySearchPlanes(someSortedPlanes, aPoint, aLow, middle - 1, anAreVerticalPlanes);
			}
		}
		else if (currentPlane.has_on_positive_side(aPoint))
		{
			if (middle == someSortedPlanes.size() - 1)
			{
				return std::make_pair(-1, -1);
			}
			else if (someSortedPlanes[middle + 1].has_on_negative_side(aPoint))
			{
				if (anAreVerticalPlanes)
				{
					return std::make_pair(-1, -1);
				}
				else
				{
					return std::make_pair(middle + 1, someSortedPlanes.size());
				}
			}
			else
			{
				return BinarySearchPlanes(someSortedPlanes, aPoint, middle + 1, anHigh, anAreVerticalPlanes);
			}
		}
		else
		{
			if (anAreVerticalPlanes)
			{
				return std::make_pair(middle, -1);
			}
			else
			{
				return std::make_pair(middle, someSortedPlanes.size());
			}
		}
	}

	/*
	Currently this function must behave correctly for the following cases:
	1) No Plane intersections found -> Parallel planes (vertical/horizontal cases too)
	2) No Plane intersections found -> Only one plane


	*/
	void BatchPointLocation(const std::vector<Plane>& somePlanes, const std::vector<Point3>& somePoints, LocationResult& anOutResult)
	{
		if (somePlanes.size() == 0 || somePoints.size() == 0)
		{
			return;
		}

		// There must be at least one plane (could be potentially removed in the future
		// CGAL_precondition(somePlanes.size() > 0);

		// There must be at least one point (could be potentially removed in the future
		// CGAL_precondition(somePoints.size() > 0);

		// There must be only unique planes in the list
		CGAL_precondition(Debug::AreItemsUnique(somePlanes));

		// Find 3D intersection lines for each pair of planes
		const auto& planeLines = FindPlaneIntersectionLines(somePlanes);
		auto& projectedLines = ProjectToXZ(planeLines);

		// Keep only unique intersection lines
		projectedLines.erase(std::unique(projectedLines.begin(), projectedLines.end()), projectedLines.end());

		// Find 2D intersection points for each pair of lines
		auto& lineIntersections = FindLineIntersections(projectedLines);

		// Keep only unique intersection points
		lineIntersections.erase(std::unique(lineIntersections.begin(), lineIntersections.end()), lineIntersections.end());

		// If intersections.size() == 0
		// -> Then we have three cases:
		// 1) Just one plane in the list -> for each point, either the plane is above or below. DONE
		// 2) All planes are vertically parallel -> if a point lies on one of these planes then return that plane, otherwise an empty list. DONE
		// 3) All planes are parallel -> sort by z each plane and use binary search for each point. DONE

		// Corner case: Parallel planes
		if (lineIntersections.size() == 0)
		{
			if (somePlanes.size() == 1)
			{
				for (const auto& queryPoint : somePoints)
				{
					anOutResult.myRanges.push_back(std::make_pair(-1, -1));

					if (somePlanes[0].has_on_negative_side(queryPoint) ||
						somePlanes[0].has_on(queryPoint))
					{
						anOutResult.myRanges.back().first = 0;
						anOutResult.myRanges.back().second = 1; // Planes count
					}
				}

				anOutResult.myType = ResultType::SINGLE_PLANE;
				return;
			}
			else
			{
				const auto orthogonalLine = somePlanes[0].perpendicular_line(somePlanes[0].point());

				// Find all intersections with each plane
				std::vector<std::pair<FT, int>> distanceIndexPairs;

				for (auto i = 0; i < somePlanes.size(); ++i)
				{
					const auto& plane = somePlanes[i];
					const auto& result = CGAL::intersection(orthogonalLine, plane);

					if (const auto* data = result.get_ptr())
					{
						Line3DIntersectionVisitor visitor{ i, orthogonalLine, distanceIndexPairs };
						data->apply_visitor(visitor);
					}
					else
					{
						// How is this possible?
						CGAL_precondition(false);
					}
				}

				// Sort planes based on the intersection point
				std::sort(distanceIndexPairs.begin(), distanceIndexPairs.end(),
					[](const auto& aFirst, const auto& aSecond) { return aFirst.first < aSecond.first; });

				// Sort indices for anOutResult
				std::transform(distanceIndexPairs.begin(), distanceIndexPairs.end(), std::back_inserter(anOutResult.mySortedPlanesIndices),
					[](const auto& anItem) {
						return anItem.second;
					});

				std::vector<Plane> sortedPlanes;

				std::transform(distanceIndexPairs.begin(), distanceIndexPairs.end(), std::back_inserter(sortedPlanes),
					[&orthogonalLine, &somePlanes](const auto& anItem) {
						// use line normal to force same direction on every plane
						return Plane{ somePlanes[anItem.second].point(), orthogonalLine.direction() };
					});

				// For each query point p, apply binary search
				for (const auto& point : somePoints)
				{
					const auto planesRange = BinarySearchPlanes(sortedPlanes, point, 0, sortedPlanes.size(), false);
					anOutResult.myRanges.push_back(planesRange);
				}

				// Return result
				anOutResult.myType = ResultType::PARALLEL_PLANES;
				return;
			}
		}
		else
		{

		}



		// sort intersection points
		//std::sort(intersections.begin(), intersections.end());

		//// All horizontal or vertical lines
		//if (sortedIntersections.size() == 0)
		//{
		//	if (projectedLines[0].is_horizontal())
		//	{
		//		// All horizontal
		//		// TODO: Compute slabs
		//	}
		//	else
		//	{
		//		// All vertical
		//		// TODO: Compute slabs
		//	}
		//}
		//else
		//{
		//	// Usual case
		//}

		//// Consider only non-vertical lines when sorting for each slab
		//std::vector<Line2> validLines;
		//std::copy_if(projectedLines.begin(), projectedLines.end(), std::back_inserter(validLines),
		//	[](const auto& aLine) {return !aLine.is_vertical(); });

		//const auto& slabs = ComputeSortedSlabs(validLines, sortedIntersections);
		//
		//
		//

		//

		//for (auto& p : slabs)
		//{
		//	std::cout << p.myValue << std::endl;
		//}

		return;
	}
}