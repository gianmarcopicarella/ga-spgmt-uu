#include <iostream>
#include <vector>
#include <algorithm>

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Line_3.h>

#include <CGAL/intersections.h>
#include <boost/range/algorithm/set_algorithm.hpp>

typedef CGAL::Exact_rational ExactReal;
typedef CGAL::Cartesian<ExactReal> Cartesian;
typedef Cartesian::Plane_3 Plane;
typedef Cartesian::Point_3 Point3;
typedef Cartesian::Direction_3 Dir3;
typedef Cartesian::Vector_3 Vec3;

typedef Cartesian::Line_3 Line3;


/*
04/09/2023

1) Setup kernel and mapping. DONE
2) Setup asserts/preconditions. DONE
3) Understand what is the default plane orientation (positive or negative). DONE (I hope)
4) Implement Plane-Plane-Plane intersection. (Not needed, CGAL has one)
5) Check if a 3d point lies below (or lies on the surface) all the other plane. // Done
*/

/* Not needed for now
// Ref, https://mathoverflow.net/questions/87688/find-an-n-dimensional-vector-orthogonal-to-a-given-vector
Dir3 ComputeOrthogonalDir3(const Dir3& aDirection)
{
	if (aDirection.dx() == ExactReal(0))
	{
		return Dir3{ 1,0,0 };
	}
	else if (aDirection.dy() == ExactReal(0))
	{
		return Dir3{ 0,1,0 };
	}
	else
	{
		return Dir3{ -aDirection.dy(), aDirection.dx(), 0 };
	}
}
*/


// Ref. https://math.stackexchange.com/questions/1755856/calculate-arbitrary-points-from-a-plane-equation
Plane GetPlaneFromEqCoeffs(const ExactReal& aA, const ExactReal& aB, const ExactReal& aC, const ExactReal& aD)
{
	const auto sqrSum = aA * aA + aB * aB + aC * aC;
	const auto tempPoint = (Vec3{ aA, aB, aC } *aD) / sqrSum;

	const Point3 point = { tempPoint.x(), tempPoint.y(), tempPoint.z() };
	const Dir3 normal = { aA, aB, aC };

	return Plane{ point, normal };
}


Plane GetPlaneFromPoint(const Point3& aPoint)
{
	// Point-Plane mapping requirement: at least one coordinate must be != 0
	CGAL_precondition(
		aPoint.x() != ExactReal(0) ||
		aPoint.y() != ExactReal(0) ||
		aPoint.z() != ExactReal(0)
	);

	// Plane Eq. -Ax -By + (1)z + d = 0 -> This way every plane normal should point upwards
	const auto a = -aPoint.x();
	const auto b = -aPoint.y();
	const ExactReal c = 1;
	const auto d = aPoint.z();

	return GetPlaneFromEqCoeffs(a, b, c, d);
}


std::vector<Plane> GetPlanesFromPoints(const std::vector<Point3>& somePoints)
{
	std::vector<Plane> planes;

	for (const auto& aPoint : somePoints)
	{
		Plane& plane = GetPlaneFromPoint(aPoint);
		planes.push_back(plane);
	}

	return planes;
}


struct VertexWrapper
{
	const Point3 myPoint;
	const std::array<int, 3> myPlanesIndices;
};


struct PlaneIntersectionVisitor
{
	typedef void result_type;

	void operator()(const Point3& aPoint)
	{
		for (int i = 0; i < myPlanes.size(); ++i)
		{
			const auto& iter = std::find(myPlanesToSkipIndices.begin(), myPlanesToSkipIndices.end(), i);
			if (iter != myPlanesToSkipIndices.end())
			{
				// The intersection point should always be lie the negative side of each plane, otherwise it is not valid
				if (myPlanes[i].has_on(aPoint) == CGAL::ON_POSITIVE_SIDE)
				{
					return;
				}
			}
		}

		myOutValidVertices.push_back(VertexWrapper{ aPoint, myPlanesToSkipIndices });
	}

	void operator()(const Line3& aLine)
	{
		CGAL_precondition(0);
	}

	void operator()(const Plane& aPlane)
	{
		CGAL_precondition(0);
	}


	const std::vector<Plane>& myPlanes;
	const std::array<int, 3> myPlanesToSkipIndices;
	std::vector<VertexWrapper>& myOutValidVertices;
};


#define ADJ_IDX(s, e, WIDTH) ((s) + (WIDTH) * (e))


std::vector<int> ComputeLowerEnvelopeBF(const std::vector<Point3>& somePoints)
{
	// 0.0) Points to Planes mapping
	const auto planes = GetPlanesFromPoints(somePoints);

	// 1.0) Compute all the valid vertices of the lower envelope
	std::vector<VertexWrapper> vertices;

	for (int i = 0; i < planes.size(); ++i)
	{
		for (int j = i + 1; j < planes.size(); ++j)
		{
			for (int k = j + 1; k < planes.size(); ++k)
			{
				PlaneIntersectionVisitor visitor{ planes, std::array<int, 3>{i, j, k}, vertices };

				if (const auto resultOpt = CGAL::intersection(planes[i], planes[j], planes[k]))
				{
					boost::apply_visitor(visitor, *resultOpt);
				}
			}
		}
	}

	// 1.0) Compute all the valid edges of the lower envelope
	std::vector<bool> adjacency;
	adjacency.resize(vertices.size() * vertices.size(), false);

	for (int i = 0; i < vertices.size(); ++i)
	{
		for (int j = i + 1; j < vertices.size(); ++j)
		{
			const auto& u = vertices[i];
			const auto& v = vertices[j];
			
			std::array<int, 3> intersection;

			std::set_intersection(
				u.myPlanesIndices.begin(), u.myPlanesIndices.end(), 
				v.myPlanesIndices.begin(), v.myPlanesIndices.end(), 
				intersection.begin()
			);
			
			const auto isSegmentValid = intersection.size() > 1;

			if (isSegmentValid)
			{
				adjacency[ADJ_IDX(i, j, vertices.size())] = true;
			}
		}
	}


	return std::vector<int>{};
}

int main()
{
	std::vector<Point3> points;
	points.push_back(Point3(0, 0, 1, 1));

	ComputeLowerEnvelopeBF(points);

	return 0;
}