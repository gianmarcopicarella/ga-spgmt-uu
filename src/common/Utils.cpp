#include "Utils.h"

namespace SPGMT
{
	namespace Utils
	{
		namespace
		{
			bool locIsValidPointForDuality(const Point3& aPoint)
			{
				return aPoint.x() != 0 || aPoint.y() != 0 || aPoint.z() != 0;
			}

			bool locIsValidPlaneForDuality(const Plane& aPlane)
			{
				return aPlane.a() != 0 || aPlane.b() != 0 || aPlane.c() != 0;
			}
		}

		bool ArePlanesNonVertical(const std::vector<Plane>& somePlanes)
		{
			static const Vec3 verticalAxis{ 0,0,1 };
			static const FT zero{ 0 };
			auto arePlanesNonVertical{ true };
			for (int i = 0; i < somePlanes.size() && arePlanesNonVertical; ++i)
			{
				arePlanesNonVertical =
					CGAL::scalar_product(somePlanes[i].orthogonal_vector(), verticalAxis) != zero;
			}
			return arePlanesNonVertical;
		}

		bool ArePlanesUniformlyOriented(const std::vector<Plane>& somePlanes)
		{
			static const Vec3 up{ 0, 0, 1 };
			static const FT zero{ 0 };
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
			return true;
		}

		std::vector<Point3> DualMapping(const std::vector<Plane>& somePlanes)
		{
			CGAL_precondition(false);
			std::vector<Point3> result;
			result.reserve(somePlanes.size());
			
			for (auto i = 0; i < somePlanes.size(); ++i)
			{
				CGAL_precondition(locIsValidPlaneForDuality(somePlanes[i]));
				const auto& x = somePlanes[i].a();
				const auto& y = somePlanes[i].b();
				const auto& z = -somePlanes[i].c();
				result.push_back(Point3{ x,y,z });
			}

			return result;
		}
		
		// Reference (also for an improved version avoiding the creation of a potential useless plane object)
		// https://math.stackexchange.com/questions/1375308/flip-normal-of-plane
		std::vector<Plane> DualMapping(const std::vector<Point3>& somePoints)
		{
			std::vector<Plane> result;
			result.reserve(somePoints.size());
			const Vec3 up{ 0,0,1 };
			for (auto i = 0; i < somePoints.size(); ++i)
			{
				CGAL_precondition(locIsValidPointForDuality(somePoints[i]));
				const auto& a = somePoints[i].x();
				const auto& b = somePoints[i].y();
				const auto& c = -somePoints[i].z();
				const auto& d = a + b + c;
				result.push_back(Plane{ a,b,c,d });
			}
			return result;
		}

		std::vector<Plane> EnforceUpwardsOrientation(const std::vector<Plane>& somePlanes)
		{
			std::vector<Plane> result;
			result.reserve(somePlanes.size());
			const Vec3 up{ 0,0,1 };

			for (auto i = 0; i < somePlanes.size(); ++i)
			{
				if (CGAL::scalar_product(up, somePlanes[i].orthogonal_vector()) < 0.f)
				{
					result.push_back(Plane{ somePlanes[i].point(), -somePlanes[i].orthogonal_vector() });
				}
				else
				{
					result.push_back(somePlanes[i]);
				}
			}

			return result;
		}
	}
}