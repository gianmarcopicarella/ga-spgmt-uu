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
				return (aPlane.a() != 0 || aPlane.b() != 0) && aPlane.c() == FT(-1);
			}
		}

		bool ArePlanesNonVertical(const std::vector<Plane>& somePlanes)
		{
			const Vec3 verticalAxis{ 0,0,1 };
			const FT zero{ 0 };
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
			const Vec3 up{ 0, 0, 1 };
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
			return true;
		}

		Point3 SingleDualMapping(const Plane& aPlane)
		{
			CGAL_precondition(locIsValidPlaneForDuality(aPlane));
			const auto& x = aPlane.a();
			const auto& y = aPlane.b();
			const auto& z = -aPlane.d();
			return Point3{ x,y,z };
		}

		Plane SingleDualMapping(const Point3& aPoint)
		{
			CGAL_precondition(locIsValidPointForDuality(aPoint));
			const auto& a = aPoint.x();
			const auto& b = aPoint.y();
			const auto& c = -1.f;
			const auto& d = -aPoint.z();
			return Plane{ a,b,c,d };
		}

		std::vector<Point3> DualMapping(const std::vector<Plane>& somePlanes)
		{
			std::vector<Point3> result;
			result.reserve(somePlanes.size());
			
			for (auto i = 0; i < somePlanes.size(); ++i)
			{
				result.emplace_back(SingleDualMapping(somePlanes[i]));
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
				result.emplace_back(SingleDualMapping(somePoints[i]));
			}
			return result;
		}

		bool IsPlaneFacingUp(const Plane& aPlane)
		{
			static const Vec3 up{ 0,0,1 };
			return CGAL::sign(CGAL::scalar_product(up, aPlane.orthogonal_vector())) != CGAL::Sign::NEGATIVE;
		}

		void FlipPlaneNormalsIfFacingDownwards(std::vector<Plane>& someOutPlanes)
		{
			for (auto i = 0; i < someOutPlanes.size(); ++i)
			{
				if (!IsPlaneFacingUp(someOutPlanes[i]))
				{
					someOutPlanes[i] = Plane{ someOutPlanes[i].point(), -someOutPlanes[i].orthogonal_direction() };
				}
			}
		}
	}
}