#include "Utils.h"

namespace SPGMT
{
	namespace Utils
	{
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
	}
}