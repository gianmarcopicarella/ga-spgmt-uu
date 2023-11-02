#pragma once

#include "Types.h"

#include <vector>

namespace SPGMT
{
	namespace Utils
	{
		template<typename T>
		bool AreItemsUnique(const std::vector<T>& someItems)
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

		bool ArePlanesNonVertical(const std::vector<Plane>& somePlanes);
		bool ArePlanesUniformlyOriented(const std::vector<Plane>& somePlanes);

		template <typename Strategy, typename T>
		std::vector<typename Strategy::payload_type> FindIntersections(const std::vector<T>& someItems, const bool aKeepOnlyUnique = false)
		{
			std::vector<typename Strategy::payload_type> result;
			for (int i = 0; i < someItems.size(); ++i)
			{
				for (int k = i + 1; k < someItems.size(); ++k)
				{
					const auto intersection = CGAL::intersection(someItems[i], someItems[k]);
					if (const auto* data = intersection.get_ptr())
					{
						Strategy visitor{ result };
						data->apply_visitor(visitor);
					}
				}
			}
			if (aKeepOnlyUnique)
			{
				result.erase(std::unique(result.begin(), result.end(), CGAL::Equal_to<typename Strategy::payload_type,
					typename Strategy::payload_type>()), result.end());
			}
			return result;
		}

		template <typename Strategy, typename T, typename K>
		std::vector<typename Strategy::payload_type> FindIntersections(const std::vector<T>& someItems, const K& anotherItem)
		{
			constexpr auto policy = Strategy::GetExecutionPolicy();
			std::vector<typename Strategy::payload_type> result;
			if constexpr (policy != ExecutionPolicy::SEQ)
			{
				result.resize(someItems.size());
			}

			BindExecutionPolicy<policy>(hpx::for_loop, 0, someItems.size(), [&](size_t anItemIdx) {
				const auto intersection = CGAL::intersection(someItems[anItemIdx], anotherItem);
				if (const auto* data = intersection.get_ptr())
				{
					CGAL_precondition(data != nullptr);
					Strategy visitor{ anotherItem, anItemIdx, result };
					data->apply_visitor(visitor);
				}
				});

			return result;
		}

		template<ExecutionPolicy E, typename T, typename F>
		bool AreItemsParallel(
			const std::vector<T>& someItems,
			const F&& aGetter)
		{
			CGAL_precondition(someItems.size() > 0);
			const auto& direction = aGetter(someItems[0]);
			const auto& oppDirection = -direction;
			const auto result = BindExecutionPolicy<E>(hpx::all_of, std::next(someItems.begin()), someItems.end(),
				[&](const auto& anItem) {
					const auto& currDirection = aGetter(anItem);
					return currDirection == direction || currDirection == oppDirection;
				});
			return result;
		}

		bool IsPlaneFacingUp(const Plane& aPlane);

		template<typename T>
		void TrimToLastValidItem(std::vector<T>& someOutItems)
		{
			// Still using a sequential version because it is fast
			auto it = someOutItems.begin();
			while (std::get<0>(*it) == STATUS::INIT && ++it != someOutItems.end());
			someOutItems.resize(std::distance(someOutItems.begin(), it));
		}

		Point3 SingleDualMapping(const Plane& aPlane);
		Plane SingleDualMapping(const Point3& aPoint);

		std::vector<Point3> DualMapping(const std::vector<Plane>& somePlanes);
		std::vector<Plane> DualMapping(const std::vector<Point3>& somePoints);
	}
}