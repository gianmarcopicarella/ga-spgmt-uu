#pragma once

#include "Types.h"

#include <vector>

namespace SPGMT
{
	namespace Utils
	{
		struct PlanePlaneVisitor
		{
			typedef Line3 payload_type;
			typedef void result_type;
			void operator()(const Line3& aLine)
			{
				myOutResult.push_back(aLine);
			}
			void operator()(const Plane& aPlane)
			{
				CGAL_precondition(false);
			}

			std::vector<payload_type>& myOutResult;
		};

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
			std::vector<typename Strategy::payload_type> result;
			for (int i = 0; i < someItems.size(); ++i)
			{
				const auto intersection = CGAL::intersection(someItems[i], anotherItem);
				if (const auto* data = intersection.get_ptr())
				{
					CGAL_precondition(data != nullptr);
					Strategy visitor{ anotherItem, i, result };
					data->apply_visitor(visitor);
				}
			}
			return result;
		}

		std::vector<Plane> EnforceUpwardsOrientation(const std::vector<Plane>& somePlanes);
		std::vector<Point3> DualMapping(const std::vector<Plane>& somePlanes);
		std::vector<Plane> DualMapping(const std::vector<Point3>& somePoints);
	}
}