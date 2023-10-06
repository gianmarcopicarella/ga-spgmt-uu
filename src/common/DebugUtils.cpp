#include "DebugUtils.h"

#include <CGAL/point_generators_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Random.h>
#include <algorithm>
#include <random>

namespace SPGMT
{
	/*const unsigned long int SEED = 23322392398;
	const unsigned long int SEED2 = 22392398232;*/
	CGAL::Random& GetDefaultRandom()
	{
		static CGAL::Random rand{ 3429112288/*3222535971*//*1306513302*//*4169948633*/ };
		return rand;
	}

	namespace
	{
		struct Visitor
		{
			typedef void result_type;
			void operator()(const Line3& aLine)
			{
				myOutResult.push_back(aLine);
			}
			void operator()(const Plane&) {}

			std::vector<Line3>& myOutResult;
		};

		struct VisitorPoint
		{
			typedef void result_type;
			void operator()(const Point3& aLine)
			{
				myOutResult.push_back(aLine);
			}
			void operator()(const Plane&) {}
			void operator()(const Line3&) {}

			std::vector<Point3>& myOutResult;
		};

		std::vector<Line3> locGetLinesIntersections(const std::vector<Plane>& somePlanes)
		{
			std::vector<Line3> result;
			for (int i = 0; i < somePlanes.size(); ++i)
			{
				for (int k = i + 1; k < somePlanes.size(); ++k)
				{
					const auto intersection = CGAL::intersection(somePlanes[i], somePlanes[k]);
					if (const auto* data = intersection.get_ptr())
					{
						Visitor visitor{ result };
						data->apply_visitor(visitor);
					}
				}
			}
			return result;
		}

		std::vector<Point3> locGetPointsIntersections(const std::vector<Plane>& somePlanes)
		{
			std::vector<Point3> result;
			for (int i = 0; i < somePlanes.size(); ++i)
			{
				for (int k = i + 1; k < somePlanes.size(); ++k)
				{
					for (int j = k + 1; j < somePlanes.size(); ++j)
					{
						const auto intersection = CGAL::intersection(somePlanes[i], somePlanes[k], somePlanes[j]);
						if (const auto* data = intersection.get_ptr())
						{
							VisitorPoint visitor{ result };
							data->apply_visitor(visitor);
						}
					}
				}
			}
			return result;
		}

		// Trivial and inefficient implementation. 
		// The numeric data structure doesn't seem to behave correctly when using std::swap
		template<typename T>
		std::vector<T> locShuffle(const std::vector<T>& someItems)
		{
			CGAL::Random rand{ GetDefaultRandom() };
			std::vector<int> indices(someItems.size());
			std::iota(indices.begin(), indices.end(), 0);
			for (int i = 0; i < someItems.size(); ++i)
			{
				std::swap(indices[i], indices[rand.get_int(i, someItems.size())]);
			}
			std::vector<T> items(someItems.size());
			for (int i = 0; i < someItems.size(); ++i)
			{
				items[indices[i]] = someItems[i];
			}
			return items;
		}
	}

	namespace Debug
	{
		Envelope_diagram_2 GetLowerEnvelopeOfPlanes(const std::vector<Plane>& somePlanes)
		{
			Envelope_diagram_2      min_diag;
			CGAL::lower_envelope_3(somePlanes.begin(), somePlanes.end(), min_diag);
			return min_diag;
		}

		std::vector<Point3> SampleTriplePlaneIntersectionPoints(const std::vector<Plane>& somePlanes, const int aSampleCount)
		{
			std::vector<Point3> samples;
			const auto intersectionPoints = locGetPointsIntersections(somePlanes);

			if (intersectionPoints.empty())
			{
				return samples;
			}

			CGAL::Random random{ GetDefaultRandom() };
			while (samples.size() < aSampleCount)
			{
				samples.push_back(intersectionPoints[random.get_int(0, intersectionPoints.size())]);
			}
			return samples;
		}

		std::vector<Point3> SamplePointsAlongPlaneIntersections(const std::vector<Plane>& somePlanes, const int aSampleCount)
		{
			std::vector<Point3> samples;
			const auto intersectionLines = locGetLinesIntersections(somePlanes);

			if (intersectionLines.empty())
			{
				return samples;
			}

			CGAL::Random random{ GetDefaultRandom() };
			while (samples.size() < aSampleCount)
			{
				const auto& line = intersectionLines[random.get_int(0, intersectionLines.size())];
				const auto alongLinePoint = line.point() + line.to_vector() * random.get_double(-100.f, 100.f);
				const auto alongVerticalPoint = alongLinePoint + Vec3{ 0,0,1 } *random.get_double(-100.f, 100.f);

				samples.push_back(alongLinePoint);
				if (samples.size() < aSampleCount)
				{
					samples.push_back(alongVerticalPoint);
				}
			}
			return samples;
		}

		std::vector<Point3> Uniform3DCubeSampling(const double anHalfSide, const int aSampleCount)
		{
			std::vector<Point3> samples;
			samples.reserve(aSampleCount);
			CGAL::Random_points_in_cube_3<Point3> generator{ anHalfSide, GetDefaultRandom() };
			std::copy_n(generator, aSampleCount, std::back_inserter(samples));
			// (NOT MANDATORY) Use a random permutation to hide the creation history of the point set.
			return samples;
		}

		SamplesSeparatedByPlane RandomPointsPartitionedByPlane(const int aSampleCount, const Plane& aPlane,
			const double aMaxPointDistance, const bool anAllowSamplesOverPlaneFlag)
		{
			constexpr auto squareSide = 50.f;

			CGAL::Random_points_on_square_2<Vec2> generator{ squareSide, GetDefaultRandom() };
			CGAL::Random random{ GetDefaultRandom() };

			SamplesSeparatedByPlane result;
			auto samplesCounter{ -1 };
			const auto upperSamplesLimit = anAllowSamplesOverPlaneFlag ? static_cast<int>(std::ceil(aSampleCount * 0.8f)) : aSampleCount;
			const auto refPlanePoint = aPlane.point();

			while (++samplesCounter < upperSamplesLimit)
			{
				const auto sign = samplesCounter % 2 == 0 ? 1 : -1;
				const auto distCoeff = sign * random.uniform_real(0.0, aMaxPointDistance);

				const auto sampleCoeffs = *generator++;
				const auto samplePoint = refPlanePoint +
					aPlane.base1() * sampleCoeffs.x() +
					aPlane.base2() * sampleCoeffs.y() +
					aPlane.orthogonal_vector() * distCoeff;

				if (distCoeff > 0.f)
				{
					result.myPositiveCount += 1;
				}
				else if (distCoeff < 0.f)
				{
					result.myNegativeCount += 1;
				}
				else
				{
					result.myOverSurfaceCount += 1;
				}

				result.mySamples.push_back(samplePoint);
			}

			const auto remainingSamplesToFind = aSampleCount - upperSamplesLimit;
			if (remainingSamplesToFind > 0)
			{
				for (samplesCounter = 0; samplesCounter < remainingSamplesToFind; ++samplesCounter, ++generator)
				{
					const auto& sampleCoeffs = *generator;
					const auto samplePoint = refPlanePoint + aPlane.base1() * sampleCoeffs.x() + aPlane.base2() * sampleCoeffs.y();

					// Points must be on the surface. Seems to work ONLY with exact_kernel, so an assert is expected with floating point arithmetic
					CGAL_precondition(aPlane.has_on(samplePoint));

					result.mySamples.push_back(samplePoint);
					result.myOverSurfaceCount += 1;
				}
			}

			CGAL_postcondition((
				result.myOverSurfaceCount +
				result.myPositiveCount +
				result.myNegativeCount) == aSampleCount);

			return result;
		}

		std::vector<Plane> RandomParallelPlanesSampling(const int aSampleCount, const double aMinPlaneDistance, const double aMaxPlaneDistance)
		{
			std::vector<Plane> parallelPlanes = RandomPlaneSampling(1);
			CGAL::Random random{ GetDefaultRandom() };

			CGAL_precondition(aMinPlaneDistance > 0.f);
			CGAL_precondition(aMinPlaneDistance < aMaxPlaneDistance);

			while (parallelPlanes.size() < aSampleCount)
			{
				const Point3 planePoint = parallelPlanes.back().point() + parallelPlanes.back().orthogonal_vector() *
					random.uniform_real(aMinPlaneDistance, aMaxPlaneDistance);

				const Plane plane{ planePoint, parallelPlanes.back().orthogonal_vector() };

				parallelPlanes.push_back(plane);
			}

			parallelPlanes = locShuffle(parallelPlanes);
			// My intuition is that random_shuffle internally uses std::swap and the Plane structure or 
			// other inner data structure into it are not swapped properly, thus causing a crash after the swap has happened.
			// The only working way seems to copy all the data at different locations as in "locShuffle"
			//CGAL::cpp98::random_shuffle(parallelPlanes.begin(), parallelPlanes.end(), GetDefaultRandom());
			return parallelPlanes;
		}


		std::vector<Plane> RandomPlaneSampling(const int aSampleCount, const double aMinPlaneHeight, const double aMaxPlaneHeight)
		{
			const Vec3 positiveVerticalAxis{ 0,0,1 };
			constexpr auto radius = 1.f;

			// Just to be sure there is room for a bit of plane height variety
			// TODO: I will need to guarantee a minimum distance between parallel planes (if i use an inaccurate numeric representation)
			CGAL_precondition((aMaxPlaneHeight - aMinPlaneHeight) >= aSampleCount);

			CGAL::Random_points_on_sphere_3<Vec3> generator{ radius, GetDefaultRandom() };
			CGAL::Random random{ GetDefaultRandom() };
			std::vector<Plane> samples;

			samples.reserve(aSampleCount);

			while (samples.size() < aSampleCount)
			{
				for (int i = 0; i < aSampleCount && samples.size() < aSampleCount; ++i, ++generator)
				{
					auto planeNormal = *generator;
					const auto dot = CGAL::scalar_product(positiveVerticalAxis, planeNormal);

					// Allow vertical planes
					/*if (dot != 0.f)
					{*/
					if (dot < 0.f)
					{
						planeNormal = Vec3{ planeNormal.x(),planeNormal.y(), CGAL::abs(planeNormal.z()) };
					}

					const auto sampleAlongNormal = planeNormal * random.get_double(aMinPlaneHeight, aMaxPlaneHeight);
					const Point3 planePoint{ sampleAlongNormal.x(), sampleAlongNormal.y(), sampleAlongNormal.z() };
					const Plane plane{ planePoint, planeNormal };

					const auto isOrientationUnique = std::find(samples.begin(), samples.end(), plane) == samples.end();

					if (isOrientationUnique)
					{
						samples.push_back(plane);
					}
					//}
				}
			}

			return samples;
		}
	}
}