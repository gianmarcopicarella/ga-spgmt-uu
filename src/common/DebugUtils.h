#pragma once

#include "Types.h"
#include <fstream>


/*
  int nb_points = 10;
  int dim =5;
  double size = 100.0;
  std::cout << "Generating "<<nb_points<<" random points on the surface of "
			<<"a sphere in "<<dim<<"D of center 0 and radius "<<size<<std::endl;
  std::vector<Point> v;
  v.reserve (nb_points);
  CGAL::Random_points_on_sphere_d<Point> gen (dim, 100.0);
  for (int i = 0; i < nb_points; ++i)  v.push_back (*gen++);
  */

namespace SPGMT
{
	namespace Debug
	{
		struct SamplesSeparatedByPlane
		{
			std::vector<Point3> mySamples;
			int myPositiveCount{ 0 }, myNegativeCount{ 0 }, myOverSurfaceCount{ 0 };
		};

		std::vector<Point3> Uniform3DCubeSampling(const double anHalfSide, const int aSampleCount);
		std::vector<Plane> RandomPlaneSampling(const int aSampleCount, const double aMinPlaneHeight = -100.f, const double aMaxPlaneHeight = 100.f);
		std::vector<Plane> RandomParallelPlanesSampling(const int aSampleCount, const double aMinPlaneDistance = -100.f, const double aMaxPlaneDistance = 100.f);
		SamplesSeparatedByPlane RandomPointsPartitionedByPlane(const int aSampleCount, const Plane& aPlane,
			const double aMaxPointDistance = 50.f, const bool anAllowSamplesOverPlaneFlag = true);

		// aPathToFile must be local to CGAL_DATA_DIR
		template<typename T>
		bool SerializeToPLY(const std::string& aPathToFile, const T& aDataStructure, const bool aBinaryFlag = false, const int aPrecisionCount = 17)
		{
			const auto& params = CGAL::parameters::stream_precision(aPrecisionCount)
				.use_binary_mode(aBinaryFlag);
			return CGAL::IO::write_PLY(CGAL::data_file_path(aPathToFile), aDataStructure, params);
		}

		// aPathToFile must be local to CGAL_DATA_DIR
		template<typename T>
		T DeserializeFromPLY(const std::string& aPathToFile)
		{
			// std::ios_base::binary is mandatory on Windows if input is binary PLY
			std::ifstream file(CGAL::data_file_path(aPathToFile));
			T data;
			// Only Points and Normals are loaded by default, any other property is ignored (requires additional code to be loaded)
			const auto isReadSucceeded = CGAL::IO::read_PLY(file, data);
			file.close();
			CGAL_precondition(isReadSucceeded);
			return data;
		}
	}
}