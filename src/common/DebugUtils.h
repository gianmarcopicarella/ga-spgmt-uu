#pragma once

#include "Types.h"
#include <fstream>

namespace SPGMT
{
	namespace Debug
	{
		PointSet3 Uniform3DCubeSampling(const double anHalfSide, const int aSampleCount);

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
			std::ifstream file (CGAL::data_file_path(aPathToFile));
			T data;
			// Only Points and Normals are loaded by default, any other property is ignored (requires additional code to be loaded)
			const auto isReadSucceeded = CGAL::IO::read_PLY(file, data);
			file.close();
			CGAL_precondition(isReadSucceeded);
			return data;
		}
	}
}