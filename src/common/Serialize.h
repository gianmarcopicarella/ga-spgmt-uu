#pragma once

#include "BruteForce.h"
#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/read_ply_points.h>

#include <utility>
#include <vector>
#include <fstream>

namespace SPGMT
{
	namespace Serialization
	{
		typedef std::tuple<Point3, std::vector<int>, std::vector<int>, int> SerializableVertex;
		typedef CGAL::Nth_of_tuple_property_map<0, SerializableVertex> PointMap;
		typedef CGAL::Nth_of_tuple_property_map<1, SerializableVertex> NeighboursMap;
		typedef CGAL::Nth_of_tuple_property_map<2, SerializableVertex> LowestLeftPlanesMap;
		typedef CGAL::Nth_of_tuple_property_map<3, SerializableVertex> TypeMap;

		//typedef std::tuple<int, std::vector<int>, int> SerializableFace;
		//typedef CGAL::Nth_of_tuple_property_map<0, SerializableFace> FaceTypeMap;
		//typedef CGAL::Nth_of_tuple_property_map<1, SerializableFace> VertexIndicesMap;
		//typedef CGAL::Nth_of_tuple_property_map<2, SerializableFace> PlaneIndexMap;
		/*
		void SerializeLowerEnvelope(const std::string& aPathToFile, const std::vector<Vertex>& someVertices)
		{
			std::vector<SerializableVertex> vertices;

			for (int i = 0; i < someVertices.size(); ++i)
			{
				vertices.push_back(std::make_tuple(
					someVertices[i].myPoint,
					someVertices[i].mySortedNeighboursIndices,
					someVertices[i].myLowestLeftPlanes,
					static_cast<int>(someVertices[i].myType)));
			}

			constexpr auto precision = 17;
			constexpr auto binary = true;
			const auto& params = CGAL::parameters::stream_precision(precision)
				.use_binary_mode(binary);

			std::ofstream file{ CGAL::data_file_path(aPathToFile), std::ios::binary };
			CGAL::IO::write_PLY_with_properties(file, vertices,
				CGAL::make_ply_point_writer(PointMap()),
				std::make_pair(NeighboursMap(), CGAL::IO::PLY_property<std::vector<int>>("sorted_neighbours_indices")),
				std::make_pair(LowestLeftPlanesMap(), CGAL::IO::PLY_property<std::vector<int>>("lowest_left_planes")),
				std::make_pair(TypeMap(), CGAL::IO::PLY_property<int>("vertex_type"))
			);
		}

		std::vector<Vertex> ReadLowerEnvelope(const std::string& aPathToFile)
		{
			std::vector<Vertex> result;
			std::vector<SerializableVertex> vertices;
			std::ifstream file{ CGAL::data_file_path(aPathToFile), std::ios::binary };

			if (!CGAL::IO::read_PLY_with_properties(file, std::back_inserter(vertices),
				CGAL::make_ply_point_reader(PointMap()),
				std::make_pair(NeighboursMap(), CGAL::IO::PLY_property<std::vector<int>>("sorted_neighbours_indices")),
				std::make_pair(LowestLeftPlanesMap(), CGAL::IO::PLY_property<std::vector<int>>("lowest_left_planes")),
				std::make_pair(TypeMap(), CGAL::IO::PLY_property<int>("vertex_type"))))
			{
				CGAL_precondition(false);
			}

			for (std::size_t i = 0; i < vertices.size(); ++i)
			{
				Vertex vertex{
					static_cast<VertexType>(std::get<3>(vertices[i])),
					std::get<0>(vertices[i]),
					std::get<1>(vertices[i]),
					std::get<2>(vertices[i])
				};
				result.push_back(vertex);
			}

			return result;
		}
		*/
	}
}