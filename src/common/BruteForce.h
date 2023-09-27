#pragma once

#include "Types.h"

namespace SPGMT
{
	enum class VertexType
	{
		FINITE,
		INFINITE
	};

	struct Vertex
	{
		VertexType myType { VertexType::FINITE };
		Point3 myPoint;
		std::vector<int> mySortedNeighboursIndices;
		std::vector<int> myLowestLeftPlanes;
	};

	std::vector<Vertex> ComputeLowerEnvelope(const std::vector<Plane>& somePlanes);
}