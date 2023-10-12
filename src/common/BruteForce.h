#pragma once

#include "Types.h"

namespace SPGMT
{
	enum class EdgeType
	{
		LINE,
		HALF_EDGE_EF,
		HALF_EDGE_SF,
		SEGMENT
	};

	template<typename T>
	struct Edge
	{
		T myStart, myEnd;
		EdgeType myType{ EdgeType::SEGMENT };
		int myLowestLeftPlane{ -1 };
	};

	using LowerEnvelope3d = std::variant<std::monostate, int, std::vector<Edge<Point3>>>;
	using Triangles3d = std::variant<std::monostate, std::vector<Point3>>;

	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes);
	Triangles3d TriangulateLowerEnvelope(const LowerEnvelope3d& aLowerEnvelope);
}