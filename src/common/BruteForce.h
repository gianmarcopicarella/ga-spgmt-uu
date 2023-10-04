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

	enum class FaceType
	{
		BOUNDED,
		UNBOUNDED,
		UNBOUNDED_ONE_EDGE
	};

	struct Face
	{
		FaceType myType{ FaceType::BOUNDED };
		std::vector<int> myVertexIndices;
		int myPlaneIndex{ -1 };
	};

	std::vector<Vertex> ComputeLowerEnvelope(const std::vector<Plane>& somePlanes);
	std::vector<Face> ExtractLowerEnvelopeFaces(const std::vector<Vertex>& someVertices);
	std::vector<int> TriangulateLowerEnvelopeFaces(const std::vector<Vertex>& someVertices, const std::vector<Face>& someFaces);
}