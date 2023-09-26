#pragma once

#include "Types.h"

namespace SPGMT
{
	/*
	Each pair of lines can intersect one time at most





		LowerEnvelope(Planes P) :
		- If P is empty then return null
		- Find all intersection points K between p_i, p_j, p_k in P | i != j != k
		- Find all intersection lines L between p_i, p_j in P | i != j
		- If K is empty and L is empty then:
		--	Return p_i in P | p_i.z is minimum in P
		- Find Bounding Box B containing all the intersection points with some padding
		- Find intersection points between l_i in L and B and set these vertices as infinity

	
	
	
	
	
	
	
	
	
	
	
	
	*/

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