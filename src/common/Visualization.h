#pragma once

#include <vector>

namespace SPGMT
{
	struct Vertex;
	namespace Visualization 
	{
		void VisualizeLowerEnvelope(const std::vector<Vertex>& someVertices);
	}
}