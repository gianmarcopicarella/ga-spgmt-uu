#pragma once

#include <vector>
#include "Types.h"
#include "BruteForce.h"

namespace SPGMT
{
	
	LowerEnvelope3d ComputeLowerEnvelopeSmart(const std::vector<Plane>& somePlanes);

}