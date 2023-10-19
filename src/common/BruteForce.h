#pragma once

#include "Types.h"

namespace SPGMT
{
	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes);

	LowerEnvelope3d ParallelComputeLowerEnvelope(const std::vector<Plane>& somePlanes);

	void TriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope);

	void ParallelTriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope);
}