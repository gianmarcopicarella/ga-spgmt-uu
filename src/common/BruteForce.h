#pragma once

#include "Types.h"

namespace SPGMT
{
	template<ExecutionPolicy E>
	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes);

	template<ExecutionPolicy E>
	void TriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope);

	template<ExecutionPolicy E>
	size_t CountVerticesInLowerEnvelope(const LowerEnvelope3d& aLowerEnvelope);
}