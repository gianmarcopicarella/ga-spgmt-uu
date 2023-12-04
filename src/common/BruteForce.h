#pragma once

#include "Types.h"

namespace SPGMT
{
	LowerEnvelope3d ComputeLowerEnvelope(const std::vector<Plane>& somePlanes);
	void TriangulateLowerEnvelope(LowerEnvelope3d& anOutLowerEnvelope);
	size_t CountVerticesInLowerEnvelope(const LowerEnvelope3d& aLowerEnvelope);
	std::vector<size_t> GetLowerEnvelopePlanesIndices(const LowerEnvelope3d& aLowerEnvelope);
}