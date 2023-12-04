#include "SmartLowerEnvelope.h"

#include "BatchPointLocation.h"
#include "Utils.h"

#include "DebugUtils.h"

#include <CGAL/Random.h>


namespace SPGMT
{
	/*
	namespace
	{
		std::vector<size_t> locSampleIndices(const std::vector<size_t>& someItemsIndices, const double aProbability)
		{
			CGAL_precondition(aProbability >= 0.f);
			static CGAL::Random random{};
			std::vector<size_t> samples;
			samples.reserve(static_cast<size_t>(someItemsIndices.size() * aProbability + 1));

			for (int i = 0; i < someItemsIndices.size(); ++i)
			{
				if (random.get_double() < aProbability)
				{
					samples.emplace_back(someItemsIndices[i]);
				}
			}
			return samples;
		}

		void locAdjustInfinityEdges(
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& somePlanesIndices,
			LowerEnvelope3d& anOutLowerEnvelope)
		{
			CGAL_precondition(std::holds_alternative<std::vector<Edge<Point3>>>(anOutLowerEnvelope));
			auto& edges = std::get<std::vector<Edge<Point3>>>(anOutLowerEnvelope);
			CGAL_precondition(edges.size() > 2);

			for (size_t i = 0; i < edges.size(); ++i)
			{
				if (edges[i].myType == EdgeType::HALF_EDGE_SF)
				{
					std::vector<std::tuple<FT, Point3>> intersections;

					for (int k = 0; k < somePlanesIndices.size(); ++k)
					{
						const auto& plane = somePlanes[somePlanesIndices[k]];

						const auto firstIntersection = CGAL::intersection(Line3{ edges[i].myStart, Dir3{0,0,1} }, plane);
						CGAL_precondition(firstIntersection.has_value());
						const Point3* firstProjectedPoint = boost::get<Point3>(&*firstIntersection);
						CGAL_precondition(firstProjectedPoint != nullptr);

						const auto secondIntersection = CGAL::intersection(Line3{ edges[i].myEnd, Dir3{0,0,1} }, plane);
						CGAL_precondition(secondIntersection.has_value());
						const Point3* secondProjectedPoint = boost::get<Point3>(&*secondIntersection);
						CGAL_precondition(secondProjectedPoint != nullptr);


						if (secondProjectedPoint->z() < firstProjectedPoint->z() &&
							secondProjectedPoint->z() >= edges[i].myEnd.z())
						{
							//std::cout << plane << "   ---   " << edges[i].myEnd << "   ---   " << Vec3{ edges[i].myStart, edges[i].myEnd } << std::endl;

							const auto rayPlaneIntersection = CGAL::intersection(Ray3{ edges[i].myEnd, Vec3{edges[i].myStart, edges[i].myEnd} }, plane);
							if (rayPlaneIntersection)
							{
								const Point3* point = boost::get<Point3>(&*rayPlaneIntersection);
								CGAL_precondition(point != nullptr);
								intersections.emplace_back(std::make_tuple(Vec3{ edges[i].myEnd, *point }.squared_length(), *point));
							}
						}
					}

					if (intersections.size() > 0)
					{
						// Pick farthest point and update points at infinity
						const auto maxItemIt = std::max_element(intersections.begin(), intersections.end(), [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
						CGAL_postcondition(maxItemIt != intersections.end());

						edges[i].myEnd = std::get<1>(*maxItemIt);

						// Search for opposite edge
						const auto oppEdgeIt = std::find_if(edges.begin(), edges.end(), [&](const auto& a) {
							return a.myType == EdgeType::HALF_EDGE_EF && a.myEnd == edges[i].myStart;
							});
						CGAL_precondition(oppEdgeIt != edges.end());

						oppEdgeIt->myStart = std::get<1>(*maxItemIt);
					}
				}
			}
		}

		std::tuple<size_t, size_t> locBinarySearch(
			const std::vector<std::tuple<Point3, size_t>>& somePoints,
			const Point3& aPoint)
		{
			CGAL_precondition(somePoints.size() > 0);
			int left = 0, right = somePoints.size();
			while (left < right)
			{
				const auto mid = left + (right - left) / 2;
				if (aPoint < std::get<0>(somePoints[mid]))
				{
					right = mid;
				}
				else if (aPoint > std::get<0>(somePoints[mid]))
				{
					left = mid + 1;
				}
				else
				{
					CGAL_postcondition(aPoint == std::get<0>(somePoints[mid]));
					return std::make_tuple(mid, std::get<1>(somePoints[mid]));
				}
			}

			// Item not found
			CGAL_precondition(false);
			return std::make_tuple(SIZE_MAX, SIZE_MAX);
		}

		struct MergeData
		{
			std::vector<std::tuple<Point3, size_t>> myUniqueVertices;
			std::vector<std::tuple<size_t, size_t, size_t>> myUniqueTriangles;

		};

		LowerEnvelope3d locMergeLowerEnvelopes(
			const std::vector<Plane>& somePlanes,
			const std::vector<LowerEnvelope3d>& someLowerEnvelopes,
			const MergeData& aMergeData)
		{
			std::vector<std::tuple<size_t, size_t, size_t>> prevEdges;
			prevEdges.reserve(aMergeData.myUniqueTriangles.size() * 3);

			for (auto i = 0; i < aMergeData.myUniqueTriangles.size(); ++i)
			{
				const auto& triangle = aMergeData.myUniqueTriangles[i];
				prevEdges.emplace_back(std::make_tuple(std::get<0>(triangle), std::get<1>(triangle), i));
				prevEdges.emplace_back(std::make_tuple(std::get<0>(triangle), std::get<2>(triangle), i));
				prevEdges.emplace_back(std::make_tuple(std::get<1>(triangle), std::get<2>(triangle), i));
			}

			// Sort by edges
			hpx::sort(hpx::execution::seq, prevEdges.begin(), prevEdges.end());

			// Remove boundary edges (edges enclosing the lower envelope that are defined by only one triangle)
			std::vector<size_t> removeIndices;

			if (std::get<0>(prevEdges[0]) != std::get<0>(prevEdges[1]) && std::get<1>(prevEdges[0]) != std::get<1>(prevEdges[1]))
			{
				removeIndices.emplace_back(0);
			}

			for (auto i = 1; i < prevEdges.size() - 1; ++i)
			{
				if (std::get<0>(prevEdges[i]) != std::get<0>(prevEdges[i - 1]) && std::get<1>(prevEdges[i]) != std::get<1>(prevEdges[i - 1]) &&
					std::get<0>(prevEdges[i]) != std::get<0>(prevEdges[i + 1]) && std::get<1>(prevEdges[i]) != std::get<1>(prevEdges[i + 1]))
				{
					removeIndices.emplace_back(i);
				}
			}

			if (std::get<0>(prevEdges.back()) != std::get<0>(prevEdges[prevEdges.size() - 2]) && std::get<1>(prevEdges.back()) != std::get<1>(prevEdges[prevEdges.size() - 2]))
			{
				removeIndices.emplace_back(prevEdges.size() - 1);
			}

			std::reverse(removeIndices.begin(), removeIndices.end());

			for (int i = prevEdges.size() - 1, k = 0; i > -1; --i)
			{
				if (i == removeIndices[k])
				{
					prevEdges.erase(prevEdges.begin() + i, prevEdges.begin() + 1 + i);
					k++;
				}
			}

			if (prevEdges.size() == 0)
			{
				CGAL_precondition(someLowerEnvelopes.size() == 1);
				return someLowerEnvelopes[0];
			}

			// Run merge operation

			std::vector<std::tuple<size_t, size_t, size_t>> nextEdges;

			std::vector<bool> isTriangleMerged(aMergeData.myUniqueTriangles.size(), false);
			
			std::vector<LowerEnvelope3d> prevLowerEnvelopes(someLowerEnvelopes);
			std::vector<LowerEnvelope3d> nextLowerEnvelopes;

			bool areThereMultipleTriangles;

			do
			{
				for (int i = 0; i < prevEdges.size() - 1; i += 2)
				{
					const auto firstTriIdx = std::get<2>(prevEdges[i]);
					const auto secondTriIdx = std::get<2>(prevEdges[i + 1]);
					if (isTriangleMerged[firstTriIdx] || isTriangleMerged[secondTriIdx])
					{
						continue;
					}
					isTriangleMerged[firstTriIdx] = true;
					isTriangleMerged[secondTriIdx] = true;

					using LowerEnvelopeEdges = std::vector<Edge<Point3>>;

					if (std::holds_alternative<std::monostate>(prevLowerEnvelopes[firstTriIdx]) &&
						std::holds_alternative<std::monostate>(prevLowerEnvelopes[secondTriIdx]))
					{
						nextLowerEnvelopes.emplace_back(LowerEnvelope3d{});
					}
					else if (std::holds_alternative<std::monostate>(prevLowerEnvelopes[firstTriIdx]))
					{
						nextLowerEnvelopes.emplace_back(prevLowerEnvelopes[secondTriIdx]);
					}
					else if (std::holds_alternative<std::monostate>(prevLowerEnvelopes[secondTriIdx]))
					{
						nextLowerEnvelopes.emplace_back(prevLowerEnvelopes[firstTriIdx]);
					}
					else if (std::holds_alternative<size_t>(prevLowerEnvelopes[firstTriIdx]) &&
						std::holds_alternative<size_t>(prevLowerEnvelopes[secondTriIdx]))
					{
						std::vector<Plane> planes(2);

						const auto firstPlaneIdx = std::get<size_t>(prevLowerEnvelopes[firstTriIdx]);
						const auto secondPlaneIdx = std::get<size_t>(prevLowerEnvelopes[secondTriIdx]);

						planes[0] = somePlanes[firstPlaneIdx];
						planes[1] = somePlanes[secondPlaneIdx];

						const auto& newEnvelope = ComputeLowerEnvelope<ExecutionPolicy::SEQ>(planes);
						nextLowerEnvelopes.emplace_back(newEnvelope);
					}
					else
					{
						// Both lower envelopes have edges
						// Merge lower envelopes along edges defined by Vertices[std::get<0>(edge)] and Vertices[std::get<1>(edge)]
						auto& le1 = std::get<LowerEnvelopeEdges>(prevLowerEnvelopes[firstTriIdx]);
						auto& le2 = std::get<LowerEnvelopeEdges>(prevLowerEnvelopes[secondTriIdx]);

						const auto v1Idx = std::get<0>(prevEdges[i]);
						const auto v2Idx = std::get<1>(prevEdges[i]);

						const Segment3 sharedBorder{ std::get<0>(aMergeData.myUniqueVertices[v1Idx]), std::get<0>(aMergeData.myUniqueVertices[v2Idx]) };

						std::vector<std::tuple<Point3, size_t, size_t>> firstEnvelopeInter, secondEnvelopeInter;

						for (auto i = 0; i < le1.size(); ++i)
						{
							if (le1[i].myType == EdgeType::HALF_EDGE_SF)
							{
								const auto intersection = CGAL::intersection(Ray3{ le1[i].myStart,le1[i].myEnd }, sharedBorder);
								if (intersection)
								{
									const Point3* point = boost::get<Point3>(&*intersection);
									CGAL_precondition(point != nullptr);

									const auto oppIter = std::find_if(le1.begin(), le1.end(),
										[&](const auto& anEdge) { return anEdge.myEnd == le1[i].myStart && anEdge.myType == EdgeType::HALF_EDGE_EF; });

									CGAL_precondition(oppIter != le1.end());

									firstEnvelopeInter.emplace_back(std::make_tuple(*point, i, std::distance(le1.begin(), oppIter)));
								}

							}
						}

						for (auto i = 0; i < le2.size(); ++i)
						{
							if (le2[i].myType == EdgeType::HALF_EDGE_SF)
							{
								const auto intersection = CGAL::intersection(Ray3{ le2[i].myStart,le2[i].myEnd }, sharedBorder);
								if (intersection)
								{
									const Point3* point = boost::get<Point3>(&*intersection);
									CGAL_precondition(point != nullptr);

									const auto oppIter = std::find_if(le2.begin(), le2.end(),
										[&](const auto& anEdge) { return anEdge.myEnd == le2[i].myStart && anEdge.myType == EdgeType::HALF_EDGE_EF; });

									CGAL_precondition(oppIter != le2.end());

									secondEnvelopeInter.emplace_back(std::make_tuple(*point, i, std::distance(le2.begin(), oppIter)));
								}

							}
						}

						std::sort(firstEnvelopeInter.begin(), firstEnvelopeInter.end());
						std::sort(secondEnvelopeInter.begin(), secondEnvelopeInter.end());

						size_t i = 0, k = 0;
						std::vector<size_t> edgesToErase;

						while (i < firstEnvelopeInter.size() && k < secondEnvelopeInter.size())
						{
							while (std::get<0>(firstEnvelopeInter[i]) < std::get<0>(secondEnvelopeInter[k]))
							{
								++i;
							}

							if (std::get<0>(firstEnvelopeInter[i]) == std::get<0>(secondEnvelopeInter[k]))
							{
								le1[std::get<1>(firstEnvelopeInter[i])].myType = EdgeType::SEGMENT;
								le1[std::get<1>(firstEnvelopeInter[i])].myEnd = le2[std::get<1>(secondEnvelopeInter[k])].myStart;

								le1[std::get<2>(firstEnvelopeInter[i])].myType = EdgeType::SEGMENT;
								le1[std::get<2>(firstEnvelopeInter[i])].myEnd = le2[std::get<1>(secondEnvelopeInter[k])].myStart;

								edgesToErase.emplace_back(std::get<1>(secondEnvelopeInter[k]));
								edgesToErase.emplace_back(std::get<2>(secondEnvelopeInter[k]));
							}

							++i;
							++k;

							while (std::get<0>(firstEnvelopeInter[k]) < std::get<0>(secondEnvelopeInter[i]))
							{
								++k;
							}
						}

						std::sort(edgesToErase.begin(), edgesToErase.end());
						std::reverse(edgesToErase.begin(), edgesToErase.end());

						for (int i = le2.size() - 1, k = 0; i > -1; --i)
						{
							if (i == edgesToErase[k])
							{
								le2.erase(le2.begin() + i, le2.begin() + 1 + i);
								k++;
							}
						}

						std::copy(le2.begin(), le2.end(), std::back_inserter(le1));
						// Sort edges

						nextLowerEnvelopes.emplace_back(le1);

						const auto newEnvelopeIndex = nextLowerEnvelopes.size() - 1;

						// Find all edges for le1 and le2, remove the merged one and store the result  with idx=newEnvelopeIndex
						std::vector<std::tuple<size_t, size_t>> le1Entries, le2Entries;

						for (int i = 0; i < prevEdges.size(); ++i)
						{
							if (std::get<2>(prevEdges[i]) == v1Idx)
							{
								le1Entries.emplace_back(std::make_tuple(std::get<0>(prevEdges[i]), std::get<1>(prevEdges[i])));
							}
							else if (std::get<2>(prevEdges[i]) == v2Idx)
							{
								le2Entries.emplace_back(std::make_tuple(std::get<0>(prevEdges[i]), std::get<1>(prevEdges[i])));
							}
						}
						
						std::vector<std::tuple<size_t, size_t>> res;
						std::set_symmetric_difference(le1Entries.begin(), le1Entries.end(), le2Entries.begin(), le2Entries.end(), std::back_inserter(res));
						// Create new edges array for next iteration using newEnvelopeIndex as new "triangle" index during merge
						for(const auto& r : res) nextEdges.emplace_back(std::make_tuple(std::get<0>(r), std::get<1>(r), newEnvelopeIndex));
					}
				}

				// Repeat process until the number of different triangles in edges is 1
				areThereMultipleTriangles = std::all_of(nextEdges.begin(), nextEdges.end(),
					[&](const auto& anEdge) { return std::get<2>(anEdge) == std::get<2>(*nextEdges.begin()); });

				prevEdges = nextEdges;
				nextEdges.clear();

				prevLowerEnvelopes = nextLowerEnvelopes;
				nextLowerEnvelopes.clear();

			} while (areThereMultipleTriangles);


			return nextLowerEnvelopes[0];
		}

		std::vector<size_t> locGetSortedPrismConflictList(
			const std::tuple<size_t, size_t, size_t>& aTriangle,
			const std::vector<std::vector<size_t>>& someConflictLists)
		{
			const auto& firstList = someConflictLists[std::get<0>(aTriangle)];
			const auto& secondList = someConflictLists[std::get<1>(aTriangle)];
			const auto& thirdList = someConflictLists[std::get<2>(aTriangle)];

			std::vector<size_t> tempUnion, prismConflictList;

			std::set_union(firstList.begin(), firstList.end(), secondList.begin(), secondList.end(), std::back_inserter(tempUnion));
			std::set_union(tempUnion.begin(), tempUnion.end(), thirdList.begin(), thirdList.end(), std::back_inserter(prismConflictList));

			CGAL_postcondition(std::is_sorted(prismConflictList.begin(), prismConflictList.end()));
			return prismConflictList;
		}

		void locAddTriangle(
			const Point3& aFirstVertex,
			const Point3& aSecondVertex,
			const Point3& aThirdVertex,
			MergeData& anOutMergeData)
		{
			// Collect triangle indices
			const auto firstIndices = locBinarySearch(anOutMergeData.myUniqueVertices, aFirstVertex);
			const auto secondIndices = locBinarySearch(anOutMergeData.myUniqueVertices, aSecondVertex);
			const auto thirdIndices = locBinarySearch(anOutMergeData.myUniqueVertices, aThirdVertex);

			std::vector<size_t> indices(3);
			indices[0] = std::get<0>(firstIndices);
			indices[1] = std::get<0>(secondIndices);
			indices[2] = std::get<0>(thirdIndices);
			hpx::sort(hpx::execution::seq, indices.begin(), indices.end());

			CGAL_postcondition(indices[0] != indices[1] && indices[0] != indices[2] && indices[1] != indices[2]);

			anOutMergeData.myUniqueTriangles.emplace_back(std::make_tuple(indices[0], indices[1], indices[2]));
		}

		void locComputeMergeData(const LowerEnvelope3d& aLowerEnvelope, MergeData& anOutMergeData)
		{
			using EdgesList = std::vector<Edge<Point3>>;
			CGAL_precondition(std::holds_alternative<EdgesList>(aLowerEnvelope));
			const auto& edges = std::get<std::vector<Edge<Point3>>>(aLowerEnvelope);
			CGAL_precondition(edges.size() > 0);

			//Debug::PrintLowerEnvelope(aLowerEnvelope);
			// Collect and sort unique vertices
			size_t i = 0;
			while (i < edges.size())
			{
				anOutMergeData.myUniqueVertices.emplace_back(std::make_tuple(edges[i].myStart, i));
				while (++i < edges.size() && std::get<0>(anOutMergeData.myUniqueVertices.back()) == edges[i].myStart);
			}

			// Collect unique triangles
			i = 0;
			while (i < edges.size())
			{
				const size_t lastStartIdx = i;
				while (++i < edges.size() - 1 && edges[lastStartIdx].myStart == edges[i].myStart)
				{
					locAddTriangle(edges[i - 1].myStart, edges[i - 1].myEnd, edges[i].myEnd, anOutMergeData);
				}

				// Check if there is an edge from point at lastStartIdx and i-1
				// If yes, then add the additional triangle
				const auto lastStartIndices = locBinarySearch(anOutMergeData.myUniqueVertices, edges[lastStartIdx].myEnd);
				size_t k = std::get<1>(lastStartIndices);
				while (k < edges.size() && edges[k].myStart == edges[lastStartIdx].myEnd)
				{
					if (edges[k].myEnd == edges[i - 1].myEnd)
					{
						locAddTriangle(edges[k].myStart, edges[k].myEnd, edges[lastStartIdx].myStart, anOutMergeData);
						break;
					}
					k++;
				}
			}

			hpx::sort(hpx::execution::seq, anOutMergeData.myUniqueTriangles.begin(), anOutMergeData.myUniqueTriangles.end());
			anOutMergeData.myUniqueTriangles.erase(
				hpx::unique(hpx::execution::seq, anOutMergeData.myUniqueTriangles.begin(), anOutMergeData.myUniqueTriangles.end()), anOutMergeData.myUniqueTriangles.end());
		}

		std::vector<Plane> locGetPlanesSubset(
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& somePlanesIndices)
		{
			std::vector<Plane> subsetPlanes(somePlanesIndices.size());
			hpx::transform(hpx::execution::seq, somePlanesIndices.begin(), somePlanesIndices.end(), subsetPlanes.begin(),
				[&](size_t aPlaneIdx) { return somePlanes[aPlaneIdx]; });
			return subsetPlanes;
		}

		std::vector<size_t> locGetPlanesIntersectingLE(
			const LowerEnvelope3d& aLowerEnvelope,
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& somePlanesIndices)
		{
			const auto& planesIndicesLE = GetLowerEnvelopePlanesIndices<ExecutionPolicy::SEQ>(aLowerEnvelope);

			std::vector<size_t> result;

			for (auto i = 0; i < somePlanesIndices.size(); ++i)
			{
				for (auto k = 0; k < planesIndicesLE.size(); ++k)
				{
					const auto& planeInLEDir = somePlanes[planesIndicesLE[k]].orthogonal_direction();
					const auto& planeTestedDir = somePlanes[somePlanesIndices[i]].orthogonal_direction();
					CGAL_precondition((-planeInLEDir) == planeTestedDir);

					if (planeInLEDir != planeTestedDir ||
						somePlanes[somePlanesIndices[i]].d() < somePlanes[planesIndicesLE[k]].d())
					{
						result.emplace_back(somePlanesIndices[i]);
						break;
					}
				}
			}

			return result;
		}


		LowerEnvelope3d locRecursiveLowerEnvelope(
			const std::vector<Plane>& somePlanes,
			const std::vector<size_t>& someAllPlanesIndices,
			const std::vector<size_t>& somePlanesIndices,
			bool aFirstEnvelopeWasComputed = false)
		{
			constexpr auto threshold = 10;
			const auto& planesSubset = locGetPlanesSubset(somePlanes, somePlanesIndices);

			if (somePlanesIndices.size() <= threshold)
			{
				std::vector<Plane> planesUpwards(planesSubset);
				Utils::FlipPlaneNormalsIfFacingDownwards(planesUpwards);

				const auto lowerEnvelope = ComputeLowerEnvelope<ExecutionPolicy::SEQ>(planesUpwards);
				CGAL_postcondition(!std::holds_alternative<std::monostate>(lowerEnvelope));
				return lowerEnvelope;
			}

			// Find sample of planes S from ALL and compute ALL - S
			CGAL_precondition(std::is_sorted(somePlanesIndices.begin(), somePlanesIndices.end()));

			std::vector<size_t> remainingIndices;
			std::set_difference(
				someAllPlanesIndices.begin(), someAllPlanesIndices.end(),
				somePlanesIndices.begin(), somePlanesIndices.end(),
				std::back_inserter(remainingIndices));

			LowerEnvelope3d lowerEnvelope;

			{
				std::vector<Plane> planesUpwards(planesSubset);
				Utils::FlipPlaneNormalsIfFacingDownwards(planesUpwards);

				lowerEnvelope = ComputeLowerEnvelope<ExecutionPolicy::SEQ>(planesUpwards);
				CGAL_postcondition(!std::holds_alternative<std::monostate>(lowerEnvelope));
				//CGAL_postcondition(Debug::IsLowerEnvelopeCorrect(lowerEnvelope, planesSubset));

				// What happens if the lower envelope is just a line with two faces or simply a plane ?
				// Compute intersecting planes in subspace and return recursive call
				if (std::holds_alternative<size_t>(lowerEnvelope) ||
					std::get<std::vector<Edge<Point3>>>(lowerEnvelope).size() == 2)
				{
					const auto recursivePlanesIndices = locGetPlanesIntersectingLE(lowerEnvelope, somePlanes, remainingIndices);
					return locRecursiveLowerEnvelope(somePlanes, someAllPlanesIndices, recursivePlanesIndices, aFirstEnvelopeWasComputed);
				}
				else
				{
					// Adjust edges at infinity needed only for the first LE computed in the main calling function
					if (!aFirstEnvelopeWasComputed)
					{
						locAdjustInfinityEdges(somePlanes, remainingIndices, lowerEnvelope);
						aFirstEnvelopeWasComputed = true;
					}

					TriangulateLowerEnvelope<ExecutionPolicy::SEQ>(lowerEnvelope);
				}
			}

			// Compute unique vertices and prisms
			MergeData data;
			locComputeMergeData(lowerEnvelope, data);

			std::vector<std::vector<size_t>> planesBelowVertex(data.myUniqueVertices.size());

			{
				// Collect unique vertices
				std::vector<Point3> uniqueVertices;
				std::transform(data.myUniqueVertices.begin(), data.myUniqueVertices.end(), std::back_inserter(uniqueVertices),
					[](const auto& aVertexPair) { return std::get<0>(aVertexPair); });

				const auto& planesMinusSubset = locGetPlanesSubset(somePlanes, remainingIndices);

				// Apply duality transform
				const auto& dualizedPlanes = Utils::DualMapping(planesMinusSubset);
				auto dualizedPoints = Utils::DualMapping(uniqueVertices);
				// TODO: force normals upwards
				Utils::FlipPlaneNormalsIfFacingDownwards(dualizedPoints);

				const auto& batchPointResult = BatchPointLocation<ExecutionPolicy::SEQ>(dualizedPoints, dualizedPlanes);
				CGAL_postcondition(batchPointResult.myRangeWrappers.size() == planesMinusSubset.size());

				for (auto i = 0; i < batchPointResult.myRangeWrappers.size(); ++i)
				{
					const auto& zoneRange = batchPointResult.myRangeWrappers[i];
					const auto& sortedPlanesIndices =
						batchPointResult.mySortedPlanesCache[zoneRange.myCacheIndex].at(zoneRange.myIndex);

					for (int k = zoneRange.myRange.first; k < zoneRange.myRange.second; ++k)
					{
						// Store Vertex-PlaneIdxBelowIt
						const auto vertexIdx = sortedPlanesIndices[k];

						planesBelowVertex[vertexIdx].emplace_back(remainingIndices[i]);
					}
				}

				// Sort planes indices in increasing order
				for (auto& planes : planesBelowVertex)
				{
					hpx::sort(hpx::execution::seq, planes.begin(), planes.end());
				}
			}

			std::vector<LowerEnvelope3d> lowerEnvelopes;

			for (int i = 0; i < data.myUniqueTriangles.size(); ++i)
			{
				const auto& conflictList = locGetSortedPrismConflictList(data.myUniqueTriangles[i], planesBelowVertex);

				const auto& recLowerEnvelope = locRecursiveLowerEnvelope(somePlanes, somePlanesIndices, conflictList, aFirstEnvelopeWasComputed);
				lowerEnvelopes.emplace_back(recLowerEnvelope);
			}

			return locMergeLowerEnvelopes(somePlanes, lowerEnvelopes, data);
		}
	}*/

	LowerEnvelope3d ComputeLowerEnvelopeSmart(const std::vector<Plane>& somePlanes)
	{
		// If the input size is within a maximum threshold 
		// then use the brute force algorithm
		constexpr auto threshold = 20;
		constexpr auto exponent = 7.f / 8.f;

		return LowerEnvelope3d{};

		//if (somePlanes.size() <= threshold)
		//{
		//	return ComputeLowerEnvelope<ExecutionPolicy::SEQ>(somePlanes);
		//}

		//// Create the set of all plane indices 
		//std::vector<size_t> allIndices(somePlanes.size());
		//std::iota(allIndices.begin(), allIndices.end(), 0);

		//// Compute sampling probability
		//const auto coefficient = std::pow(somePlanes.size(), exponent);
		//const auto probability = coefficient / somePlanes.size();

		//const auto& sampledPlaneIndices = locSampleIndices(allIndices, probability);

		//return locRecursiveLowerEnvelope(somePlanes, allIndices, sampledPlaneIndices, false);
	}
}