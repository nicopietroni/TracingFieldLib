#ifndef DUAL_PATCHES
#define DUAL_PATCHES

#include <vector>
#include <UV_param.h>
#include <UV_utils.h>
#include "mesh_functions.h"
#include "mesh_refine.h"
#include "patch_managing.h"
#include "pos3.h"
#include "polyline2.h"
#include "patch_adjacency.h"
#include <Space_Query/segment_grid_2D.h>
#include <chrono>

#define MIN_ALPHA 0.1
namespace Geo
{

    template <class ScalarType>
    class DualPatch
    {
        // static bool ExtractDualOLD(const std::vector<Geo::Point3<ScalarType>> &VertPos,
        //                         const std::vector<std::vector<int>> &Faces,
        //                         PatchManaging<ScalarType> &PM,
        //                         DualParameter &DualParam,
        //                         std::vector<Geo::Point3<ScalarType>> &DualVertPos,
        //                         std::vector<std::vector<int>> &DualFaces)
        // {
        //     DualParam.SidesMidP.clear();

        //     // first get all corners
        //     std::set<int> CornerSet;
        //     for (int i = 0; i < PM.NumPatches(); i++)
        //     {
        //         std::vector<int> CornerVert;
        //         std::vector<CornerType> CurrType;
        //         PM.GetAllCorners(i, CornerVert, CurrType);
        //         PM.MapToOriginal(i, CornerVert);
        //         CornerSet.insert(CornerVert.begin(), CornerVert.end());
        //     }

        //     // compute mid point of edge
        //     std::set<std::pair<int, int>> AddedMidP;

        //     // std::set<EdgeInt> InterPEdge;

        //     std::map<std::pair<int, int>, ScalarType> EdgeAlpha;

        //     assert(PM.PData.SubPatchFaces.size() == PM.PData.SubPatchVertPos.size());
        //     for (size_t i = 0; i < PM.PData.SubPatchFaces.size(); i++)
        //     {
        //         // avoid empty patches
        //         if (PM.PData.SubPatchFaces[i].size() == 0)
        //             continue;

        //         std::vector<std::pair<int, int>> CornersExtremes;
        //         // std::vector<EdgeInt> CurrMidP;

        //         std::vector<std::pair<int, int>> EdgeInt;
        //         std::vector<ScalarType> Alpha;
        //         bool foundMidP = PM.GetSideMidpoints(i, CornerSet, CornersExtremes, EdgeInt, Alpha);
        //         if (!foundMidP)
        //             return false;

        //         assert(CornersExtremes.size() == EdgeInt.size());
        //         assert(CornersExtremes.size() == Alpha.size());

        //         for (size_t j = 0; j < CornersExtremes.size(); j++)
        //         {
        //             int IndexC0 = CornersExtremes[j].first;
        //             int IndexC1 = CornersExtremes[j].second;
        //             std::pair<int, int> Entry(std::min(IndexC0, IndexC1), std::max(IndexC0, IndexC1));
        //             if (AddedMidP.count(Entry) > 0)
        //                 continue;
        //             AddedMidP.insert(Entry);

        //             int IndexV0 = EdgeInt[j].first;
        //             int IndexV1 = EdgeInt[j].second;
        //             assert(IndexV0 < IndexV1);
        //             ScalarType CurrAlpha = Alpha[j];
        //             assert(IndexV0 >= 0);
        //             assert(IndexV1 >= 0);
        //             assert(IndexV0 < VertPos.size());
        //             assert(IndexV1 < VertPos.size());
        //             EdgeAlpha[EdgeInt[j]] = CurrAlpha;

        //             Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        //             Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        //             Geo::Point3<ScalarType> NewPos = (P0 * CurrAlpha + P1 * (1 - CurrAlpha));

        //             DualParam.SidesMidP.push_back(NewPos);
        //         }
        //     }

        //     DualVertPos.clear();
        //     DualFaces.clear();
        //     // std::vector<Geo::Point3<ScalarType>> SaveMPos;
        //     // std::vector<std::vector<int>> SaveMFaces;

        //     // compute tutte embedding and split faces
        //     for (size_t i = 0; i < PM.PData.SubPatchFaces.size(); i++)
        //     {
        //         // avoid empty patches
        //         if (PM.PData.SubPatchFaces[i].size() == 0)
        //             continue;
        //         // avoid patches with multiple boundaries
        //         if (PM.PData.PatchCorners[i].size() != 1)
        //             continue;
        //         ;

        //         std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
        //         std::vector<std::vector<int>> PatchFaces = PM.PData.SubPatchFaces[i];

        //         if ((DualParam.OnlyNonQuad) && (PM.NumCorners(i) == 4))
        //         {
        //             AppendMesh(DualVertPos, DualFaces, PatchVertPos, PatchFaces);
        //             continue;
        //         }
        //         // compute tutte's embedding
        //         std::vector<Geo::Point2<ScalarType>> VertUV;
        //         bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(PatchVertPos, PatchFaces, VertUV);
        //         assert(done_param);

        //         // get the split segments
        //         // std::vector<Geo::Segment2<ScalarType> > SplitSeg;
        //         std::vector<std::vector<Geo::Point2<ScalarType>>> SplitPoly;
        //         for (size_t j = 0; j < PatchFaces.size(); j++)
        //         {
        //             for (size_t k = 0; k < 3; k++)
        //             {
        //                 int IndexV0 = PatchFaces[j][k];
        //                 int IndexV1 = PatchFaces[j][(k + 1) % 3];

        //                 Geo::Point2<ScalarType> UV0 = VertUV[IndexV0];
        //                 Geo::Point2<ScalarType> UV1 = VertUV[IndexV1];

        //                 IndexV0 = PM.MapToOriginal(i, IndexV0);
        //                 IndexV1 = PM.MapToOriginal(i, IndexV1);

        //                 std::pair<int, int> Key(std::min(IndexV0, IndexV1),
        //                                         std::max(IndexV0, IndexV1));

        //                 if (EdgeAlpha.count(Key) == 0)
        //                     continue;

        //                 ScalarType alpha = EdgeAlpha[Key];

        //                 Geo::Point2<ScalarType> InterpPUV = UV0 * alpha + UV1 * (1 - alpha);
        //                 // inflate a bit
        //                 InterpPUV = InterpPUV * 1.01;
        //                 SplitPoly.resize(SplitPoly.size() + 1);
        //                 SplitPoly.back().push_back(Geo::Point2<ScalarType>(0, 0));
        //                 SplitPoly.back().push_back(InterpPUV);
        //             }
        //         }

        //         std::vector<std::pair<int, int>> NewEdges;
        //         // SplitMeshWithUVPolyline(TestPos,VertUV,SubPatchFaces[i],SplitSeg,NewEdges);
        //         RefineMesh<ScalarType>::SplitMeshWithUVPolylines(PatchVertPos, VertUV,
        //                                                          PatchFaces, SplitPoly, NewEdges);

        //         std::vector<std::vector<int>> NextFP;
        //         std::vector<std::vector<int>> NextEP;
        //         getFFAdjacency(PatchFaces, NextFP, NextEP);

        //         SplitMeshFromEdges(PatchVertPos, PatchFaces, NextFP, NextEP, NewEdges);

        //         AppendMesh(DualVertPos, DualFaces, PatchVertPos, PatchFaces);
        //     }
        //     // WriteOBJ("./dualized.obj", SaveMPos, SaveMFaces);
        //     return true;
        // }

        //**** NEW VERSION ****

        typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;

        static bool BestEdgeSide(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                                 const std::vector<std::pair<int, int>> &PatchEdges,
                                 std::map<PosPair, ScalarType> &DistMap,
                                 PosPair &SplEdge, ScalarType &MinD)
        {
            MinD = std::numeric_limits<ScalarType>::max();
            bool found = false;
            for (size_t j = 0; j < PatchEdges.size(); j++)
            {
                int IndexV0 = PatchEdges[j].first;
                int IndexV1 = PatchEdges[j].second;
                assert(IndexV0 >= 0);
                assert(IndexV1 >= 0);
                assert(IndexV0 < PatchVertPos.size());
                assert(IndexV1 < PatchVertPos.size());

                Geo::Point3<ScalarType> SegP0 = PatchVertPos[IndexV0];
                Geo::Point3<ScalarType> SegP1 = PatchVertPos[IndexV1];
                PosPair PP(std::min(SegP0, SegP1), std::max(SegP0, SegP1));
                assert(DistMap.count(PP) > 0);

                if (DistMap.count(PP) == 0)
                    exit(0);

                ScalarType DistTest = DistMap[PP];
                if (DistTest < 0)
                {
                    exit(0);
                    continue;
                }
                if (DistTest >= MinD)
                    continue;
                MinD = DistTest;
                SplEdge = PP; // PosPair(std::min(SegP0, SegP1),
                              //      std::max(SegP0, SegP1));
                found = true;
                // bestIdx=j;
            }
            if (!found)
            {
                exit(0);
                MinD = -1;
            }
            return found;
        }

        // static void BestSplitEdgeSide(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
        //                               const std::vector<std::pair<int, int>> &PatchEdges,
        //                               std::map<PosPair, ScalarType> &DistMap,
        //                               PosPair &SplEdge, ScalarType &MinD)
        // {
        //     MinD = std::numeric_limits<ScalarType>::max();
        //     bool found = false;
        //     for (size_t j = 0; j < PatchEdges.size(); j++)
        //     {
        //         int IndexV0 = PatchEdges[j].first;
        //         int IndexV1 = PatchEdges[j].second;
        //         assert(IndexV0 >= 0);
        //         assert(IndexV1 >= 0);
        //         assert(IndexV0 < PatchVertPos.size());
        //         assert(IndexV1 < PatchVertPos.size());

        //         Geo::Point3<ScalarType> SegP0 = PatchVertPos[IndexV0];
        //         Geo::Point3<ScalarType> SegP1 = PatchVertPos[IndexV1];
        //         PosPair PP(std::min(SegP0, SegP1), std::max(SegP0, SegP1));
        //         assert(DistMap.count(PP) > 0);
        //         // Geo::Segment3<ScalarType> SegTest(SegP0,SegP1);
        //         // ScalarType DistTest=Geo::Distance(SegTest,CentralPos);
        //         ScalarType DistTest = DistMap[PP];
        //         if (DistTest > MinD)
        //             continue;
        //         MinD = DistTest;
        //         SplEdge = PosPair(std::min(SegP0, SegP1),
        //                           std::max(SegP0, SegP1));
        //         found = true;
        //         // bestIdx=j;
        //     }
        //     if (!found)
        //         MinD = -1;
        //     assert(found);
        //     // SplEdges.push_back(BestPair);
        // }

        static bool BestJoinSplitEdgeSide(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                                          const std::vector<std::pair<int, int>> &PatchEdges,
                                          std::map<PosPair, ScalarType> &JointDistMap,
                                          PosPair &SplEdge, ScalarType &MinD)
        {
            MinD = std::numeric_limits<ScalarType>::max();
            bool found = false;
            for (size_t j = 0; j < PatchEdges.size(); j++)
            {
                int IndexV0 = PatchEdges[j].first;
                int IndexV1 = PatchEdges[j].second;
                assert(IndexV0 >= 0);
                assert(IndexV1 >= 0);
                assert(IndexV0 < PatchVertPos.size());
                assert(IndexV1 < PatchVertPos.size());

                Geo::Point3<ScalarType> SegP0 = PatchVertPos[IndexV0];
                Geo::Point3<ScalarType> SegP1 = PatchVertPos[IndexV1];
                PosPair PP(std::min(SegP0, SegP1), std::max(SegP0, SegP1));
                assert(JointDistMap.count(PP) > 0);

                ScalarType DistTest = JointDistMap[PP];
                if (DistTest < 0)
                    continue;

                if (DistTest > MinD)
                    continue;
                MinD = DistTest;
                SplEdge = PosPair(std::min(SegP0, SegP1),
                                  std::max(SegP0, SegP1));
                found = true;
                // bestIdx=j;
            }
            if (!found)
                MinD = -1;
            return found;
            // SplEdges.push_back(BestPair);
        }

        static void InitPatchBestCandidates(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                                            const std::vector<std::vector<std::pair<int, int>>> &PatchEdges,
                                            std::map<PosPair, ScalarType> &DistMap,
                                            std::vector<PosPair> &BestSplEdges,
                                            std::vector<ScalarType> &BestEdgeDist)
        {
            BestSplEdges.clear();
            PosPair DefPos(Geo::Point3<ScalarType>(0, 0, 0), Geo::Point3<ScalarType>(0, 0, 0));
            BestSplEdges.resize(PatchEdges.size(), DefPos);

            BestEdgeDist.clear();
            BestEdgeDist.resize(PatchEdges.size(), -1);

            for (size_t i = 0; i < PatchEdges.size(); i++)
                BestEdgeSide(PatchVertPos, PatchEdges[i], DistMap, BestSplEdges[i], BestEdgeDist[i]);
            // BestSplitEdgeSide(PatchVertPos, PatchEdges[i], DistMap, BestSplEdges[i], BestEdgeDist[i]);
        }

        static void InitPatchBestJointCandidates(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                                                 const std::vector<std::vector<std::pair<int, int>>> &PatchEdges,
                                                 std::map<PosPair, ScalarType> &DistMap,
                                                 std::vector<PosPair> &BestSplEdges,
                                                 std::vector<ScalarType> &BestEdgeDist)
        {
            PosPair DefPos(Geo::Point3<ScalarType>(0, 0, 0), Geo::Point3<ScalarType>(0, 0, 0));
            BestSplEdges.resize(PatchEdges.size(), DefPos);

            BestEdgeDist.clear();
            BestEdgeDist.resize(PatchEdges.size(), -1);

            for (size_t i = 0; i < PatchEdges.size(); i++)
            {
                bool found = BestEdgeSide(PatchVertPos, PatchEdges[i],
                                          DistMap, BestSplEdges[i],
                                          BestEdgeDist[i]);
                if (!found)
                {
                    exit(0);
                    BestEdgeDist[i] = -1;
                }
            }
        }

        static void ComputeSplitEdges(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                                      const std::vector<PosPair> &BestSplEdges,
                                      const std::vector<PosPair> &BestSplEdgesJoint,
                                      const std::vector<ScalarType> &BestSplEdgesJointDist,
                                      std::vector<PosPair> &SplEdges,
                                      bool UseJointBorder)
        {
            SplEdges.clear();
            assert(BestSplEdges.size() == BestSplEdgesJoint.size());
            assert(BestSplEdges.size() == BestSplEdgesJointDist.size());
            for (size_t i = 0; i < BestSplEdges.size(); i++)
            {
                PosPair BestPair = BestSplEdges[i];
                PosPair BestPairJ = BestSplEdgesJoint[i];
                ScalarType DistJoint = BestSplEdgesJointDist[i];

                if ((DistJoint < 0) || (!UseJointBorder))
                {
                    // exit(0);
                    SplEdges.push_back(BestPair);
                }
                else
                    SplEdges.push_back(BestPairJ);
            }
        }

        static void PosPairToIdx(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                 const std::vector<std::vector<int>> &PatchFaces,
                                 const std::vector<PosPair> &SplEdgesPos,
                                 std::vector<std::pair<int, int>> &SplEdgesIdx)
        {
            std::map<PosPair, std::pair<int, int>> EdgeMap;
            for (size_t i = 0; i < PatchFaces.size(); i++)
            {
                assert(PatchFaces[i].size() == 3);
                for (size_t j = 0; j < 3; j++)
                {
                    int IndexV0 = PatchFaces[i][j];
                    int IndexV1 = PatchFaces[i][(j + 1) % 3];
                    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                    PosPair PP(std::min(P0, P1), std::max(P0, P1));
                    std::pair<int, int> EdgeI(std::min(IndexV0, IndexV1),
                                              std::max(IndexV0, IndexV1));
                    EdgeMap[PP] = EdgeI;
                }
            }
            for (size_t i = 0; i < SplEdgesPos.size(); i++)
            {
                PosPair PP = SplEdgesPos[i];
                assert(EdgeMap.count(PP) > 0);
                if (EdgeMap.count(PP) == 0)
                    exit(0);
                SplEdgesIdx.push_back(EdgeMap[PP]);
            }
        }

        // static void SplitPatchWith(std::vector<Geo::Point3<ScalarType>> &VertPos,
        //                            std::vector<std::vector<int>> &PatchFaces,
        //                            std::vector<Geo::Point2<ScalarType>> &PatchUV,
        //                            const int CentralV,
        //                            const std::vector<PosPair> &SplEdges)
        // {
        //     std::vector<std::pair<int, int>> SplEdgesIdx;

        //     std::map<PosPair, std::pair<int, int>> EdgeMap;
        //     for (size_t i = 0; i < PatchFaces.size(); i++)
        //     {
        //         assert(PatchFaces[i].size() == 3);
        //         for (size_t j = 0; j < 3; j++)
        //         {
        //             int IndexV0 = PatchFaces[i][j];
        //             int IndexV1 = PatchFaces[i][(j + 1) % 3];
        //             Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        //             Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        //             PosPair PP(std::min(P0, P1), std::max(P0, P1));
        //             std::pair<int, int> EdgeI(std::min(IndexV0, IndexV1),
        //                                       std::max(IndexV0, IndexV1));
        //             EdgeMap[PP] = EdgeI;
        //         }
        //     }

        //     for (size_t i = 0; i < SplEdges.size(); i++)
        //     {
        //         PosPair PP = SplEdges[i];
        //         assert(EdgeMap.count(PP) > 0);
        //         SplEdgesIdx.push_back(EdgeMap[PP]);
        //     }

        //     std::vector<std::vector<Geo::Point2<ScalarType>>> SplitPoly;
        //     assert(CentralV != -1);
        //     Geo::Point2<ScalarType> CentralUV = PatchUV[CentralV];
        //     for (size_t i = 0; i < SplEdgesIdx.size(); i++)
        //     {
        //         Geo::Point2<ScalarType> P0 = CentralUV;
        //         int IndexV0 = SplEdgesIdx[i].first;
        //         int IndexV1 = SplEdgesIdx[i].second;
        //         Geo::Point2<ScalarType> P1 = (PatchUV[IndexV0] + PatchUV[IndexV1]) / 2;
        //         // inflate a bit
        //         P1 = P0 + (P1 - P0) * 1.01;
        //         SplitPoly.resize(SplitPoly.size() + 1);
        //         SplitPoly.back().push_back(P0);
        //         SplitPoly.back().push_back(P1);
        //     }

        //     std::vector<std::pair<int, int>> NewEdges;
        //     RefineMesh<ScalarType>::SplitMeshWithUVPolylines(VertPos, PatchUV, PatchFaces, SplitPoly, NewEdges);

        //     std::vector<std::vector<int>> NextFP;
        //     std::vector<std::vector<int>> NextEP;
        //     getFFAdjacency(PatchFaces, NextFP, NextEP);
        //     SplitMeshFromEdges(VertPos, PatchFaces, NextFP, NextEP, NewEdges);
        // }

        static void SplitPatchWith(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                   std::vector<std::vector<int>> &PatchFaces,
                                   std::vector<Geo::Point2<ScalarType>> &PatchUV,
                                   const int CentralV,
                                   const std::vector<std::pair<int, int>> &SplEdgesIdx,
                                   std::vector<int> &NewCorners)
        {
            // std::vector<std::pair<int, int>> SplEdgesIdx;

            // std::map<PosPair, std::pair<int, int>> EdgeMap;
            // for (size_t i = 0; i < PatchFaces.size(); i++)
            // {
            //     assert(PatchFaces[i].size() == 3);
            //     for (size_t j = 0; j < 3; j++)
            //     {
            //         int IndexV0 = PatchFaces[i][j];
            //         int IndexV1 = PatchFaces[i][(j + 1) % 3];
            //         Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
            //         Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
            //         PosPair PP(std::min(P0, P1), std::max(P0, P1));
            //         std::pair<int, int> EdgeI(std::min(IndexV0, IndexV1),
            //                                   std::max(IndexV0, IndexV1));
            //         EdgeMap[PP] = EdgeI;
            //     }
            // }

            // for (size_t i = 0; i < SplEdges.size(); i++)
            // {
            //     PosPair PP = SplEdges[i];
            //     assert(EdgeMap.count(PP) > 0);
            //     SplEdgesIdx.push_back(EdgeMap[PP]);
            // }
            NewCorners.clear();

            // get all border vert and save pos
            std::vector<size_t> OldBorderVIDX;
            getBorderVert(PatchFaces, OldBorderVIDX);
            std::set<Geo::Point3<ScalarType>> OldBorderPos;
            for (size_t i = 0; i < OldBorderVIDX.size(); i++)
                OldBorderPos.insert(VertPos[OldBorderVIDX[i]]);

            std::vector<std::vector<Geo::Point2<ScalarType>>> SplitPoly;
            assert(CentralV != -1);
            Geo::Point2<ScalarType> CentralUV = PatchUV[CentralV];
            for (size_t i = 0; i < SplEdgesIdx.size(); i++)
            {
                Geo::Point2<ScalarType> P0 = CentralUV;
                int IndexV0 = SplEdgesIdx[i].first;
                int IndexV1 = SplEdgesIdx[i].second;
                Geo::Point2<ScalarType> P1 = (PatchUV[IndexV0] + PatchUV[IndexV1]) / 2;
                // inflate a bit
                P1 = P0 + (P1 - P0) * 1.01;
                SplitPoly.resize(SplitPoly.size() + 1);
                SplitPoly.back().push_back(P0);
                SplitPoly.back().push_back(P1);
            }

            std::vector<std::pair<int, int>> NewEdges;
            RefineMesh<ScalarType>::SplitMeshWithUVPolylines(VertPos, PatchUV, PatchFaces, SplitPoly, NewEdges);
            std::vector<size_t> NewBorderVIDX;
            getBorderVert(PatchFaces, NewBorderVIDX);
            for (size_t i = 0; i < NewBorderVIDX.size(); i++)
            {
                int IndexV = NewBorderVIDX[i];
                Geo::Point3<ScalarType> PosTest = VertPos[IndexV];
                if (OldBorderPos.count(PosTest) > 0)
                    continue;
                NewCorners.push_back(IndexV);
            }

            std::vector<std::vector<int>> NextFP;
            std::vector<std::vector<int>> NextEP;
            getFFAdjacency(PatchFaces, NextFP, NextEP);
            SplitMeshFromEdges(VertPos, PatchFaces, NextFP, NextEP, NewEdges);
        }

        // static void InitPatchDistanceMap(const PatchManaging<ScalarType> &PM,
        //                                  std::vector<int> &CenterV,
        //                                  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUV,
        //                                  std::vector<std::map<PosPair, ScalarType>> &PatchDistMap)
        // {
        //     PatchDistMap.clear();
        //     assert(CenterV.size() == PM.NumPatches());
        //     PatchDistMap.resize(CenterV.size());
        //     for (size_t i = 0; i < CenterV.size(); i++)
        //     {
        //         if (CenterV[i] == -1)
        //             continue;
        //         if (PM.isEmpty(i))
        //             continue;
        //         int IndexC = CenterV[i];
        //         Geo::Point3<ScalarType> CentralPos = PM.PData.SubPatchVertPos[i][IndexC];
        //         assert(PM.NumBorders(i) == 1);
        //         for (int s = 0; s < PM.NumSides(i); s++)
        //         {
        //             std::vector<int> SideV;
        //             PM.GetAllSideVert(i, s, SideV);
        //             assert(SideV.size() >= 2);
        //             for (int e = 0; e < SideV.size() - 1; e++)
        //             {
        //                 int IndexV0 = SideV[e];
        //                 int IndexV1 = SideV[e + 1];
        //                 Geo::Point3<ScalarType> SegP0 = PM.PData.SubPatchVertPos[i][IndexV0];
        //                 Geo::Point3<ScalarType> SegP1 = PM.PData.SubPatchVertPos[i][IndexV1];
        //                 Geo::Point3<ScalarType> AvgP=(SegP0+SegP1)/(ScalarType)2;
        //                 // Geo::Segment3<ScalarType> SegTest(SegP0, SegP1);
        //                 // ScalarType DistTest = Geo::Distance(SegTest, CentralPos);
        //                 ScalarType DistTest = (CentralPos-AvgP).Norm();
        //                 if (DistTest<0)exit(0);
        //                 PosPair Key(std::min(SegP0, SegP1),
        //                             std::max(SegP0, SegP1));
        //                 PatchDistMap[i][Key] = DistTest;
        //             }
        //         }
        //     }
        // }

        static ScalarType DistanceFunct(const Geo::Point2<ScalarType> &SegP0UV,
                                        const Geo::Point2<ScalarType> &SegP1UV,
                                        const Geo::Point2<ScalarType> &CornP0UV,
                                        const Geo::Point2<ScalarType> &CornP1UV,
                                        const Geo::Point2<ScalarType> &CentralP)
        {
            ScalarType DCentr_Corn0 = (CentralP - CornP0UV).Norm();
            ScalarType DCentr_Corn1 = (CentralP - CornP1UV).Norm();
            Geo::Point2<ScalarType> AvgP = (SegP0UV + SegP1UV) / 2;
            ScalarType DAvg_Corn0 = (AvgP - CornP0UV).Norm();
            ScalarType DAvg_Corn1 = (AvgP - CornP1UV).Norm();
            ScalarType Ratio0 = DCentr_Corn0 / DCentr_Corn1;
            ScalarType Ratio1 = DAvg_Corn0 / DAvg_Corn1;
            return (pow(Ratio0 - Ratio1, 2));
        }

        static void InitPatchDistanceMap(const PatchManaging<ScalarType> &PM,
                                         const std::vector<int> &CenterV,
                                         const std::vector<std::vector<Geo::Point2<ScalarType>>> &PatchVertUV,
                                         std::vector<std::map<PosPair, ScalarType>> &PatchDistMap)
        {
            PatchDistMap.clear();
            assert(CenterV.size() == PM.NumPatches());
            PatchDistMap.resize(CenterV.size());
            for (size_t i = 0; i < CenterV.size(); i++)
            {
                if (CenterV[i] == -1)
                    continue;
                if (PM.isEmpty(i))
                    continue;
                int IndexC = CenterV[i];
                Geo::Point2<ScalarType> CentralPos = PatchVertUV[i][IndexC];
                assert(PM.NumBorders(i) == 1);
                for (int s = 0; s < PM.NumSides(i); s++)
                {
                    std::vector<int> SideV;
                    PM.GetAllSideVert(i, s, SideV);
                    assert(SideV.size() >= 2);
                    int IndexC0 = SideV[0];
                    int IndexC1 = SideV.back();
                    Geo::Point2<ScalarType> UVC0 = PatchVertUV[i][IndexC0];
                    Geo::Point2<ScalarType> UVC1 = PatchVertUV[i][IndexC1];
                    for (int e = 0; e < SideV.size() - 1; e++)
                    {
                        int IndexV0 = SideV[e];
                        int IndexV1 = SideV[e + 1];
                        Geo::Point2<ScalarType> SegUVP0 = PatchVertUV[i][IndexV0];
                        Geo::Point2<ScalarType> SegUVP1 = PatchVertUV[i][IndexV1];
                        // Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
                        // // Geo::Segment3<ScalarType> SegTest(SegP0, SegP1);
                        // // ScalarType DistTest = Geo::Distance(SegTest, CentralPos);
                        // ScalarType DistTest = (CentralPos - AvgP).Norm();
                        ScalarType DistTest = DistanceFunct(SegUVP0, SegUVP1, UVC0, UVC1, CentralPos);
                        if (DistTest < 0)
                            exit(0);
                        Geo::Point3<ScalarType> SegP0_3D = PM.PData.SubPatchVertPos[i][IndexV0];
                        Geo::Point3<ScalarType> SegP1_3D = PM.PData.SubPatchVertPos[i][IndexV1];
                        PosPair Key(std::min(SegP0_3D, SegP1_3D),
                                    std::max(SegP0_3D, SegP1_3D));
                        PatchDistMap[i][Key] = DistTest;
                    }
                }
            }
        }

        // static void InitPatchDistanceMap(const PatchManaging<ScalarType> &PM,
        //                                  const std::vector<int> &CenterV,
        //                                  const std::vector<std::vector<Geo::Point2<ScalarType>>> &PatchVertUV,
        //                                  std::vector<std::map<PosPair, ScalarType>> &PatchDistMap)
        // {
        //     PatchDistMap.clear();
        //     assert(CenterV.size() == PM.NumPatches());
        //     PatchDistMap.resize(CenterV.size());
        //     for (size_t i = 0; i < CenterV.size(); i++)
        //     {
        //         if (CenterV[i] == -1)
        //             continue;
        //         if (PM.isEmpty(i))
        //             continue;
        //         int IndexC = CenterV[i];
        //         Geo::Point2<ScalarType> CentralPos = PatchVertUV[i][IndexC];
        //         assert(PM.NumBorders(i) == 1);
        //         for (int s = 0; s < PM.NumSides(i); s++)
        //         {
        //             std::vector<int> SideV;
        //             PM.GetAllSideVert(i, s, SideV);
        //             assert(SideV.size() >= 2);
        //             for (int e = 0; e < SideV.size() - 1; e++)
        //             {
        //                 int IndexV0 = SideV[e];
        //                 int IndexV1 = SideV[e + 1];
        //                 Geo::Point2<ScalarType> SegP0 = PatchVertUV[i][IndexV0];
        //                 Geo::Point2<ScalarType> SegP1 = PatchVertUV[i][IndexV1];
        //                 Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
        //                 // Geo::Segment3<ScalarType> SegTest(SegP0, SegP1);
        //                 // ScalarType DistTest = Geo::Distance(SegTest, CentralPos);
        //                 ScalarType DistTest = (CentralPos - AvgP).Norm();
        //                 if (DistTest < 0)
        //                     exit(0);
        //                 Geo::Point3<ScalarType> SegP0_3D = PM.PData.SubPatchVertPos[i][IndexV0];
        //                 Geo::Point3<ScalarType> SegP1_3D = PM.PData.SubPatchVertPos[i][IndexV1];
        //                 PosPair Key(std::min(SegP0_3D, SegP1_3D),
        //                             std::max(SegP0_3D, SegP1_3D));
        //                 PatchDistMap[i][Key] = DistTest;
        //             }
        //         }
        //     }
        // }

        static void UpdateJointAlignDistanceMap(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                                const std::vector<std::vector<std::pair<int, int>>> &PatchEdges,
                                                const std::vector<PosPair> &BestEdge,
                                                const std::vector<ScalarType> &BestDist,
                                                std::map<PosPair, ScalarType> &PatchDistMap,
                                                std::map<PosPair, ScalarType> &JointDistMap)
        {
            // for each side
            for (size_t i = 0; i < PatchEdges.size(); i++)
            {
                // for each edge of the side
                for (size_t j = 0; j < PatchEdges[i].size(); j++)
                {
                    std::pair<int, int> EdgeV = PatchEdges[i][j];
                    Geo::Point3<ScalarType> VPos0 = VertPos[EdgeV.first];
                    Geo::Point3<ScalarType> VPos1 = VertPos[EdgeV.second];
                    PosPair Key(std::min(VPos0, VPos1), std::max(VPos0, VPos1));
                    assert(PatchDistMap.count(Key) > 0);
                    if (PatchDistMap.count(Key) == 0)
                        exit(0);

                    // get the distance of the patch
                    ScalarType PatchD = PatchDistMap[Key];

                    bool HasInserted = (JointDistMap.count(Key) > 0);
                    if (HasInserted)
                    {
                        // //means cannot join from previous side
                        // if (JointDistMap[Key]<0)continue;
                        JointDistMap[Key] += PatchD;
                    }
                    else
                        JointDistMap[Key] = PatchD;
                    // VertPos
                }
            }
        }

        // static ScalarType SingPosEnegy(Geo::Point2<ScalarType> &CandidatePos,
        //                                std::vector<Geo::Point2<ScalarType> > &SidePos)
        // {
        //     ScalarType DotSum=0;
        //     ScalarType IdealA=2*M_PI/(SidePos.size());
        //     for (size_t s=0;s<SidePos.size();s++)
        //     {
        //         Geo::Point2<ScalarType> SideP0=SidePos[s];
        //         Geo::Point2<ScalarType> SideP1=SidePos[(s+1)%SidePos.size()];
        //         Geo::Point2<ScalarType> Dir0=SideP0-CandidatePos;
        //         Geo::Point2<ScalarType> Dir1=SideP1-CandidatePos;
        //         Dir0.Normalize();
        //         Dir1.Normalize();
        //         ScalarType AngleVal=Geo::AngleRad(Dir0,Dir1);
        //         DotSum+=pow(AngleVal-IdealA,2);
        //         //DotSum=std::max(DotSum,pow(AngleVal-IdealA,2));
        //     }
        //     return DotSum;
        // }

        static ScalarType SingPosEnegy(const Geo::Point2<ScalarType> &CandidatePos,
                                       const std::vector<Geo::Point2<ScalarType>> &CornerPos,
                                       const std::vector<Geo::Point2<ScalarType>> &SidePos)
        {
            ScalarType DotSum = 0;
            ScalarType IdealA = 2 * M_PI / (SidePos.size());
            for (size_t s = 0; s < SidePos.size(); s++)
            {
                // Geo::Point2<ScalarType> SideP0=SidePos[s];
                // Geo::Point2<ScalarType> SideP1=SidePos[(s+1)%SidePos.size()];
                // Geo::Point2<ScalarType> Dir0=SideP0-CandidatePos;
                // Geo::Point2<ScalarType> Dir1=SideP1-CandidatePos;
                // Dir0.Normalize();
                // Dir1.Normalize();
                // ScalarType AngleVal=Geo::AngleRad(Dir0,Dir1);
                // DotSum+=pow(AngleVal-IdealA,2);
                // DotSum=std::max(DotSum,pow(AngleVal-IdealA,2));
                // get the quad four points
                Geo::Point2<ScalarType> C0 = SidePos[s];
                Geo::Point2<ScalarType> C1 = CornerPos[(s + 1) % CornerPos.size()];
                Geo::Point2<ScalarType> C2 = SidePos[(s + 1) % SidePos.size()];
                Geo::Point2<ScalarType> C3 = CandidatePos;
                Geo::Point2<ScalarType> Dir0 = C1 - C0;
                Geo::Point2<ScalarType> Dir1 = C2 - C1;
                Geo::Point2<ScalarType> Dir2 = C3 - C2;
                Geo::Point2<ScalarType> Dir3 = C0 - C3;
                Dir0.Normalize();
                Dir1.Normalize();
                Dir2.Normalize();
                Dir3.Normalize();
                // ScalarType AngleVal0=Geo::AngleRad(Dir0,Dir1);
                ScalarType AngleVal1 = Geo::AngleRad(Dir1, Dir2);
                ScalarType AngleVal2 = Geo::AngleRad(Dir2, Dir3);
                ScalarType AngleVal3 = Geo::AngleRad(Dir3, Dir0);
                // ScalarType Dev0=pow(AngleVal0-M_PI/2,2);
                ScalarType Dev1 = fabs(AngleVal1 - M_PI / 2); //,2);
                ScalarType Dev2 = fabs(AngleVal2 - IdealA);   //,2);
                ScalarType Dev3 = fabs(AngleVal3 - M_PI / 2); //,2);
                // DotSum+=Dev0;
                DotSum += Dev1;
                DotSum += Dev2;
                DotSum += Dev3;
            }
            return DotSum;
        }

        static int OptimizeSingularPos(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                       const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
                                       const std::vector<std::pair<int, int>> &SplEdgesIdx,
                                       const std::vector<int> &CornerV,
                                       const std::vector<bool> &PerPatchIsB)
        {
            int BestIdx = -1;
            ScalarType BestD = std::numeric_limits<ScalarType>::max();
            assert(PerPatchIsB.size() == VertPos.size());

            std::vector<Geo::Point2<ScalarType>> SidePos;
            for (size_t e = 0; e < SplEdgesIdx.size(); e++)
            {
                int IndexV0 = SplEdgesIdx[e].first;
                int IndexV1 = SplEdgesIdx[e].second;
                Geo::Point2<ScalarType> SegP0 = PatchVertUV[IndexV0];
                Geo::Point2<ScalarType> SegP1 = PatchVertUV[IndexV1];
                Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
                SidePos.push_back(AvgP);
            }

            std::vector<Geo::Point2<ScalarType>> CornerPos;
            for (size_t e = 0; e < CornerV.size(); e++)
                CornerPos.push_back(PatchVertUV[CornerV[e]]);

            for (size_t i = 0; i < VertPos.size(); i++)
            {
                if (PerPatchIsB[i])
                    continue;

                Geo::Point2<ScalarType> CurrUV = PatchVertUV[i];
                ScalarType DistSqTest = SingPosEnegy(CurrUV, CornerPos, SidePos);

                if (DistSqTest > BestD)
                    continue;
                BestD = DistSqTest;
                BestIdx = i;
            }
            assert(BestIdx != -1);
            return BestIdx;
        }

        // static int OptimizeSingularPos(const std::vector<Geo::Point3<ScalarType>> &VertPos,
        //                                const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
        //                                const std::vector<std::pair<int, int>> &SplEdgesIdx,
        //                                const std::vector<bool> &PerPatchIsB)
        // {
        //     int BestIdx = -1;
        //     ScalarType BestD = std::numeric_limits<ScalarType>::max();
        //     assert(PerPatchIsB.size() == VertPos.size());

        //     for (size_t i = 0; i < VertPos.size(); i++)
        //     {
        //         if (PerPatchIsB[i])
        //             continue;
        //         ScalarType DistSqTest = 0;
        //         Geo::Point2<ScalarType> CurrUV = PatchVertUV[i];
        //         for (size_t e = 0; e < SplEdgesIdx.size(); e++)
        //         {
        //             int IndexV0 = SplEdgesIdx[e].first;
        //             int IndexV1 = SplEdgesIdx[e].second;
        //             Geo::Point2<ScalarType> SegP0 = PatchVertUV[IndexV0];
        //             Geo::Point2<ScalarType> SegP1 = PatchVertUV[IndexV1];
        //             Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
        //             DistSqTest += pow((CurrUV - AvgP).Norm(), 2);
        //         }
        //         if (DistSqTest > BestD)
        //             continue;
        //         BestD = DistSqTest;
        //         BestIdx = i;
        //     }
        //     assert(BestIdx != -1);
        //     return BestIdx;
        // }

        static void MakeConsistent(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                   std::vector<std::vector<int>> &Faces,
                                   const std::vector<PosPair> &SplEdges)
        {
            // SAVE INITIAL BORDER (REFINE unfortunately CANCEL IT)
            std::vector<std::pair<int, int>> BorderE;
            getBorderEdges(Faces, BorderE);

            std::set<PosPair> BorderPosSet;
            for (size_t i = 0; i < BorderE.size(); i++)
            {
                int IndexV0 = BorderE[i].first;
                int IndexV1 = BorderE[i].second;
                Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                PosPair PP(std::min(P0, P1), std::max(P0, P1));
                BorderPosSet.insert(PP);
            }

            // SPLIT ON BORDER
            std::set<PosPair> SplEdgesSet(SplEdges.begin(), SplEdges.end());
            std::vector<std::pair<int, int>> ToSplit;
            for (size_t i = 0; i < Faces.size(); i++)
            {
                assert(Faces[i].size() == 3);
                for (size_t j = 0; j < 3; j++)
                {
                    int IndexV0 = Faces[i][j];
                    int IndexV1 = Faces[i][(j + 1) % 3];
                    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                    PosPair PP(std::min(P0, P1), std::max(P0, P1));

                    if (SplEdgesSet.count(PP) == 0)
                        continue;

                    std::pair<int, int> KeyVert(std::min(IndexV0, IndexV1),
                                                std::max(IndexV0, IndexV1));
                    ToSplit.push_back(KeyVert);
                }
            }

            std::vector<Geo::Point2<ScalarType>> PosUV(VertPos.size(), Geo::Point2<ScalarType>(0, 0));
            std::vector<ScalarType> Alpha(ToSplit.size(), 0.5);
            std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> MarkedE;
            RefineMesh<ScalarType>::SplitTriangleFaces(VertPos, PosUV, Faces, ToSplit, Alpha, MarkedE);

            // REINTRODUCE INITIAL BORDER (REFINE unfortunately CANCEL IT)
            BorderE.clear();
            for (size_t i = 0; i < Faces.size(); i++)
            {
                assert(Faces[i].size() == 3);
                for (size_t j = 0; j < 3; j++)
                {
                    int IndexV0 = Faces[i][j];
                    int IndexV1 = Faces[i][(j + 1) % 3];
                    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                    PosPair PP(std::min(P0, P1), std::max(P0, P1));
                    if (BorderPosSet.count(PP) == 0)
                        continue;

                    std::pair<int, int> KeyVert(std::min(IndexV0, IndexV1),
                                                std::max(IndexV0, IndexV1));
                    BorderE.push_back(KeyVert);
                }
            }
            std::vector<std::vector<int>> NextF;
            std::vector<std::vector<int>> NextE;
            getFFAdjacency(Faces, NextF, NextE);
            SplitMeshFromEdges(VertPos, Faces, NextF, NextE, BorderE);
        }

        static void FindSplitPos(const PatchManaging<ScalarType> &PM,
                                 const std::vector<int> &CentralV,
                                 const std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVArap,
                                 std::vector<PosPair> &AllSplEdgesPos,
                                 std::vector<std::vector<std::pair<int, int>>> &SplEdgesIDX,
                                 bool UseJointBorder)
        {
            std::vector<std::vector<PosPair>> SplEdgesPos;
            SplEdgesPos.resize(PM.NumPatches());

            std::vector<std::map<PosPair, ScalarType>> DistMap;
            InitPatchDistanceMap(PM, CentralV, PatchVertUVArap, DistMap);

            std::vector<std::vector<PosPair>> BestSplEdges;
            BestSplEdges.resize(PM.NumPatches());

            std::vector<std::vector<ScalarType>> BestSplEdgesDist;
            BestSplEdgesDist.resize(PM.NumPatches());
            std::map<PosPair, ScalarType> JointDistMap;
            for (size_t i = 0; i < PM.NumPatches(); i++)
            {
                if (CentralV[i] == -1)
                    continue;

                std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
                std::vector<std::vector<std::pair<int, int>>> PatchEdges = PM.PData.EdgePatchSides[i][0];

                InitPatchBestCandidates(PatchVertPos, PatchEdges, DistMap[i], BestSplEdges[i], BestSplEdgesDist[i]);

                UpdateJointAlignDistanceMap(PatchVertPos, PatchEdges, BestSplEdges[i],
                                            BestSplEdgesDist[i], DistMap[i], JointDistMap);
            }

            std::vector<std::vector<PosPair>> BestSplJointEdges;
            BestSplJointEdges.resize(PM.NumPatches());

            std::vector<std::vector<ScalarType>> BestSplJointEdgesDist;
            BestSplJointEdgesDist.resize(PM.NumPatches());
            for (size_t i = 0; i < PM.NumPatches(); i++)
            {
                if (CentralV[i] == -1)
                    continue;

                std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
                std::vector<std::vector<std::pair<int, int>>> PatchEdges = PM.PData.EdgePatchSides[i][0];

                InitPatchBestJointCandidates(PatchVertPos, PatchEdges, JointDistMap,
                                             BestSplJointEdges[i], BestSplJointEdgesDist[i]);
            }

            for (size_t i = 0; i < PM.NumPatches(); i++)
            {
                if (CentralV[i] == -1)
                    continue;

                std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
                std::vector<std::vector<int>> PatchFaces = PM.PData.SubPatchFaces[i];
                std::vector<std::vector<std::pair<int, int>>> PatchEdges = PM.PData.EdgePatchSides[i][0];
                // ComputeSplitEdges(PatchVertPos, PatchEdges, DistMap[i], SplEdgesPos[i]);
                ComputeSplitEdges(PatchVertPos, BestSplEdges[i],
                                  BestSplJointEdges[i], BestSplJointEdgesDist[i],
                                  SplEdgesPos[i], UseJointBorder);
            }

            // then propagate on the non done so far
            // std::vector<PosPair> AllSplEdgesPos;
            AllSplEdgesPos.clear();
            for (size_t i = 0; i < SplEdgesPos.size(); i++)
                AllSplEdgesPos.insert(AllSplEdgesPos.end(), SplEdgesPos[i].begin(), SplEdgesPos[i].end());

            // map pair of pos to pair of indexes
            // std::vector<std::vector<std::pair<int, int>>> SplEdgesIDX;
            SplEdgesIDX.clear();
            SplEdgesIDX.resize(PM.NumPatches());
            for (size_t i = 0; i < PM.NumPatches(); i++)
                PosPairToIdx(PM.PData.SubPatchVertPos[i], PM.PData.SubPatchFaces[i],
                             SplEdgesPos[i], SplEdgesIDX[i]);
        }

        static int FindInitialCentralPos(const std::vector<Geo::Point2<ScalarType>> &VertUV,
                                         const std::vector<bool> &IsB)

        {
            ScalarType MinD = std::numeric_limits<ScalarType>::max();
            int IndexMin = -1;
            for (size_t i = 0; i < VertUV.size(); i++)
            {
                if (IsB[i])
                    continue;
                // get closest to (0,0)
                if (VertUV[i].Norm() > MinD)
                    continue;
                MinD = VertUV[i].Norm();
                IndexMin = i;
            }
            assert(IndexMin >= 0);
            return IndexMin;
        }

        static void ExtractCornerVisPos(const std::vector<Geo::Point3<ScalarType>> &DualPos,
                                        const std::vector<std::vector<int>> &DualFaces,
                                        const std::vector<int> &PatchFaceIDX,
                                        const std::vector<int> &PatchCornerIDX,
                                        std::vector<Geo::Point3<ScalarType>> &VisCornerPos)
        {
            // for (size_t i=0;i<PatchCornerIDX.size();i++)
            // {
            //     int IndexV=PatchCornerIDX[i];
            //     Geo::Point3<ScalarType> Pos=DualPos[IndexV];
            //     VisCornerPos.push_back(Pos);
            // }
            std::set<int> PatchCornerIDXSet(PatchCornerIDX.begin(), PatchCornerIDX.end());

            typedef std::pair<Geo::Point3<ScalarType>, int> AvgPos;
            VisCornerPos.clear();
            std::map<int, AvgPos> AvgCornMap;
            for (size_t i = 0; i < PatchFaceIDX.size(); i++)
            {
                int IndexF = PatchFaceIDX[i];
                assert(DualFaces[IndexF].size() == 3);

                // compute barycenter
                int IndexV0 = DualFaces[IndexF][0];
                int IndexV1 = DualFaces[IndexF][1];
                int IndexV2 = DualFaces[IndexF][2];
                Geo::Point3<ScalarType> Pos0 = DualPos[IndexV0];
                Geo::Point3<ScalarType> Pos1 = DualPos[IndexV1];
                Geo::Point3<ScalarType> Pos2 = DualPos[IndexV2];
                Geo::Point3<ScalarType> BaryFace = (Pos0 + Pos1 + Pos2) / 3;

                for (size_t j = 0; j < 3; j++)
                {
                    int VertIdx = DualFaces[IndexF][j];
                    if (PatchCornerIDXSet.count(VertIdx) == 0)
                        continue;
                    if (AvgCornMap.count(VertIdx) == 0)
                    {
                        AvgCornMap[VertIdx] = AvgPos(BaryFace, 1);
                    }
                    else
                    {
                        AvgCornMap[VertIdx].first += BaryFace;
                        AvgCornMap[VertIdx].second++;
                    }
                }
            }

            for (size_t i = 0; i < PatchCornerIDX.size(); i++)
            {
                int VertIdx = PatchCornerIDX[i];
                assert(AvgCornMap.count(VertIdx) > 0);
                Geo::Point3<ScalarType> PosSum = AvgCornMap[VertIdx].first;
                int Num = AvgCornMap[VertIdx].second;
                VisCornerPos.push_back(PosSum / Num);
            }
        }

        static void ExtractPatchCorners(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
                                        const std::vector<std::vector<int>> &PatchFaces,
                                        const int BaseNum,
                                        const std::set<Geo::Point3<ScalarType>> &PatchCornerPos,
                                        std::vector<int> &CornerIDX)
        {
            CornerIDX.clear();

            std::vector<std::vector<int>> BorderSeq;
            getBorderSequncesVert(PatchFaces, BorderSeq);
            // only one loop
            // std::cout<<"Size B:"<<BorderSeq.size()<<std::endl;
            assert(BorderSeq.size() == 1);
            for (int i = 0; i < BorderSeq[0].size(); i++)
            {
                int VIndex = BorderSeq[0][i];
                Geo::Point3<ScalarType> Pos = PatchPos[VIndex];
                if (PatchCornerPos.count(Pos) == 0)
                    continue;
                CornerIDX.push_back(VIndex + BaseNum);
            }
        }

    public:
        enum RegType
        {
            RTypeEven,
            RTypeOrtho
        };

        struct DualParameter
        {
            // input
            std::vector<int> OnlyPatch;
            bool OnlyNonQuad;
            bool InitSingPos;
            bool FixSingPos;
            bool JointBorder;
            int MaxIte;
            RegType RType;
            int SmoothSteps;

            DualParameter()
            {
                OnlyNonQuad = false;
                InitSingPos = false;
                FixSingPos = false;
                JointBorder = true;
                MaxIte = 5;
                RType = RTypeEven;
                //RType = RTypeOrtho;
                SmoothSteps = 5;
            }
        };

        static void ExtractCornerVisPos(const std::vector<Geo::Point3<ScalarType>> &DualPos,
                                        const std::vector<std::vector<int>> &DualFaces,
                                        const std::vector<std::vector<int>> &PerPatchFaceIDX,
                                        const std::vector<std::vector<int>> &PerPatchCornerIDX,
                                        std::vector<Geo::Point3<ScalarType>> &VisCornerPos)
        {
            VisCornerPos.clear();
            for (size_t i = 0; i < PerPatchFaceIDX.size(); i++)
            {
                std::vector<Geo::Point3<ScalarType>> PatchVisCornerPos;
                ExtractCornerVisPos(DualPos, DualFaces, PerPatchFaceIDX[i],
                                    PerPatchCornerIDX[i], PatchVisCornerPos);
                VisCornerPos.insert(VisCornerPos.end(),
                                    PatchVisCornerPos.begin(),
                                    PatchVisCornerPos.end());
            }
        }

        static void DisambiguateForDual(PatchManaging<ScalarType> &PM)
        {
             // INITIALIZE
            for (size_t i = 0; i < PM.NumPatches(); i++)
            {
                if (PM.NumBorders(i) != 1)
                    continue;
                if (PM.isEmpty(i))
                    continue;
                //ScalarType delta= AvgEdgeLen(PM.PData.SubPatchVertPos[i],PM.PData.SubPatchFaces[i])*0.0001;
                ScalarType delta = std::numeric_limits<float>::epsilon() * (ScalarType)10;
                // DesambiguateSamePosVertices<ScalarType>(PM.PData.SubPatchVertPos[i],
                //                                         PM.PData.SubPatchFaces[i],delta);
                DesambiguateSamePosVertices<ScalarType>(PM.PData.SubPatchVertPos[i],
                                                        PM.PData.SubPatchFaces[i], delta);
                // return (IsCatmullClarkable(PatchM.PData.SideLenghts[IndexPatch][0],SideThr,true));
            }

        }

        static bool ExtractDual(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                const std::vector<std::vector<int>> &Faces,
                                const PatchManaging<ScalarType> &PM,
                                DualParameter &DualParam,
                                std::vector<Geo::Point3<ScalarType>> &DualVertPos,
                                std::vector<std::vector<int>> &DualFaces,
                                std::vector<std::vector<int>> &PerPatchVertIDX,
                                std::vector<std::vector<int>> &PerPatchFaceIDX,
                                std::vector<std::vector<int>> &PerPatchCornerIDX,
                                bool find_corners = false,
                                bool write_timing = true)
        {
            DualVertPos.clear();
            DualFaces.clear();
            PerPatchFaceIDX.clear();
            PerPatchVertIDX.clear();
            PerPatchCornerIDX.clear();
            // first round compute the splitting points
            // std::vector<std::vector<PosPair>> SplEdgesPos;
            // SplEdgesPos.resize(PM.NumPatches());
            std::vector<int> CentralV(PM.NumPatches(), -1);
            std::vector<bool> HasSing(PM.NumPatches(), false);

            std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUV;
            std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVArap;
            PatchVertUV.resize(PM.NumPatches());
            PatchVertUVArap.resize(PM.NumPatches());

            std::vector<std::vector<bool>> PerPatchIsB;
            PerPatchIsB.resize(PM.NumPatches());

            auto t0 = std::chrono::high_resolution_clock::now();

            // // INITIALIZE
            // for (size_t i = 0; i < PM.NumPatches(); i++)
            // {
            //     assert(IndexPatch < PatchM.PData.SideLenghts.size());
            //     if (PM.NumBorders(i) != 1)
            //         continue;
            //     if (PM.isEmpty(i))
            //         continue;
            //     // ScalarType delta= AvgEdgeLen(Pos,Faces)*0.0001;
            //     ScalarType delta = std::numeric_limits<float>::epsilon() * (ScalarType)2;
            //     // DesambiguateSamePosVertices<ScalarType>(PM.PData.SubPatchVertPos[i],
            //     //                                         PM.PData.SubPatchFaces[i],delta);
            //     DesambiguateSamePosVertices<ScalarType>(PM.PData.SubPatchVertPos[i],
            //                                             PM.PData.SubPatchFaces[i], delta);
            //     // return (IsCatmullClarkable(PatchM.PData.SideLenghts[IndexPatch][0],SideThr,true));
            // }

            for (size_t i = 0; i < PM.NumPatches(); i++)
            {
                // avoid empty patches
                if (PM.PData.SubPatchFaces[i].size() == 0)
                {
                    std::cout << "WARNING PATCH ZERO SIDES" << std::endl;
                    continue;
                }
                // avoid patches with multiple boundaries
                if (PM.PData.PatchCorners[i].size() != 1)
                {
                    std::cout << "WARNING PATCH MULTIPL SIDES" << std::endl;
                    continue;
                }
                // avoid valence 2 (or smaller) paches
                if (PM.PData.PatchCorners[i][0].size() < 3)
                {
                    std::cout << "WARNING PATCH SIDES < 3" << std::endl;
                    continue;
                }
                int NumSides = PM.NumSides(i);

                if ((NumSides == 4) && (DualParam.OnlyNonQuad))
                    continue;

                std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
                std::vector<std::vector<int>> PatchFaces = PM.PData.SubPatchFaces[i];
                getBorderVert(PatchFaces, PerPatchIsB[i]);

                bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(PatchVertPos, PatchFaces,
                                                                           PM.PData.PatchCorners[i][0],
                                                                           PatchVertUV[i], true);
                // bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(PatchVertPos, PatchFaces,PatchVertUV[i]);
                assert(done_param);
                bool done_param_1 = Geo::UVParam<ScalarType>::LeastSquaresConformal(PatchVertPos, PatchFaces, PatchVertUVArap[i]);
                // bool done_param_1 = Geo::UVParam<ScalarType>::Arap(PatchVertPos, PatchFaces, PatchVertUVArap[i]);
                assert(done_param_1);

                // see if it has singularities
                std::vector<int> PatchSing;
                PM.PatchSingularities(i, PatchSing);
                if ((PatchSing.size() == 1) && (DualParam.InitSingPos))
                {
                    HasSing[i] = true;
                    CentralV[i] = PatchSing[0];
                }
                else
                {
                    CentralV[i] = FindInitialCentralPos(PatchVertUV[i], PerPatchIsB[i]);
                    // assert(PM.PData.PatchCorners[i][0].size() == 4);
                }
            }

            // // then propagate on the non done so far
            std::vector<PosPair> AllSplEdgesPos;

            // // map pair of pos to pair of indexes
            std::vector<std::vector<std::pair<int, int>>> SplEdgesIDX;

            for (size_t s = 0; s < DualParam.MaxIte; s++)
            {
                if (DualParam.RType == RTypeEven)
                    FindSplitPos(PM, CentralV, PatchVertUV, AllSplEdgesPos, SplEdgesIDX, DualParam.JointBorder);
                else
                    FindSplitPos(PM, CentralV, PatchVertUVArap, AllSplEdgesPos, SplEdgesIDX, DualParam.JointBorder);
                // optimize singularities
                for (size_t i = 0; i < PM.NumPatches(); i++)
                {
                    if (CentralV[i] == -1)
                        continue;

                    std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
                    std::vector<std::vector<int>> PatchFaces = PM.PData.SubPatchFaces[i];
                    std::vector<std::vector<std::pair<int, int>>> PatchEdges = PM.PData.EdgePatchSides[i][0];

                    if (HasSing[i] && DualParam.FixSingPos)
                        continue;

                    if (DualParam.RType == RTypeEven)
                    {
                        CentralV[i] = OptimizeSingularPos(PatchVertPos, PatchVertUV[i],
                                                          SplEdgesIDX[i], PM.PData.PatchCorners[i][0],
                                                          PerPatchIsB[i]);
                    }
                    else
                    {
                        CentralV[i] = OptimizeSingularPos(PatchVertPos, PatchVertUVArap[i],
                                                          SplEdgesIDX[i], PM.PData.PatchCorners[i][0],
                                                          PerPatchIsB[i]);
                    }
                    // CentralV[i] = OptimizeSingularPos(PatchVertPos, PatchVertUV[i], SplEdgesIDX[i], PerPatchIsB[i]);
                }
            }

            FindSplitPos(PM, CentralV, PatchVertUVArap, AllSplEdgesPos, SplEdgesIDX, DualParam.JointBorder);
            auto t1 = std::chrono::high_resolution_clock::now();

            // SAVE CORNER POS
            //  std::vector<Geo::Point3<ScalarType>> NewCornerPos;
            //  for (size_t i = 0; i < AllSplEdgesPos.size(); i++)
            //  {
            //      Geo::Point3<ScalarType> P0 = AllSplEdgesPos[i].first;
            //      Geo::Point3<ScalarType> P1 = AllSplEdgesPos[i].second;
            //      NewCornerPos.push_back(P0 * 0.5 + P1 *(1- 0.5));
            //      NewCornerPos.push_back(P1 * 0.5 + P0 * (1- 0.5));
            //  }
            //  //SAVE CENTRAL POS
            //  std::vector<Geo::Point3<ScalarType>> CentralPos;
            //  for (size_t i = 0; i < CentralV.size(); i++)
            //  {
            //      if (CentralV[i] == -1)
            //          continue;
            //      int IndexV = CentralV[i];
            //      CentralPos.push_back(PM.PData.SubPatchVertPos[i][IndexV]);
            //  }
            //  std::set<Geo::Point3<ScalarType>> NewCornerSet(NewCornerPos.begin(),
            //                                                 NewCornerPos.end());
            //  NewCornerSet.insert(CentralPos.begin(),CentralPos.end());

            std::set<int> SplitPatchSet(DualParam.OnlyPatch.begin(), DualParam.OnlyPatch.end());
            for (size_t i = 0; i < PM.NumPatches(); i++)
            {
                // avoid empty patches
                if (PM.PData.SubPatchFaces[i].size() == 0)
                    continue;
                // avoid patches with multiple boundaries
                if (PM.PData.PatchCorners[i].size() != 1)
                    continue;

                // avoid valence 2 (or smaller) paches
                if (PM.PData.PatchCorners[i][0].size() < 3)
                    continue;

                // tutte embedding
                std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
                std::vector<std::vector<int>> PatchFaces = PM.PData.SubPatchFaces[i];

                // if (PM.NumSides(i) == 4)
                if ((SplitPatchSet.size() != 0) && (SplitPatchSet.count(i) == 0))
                    continue;

                std::vector<int> NewCornerIdX;
                if (CentralV[i] != -1)
                {
                    SplitPatchWith(PatchVertPos, PatchFaces,
                                   PatchVertUV[i], CentralV[i],
                                   SplEdgesIDX[i], NewCornerIdX);
                    NewCornerIdX.push_back(CentralV[i]);
                }

                std::vector<Geo::Point3<ScalarType>> NewCornerPos;
                for (size_t i = 0; i < NewCornerIdX.size(); i++)
                {
                    int IndexV = NewCornerIdX[i];
                    Geo::Point3<ScalarType> CornerPos = PatchVertPos[IndexV];
                    NewCornerPos.push_back(CornerPos);
                }

                std::vector<Geo::Point3<ScalarType>> PatchCornerPos;
                PM.GetAllCornerPos(i, PatchCornerPos);
                std::set<Geo::Point3<ScalarType>> PatchCornerSet(NewCornerPos.begin(),
                                                                 NewCornerPos.end());
                PatchCornerSet.insert(PatchCornerPos.begin(), PatchCornerPos.end());

                // MAP VERTICES AND FACES TO COMPONENTS AND CORNERS
                //  then get connected components
                // std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshesV;
                // std::vector<std::vector<std::vector<int>>> SubMeshesF;
                // std::vector<std::vector<int>> OriginalV, OriginalF;
                // SplitByFFConnectedComponenent(PatchVertPos, PatchFaces,
                //                               SubMeshesV, SubMeshesF,
                //                               OriginalV, OriginalF);

                // if (find_corners)
                // {
                //     for (size_t c = 0; c < SubMeshesF.size(); c++)
                //     {
                //         int InitialVNum = DualVertPos.size();
                //         int InitialFNum = DualFaces.size();

                //         PerPatchCornerIDX.resize(PerPatchCornerIDX.size() + 1);

                //         // std::cout<<"A"<<std::endl;

                //         ExtractPatchCorners(SubMeshesV[c], SubMeshesF[c],
                //                             InitialVNum, PatchCornerSet,
                //                             PerPatchCornerIDX.back());

                //         // std::cout<<"B"<<std::endl;

                //         // MakeConsistent(SubMeshesV[c], SubMeshesF[c], AllSplEdgesPos);

                //         // std::vector<int> NewVidx(SubMeshesV[c].size());
                //         // std::vector<int> NewFidx(SubMeshesF[c].size());

                //         // // Fill the vector with sequential integers starting from previous idx
                //         // std::iota(NewVidx.begin(), NewVidx.end(), InitialVNum);
                //         // std::iota(NewFidx.begin(), NewFidx.end(), InitialFNum);

                //         // PerPatchVertIDX.push_back(NewVidx);
                //         // PerPatchFaceIDX.push_back(NewFidx);

                //         // AppendMesh(DualVertPos, DualFaces, SubMeshesV[c], SubMeshesF[c]);
                //     }
                // }
                // std::cout<<"A"<<std::endl;
                // for (size_t c = 0; c < SubMeshesF.size(); c++)
                // {
                //     // int InitialVNum = DualVertPos.size();
                //     // int InitialFNum = DualFaces.size();

                //     //  std::cout << "De" << std::endl;
                //     // WriteOBJ("testMesh.obj", SubMeshesV[c], SubMeshesF[c]);
                //     // MakeConsistent(SubMeshesV[c], SubMeshesF[c], AllSplEdgesPos);

                //     // std::vector<int> NewVidx(SubMeshesV[c].size());
                //     // std::vector<int> NewFidx(SubMeshesF[c].size());

                //     // // Fill the vector with sequential integers starting from previous idx
                //     // std::iota(NewVidx.begin(), NewVidx.end(), InitialVNum);
                //     // std::iota(NewFidx.begin(), NewFidx.end(), InitialFNum);

                //     // PerPatchVertIDX.push_back(NewVidx);
                //     // PerPatchFaceIDX.push_back(NewFidx);

                //     AppendMesh(DualVertPos, DualFaces, SubMeshesV[c], SubMeshesF[c]);
                // }
                // std::cout<<"B"<<std::endl;
                // std::vector<std::vector<int>> PerPatchFaceIDX;
                // std::vector<std::vector<int>> PerPatchFaceVert;

                // MakeConsistent(PatchVertPos, PatchFaces, AllSplEdgesPos);
                AppendMesh(DualVertPos, DualFaces, PatchVertPos, PatchFaces);
            }

            if (DualFaces.size() > 0)
            {
                MakeConsistent(DualVertPos, DualFaces, AllSplEdgesPos);
                RemoveUnfererencedVert(DualVertPos, DualFaces);
            }
            auto t2 = std::chrono::high_resolution_clock::now();

            typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam SmParam;
            SmParam.Steps = DualParam.SmoothSteps;
            if ((DualFaces.size() > 0) && (DualParam.SmoothSteps > 0))
                Geo::Local_Param_Smooth<ScalarType>::SmoothWithInternalBorders(DualVertPos, DualFaces, true, SmParam);

            auto t3 = std::chrono::high_resolution_clock::now();

            if (write_timing)
            {
                auto duration0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
                auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
                auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);

                std::cout << "Time Deriving Pos: " << duration0.count() << " milliseconds" << std::endl;
                std::cout << "Time Cutting Mesh: " << duration1.count() << " milliseconds" << std::endl;
                std::cout << "Time Smoothing: " << duration3.count() << " milliseconds" << std::endl;
            }
            //
            // std::cout << "DONE" << std::endl;
            return true;
        }
    };

}
#endif
