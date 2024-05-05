#ifndef DUAL_PATCHES
#define DUAL_PATCHES

#include <vector>
#include <UV_param.h>
#include <UV_utils.h>
#include "mesh_functions.h"
#include "mesh_refine.h"
#include "patch_managing.h"

namespace Geo
{

    template <class ScalarType>
    class DualPatch
    {

        static bool GetCornersSplitAlpha(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                         std::vector<std::vector<int>> &Faces,
                                         const std::vector<Geo::Point2<ScalarType>> &VertUV,
                                         const std::vector<int> &Corners,
                                         const std::vector<ScalarType> &SideAlpha,
                                         const int &CentralV)
        {
            assert(CentralV >= 0);
            assert(CentralV < VertUV.size());
            assert(CentralV < VertPos.size());
        }

    public:
        //     static bool ExtractDual(std::vector<Geo::Point3<ScalarType > > &VertPos,
        //                             std::vector<std::vector<int> > &Connectivity,
        //                             PatchManaging<ScalarType> &PM,
        //                             std::vector<Geo::Point3<ScalarType> > &SidesMidP)
        //     {
        //         SidesMidP.clear();

        //         //first get all corners
        //         std::set<int> CornerSet;
        //         for (int i=0;i<PM.NumPatches();i++)
        //         {
        //             std::vector<int> CornerVert;
        //             std::vector<CornerType> CurrType;
        //             PM.GetAllCorners(i,CornerVert,CurrType);
        //             PM.MapToOriginal(i,CornerVert);
        //             CornerSet.insert(CornerVert.begin(),CornerVert.end());
        //         }

        //         //compute mid point of edge
        //         std::set<std::pair<int,int> > AddedMidP;

        //         //std::set<EdgeInt> InterPEdge;

        //         std::map<std::pair<int,int> , ScalarType> EdgeAlpha;

        //         assert(PM.PData.SubPatchFaces.size()==PM.PData.SubPatchVertPos.size());
        //         for (size_t i=0;i<PM.PData.SubPatchFaces.size();i++)
        //         {
        //             //avoid empty patches
        //             if (PM.PData.SubPatchFaces[i].size()==0)continue;

        //             std::vector<std::pair<int,int> > CornersExtremes;
        //             //std::vector<EdgeInt> CurrMidP;

        //             std::vector<std::pair<int,int> > EdgeInt;
        //             std::vector<ScalarType> Alpha;
        //             bool foundMidP=PM.GetSideMidpoints(i,CornerSet,CornersExtremes,EdgeInt,Alpha);
        //             if (!foundMidP)return false;

        //             assert(CornersExtremes.size()==EdgeInt.size());
        //             assert(CornersExtremes.size()==Alpha.size());

        //             for (size_t j=0;j<CornersExtremes.size();j++)
        //             {
        //                 int IndexC0=CornersExtremes[j].first;
        //                 int IndexC1=CornersExtremes[j].second;
        //                 std::pair<int,int> Entry(std::min(IndexC0,IndexC1),std::max(IndexC0,IndexC1));
        //                 if (AddedMidP.count(Entry)>0)continue;
        //                 AddedMidP.insert(Entry);

        //                 int IndexV0=EdgeInt[j].first;
        //                 int IndexV1=EdgeInt[j].second;
        //                 assert(IndexV0<IndexV1);
        //                 ScalarType CurrAlpha=Alpha[j];
        //                 assert(IndexV0>=0);
        //                 assert(IndexV1>=0);
        //                 assert(IndexV0<VertPos.size());
        //                 assert(IndexV1<VertPos.size());
        //                 EdgeAlpha[EdgeInt[j]]=CurrAlpha;

        //                 Geo::Point3<ScalarType> P0=VertPos[IndexV0];
        //                 Geo::Point3<ScalarType> P1=VertPos[IndexV1];
        //                 Geo::Point3<ScalarType> NewPos= (P0*CurrAlpha+P1*(1-CurrAlpha));

        //                 SidesMidP.push_back(NewPos);
        //             }
        //         }

        //         std::vector<Geo::Point3<ScalarType> > SaveMPos;
        //         std::vector<std::vector<int> > SaveMFaces;

        //         //compute tutte embedding and split faces
        //         for (size_t i=0;i<PM.PData.SubPatchFaces.size();i++)
        //         {
        //             //avoid empty patches
        //             if (PM.PData.SubPatchFaces[i].size()==0)continue;

        //             //compute tutte's embedding
        //             std::vector<Geo::Point2<ScalarType> > VertUV;
        //             bool done_param=Geo::UVParam<ScalarType>::TutteEmbedding(PM.PData.SubPatchVertPos[i],
        //                                                                      PM.PData.SubPatchFaces[i],VertUV);
        //             assert(done_param);

        //             //get the split segments
        //             //std::vector<Geo::Segment2<ScalarType> > SplitSeg;
        //             std::vector<std::vector<Geo::Point2<ScalarType> > > SplitPoly;
        //             for (size_t j=0;j<PM.PData.SubPatchFaces[i].size();j++)
        //             {
        //                 for (size_t k=0;k<3;k++)
        //                 {
        //                     int IndexV0=PM.PData.SubPatchFaces[i][j][k];
        //                     int IndexV1=PM.PData.SubPatchFaces[i][j][(k+1)%3];

        //                     Geo::Point2<ScalarType> UV0=VertUV[IndexV0];
        //                     Geo::Point2<ScalarType> UV1=VertUV[IndexV1];

        //                     IndexV0=PM.MapToOriginal(i,IndexV0);
        //                     IndexV1=PM.MapToOriginal(i,IndexV1);

        //                     std::pair<int,int> Key(std::min(IndexV0,IndexV1),
        //                                            std::max(IndexV0,IndexV1));

        //                     if (EdgeAlpha.count(Key)==0)continue;
        //                     ScalarType alpha=EdgeAlpha[Key];

        //                     Geo::Point2<ScalarType> InterpPUV=UV0*alpha+UV1*(1-alpha);
        //                     //inflate a bit
        //                     InterpPUV=InterpPUV*1.01;
        //                     SplitPoly.resize(SplitPoly.size()+1);
        //                     SplitPoly.back().push_back(Geo::Point2<ScalarType>(0,0));
        //                     SplitPoly.back().push_back(InterpPUV);
        // //                    Geo::Segment2<ScalarType> CurrSeg(Geo::Point2<ScalarType>(0,0),InterpPUV);
        // //                    SplitSeg.push_back(CurrSeg);
        //                 }
        //             }
        //             //std::cout<<"Size split:"<<SplitP.size()<<std::endl;

        // //            std::vector<Geo::Point3<ScalarType> > TestPos;
        // //            for (size_t i=0;i<VertUV.size();i++)
        // //            {
        // //                Geo::Point2<ScalarType> uvP=VertUV[i];
        // //                TestPos.push_back(Geo::Point3<ScalarType>(uvP.X,uvP.Y,0));
        // //            }
        //             std::vector<std::pair<int,int> > NewEdges;
        //             //SplitMeshWithUVPolyline(TestPos,VertUV,SubPatchFaces[i],SplitSeg,NewEdges);
        //             RefineMesh<ScalarType>::SplitMeshWithUVPolylines(PM.PData.SubPatchVertPos[i],VertUV,
        //                                                             PM.PData.SubPatchFaces[i],SplitPoly,NewEdges);

        //             std::vector<std::vector<int> > NextFP;
        //             std::vector<std::vector<int> > NextEP;
        //             getFFAdjacency(PM.PData.SubPatchFaces[i],NextFP,NextEP);

        //             SplitMeshFromEdges(PM.PData.SubPatchVertPos[i],PM.PData.SubPatchFaces[i],NextFP,NextEP,NewEdges);

        //             AppendMesh(SaveMPos,SaveMFaces,PM.PData.SubPatchVertPos[i],PM.PData.SubPatchFaces[i]);

        // //            std::string NameF= "Test_split" + std::to_string(i) + ".obj";
        // //            WriteOBJ(NameF,SubPatchVertPos[i],SubPatchFaces[i]);

        // //            SplitMeshFromEdges(TestPos,SubPatchFaces[i],NextFP,NextEP,NewEdges);
        // //            std::string NameF= "Test_split" + std::to_string(i) + ".obj";
        // //            WriteOBJ(NameF,TestPos,SubPatchFaces[i]);

        // //            std::vector<std::pair<int,int> > NewEdges;
        // //            SplitMeshWithUVPolyline(SubPatchVertPos[i],VertUV,SubPatchFaces[i],
        // //                                    SplitSeg,NewEdges);

        // //            std::vector<std::vector<int> > NextFP;
        // //            std::vector<std::vector<int> > NextEP;
        // //            getFFAdjacency(SubPatchFaces[i],NextFP,NextEP);

        // //            SplitMeshFromEdges(SubPatchVertPos[i],SubPatchFaces[i],NextFP,NextEP,NewEdges);
        // //            std::string NameF= "Test_split" + std::to_string(i) + ".obj";
        // //            WriteOBJ(NameF,SubPatchVertPos[i],SubPatchFaces[i]);
        //         }
        //         WriteOBJ("./dualized.obj",SaveMPos,SaveMFaces);
        //         return true;
        //     }

        static bool SubdividePatch(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                   const std::vector<std::vector<int>> &Faces,
                                   const std::vector<int> &Corners,
                                   const std::vector<int> &SideAlpha,
                                   const int &CentralV,
                                   std::vector<Geo::Point3<ScalarType>> &SplitVertPos,
                                   std::vector<std::vector<int>> &SplitFaces,
                                   std::vector<std::pair<int, int>> &SplitEdges,
                                   std::vector<ScalarType> &SplitAlpha)
        {
            SplitVertPos.clear();
            SplitFaces.clear();
            SplitEdges.clear();
            SplitAlpha.clear();

            // compute tutte's embedding
            std::vector<Geo::Point2<ScalarType>> VertUV;
            bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(VertPos, Faces, VertUV);
            assert(done_param);
            if (CentralV==-1)
                CentralV=Geo::UVUtils<ScalarType>::ClosestTo(VertUV,Geo::Point2<ScalarType>(0,0));
            
            assert(CentralV >= 0);
            assert(CentralV < VertUV.size());
            assert(CentralV < VertPos.size());
            Geo::Point2<ScalarType> CentralUV = VertUV[CentralV];
            // get the split segments
            // std::vector<Geo::Segment2<ScalarType> > SplitSeg;
            std::vector<std::vector<Geo::Point2<ScalarType>>> SplitPoly;
            for (size_t i = 0; i < Corners.size(); i++)
            {
                int IndexV0 = Corners[i] int IndexV1 = Corners[(i + 1) * Corners.size()];
                assert(IndexV0 >= 0);
                assert(IndexV1 >= 0);
                assert(IndexV0 < VertUV.size());
                assert(IndexV0 < VertPos.size());
                assert(IndexV1 < VertUV.size());
                assert(IndexV1 < VertPos.size());

                Geo::Point2<ScalarType> UV0 = VertUV[IndexV0];
                Geo::Point2<ScalarType> UV1 = VertUV[IndexV1];

                ScalarType Alpha = SplitAlpha[i]

                    Geo::Point2<ScalarType>
                        InterpPUV = UV0 * alpha + UV1 * (1 - alpha);

                // inflate a bit
                Geo::Point2<ScalarType> DirCut = InterpPUV - CentralUV;
                ScalarType NormDist = (DirCut).Norm();
                DirCut.Normalize();
                NormDist *= 1.01;
                InterpPUV = CentralUV + DirCut * NormDist;
                SplitPoly.resize(SplitPoly.size() + 1);
                SplitPoly.back().push_back(CentralUV);
                SplitPoly.back().push_back(InterpPUV);
            }

            std::vector<std::pair<int, int>> NewEdges;
            SplitVertPos=VertPos;
            SplitFaces=Faces;
            RefineMesh<ScalarType>::SplitMeshWithUVPolylines(SplitVertPos, VertUV,
                                                             SplitFaces, SplitPoly, 
                                                             NewEdges);

            std::vector<std::vector<int>> NextFP;
            std::vector<std::vector<int>> NextEP;
            getFFAdjacency(PatchFaces, NextFP, NextEP);

            SplitMeshFromEdges(SplitVertPos, SplitVertPos, NextFP, NextEP, NewEdges);
        }

        struct DualParameter
        {
            // input
            std::vector<int> OnlyPatch;
            bool OnlyNonQuad;
            bool UseSing;
            std::vector<int> SingIndex;

            // output
            std::vector<Geo::Point3<ScalarType>> SidesMidP;

            DualParameter()
            {
                OnlyNonQuad = false;
                UseSing = false;
            }
        };

        // static bool ExtractDual(const std::vector<Geo::Point3<ScalarType>> &VertPos,
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
        //     DualConnectivity.clear();
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
        //             AppendMesh(DualVertPos, DualConnectivity, PatchVertPos, PatchFaces);
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

        //         AppendMesh(DualVertPos, DualConnectivity, PatchVertPos, PatchFaces);
        //     }
        //     // WriteOBJ("./dualized.obj", SaveMPos, SaveMFaces);
        //     return true;
        // }

       

        static bool ExtractDual(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                const std::vector<std::vector<int>> &Faces,
                                PatchManaging<ScalarType> &PM,
                                DualParameter &DualParam,
                                std::vector<Geo::Point3<ScalarType>> &DualVertPos,
                                std::vector<std::vector<int>> &DualFaces)
        {
            DualParam.SidesMidP.clear();

            // first get all corners
            std::set<int> CornerSet;
            for (int i = 0; i < PM.NumPatches(); i++)
            {
                std::vector<int> CornerVert;
                std::vector<CornerType> CurrType;
                PM.GetAllCorners(i, CornerVert, CurrType);
                PM.MapToOriginal(i, CornerVert);
                CornerSet.insert(CornerVert.begin(), CornerVert.end());
            }

            // compute mid point of edge
            std::set<std::pair<int, int>> AddedMidP;

            // std::set<EdgeInt> InterPEdge;

            std::map<std::pair<int, int>, ScalarType> EdgeAlpha;

            assert(PM.PData.SubPatchFaces.size() == PM.PData.SubPatchVertPos.size());
            for (size_t i = 0; i < PM.PData.SubPatchFaces.size(); i++)
            {
                // avoid empty patches
                if (PM.PData.SubPatchFaces[i].size() == 0)
                    continue;

                std::vector<std::pair<int, int>> CornersExtremes;
                // std::vector<EdgeInt> CurrMidP;

                std::vector<std::pair<int, int>> EdgeInt;
                std::vector<ScalarType> Alpha;
                bool foundMidP = PM.GetSideMidpoints(i, CornerSet, CornersExtremes, EdgeInt, Alpha);
                if (!foundMidP)
                    return false;

                assert(CornersExtremes.size() == EdgeInt.size());
                assert(CornersExtremes.size() == Alpha.size());

                for (size_t j = 0; j < CornersExtremes.size(); j++)
                {
                    int IndexC0 = CornersExtremes[j].first;
                    int IndexC1 = CornersExtremes[j].second;
                    std::pair<int, int> Entry(std::min(IndexC0, IndexC1), std::max(IndexC0, IndexC1));
                    if (AddedMidP.count(Entry) > 0)
                        continue;
                    AddedMidP.insert(Entry);

                    int IndexV0 = EdgeInt[j].first;
                    int IndexV1 = EdgeInt[j].second;
                    assert(IndexV0 < IndexV1);
                    ScalarType CurrAlpha = Alpha[j];
                    assert(IndexV0 >= 0);
                    assert(IndexV1 >= 0);
                    assert(IndexV0 < VertPos.size());
                    assert(IndexV1 < VertPos.size());
                    EdgeAlpha[EdgeInt[j]] = CurrAlpha;

                    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                    Geo::Point3<ScalarType> NewPos = (P0 * CurrAlpha + P1 * (1 - CurrAlpha));

                    DualParam.SidesMidP.push_back(NewPos);
                }
            }

            DualVertPos.clear();
            DualConnectivity.clear();
            // std::vector<Geo::Point3<ScalarType>> SaveMPos;
            // std::vector<std::vector<int>> SaveMFaces;

            // compute tutte embedding and split faces
            for (size_t i = 0; i < PM.PData.SubPatchFaces.size(); i++)
            {
                // avoid empty patches
                if (PM.PData.SubPatchFaces[i].size() == 0)
                    continue;
                // avoid patches with multiple boundaries
                if (PM.PData.PatchCorners[i].size() != 1)
                    continue;
                ;

                std::vector<Geo::Point3<ScalarType>> PatchVertPos = PM.PData.SubPatchVertPos[i];
                std::vector<std::vector<int>> PatchFaces = PM.PData.SubPatchFaces[i];

                if ((DualParam.OnlyNonQuad) && (PM.NumCorners(i) == 4))
                {
                    AppendMesh(DualVertPos, DualConnectivity, PatchVertPos, PatchFaces);
                    continue;
                }
                // compute tutte's embedding
                std::vector<Geo::Point2<ScalarType>> VertUV;
                bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(PatchVertPos, PatchFaces, VertUV);
                assert(done_param);

                // get the split segments
                // std::vector<Geo::Segment2<ScalarType> > SplitSeg;
                std::vector<std::vector<Geo::Point2<ScalarType>>> SplitPoly;
                for (size_t j = 0; j < PatchFaces.size(); j++)
                {
                    for (size_t k = 0; k < 3; k++)
                    {
                        int IndexV0 = PatchFaces[j][k];
                        int IndexV1 = PatchFaces[j][(k + 1) % 3];

                        Geo::Point2<ScalarType> UV0 = VertUV[IndexV0];
                        Geo::Point2<ScalarType> UV1 = VertUV[IndexV1];

                        IndexV0 = PM.MapToOriginal(i, IndexV0);
                        IndexV1 = PM.MapToOriginal(i, IndexV1);

                        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                                std::max(IndexV0, IndexV1));

                        if (EdgeAlpha.count(Key) == 0)
                            continue;

                        ScalarType alpha = EdgeAlpha[Key];

                        Geo::Point2<ScalarType> InterpPUV = UV0 * alpha + UV1 * (1 - alpha);
                        // inflate a bit
                        InterpPUV = InterpPUV * 1.01;
                        SplitPoly.resize(SplitPoly.size() + 1);
                        SplitPoly.back().push_back(Geo::Point2<ScalarType>(0, 0));
                        SplitPoly.back().push_back(InterpPUV);
                    }
                }

                std::vector<std::pair<int, int>> NewEdges;
                // SplitMeshWithUVPolyline(TestPos,VertUV,SubPatchFaces[i],SplitSeg,NewEdges);
                RefineMesh<ScalarType>::SplitMeshWithUVPolylines(PatchVertPos, VertUV,
                                                                 PatchFaces, SplitPoly, NewEdges);

                std::vector<std::vector<int>> NextFP;
                std::vector<std::vector<int>> NextEP;
                getFFAdjacency(PatchFaces, NextFP, NextEP);

                SplitMeshFromEdges(PatchVertPos, PatchFaces, NextFP, NextEP, NewEdges);

                AppendMesh(DualVertPos, DualConnectivity, PatchVertPos, PatchFaces);
            }
            // WriteOBJ("./dualized.obj", SaveMPos, SaveMFaces);
            return true;
        }
    };

}
#endif
