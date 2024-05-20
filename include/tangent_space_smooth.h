#ifndef TANGENT_SPACE_SMOOTH
#define TANGENT_SPACE_SMOOTH

#include "triangle2.h"
#include "connectivity.h"
#include "point3.h"
#include "pos3.h"
#include "cotangent.h"
#include "triangle3.h"
#include "plane3.h"

namespace Geo
{

    template <class ScalarType>
    class Local_Param_Smooth
    {

    public:
        struct UVSmoothParam
        {
            std::vector<int> Fixed;
            std::vector<int> SmoothOnly;
            std::vector<std::pair<int, int>> FeatureLines;
            ScalarType Damp;
            size_t Steps;

            UVSmoothParam()
            {
                Damp = 0.2;
                Steps = 5;
            }
        };

    private:
        //    static Geo::Point3<ScalarType > Interpolate3DPos(const std::vector<Geo::Traingle2<ScalarType> > &TrisUV,
        //                                                     const std::vector<Geo::Traingle3<ScalarType> > &Tris3D,
        //                                                     const Geo::Point2<ScalarType> &To_Sample)
        //    {
        //        assert(TrisUV.size()>=2);

        //        //try to interpolate inside
        //        for (size_t i=0;i<TrisUV.size();i++)
        //        {
        //            if (!Geo::IsInside(TrisUV[i],To_Sample))continue;
        //            ScalarType alpha,beta,gamma;
        //            Geo::BarycentricCoords(TrisUV[i],To_Sample,alpha,beta,gamma);

        //            Geo::Point3<ScalarType > Ret=Tris3D[i].P(0)*alpha+
        //                    Tris3D[i].P(1)*beta+
        //                    Tris3D[i].P(2)*gamma;
        //            return Ret;
        //        }

        //        //otherwise get the closest one
        //        int closestT=-1;
        //        ScalarType minD=std::numeric_limits<ScalarType>::max();
        //        for (int i=0;i<(int)TrisUV.size();i++)
        //        {
        //            ScalarType currD=Distance(TrisUV[i],To_Sample);
        //            if (currD>minD)continue;
        //            minD=currD;
        //            closestT=i;
        //        }
        //        assert(closestT!=-1);
        //        ScalarType alpha,beta,gamma;
        //        Geo::BarycentricCoords(TrisUV[closestT],To_Sample,alpha,beta,gamma);

        //    }

        static Geo::Point3<ScalarType> Interpolate3DPos(const std::vector<Geo::Traingle2<ScalarType>> &TrisUV,
                                                        const std::vector<Geo::Traingle3<ScalarType>> &Tris3D,
                                                        const Geo::Point2<ScalarType> &To_Sample)
        {
            assert(TrisUV.size() >= 2);

            // try to interpolate inside
            int InterpT = -1;
            for (size_t i = 0; i < TrisUV.size(); i++)
            {
                if (!Geo::IsInside(TrisUV[i], To_Sample))
                    continue;
                InterpT = i;
            }

            // otherwise get the closest one
            if (InterpT == -1)
            {
                ScalarType minD = std::numeric_limits<ScalarType>::max();
                for (int i = 0; i < (int)TrisUV.size(); i++)
                {
                    ScalarType currD = Distance(TrisUV[i], To_Sample);
                    if (currD > minD)
                        continue;
                    minD = currD;
                    InterpT = i;
                }
            }

            assert(InterpT != -1);
            ScalarType alpha, beta, gamma;
            Geo::BarycentricCoords(TrisUV[InterpT], To_Sample, alpha, beta, gamma);
            Geo::Point3<ScalarType> Ret = Tris3D[InterpT].P(0) * alpha +
                                          Tris3D[InterpT].P(1) * beta +
                                          Tris3D[InterpT].P(2) * gamma;
            return Ret;
        }

        static Geo::Point3<ScalarType> ReprojectOnPoliline(const std::vector<Geo::Point3<ScalarType>> &PolyL,
                                                           const Geo::Point3<ScalarType> &TestP)
        {
            ScalarType MinDist = std::numeric_limits<ScalarType>::max();
            Geo::Point3<ScalarType> minP(0, 0, 0);
            for (size_t i = 0; i < PolyL.size(); i++)
            {
                Geo::Point3<ScalarType> P0 = PolyL[i];
                Geo::Point3<ScalarType> P1 = PolyL[(i + 1) % PolyL.size()];
                Geo::Segment3<ScalarType> Seg(P0, P1);
                Geo::Point3<ScalarType> ClosP = Geo::ClosestPoint(Seg, TestP);
                ScalarType Dist = (ClosP - TestP).Norm();
                if (Dist > MinDist)
                    continue;
                MinDist = Dist;
                minP = ClosP;
            }
            return minP;
        }

        //    static void ParametrizeStar(std::vector<Geo::Point3<ScalarType > > &VertPos,
        //                                const std::vector<std::vector<int> > &Faces,
        //                                const std::vector<std::vector<int> > &NextF,
        //                                const std::vector<std::vector<int> > &NextE,
        //                                const std::vector<Pos3> &StarMesh,
        //                                std::vector<Geo::Point2<ScalarType> > &StarUV,
        //                                std::vector<Geo::Traingle2<ScalarType> > &TrisUV,
        //                                std::vector<Geo::Traingle3<ScalarType> > &Tris3D)
        //    {
        //        std::map<int,Geo::Point2<ScalarType> > VertUVsMap;

        //        TrisUV.clear();
        //        Tris3D.clear();
        //        StarUV.clear();

        //        int CentralV=StarMesh[0].V;

        //        ScalarType AngleInterval=2.0*M_PI/(ScalarType)StarMesh.size();
        //        ScalarType CurrAngle=0;

        //        Geo::Point2<ScalarType> CentralUV(0,0);
        //        ScalarType SumW=0;
        //        for (size_t i=0;i<StarMesh.size();i++)
        //        {
        //            int VIdx0=StarMesh[i].V;
        //            assert(VIdx0==CentralV);

        //            //then get the one outside on the star
        //            Pos3 CurrP=StarMesh[i];
        //            CurrP.FlipV(Faces);

        //            int VIdx1=CurrP.V;
        //            assert(VIdx0!=VIdx1);

        //            Geo::Point2<ScalarType> CurrUV(cos(CurrAngle),sin(CurrAngle));
        //            CurrAngle+=AngleInterval;

        //            VertUVsMap[VIdx1]=CurrUV;

        //            StarUV.push_back(CurrUV);

        //            ScalarType CurrW=CotangentWeightTriFace(VertPos,Faces,NextF,NextE,StarMesh[i].F,StarMesh[i].E);
        //            CentralUV+=(CurrUV*CurrW);
        //            SumW+=CurrW;
        //        }

        //        CentralUV/=SumW;
        //        assert(VertUVsMap.count(CentralV)==0);
        //        VertUVsMap[CentralV]=CentralUV;

        //        //then return the parametrization
        //        for (size_t i=0;i<StarMesh.size();i++)
        //        {
        //            int IndexF=StarMesh[i].F;
        //            assert(IndexF>=0);
        //            assert(IndexF<Faces.size());
        //            assert(Faces[IndexF].size()==3);

        //            Geo::Point2<ScalarType> P_UV[3];
        //            Geo::Point3<ScalarType> P_3D[3];

        //            for (size_t j=0;j<Faces[IndexF].size();j++)
        //            {
        //                int IndexV=Faces[IndexF][j];
        //                assert(VertUVsMap.count(IndexV)>0);
        //                P_UV[j]=VertUVsMap[IndexV];
        //                P_3D[j]=VertPos[IndexV];
        //            }
        //            Geo::Traingle2<ScalarType> T_UV(P_UV[0],P_UV[1],P_UV[2]);
        //            Geo::Traingle3<ScalarType> T_3D(P_3D[0],P_3D[1],P_3D[2]);

        //            TrisUV.push_back(T_UV);
        //            Tris3D.push_back(T_3D);
        //        }
        //    }

        static ScalarType AngleRad(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                   const std::vector<std::vector<int>> &Faces,
                                   const Pos3 &Pos0, const Pos3 &Pos1)
        {
            int VIdx0 = Pos0.V;
            Pos3 CurrP = Pos0;
            CurrP.FlipV(Faces);
            int VIdx1 = CurrP.V;
            Geo::Point3<ScalarType> Dir0 = VertPos[VIdx1] - VertPos[VIdx0];
            Dir0.Normalize();

            int VIdx2 = Pos1.V;
            CurrP = Pos1;
            CurrP.FlipV(Faces);
            int VIdx3 = CurrP.V;
            Geo::Point3<ScalarType> Dir1 = VertPos[VIdx3] - VertPos[VIdx2];
            Dir1.Normalize();

            return (Geo::AngleRad(Dir0, Dir1));
        }

        static ScalarType LenghtPos(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                    const std::vector<std::vector<int>> &Faces, const Pos3 &Pos)
        {
            int VIdx0 = Pos.V;
            Pos3 CurrP = Pos;
            CurrP.FlipV(Faces);
            int VIdx1 = CurrP.V;
            return (VertPos[VIdx0] - VertPos[VIdx1]).Norm();
        }

        static void ParametrizeStar(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                    const std::vector<std::vector<int>> &Faces,
                                    const std::vector<std::vector<int>> &NextF,
                                    const std::vector<std::vector<int>> &NextE,
                                    const std::vector<Pos3> &StarMesh,
                                    std::vector<Geo::Point2<ScalarType>> &StarUV,
                                    std::vector<Geo::Traingle2<ScalarType>> &TrisUV,
                                    std::vector<Geo::Traingle3<ScalarType>> &Tris3D)
        {
            std::map<int, Geo::Point2<ScalarType>> VertUVsMap;

            TrisUV.clear();
            Tris3D.clear();
            StarUV.clear();

            int CentralV = StarMesh[0].V;

            ScalarType ToTAngle = 0;
            ScalarType ToTLenght = 0;
            std::vector<ScalarType> AnglePosRad;
            std::vector<ScalarType> LenPos;
            std::vector<int> StarV;
            int SizeStar = StarMesh.size();
            for (size_t i = 0; i < StarMesh.size(); i++)
            {
                Pos3 CurrP = StarMesh[i];
                Pos3 NextP = StarMesh[(i + 1) % SizeStar];

                ScalarType Angle = AngleRad(VertPos, Faces, CurrP, NextP);
                ScalarType Lenght = LenghtPos(VertPos, Faces, CurrP);

                ToTAngle += Angle;
                ToTLenght += Lenght;

                AnglePosRad.push_back(Angle);
                LenPos.push_back(Lenght);

                assert(CurrP.V == CentralV);
                CurrP.FlipV(Faces);
                StarV.push_back(CurrP.V);
            }

            ScalarType AngleScale = 2.0 * M_PI / ToTAngle;
            ScalarType CurrAngle = 0;

            Geo::Point2<ScalarType> CentralUV(0, 0);
            ScalarType SumW = 0;
            for (size_t i = 0; i < StarMesh.size(); i++)
            {
                Geo::Point2<ScalarType> CurrUV(cos(CurrAngle), sin(CurrAngle));
                CurrUV = CurrUV * LenPos[i];
                CurrAngle += (AnglePosRad[i] * AngleScale);

                VertUVsMap[StarV[i]] = CurrUV;

                StarUV.push_back(CurrUV);

                ScalarType CurrW = CotangentWeightTriFace(VertPos, Faces, NextF, NextE, StarMesh[i].F, StarMesh[i].E);
                CentralUV += (CurrUV * CurrW);
                SumW += CurrW;
            }
            CentralUV /= SumW;
            // std::cout<<"Central UV "<<CentralUV.X<<","<<CentralUV.Y<<std::endl;
            assert(VertUVsMap.count(CentralV) == 0);
            VertUVsMap[CentralV] = CentralUV;

            // then return the parametrization
            for (size_t i = 0; i < StarMesh.size(); i++)
            {
                int IndexF = StarMesh[i].F;
                assert(IndexF >= 0);
                assert(IndexF < Faces.size());
                assert(Faces[IndexF].size() == 3);

                Geo::Point2<ScalarType> P_UV[3];
                Geo::Point3<ScalarType> P_3D[3];

                for (size_t j = 0; j < Faces[IndexF].size(); j++)
                {
                    int IndexV = Faces[IndexF][j];
                    assert(VertUVsMap.count(IndexV) > 0);
                    P_UV[j] = VertUVsMap[IndexV];
                    P_3D[j] = VertPos[IndexV];
                }
                Geo::Traingle2<ScalarType> T_UV(P_UV[0], P_UV[1], P_UV[2]);
                Geo::Traingle3<ScalarType> T_3D(P_3D[0], P_3D[1], P_3D[2]);

                TrisUV.push_back(T_UV);
                Tris3D.push_back(T_3D);
            }
        }

        static void SmoothStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                               const std::vector<std::vector<int>> &Faces,
                               const std::vector<std::vector<int>> &NextF,
                               const std::vector<std::vector<int>> &NextE,
                               const std::vector<std::vector<Pos3>> &StarMeshes,
                               const std::vector<bool> &IsB,
                               std::vector<std::vector<int>> BorderVertSeq,
                               UVSmoothParam UVP,
                               bool do_features,
                               bool do_nonfeatures)
        {
            std::vector<bool> IsFeature(VertPos.size(), false);
            for (size_t i = 0; i < BorderVertSeq.size(); i++)
                for (size_t j = 0; j < BorderVertSeq[i].size(); j++)
                    IsFeature[BorderVertSeq[i][j]] = true;

            for (size_t i = 0; i < UVP.FeatureLines.size(); i++)
            {
                int IndexV0 = UVP.FeatureLines[i].first;
                int IndexV1 = UVP.FeatureLines[i].second;
                IsFeature[IndexV0] = true;
                IsFeature[IndexV1] = true;
            }

            std::set<int> SmoothOnlySet(UVP.SmoothOnly.begin(), UVP.SmoothOnly.end());
            std::set<int> FixedSet(UVP.Fixed.begin(), UVP.Fixed.end());
            std::set<std::pair<int, int>> FeatureSet(UVP.FeatureLines.begin(),
                                                     UVP.FeatureLines.end());

            std::vector<Geo::Point3<ScalarType>> TargetP(VertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));
            std::vector<size_t> TargetN(VertPos.size(), 0);

            std::vector<std::vector<Geo::Point3<ScalarType>>> ReprojectPolys;

            if (do_features)
            {
                ReprojectPolys.resize(BorderVertSeq.size());
                for (size_t i = 0; i < BorderVertSeq.size(); i++)
                {
                    for (size_t j = 0; j < BorderVertSeq[i].size(); j++)
                    {
                        Geo::Point3<ScalarType> P = VertPos[BorderVertSeq[i][j]];
                        ReprojectPolys[i].push_back(P);
                    }
                }

                // calculate target positions for border
                for (size_t i = 0; i < BorderVertSeq.size(); i++)
                {
                    int SizeB = BorderVertSeq[i].size();
                    assert(SizeB > 2);

                    if (SizeB <= 3)
                        continue; // too small border sequence to be smoothed

                    for (size_t j = 0; j < BorderVertSeq[i].size(); j++)
                    {
                        int LocIndex = j;
                        int LocIndex0 = (j + (SizeB - 1)) % SizeB;
                        int LocIndex1 = (j + 1) % SizeB;

                        int IndexV = BorderVertSeq[i][LocIndex];
                        int IndexV0 = BorderVertSeq[i][LocIndex0];
                        int IndexV1 = BorderVertSeq[i][LocIndex1];

                        TargetP[IndexV] += VertPos[IndexV0];
                        TargetP[IndexV] += VertPos[IndexV1];
                        TargetN[IndexV]++;
                        TargetN[IndexV]++;
                    }
                }
            }

            for (size_t i = 0; i < StarMeshes.size(); i++)
            {
                if (IsB[i])
                    continue;

                assert(StarMeshes[i].size() >= 2);

                if ((!do_features) && (IsFeature[i]))
                    continue;
                if ((!do_nonfeatures) && (!IsFeature[i]))
                    continue;
                if (FixedSet.count(i) > 0)
                    continue;
                if ((SmoothOnlySet.size() > 0) && (SmoothOnlySet.count(i) == 0))
                    continue;

                std::vector<Geo::Traingle2<ScalarType>> TrisUV;
                std::vector<Geo::Traingle3<ScalarType>> Tris3D;
                std::vector<Geo::Point2<ScalarType>> StarUV;
                ParametrizeStar(VertPos, Faces, NextF,
                                NextE, StarMeshes[i],
                                StarUV, TrisUV, Tris3D);

                // check the point to interpolata, might be the center
                Geo::Point2<ScalarType> CenterP(0, 0);
                int CenterNum = 0;
                // check if there are features
                for (size_t j = 0; j < StarMeshes[i].size(); j++)
                {
                    Pos3 CurrPos = StarMeshes[i][j];
                    int IndexV0 = CurrPos.V;
                    assert(IndexV0 == i); // check star centered on V
                    CurrPos.FlipV(Faces);
                    int IndexV1 = CurrPos.V;
                    std::pair<int, int> Key(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
                    if ((IsFeature[i]) && FeatureSet.count(Key) == 0)
                        continue;
                    CenterP += StarUV[j];
                    CenterNum++;
                }
                if (CenterNum == 1) // only one feature should fix it
                {
                    FixedSet.insert(i);
                    CenterP = Geo::Point2<ScalarType>(0, 0);
                }
                if (CenterNum > 1)
                {
                    CenterP /= CenterNum;
                }

                TargetP[i] = Interpolate3DPos(TrisUV, Tris3D, CenterP);
                TargetN[i]++;
            }

            for (size_t i = 0; i < TargetP.size(); i++)
            {
                if (TargetN[i] > 0)
                    TargetP[i] /= (ScalarType)TargetN[i];
            }

            // reproject borders
            if (do_features)
            {
                for (size_t i = 0; i < BorderVertSeq.size(); i++)
                {
                    for (size_t j = 0; j < BorderVertSeq[i].size(); j++)
                    {
                        int indexV = BorderVertSeq[i][j];
                        Geo::Point3<ScalarType> P = VertPos[indexV];
                        TargetP[indexV] = ReprojectOnPoliline(ReprojectPolys[i], P);
                    }
                }
            }

            // update mesh position
            for (size_t i = 0; i < VertPos.size(); i++)
            {
                if ((!do_features) && (IsFeature[i]))
                    continue;
                if ((!do_nonfeatures) && (!IsFeature[i]))
                    continue;
                if (FixedSet.count(i))
                    continue;
                if ((SmoothOnlySet.size() > 0) && (SmoothOnlySet.count(i) == 0))
                    continue;

                VertPos[i] = VertPos[i] * UVP.Damp + TargetP[i] * (1 - UVP.Damp);
            }
        }

        static void GetStarMeshes(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                  const std::vector<std::vector<int>> &Faces,
                                  const std::vector<std::vector<int>> &NextF,
                                  const std::vector<std::vector<int>> &NextE,
                                  std::vector<bool> IsBorder,
                                  std::vector<std::vector<Pos3>> &StarMeshes)
        {

            StarMeshes.clear();
            StarMeshes.resize(VertPos.size());
            for (size_t i = 0; i < Faces.size(); i++)
                for (size_t j = 0; j < 3; j++)
                {

                    // then get the star of vertices
                    Pos3 TestPos(i, j, Faces[i][j]);
                    int CenterV = TestPos.V;

                    // check if already initialized or not
                    if (StarMeshes[CenterV].size() > 0)
                        continue;
                    // border do nothing
                    if (IsBorder[CenterV])
                        continue;

                    std::vector<Pos3> VertFan;
                    bool HasBorder;
                    GetSortedFacesFaceFun(TestPos, Faces, NextF, NextE, VertFan, HasBorder);
                    assert(!HasBorder);

                    if (VertFan.size() == 0)
                        continue;
                    StarMeshes[CenterV] = VertFan;
                }
        }

    public:
        static void Smooth(std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Faces,
                           UVSmoothParam UVP = UVSmoothParam())
        {

            // initialize constraints
            std::set<std::pair<int, int>> AllFeatures;
            for (size_t i = 0; i < UVP.FeatureLines.size(); i++)
            {
                int IndexV0 = UVP.FeatureLines[i].first;
                int IndexV1 = UVP.FeatureLines[i].second;
                std::pair<int, int> entry(std::min(IndexV0, IndexV1),
                                          std::max(IndexV0, IndexV1));
                AllFeatures.insert(entry);
            }

            // then add borders
            std::vector<std::vector<int>> NextF, NextE;
            Geo::getFFAdjacency(Faces, NextF, NextE);

            for (size_t i = 0; i < NextF.size(); i++)
            {
                assert(NextF[i].size() == 3); // only triangles
                for (size_t j = 0; j < NextF[i].size(); j++)
                {
                    if (NextF[i][j] != -1)
                        continue;
                    int IndexV0 = Faces[i][j];
                    int IndexV1 = Faces[i][(j + 1) % NextF[i].size()];
                    std::pair<int, int> entry(std::min(IndexV0, IndexV1),
                                              std::max(IndexV0, IndexV1));
                    AllFeatures.insert(entry);
                }
            }

            // std::vector<std::vector<PosType> > ProjBasis;
            // GetStarMeshes(mesh,StarMeshes,ProjBasis,OrigIndex);

            std::vector<bool> IsBorder;
            Geo::getBorderVert(Faces, NextF, VertPos.size(), IsBorder);
            std::vector<std::vector<int>> BorderVertSeq;
            getBorderSequncesVert(Faces, NextF, NextE, BorderVertSeq);

            // create star meshes
            std::vector<std::vector<Pos3>> StarMeshes;
            GetStarMeshes(VertPos, Faces, NextF, NextE, IsBorder, StarMeshes);

            // and finally smooth
            for (size_t i = 0; i < UVP.Steps; i++)
            {
                // first smooth features
                SmoothStep(VertPos, Faces, NextF, NextE, StarMeshes, IsBorder, BorderVertSeq, UVP, true, true);
                // then the rest
                //             SmoothStep(VertPos,Faces,NextF,NextE,StarMeshes,IsBorder,BorderVertSeq,UVP,false,true);
            }
        }

        static void SmoothWithInternalBorders(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                              const std::vector<std::vector<int>> &Faces,
                                              bool UseCloseDelta,
                                              UVSmoothParam UVP = UVSmoothParam())
        {
            //auto t0 = std::chrono::high_resolution_clock::now();

            std::vector<int> IndexF;
            std::vector<int> IndexE;
            GetInternalBorders(VertPos, Faces, IndexF, IndexE, UseCloseDelta);

            std::vector<Geo::Point3<ScalarType>> VertMerged = VertPos;
            std::vector<std::vector<int>> FacesMerged = Faces;


            ScalarType delta = std::numeric_limits<ScalarType>::epsilon() * 2;
            if (UseCloseDelta)
                MergeCloseVert(VertMerged, FacesMerged, delta);
            else
                RemoveDuplicatedVert(VertMerged, FacesMerged);

            assert(FacesMerged.size() == Faces.size());

            //auto t1 = std::chrono::high_resolution_clock::now();

            std::vector<std::pair<int, int>> FeatureE;
            for (size_t i = 0; i < IndexF.size(); i++)
            {
                int CurrF = IndexF[i];
                int CurrE = IndexE[i];
                int IndexV0 = FacesMerged[CurrF][CurrE];
                int IndexV1 = FacesMerged[CurrF][(CurrE + 1) % 3];
                std::pair<int, int> KeyF(std::min(IndexV0, IndexV1),
                                         std::max(IndexV0, IndexV1));
                FeatureE.push_back(KeyF);
            }
            std::sort(FeatureE.begin(), FeatureE.end());
            FeatureE.erase(std::unique(FeatureE.begin(), FeatureE.end()), FeatureE.end());
            //UVP.FeatureLines = FeatureE;
            UVP.FeatureLines.insert(UVP.FeatureLines.end(),FeatureE.begin(),FeatureE.end());
            Smooth(VertMerged, FacesMerged, UVP);

            // auto t2 = std::chrono::high_resolution_clock::now();

            // auto duration0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
            // auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            // std::cout<<"Time 0 "<<duration0.count()<<std::endl;
            // std::cout<<"Time 1 "<<duration1.count()<<std::endl;
            // copy back values
            for (size_t i = 0; i < Faces.size(); i++)
                for (size_t j = 0; j < Faces[i].size(); j++)
                {
                    int IndexVOrig = Faces[i][j];
                    int IndexVMerged = FacesMerged[i][j];
                    Geo::Point3<ScalarType> SmoothedPos = VertMerged[IndexVMerged];
                    VertPos[IndexVOrig] = SmoothedPos;
                }
        }
    };

}

#endif // LOCAL_PARAM_SMOOTH
