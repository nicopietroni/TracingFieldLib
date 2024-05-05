#ifndef UV_UTILS
#define UV_UTILS

#include "Space_Query/triangle_grid_2D.h"
#include <point2.h>
#include <cross_field.h>
#include <Space_Query/segment_grid_2D.h>

namespace Geo
{

    template <class ScalarType>
    class UVUtils
    {

        static void MapToAngle(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                               const std::vector<int> &VertSeq,
                               const ScalarType &Radius,
                               const ScalarType &Angle0,
                               const ScalarType &Angle1,
                               std::vector<Geo::Point2<ScalarType>> &vert_UV)
        {

            assert(Angle1 > Angle0);

            assert(Angle0 >= 0);
            assert(Angle0 <= 2 * M_PI);

            assert(Angle1 >= 0);
            assert(Angle1 <= 2 * M_PI);

            vert_UV.clear();
            ScalarType TotL = 0;
            bool IsLoop = ((Angle0 == 0) && (Angle1 == 2 * M_PI));
            int LenghtSeq = VertSeq.size();
            if (!IsLoop)
                LenghtSeq--;

            for (size_t i = 0; i < LenghtSeq; i++)
            {
                int IndexV0 = VertSeq[i];
                int IndexV1 = VertSeq[(i + 1) % VertSeq.size()];

                assert(IndexV0 < (int)VertPos.size());
                assert(IndexV1 < (int)VertPos.size());
                assert(IndexV0 >= 0);
                assert(IndexV1 >= 0);

                Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                ScalarType Dist = (P0 - P1).Norm();
                TotL += Dist;
            }

            ScalarType CurrL = 0;
            vert_UV.clear();
            ScalarType AngleInt = Angle1 - Angle0;
            for (size_t i = 0; i < VertSeq.size(); i++)
            {
                // ScalarType frac = 2 * M_PI - CurrL * (ScalarType)2 * M_PI / TotL;
                ScalarType frac = 2 * M_PI - Angle0 + CurrL * (ScalarType)AngleInt / TotL;
                vert_UV.push_back(Geo::Point2<ScalarType>(cos(frac), sin(frac)));

                int IndexV0 = VertSeq[i];
                int IndexV1 = VertSeq[(i + 1) % VertSeq.size()];
                Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                ScalarType Dist = (P0 - P1).Norm();

                CurrL += Dist;
            }

            for (size_t i = 0; i < vert_UV.size(); i++)
                vert_UV[i] = vert_UV[i] * Radius;
        }

        static int SortFromCorner(std::vector<int> &VertSeq,
                                  const int &CornerIdx)
        {
            // first get the boundery
            int firstCorner = -1;
            for (size_t i = 0; i < VertSeq.size(); i++)
                if (VertSeq[i] == CornerIdx)
                {
                    firstCorner = i;
                    break;
                }
            if (firstCorner == -1)
                return -1;
            std::rotate(VertSeq.begin(),
                        VertSeq.begin() + firstCorner,
                        VertSeq.end());
            return firstCorner;
        }

        // ScalarType LenghtSequence(const std::vector<Geo::Point3<ScalarType>> &VertPos,
        //                           const std::vector<int> &BorderSeq,
        //                           const int VIdx0, const int VIdx1)
        // {

        //     ScalarType CurrL = 0;
        //     int currentIndex = VIdx0;
        //     while (currentIndex != VIdx1)
        //     {
        //         int CurrV0 = currentIndex;
        //         int CurrV1 = (CurrV0 + 1) % BorderSeq.size();
        //         CurrV0 = BorderSeq[CurrV0];
        //         CurrV1 = BorderSeq[CurrV1];
        //         Geo::Point3<ScalarType> P0 = VertPos[CurrV0];
        //         Geo::Point3<ScalarType> P1 = VertPos[CurrV1];
        //         CurrL += (P0 - P1).Norm();
        //         currentIndex = (currentIndex + 1) % BorderSeq.size();
        //     }
        //     return CurrL;
        // }

        static ScalarType LenghtSequence(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                         const std::vector<int> &VertSeq, bool is_cycle)
        {
            ScalarType TotL = 0;
            int LastIdx = VertSeq.size();
            if (!is_cycle)
                LastIdx--;
            for (size_t i = 0; i < LastIdx; i++)
            {
                int IndexV0 = VertSeq[i];
                int IndexV1 = VertSeq[(i + 1) % VertSeq.size()];

                assert(IndexV0 < (int)VertPos.size());
                assert(IndexV1 < (int)VertPos.size());
                assert(IndexV0 >= 0);
                assert(IndexV1 >= 0);

                Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                ScalarType Dist = (P0 - P1).Norm();
                TotL += Dist;
            }
            return (TotL);
        }

        static bool FindSubSequences(const std::vector<int> &BorderSeq,
                                     const std::vector<int> &CornerIDX,
                                     std::vector<std::vector<int>> &VertSubSeq)
        {
            VertSubSeq.clear();
            assert(CornerIDX.size() > 0);
            assert(BorderSeq.size() > 1);
            if (BorderSeq[0] != CornerIDX[0])
            {
                assert(0);
                return false;
            }

            // then compute subsequences
            std::set<int> CornerIDXSet(CornerIDX.begin(), CornerIDX.end());
            VertSubSeq.resize(1);
            VertSubSeq.back().push_back(BorderSeq[0]);

            for (size_t i = 1; i < BorderSeq.size(); i++)
            {
                int CurrV = BorderSeq[i];
                VertSubSeq.back().push_back(CurrV);

                bool IsCorner = (CornerIDXSet.count(CurrV) > 0);

                if (IsCorner)
                {
                    VertSubSeq.resize(VertSubSeq.size() + 1);
                    VertSubSeq.back().push_back(CurrV);
                }
            }
            // add again the first one
            VertSubSeq.back().push_back(BorderSeq[0]);

            // final check
            for (size_t i = 0; i < VertSubSeq.size(); i++)
            {
                assert(CornerIDXSet.count(VertSubSeq[i][0]) > 0);
                assert(CornerIDXSet.count(VertSubSeq[i].back()) > 0);
                for (size_t j = 1; j < VertSubSeq[i].size() - 1; j++)
                    assert(CornerIDXSet.count(VertSubSeq[i][j]) == 0);
            }
            return true;
        }

        static Geo::Point2<ScalarType> AngleToUV(ScalarType alpha)
        {
            return (Geo::Point2<ScalarType>(cos(alpha), sin(alpha)));
        }

        static void InterpolateLinearUV(const Geo::Point2<ScalarType> &UV0,
                                        const Geo::Point2<ScalarType> &UV1,
                                        const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                        std::vector<int> &VertSeq,
                                        std::vector<Geo::Point2<ScalarType>> &UVSeq)
        {
            UVSeq.clear();
            ScalarType LSeq = LenghtSequence(VertPos, VertSeq, false);
            ScalarType CurrL = 0;
            for (size_t i = 0; i < VertSeq.size() - 1; i++)
            {
                ScalarType ratio = 1 - CurrL / LSeq; // 1-(ScalarType)i/(ScalarType)(VertSeq.size()-1);///
                assert(ratio >= 0);
                assert(ratio <= 1);
                Geo::Point2<ScalarType> CurrUV = UV0 * ratio + UV1 * (1 - ratio);
                UVSeq.push_back(CurrUV);

                int IndexV0 = VertSeq[i];
                int IndexV1 = VertSeq[(i + 1) % VertSeq.size()];

                assert(IndexV0 < (int)VertPos.size());
                assert(IndexV1 < (int)VertPos.size());
                assert(IndexV0 >= 0);
                assert(IndexV1 >= 0);

                Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
                Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
                ScalarType Dist = (P0 - P1).Norm();
                CurrL += Dist;
            }
            // make sure their precise at begin and add the one at the end
            UVSeq[0] = UV0;
            UVSeq.push_back(UV1);
        }

    public:
        enum UVType
        {
            UVArap,
            UVLSConformal,
            UVTutte,
            UVPoisson
        };

        static void MapToACircle(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                 const std::vector<int> &CircleSeq,
                                 const ScalarType &Radius,
                                 std::vector<Geo::Point2<ScalarType>> &vert_UV)
        {
            MapToAngle(VertPos, CircleSeq, Radius, 0, (ScalarType)2 * M_PI, vert_UV);
        }

        static bool MapToAPolygon(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                  const std::vector<int> &BorderSeq,
                                  const ScalarType &Radius,
                                  const std::vector<int> &CornerIDX,
                                  std::vector<Geo::Point2<ScalarType>> &vert_UV,
                                  bool uniform = false)
        {
            // for (size_t i=0;i<CornerIDX.size();i++)
            //     std::cout<<CornerIDX[i]<<std::endl;
            if (CornerIDX.size() < 3)
            {
                assert(0);
                return false;
            }
            std::vector<int> SortedSeq = BorderSeq;
            int IndexFirst = SortFromCorner(SortedSeq, CornerIDX[0]);
            if (IndexFirst == -1)
            {
                assert(0);
                return false;
            }

            std::vector<std::vector<int>> VertSubSeq;
            bool found_sub_seq = FindSubSequences(SortedSeq, CornerIDX, VertSubSeq);
            if (!found_sub_seq)
            {
                assert(0);
                return false;
            }
            assert(CornerIDX.size() == VertSubSeq.size());

            std::vector<Geo::Point2<ScalarType>> CornerUV;
            ScalarType TotL = LenghtSequence(VertPos, SortedSeq, true);

            // map corners
            if (uniform)
            {
                ScalarType curr_alpha = 0;
                ScalarType interval_alpha = ((ScalarType)(2 * M_PI)) /
                                            (ScalarType)CornerIDX.size();
                for (size_t i = 0; i < CornerIDX.size(); i++)
                {
                    Geo::Point2<ScalarType> CurrUV = AngleToUV(curr_alpha) * Radius;
                    CornerUV.push_back(CurrUV);
                    curr_alpha += interval_alpha;
                }
            }
            else
            {
                ScalarType curr_alpha = 0;
                for (size_t i = 0; i < CornerIDX.size(); i++)
                {
                    ScalarType SubSeqL = LenghtSequence(VertPos, VertSubSeq[i], false);
                    ScalarType Ratio = SubSeqL / TotL; // 1/(ScalarType)CornerIDX.size();
                    Geo::Point2<ScalarType> CurrUV = AngleToUV(curr_alpha) * Radius;
                    CornerUV.push_back(CurrUV);
                    curr_alpha += Ratio * (ScalarType)(2 * M_PI);
                }
            }

            assert(CornerUV.size() == CornerIDX.size());

            // then interpolate the ones in between
            for (size_t i = 0; i < CornerIDX.size(); i++)
            {
                Geo::Point2<ScalarType> CurrUV0 = CornerUV[i];
                Geo::Point2<ScalarType> CurrUV1 = CornerUV[(i + 1) % CornerUV.size()];
                std::vector<Geo::Point2<ScalarType>> UVSeq;
                InterpolateLinearUV(CurrUV0, CurrUV1, VertPos, VertSubSeq[i], UVSeq);
                assert(UVSeq.size() == VertSubSeq[i].size());
                // Insert elements from the beginning but the last because it is repeated in the next
                vert_UV.insert(vert_UV.end(), UVSeq.begin(), UVSeq.end() - 1);
            }
            assert(vert_UV.size() == BorderSeq.size());
            // FINAL CHECK
            // std::cout<<"First UV 0 :"<<vert_UV[0].X<<","<<vert_UV[0].Y<<std::endl;
            std::rotate(vert_UV.rbegin(),
                        vert_UV.rbegin() + IndexFirst,
                        vert_UV.rend());
            // int Corner0=CornerIDX[0];
            // for (size_t i=0;i<BorderSeq.size();i++)
            // {
            //     if (BorderSeq[i]!=Corner0)continue;
            //     std::cout<<"First UV 1 : "<<vert_UV[i].X<<","<<vert_UV[i].Y<<std::endl;
            //     break;
            // }
            // std::cout<<"First UV 1 : "<<vert_UV[CornerIDX[0]].X<<","<<vert_UV[CornerIDX[0]].Y<<std::endl;
            return true;
        }

        static bool MapBetween(const std::vector<Geo::Point2<ScalarType>> &vert0,
                               const std::vector<Geo::Point2<ScalarType>> &vert1,
                               const std::vector<std::vector<int>> &face1,
                               std::vector<int> &faceIDx, std::vector<Geo::Point3<ScalarType>> &bary)
        {
            Geo::TriangleSpatialIndex2D<ScalarType> Grid_2D;
            Grid_2D.Init(vert1, face1);
            faceIDx.resize(vert0.size(), -1);
            bary.resize(vert0.size(), Geo::Point3<ScalarType>(0, 0, 0));
            for (size_t i = 0; i < vert0.size(); i++)
            {
                bool found = Grid_2D.FirstTrianglesOverPoint(vert1, face1, vert0[i], faceIDx[i], bary[i]);
                if (!found)
                    return false;
            }
            return true;
        }

        // static bool Interpolate3DPosFromUV(const std::vector<Geo::Point2<ScalarType> > &vert0_2D,
        //                                    const std::vector<Geo::Point2<ScalarType> > &vert1_2D,
        //                                    const std::vector<Geo::Point3<ScalarType> > &vert1_3D,
        //                                    const std::vector<std::vector<int> > &face1,
        //                                    std::vector<Geo::Point3<ScalarType> > &interp_vert0_3D)
        // {
        //     std::vector<std::vector<int> > ConnectivityTris;
        //     if (!IsAllTris(face1))
        //     {
        //         assert(IsAllQuad(face1));
        //         SplitQuadIntoTris(vert1_3D,face1,ConnectivityTris);
        //     }
        //     else
        //         ConnectivityTris=face1;

        //     std::vector<int> faceIDx;
        //     std::vector<Geo::Point3<ScalarType>> bary;
        //     bool has_mapped=MapBetween(vert0_2D,vert1_2D,ConnectivityTris,faceIDx,bary);
        //     if (!has_mapped)return false;

        //     interp_vert0_3D.resize(vert0_2D.size());
        //     for (size_t i=0;i<faceIDx.size();i++)
        //     {
        //         int FaceInterp=faceIDx[i];
        //         assert(FaceInterp>=0);
        //         assert(FaceInterp<(int)ConnectivityTris.size());

        //         int VIdx0=ConnectivityTris[FaceInterp][0];
        //         int VIdx1=ConnectivityTris[FaceInterp][1];
        //         int VIdx2=ConnectivityTris[FaceInterp][2];

        //         Geo::Point3<ScalarType> P0=vert1_3D[VIdx0];
        //         Geo::Point3<ScalarType> P1=vert1_3D[VIdx1];
        //         Geo::Point3<ScalarType> P2=vert1_3D[VIdx2];

        //         Geo::Point3<ScalarType> PInterp=P0*bary[i].X+
        //                 P1*bary[i].Y+
        //                 P2*bary[i].Z;
        //         interp_vert0_3D[i]=PInterp;
        //     }
        //     return true;
        // }

        // static void ScaleToMatchArea(const std::vector<Geo::Point3<ScalarType> > &vert_3D,
        //                              const std::vector<std::vector<int> > &Faces,
        //                              std::vector<Geo::Point2<ScalarType> > &vert_UV)
        // {
        //     ScalarType Area3D=Area(vert_3D,Faces);
        //     ScalarType Area2D=Area(vert_UV,Faces);
        //     ScalarType ScaleFact=Area3D/Area2D;
        //     for (size_t i=0;i<vert_UV.size();i++)
        //         vert_UV[i]=vert_UV[i]*ScaleFact;
        // }

        static void ScaleToMatchArea(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                     const std::vector<std::vector<int>> &Faces,
                                     std::vector<Geo::Point2<ScalarType>> &VertUV)
        {
            // Rescale the parametrization to match the 3D area
            ScalarType meshArea2D = Area(VertUV, Faces);
            ScalarType meshArea3D = Area(VertPos, Faces);

            // std::cout<<"Area 2D:"<<meshArea2D<<std::endl;
            // std::cout<<"Area 3D:"<<meshArea3D<<std::endl;
            ScalarType scaleFact = std::sqrt(meshArea3D / meshArea2D);

            for (size_t i = 0; i < VertUV.size(); i++)
                VertUV[i] = VertUV[i] * scaleFact;
        }

        static void Mirror(std::vector<Geo::Point2<ScalarType>> &vert_UV, int Dir)
        {
            assert((Dir >= 0) && (Dir < 2));
            Geo::Box2<ScalarType> BB;
            BB.Init(vert_UV);
            ScalarType CenterVal = 0;
            if (Dir == 0)
                CenterVal = BB.Center().X;
            else
                CenterVal = BB.Center().Y;

            for (size_t i = 0; i < vert_UV.size(); i++)
            {
                ScalarType Diff = 0;
                if (Dir == 0)
                    Diff = vert_UV[i].X - CenterVal;
                else
                    Diff = vert_UV[i].Y - CenterVal;

                if (Dir == 0)
                    vert_UV[i].X = CenterVal - Diff;
                else
                    vert_UV[i].Y = CenterVal - Diff;
            }
        }

        static bool Bijective(const std::vector<Geo::Point2<ScalarType>> &VertUV,
                              const std::vector<std::vector<int>> &Faces)
        {
            const ScalarType eps = 0.000001;
            // first test, simply check if a coordinate appear twice
            std::set<Geo::Point2<ScalarType>> TestVertUV(VertUV.begin(), VertUV.end());
            if (TestVertUV.size() != VertUV.size())
                return false;

            // check if there are folded tris
            int PosA = 0;
            int NegA = 0;

            for (int i = 0; i < Faces.size(); i++)
            {
                Geo::Point2<ScalarType> P0 = VertUV[Faces[i][0]];
                Geo::Point2<ScalarType> P1 = VertUV[Faces[i][1]];
                Geo::Point2<ScalarType> P2 = VertUV[Faces[i][2]];
                ScalarType FaceArea = (Geo::Cross((P1 - P0), (P2 - P0)) / 2);
                if (FaceArea == 0)
                    return false;
                if (FaceArea < 0)
                    NegA++;
                else
                    PosA++;
            }
            // if not all positive or all negative has fold
            if ((NegA > 0) && (PosA > 0))
                return false;

            // check if the border intersect itself
            std::vector<Geo::Segment2<ScalarType>> BorderSeg;
            Geo::SegmentSpatialIndex2D<ScalarType> SegIndex2D;
            std::vector<std::pair<int, int>> BorderE;
            Geo::getBorderEdges(Faces, BorderE);

            for (size_t i = 0; i < BorderE.size(); i++)
            {
                int IndexV0 = BorderE[i].first;
                int IndexV1 = BorderE[i].second;
                assert(IndexV0 >= 0);
                assert(IndexV1 >= 0);
                assert(IndexV0 < VertUV.size());
                assert(IndexV1 < VertUV.size());
                Geo::Point2<ScalarType> P0 = VertUV[IndexV0];
                Geo::Point2<ScalarType> P1 = VertUV[IndexV1];
                Geo::Segment2<ScalarType> s(P0, P1);

                // scale a bit to not intersect the adjacent
                Scale(s, 1 - eps);
                BorderSeg.push_back(Geo::Segment2<ScalarType>(P0, P1));
            }
            SegIndex2D.Init(BorderSeg);

            for (size_t i = 0; i < BorderSeg.size(); i++)
            {
                std::vector<int> IntSegIdx;
                std::vector<Geo::Point2<ScalarType>> IntSegPos;
                SegIndex2D.InteresectSeg(BorderSeg, BorderSeg[i], IntSegIdx, IntSegPos);
                if (IntSegIdx.size() == 0)
                    continue;
                for (size_t j = 0; j < IntSegIdx.size(); j++)
                {
                    int SegIdx = IntSegIdx[j];
                    assert(SegIdx >= 0);
                    assert(SegIdx < BorderSeg.size());

                    // check same segment
                    if (IntSegIdx[j] == (int)i)
                        continue;

                    // check if adjacent one
                    if (BorderSeg[SegIdx].P(0) == BorderSeg[i].P(0))
                        continue;
                    if (BorderSeg[SegIdx].P(0) == BorderSeg[i].P(1))
                        continue;
                    if (BorderSeg[SegIdx].P(1) == BorderSeg[i].P(0))
                        continue;
                    if (BorderSeg[SegIdx].P(1) == BorderSeg[i].P(1))
                        continue;
                    return false;
                }
            }
            return true;
        }

        static int ClosestTo(const std::vector<Geo::Point2<ScalarType>> &VertUV,
                             const Geo::Point2<ScalarType> &TestPos)
        {
            ScalarType minD = std::numeric_limits<ScalarType>::max();
            int ClosestIDX = -1;
            for (size_t i = 0; i < VertUV.size(); i++)
            {
                ScalarType TestD = (VertUV[i] - TestPos).Norm();
                if (TestD > minD)
                    continue;
                minD = TestD;
                ClosestIDX = i;
            }
            assert(ClosestIDX != -1);
            return ClosestIDX;
        }

        /// return the number of folded faces
        static bool IsFolded(const std::vector<Geo::Point2<ScalarType>> &VertUV,
                             const std::vector<std::vector<int>> &Faces,
                             const int IndexF)
        {
            ScalarType areaUV = Area(VertUV, Faces, IndexF, false);
            return (areaUV < 0);
        }

        static int FoldedNum(const std::vector<Geo::Point2<ScalarType>> &VertUV,
                             const std::vector<std::vector<int>> &Faces)
        {
            int folded = 0;

            for (size_t i = 0; i < Faces.size(); i++)
                if (IsFolded(VertUV, Faces, i))
                    folded++;

            return folded;
        }

        static bool GloballyUnFolded(const std::vector<Geo::Point2<ScalarType>> &VertUV,
                                     const std::vector<std::vector<int>> &Faces)
        {
            int num = FoldedNum(VertUV,Faces);
            return (num > (Faces.size()) / 2);
        }
    };

}

#endif
