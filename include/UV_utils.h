#ifndef UV_UTILS
#define UV_UTILS

#include "Space_Query/triangle_grid_2D.h"
#include <point2.h>
#include <cross_field.h>
#include <Space_Query/segment_grid_2D.h>

namespace Geo{

template <class ScalarType >
class UVUtils{

public:

    static void MapToACircle(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                             const std::vector<int> &CircleSeq,
                             const ScalarType &Radius,
                             std::vector<Geo::Point2<ScalarType> > &vert_UV)
    {
        vert_UV.clear();
        //compute total lenght
        ScalarType TotL=0;

        for (size_t i=0;i<CircleSeq.size();i++)
        {
            int IndexV0=CircleSeq[i];
            int IndexV1=CircleSeq[(i+1)%CircleSeq.size()];

            assert(IndexV0<(int)VertPos.size());
            assert(IndexV1<(int)VertPos.size());
            assert(IndexV0>=0);
            assert(IndexV1>=0);

            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
            ScalarType Dist=(P0-P1).Norm();
            TotL+=Dist;
        }

        ScalarType CurrL=0;
        vert_UV.clear();
        for (size_t i=0;i<CircleSeq.size();i++)
        {
            ScalarType frac = 2*M_PI-CurrL * (ScalarType) 2 * M_PI / TotL;
            vert_UV.push_back(Geo::Point2<ScalarType>(cos(frac), sin(frac)));

            int IndexV0=CircleSeq[i];
            int IndexV1=CircleSeq[(i+1)%CircleSeq.size()];
            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
            ScalarType Dist=(P0-P1).Norm();

            CurrL+=Dist;
        }

        for (size_t i=0;i<vert_UV.size();i++)
            vert_UV[i]=vert_UV[i]*Radius;
    }

    static void MapToAPolygon(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                              const std::vector<int> &CircleSeq,
                              const ScalarType &Radius,
                              const std::vector<int> &CornerIDX,
                              std::vector<Geo::Point2<ScalarType> > &vert_UV)
    {
        MapToACircle(VertPos,CircleSeq,Radius,vert_UV);
        if (CornerIDX.size()<=0)return;

        //        std::set<int> CornerIDXSet(CornerIDX.begin(),CornerIDX.end());
        //        std::vector<std::vector<int> > BorderSeq;
        //        int startingIdx=-1;
        //        for (size_t i=0;i<CircleSeq.size();i++)
        //        {
        //            int IndexV=CircleSeq[i];
        //            if (CornerIDXSet.count(IndexV)==0)continue;
        //            startingIdx=i;
        //            break;
        //        }
        //        assert(startingIdx!=-1);

        //        for (size_t i=0;i<CircleSeq.size();i++)
        //        {
        //            int CurrIdx=(i+startingIdx)%CircleSeq.size();
        //            int IndexV=CircleSeq[CurrIdx];

        //            if ((CornerIDXSet.count(IndexV)>0)&&(i!=(CircleSeq.size()-1)))
        //            {
        //                if (BorderSeq.size()>0)
        //                    BorderSeq.back().push_back(IndexV);

        //                BorderSeq.resize(BorderSeq.size()+1);
        //            }

        //            BorderSeq.back().push_back(IndexV);
        //        }
        //        BorderSeq.back().push_back(CircleSeq[startingIdx]);

        //        std::cout<<"corner size:"<<CornerIDX.size()<<std::endl;
        //        std::cout<<"boudanry size:"<<BorderSeq.size()<<std::endl;
        //        //check
        //        for (size_t i=0;i<BorderSeq.size();i++)
        //        {
        //            int next_i=(i+1)%BorderSeq.size();
        //            int IndexVEnd=BorderSeq[i].back();
        //            int IndexVNext=BorderSeq[next_i][0];
        //            std::cout<<"Side:"<<i<<std::endl;
        //            std::cout<<IndexVEnd<<std::endl;
        //            std::cout<<IndexVNext<<std::endl;
        //            std::cout<<std::endl<<std::endl;
        //            assert(IndexVEnd==IndexVNext);
        //        }

        //        vert_UV.clear();
        //        //compute total lenght
        //        ScalarType TotL=0;

        //        for (size_t i=0;i<CircleSeq.size();i++)
        //        {
        //            int IndexV0=CircleSeq[i];
        //            int IndexV1=CircleSeq[(i+1)%CircleSeq.size()];

        //            assert(IndexV0<(int)VertPos.size());
        //            assert(IndexV1<(int)VertPos.size());
        //            assert(IndexV0>=0);
        //            assert(IndexV1>=0);

        //            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
        //            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
        //            ScalarType Dist=(P0-P1).Norm();
        //            TotL+=Dist;
        //        }

        //        ScalarType CurrL=0;
        //        vert_UV.clear();
        //        for (size_t i=0;i<CircleSeq.size();i++)
        //        {
        //            ScalarType frac = CurrL * (ScalarType) 2 * M_PI / TotL;
        //            vert_UV.push_back(Geo::Point2<ScalarType>(cos(frac), sin(frac)));

        //            int IndexV0=CircleSeq[i];
        //            int IndexV1=CircleSeq[(i+1)%CircleSeq.size()];
        //            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
        //            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
        //            ScalarType Dist=(P0-P1).Norm();

        //            CurrL+=Dist;
        //        }

        //        for (size_t i=0;i<vert_UV.size();i++)
        //            vert_UV[i]=vert_UV[i]*Radius;
    }

    static bool MapBetween(const std::vector<Geo::Point2<ScalarType> > &vert0,
                           const std::vector<Geo::Point2<ScalarType> > &vert1,
                           const std::vector<std::vector<int> > &face1,
                           std::vector<int> &faceIDx,std::vector<Geo::Point3<ScalarType>> &bary)
    {
        Geo::TriangleSpatialIndex2D<ScalarType> Grid_2D;
        Grid_2D.Init(vert1,face1);
        faceIDx.resize(vert0.size(),-1);
        bary.resize(vert0.size(),Geo::Point3<ScalarType>(0,0,0));
        for (size_t i=0;i<vert0.size();i++)
        {
            bool found=Grid_2D.FirstTrianglesOverPoint(vert1,face1,vert0[i],faceIDx[i],bary[i]);
            if (!found)return false;
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

    static void ScaleToMatchArea(const std::vector<Geo::Point3<ScalarType> > &vert_3D,
                                 const std::vector<std::vector<int> > &Faces,
                                 std::vector<Geo::Point2<ScalarType> > &vert_UV)
    {
        ScalarType Area3D=Area(vert_3D,Faces);
        ScalarType Area2D=Area(vert_UV,Faces);
        ScalarType ScaleFact=Area3D/Area2D;
        for (size_t i=0;i<vert_UV.size();i++)
            vert_UV[i]=vert_UV[i]*ScaleFact;
    }


    static void Mirror(std::vector<Geo::Point2<ScalarType> > &vert_UV,int Dir)
    {
        assert((Dir>=0)&&(Dir<2));
        Geo::Box2<ScalarType> BB;
        BB.Init(vert_UV);
        ScalarType CenterVal=0;
        if (Dir==0)
           CenterVal=BB.Center().X;
        else
           CenterVal=BB.Center().Y;

        for (size_t i=0;i<vert_UV.size();i++)
        {
            ScalarType Diff=0;
            if (Dir==0)
                Diff=vert_UV[i].X-CenterVal;
            else
                Diff=vert_UV[i].Y-CenterVal;
            
            if (Dir==0)
                vert_UV[i].X=CenterVal-Diff;
            else
                vert_UV[i].Y=CenterVal-Diff;
        }
    }

    static bool Bijective(const std::vector<Geo::Point2<ScalarType> > &VertUV,
                          const std::vector<std::vector<int> > &Faces)
    {
        const ScalarType eps=0.000001;
        //first test, simply check if a coordinate appear twice
        std::set<Geo::Point2<ScalarType> > TestVertUV(VertUV.begin(),VertUV.end());
        if (TestVertUV.size()!=VertUV.size())return false;

        //check if there are folded tris
        int PosA=0;
        int NegA=0;

        for (int i=0;i<Faces.size();i++)
        {
            Geo::Point2<ScalarType> P0=VertUV[Faces[i][0]];
            Geo::Point2<ScalarType> P1=VertUV[Faces[i][1]];
            Geo::Point2<ScalarType> P2=VertUV[Faces[i][2]];
            ScalarType FaceArea=(Geo::Cross((P1-P0),(P2-P0))/2);
            if (FaceArea==0)return false;
            if (FaceArea<0)
                NegA++;
            else
                PosA++;
        }
        //if not all positive or all negative has fold
        if ((NegA>0)&&(PosA>0))return false;

        //check if the border intersect itself
        std::vector<Geo::Segment2<ScalarType> > BorderSeg;
        Geo::SegmentSpatialIndex2D<ScalarType> SegIndex2D;
        std::vector<std::pair<int,int> > BorderE;
        Geo::getBorderEdges(Faces,BorderE);

        for (size_t i=0;i<BorderE.size();i++)
        {
            int IndexV0=BorderE[i].first;
            int IndexV1=BorderE[i].second;
            assert(IndexV0>=0);
            assert(IndexV1>=0);
            assert(IndexV0<VertUV.size());
            assert(IndexV1<VertUV.size());
            Geo::Point2<ScalarType> P0=VertUV[IndexV0];
            Geo::Point2<ScalarType> P1=VertUV[IndexV1];
            Geo::Segment2<ScalarType> s(P0,P1);

            //scale a bit to not intersect the adjacent
            Scale(s,1-eps);
            BorderSeg.push_back(Geo::Segment2<ScalarType>(P0,P1));
        }
        SegIndex2D.Init(BorderSeg);

        for (size_t i=0;i<BorderSeg.size();i++)
        {
            std::vector<int> IntSegIdx;
            std::vector<Geo::Point2<ScalarType> > IntSegPos;
            SegIndex2D.InteresectSeg(BorderSeg,BorderSeg[i],IntSegIdx,IntSegPos);
            if (IntSegIdx.size()==0)continue;
            for (size_t j=0;j<IntSegIdx.size();j++)
            {
                int SegIdx=IntSegIdx[j];
                assert(SegIdx>=0);
                assert(SegIdx<BorderSeg.size());

                //check same segment
                if (IntSegIdx[j]==(int)i)continue;

                //check if adjacent one
                if (BorderSeg[SegIdx].P(0)==BorderSeg[i].P(0))continue;
                if (BorderSeg[SegIdx].P(0)==BorderSeg[i].P(1))continue;
                if (BorderSeg[SegIdx].P(1)==BorderSeg[i].P(0))continue;
                if (BorderSeg[SegIdx].P(1)==BorderSeg[i].P(1))continue;
                return false;
            }
        }
        return true;
    }
};

}

#endif
