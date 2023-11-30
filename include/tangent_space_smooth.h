#ifndef TANGENT_SPACE_SMOOTH
#define TANGENT_SPACE_SMOOTH

#include "triangle2.h"
#include "connectivity.h"
#include "point3.h"
#include "pos3.h"
#include "cotangent.h"
#include "triangle3.h"

namespace Geo{

template <class ScalarType>
class Local_Param_Smooth
{

public:

    struct UVSmoothParam
    {
        std::vector<int> Fixed;
        std::vector<std::pair<int,int> > FeatureLines;
        ScalarType Damp;
        size_t Steps;

        UVSmoothParam()
        {
            Damp=0.2;
            Steps=5;
        }
    };

private:


    static Geo::Point3<ScalarType > Interpolate3DPos(const std::vector<Geo::Traingle2<ScalarType> > &TrisUV,
                                                     const std::vector<Geo::Traingle3<ScalarType> > &Tris3D,
                                                     const Geo::Point2<ScalarType> &To_Sample)
    {
        assert(TrisUV.size()>=2);
        for (size_t i=0;i<TrisUV.size();i++)
        {
            if (!Geo::IsInside(TrisUV[i],To_Sample))continue;
            ScalarType alpha,beta,gamma;
            Geo::BarycentricCoords(TrisUV[i],To_Sample,alpha,beta,gamma);

            Geo::Point3<ScalarType > Ret=Tris3D[i].P(0)*alpha+
                    Tris3D[i].P(1)*beta+
                    Tris3D[i].P(2)*gamma;
            return Ret;
        }
        assert(0);
    }

    static Geo::Point3<ScalarType > ReprojectOnPoliline(const std::vector<Geo::Point3<ScalarType> > &PolyL,
                                                        const Geo::Point3<ScalarType> &TestP)
    {
        ScalarType MinDist=std::numeric_limits<ScalarType>::max();
        Geo::Point3<ScalarType > minP(0,0,0);
        for (size_t i=0;i<PolyL.size();i++)
        {
            Geo::Point3<ScalarType> P0=PolyL[i];
            Geo::Point3<ScalarType> P1=PolyL[(i+1)%PolyL.size()];
            Geo::Segment3<ScalarType> Seg(P0,P1);
            Geo::Point3<ScalarType> ClosP=Geo::ClosestPoint(Seg,TestP);
            ScalarType Dist=(ClosP-TestP).Norm();
            if (Dist>MinDist)continue;
            Dist=MinDist;
            minP=ClosP;
        }
        return minP;
    }


    static void ParametrizeStar(std::vector<Geo::Point3<ScalarType > > &VertPos,
                                const std::vector<std::vector<int> > &Faces,
                                const std::vector<std::vector<int> > &NextF,
                                const std::vector<std::vector<int> > &NextE,
                                const std::vector<Pos3> &StarMesh,
                                std::vector<Geo::Traingle2<ScalarType> > &TrisUV,
                                std::vector<Geo::Traingle3<ScalarType> > &Tris3D)
    {
        std::map<int,Geo::Point2<ScalarType> > VertUVsMap;

        TrisUV.clear();
        Tris3D.clear();

        int CentralV=StarMesh[0].V;

        ScalarType AngleInterval=2.0*M_PI/(ScalarType)StarMesh.size();
        ScalarType CurrAngle=0;

        Geo::Point2<ScalarType> CentralUV(0,0);
        ScalarType SumW=0;
        for (size_t i=0;i<StarMesh.size();i++)
        {
            int VIdx0=StarMesh[i].V;
            assert(VIdx0==CentralV);

            //then get the one outside on the star
            Pos3 CurrP=StarMesh[i];
            CurrP.FlipV(Faces);

            int VIdx1=CurrP.V;
            assert(VIdx0!=VIdx1);

            Geo::Point2<ScalarType> CurrUV(cos(CurrAngle),sin(CurrAngle));
            CurrAngle+=AngleInterval;

            VertUVsMap[VIdx1]=CurrUV;

            ScalarType CurrW=CotangentWeightTriFace(VertPos,Faces,NextF,NextE,StarMesh[i].F,StarMesh[i].E);
            CentralUV+=(CurrUV*CurrW);
            SumW+=CurrW;
        }
        CentralUV/=SumW;
        assert(VertUVsMap.count(CentralV)==0);
        VertUVsMap[CentralV]=CentralUV;

        //then return the parametrization
        for (size_t i=0;i<StarMesh.size();i++)
        {
            int IndexF=StarMesh[i].F;
            assert(IndexF>=0);
            assert(IndexF<Faces.size());
            assert(Faces[IndexF].size()==3);

            Geo::Point2<ScalarType> P_UV[3];
            Geo::Point3<ScalarType> P_3D[3];

            for (size_t j=0;j<Faces[IndexF].size();j++)
            {
                int IndexV=Faces[IndexF][j];
                assert(VertUVsMap.count(IndexV)>0);
                P_UV[j]=VertUVsMap[IndexV];
                P_3D[j]=VertPos[IndexV];
            }
            Geo::Traingle2<ScalarType> T_UV(P_UV[0],P_UV[1],P_UV[2]);
            Geo::Traingle3<ScalarType> T_3D(P_3D[0],P_3D[1],P_3D[2]);

            TrisUV.push_back(T_UV);
            Tris3D.push_back(T_3D);
        }
    }

    static void SmoothStep(std::vector<Geo::Point3<ScalarType > > &VertPos,
                           const std::vector<std::vector<int> > &Faces,
                           const std::vector<std::vector<int> > &NextF,
                           const std::vector<std::vector<int> > &NextE,
                           const std::vector<std::vector<Pos3> > &StarMeshes,
                           const std::vector<bool> &IsB,
                           std::vector<std::vector<int > > BorderVertSeq,
                           UVSmoothParam UVP=UVSmoothParam())
    {
        std::vector<Geo::Point3<ScalarType > > TargetP(VertPos.size(),Geo::Point3<ScalarType >(0,0,0));
        std::vector<size_t> TargetN(VertPos.size(),0);

        std::vector<Geo::Point3<ScalarType> > ReprojectPolys;
        ReprojectPolys.resize(BorderVertSeq.size());
        for (size_t i=0;i<BorderVertSeq.size();i++)
        {
            for (size_t j=0;j<BorderVertSeq[i].size();j++)
            {
                Geo::Point3<ScalarType > P=VertPos[BorderVertSeq[i][j]];
                ReprojectPolys[i].push_back(P);
            }
        }

        //calculate target positions for border
        for (size_t i=0;i<BorderVertSeq.size();i++)
        {
            int SizeB=BorderVertSeq[i].size();
            assert(SizeB>2);

            if (SizeB<=3)continue;//too small border sequence to be smoothed

            for (size_t j=0;j<BorderVertSeq[i].size();j++)
            {
                int LocIndex=j;
                int LocIndex0=(j+(SizeB-1))%SizeB;
                int LocIndex1=(j+1)%SizeB;

                int IndexV=BorderVertSeq[i][LocIndex];
                int IndexV0=BorderVertSeq[i][LocIndex0];
                int IndexV1=BorderVertSeq[i][LocIndex1];

                TargetP[IndexV]+=VertPos[IndexV0];
                TargetP[IndexV]+=VertPos[IndexV1];
                TargetN[IndexV]++;
                TargetN[IndexV]++;
            }
        }

        for (size_t i=0;i<StarMeshes.size();i++)
        {
            if (IsB[i])continue;

            assert(StarMeshes[i].size()>=2);

            std::vector<Geo::Traingle2<ScalarType> > TrisUV;
            std::vector<Geo::Traingle3<ScalarType> > Tris3D;
            ParametrizeStar(VertPos,Faces,NextF,NextE,StarMeshes[i],TrisUV,Tris3D);


            TargetP[i]=Interpolate3DPos(TrisUV,Tris3D,Geo::Point2<ScalarType>(0,0));
            TargetN[i]++;
        }

        for (size_t i=0;i<TargetP.size();i++)
        {
            assert(TargetN[i]>0);
            TargetP[i]/=(ScalarType)TargetN[i];
        }

        //reproject borders
        for (size_t i=0;i<BorderVertSeq.size();i++)
        {
            for (size_t j=0;j<BorderVertSeq[i].size();j++)
            {
                int indexV=VertPos[BorderVertSeq[i][j]];
                Geo::Point3<ScalarType > P=VertPos[indexV];
                TargetP[i]=ReprojectOnPoliline(ReprojectPolys[i],P);
            }
        }

        //update mesh position
        for (size_t i=0;i<VertPos.size();i++)
        {
            VertPos[i]=VertPos[i]*UVP.Damp+TargetP[i]*(1-UVP.Damp);
        }

    }

    static void GetStarMeshes(const std::vector<Geo::Point3<ScalarType > > &VertPos,
                              const std::vector<std::vector<int> > &Faces,
                              const std::vector<std::vector<int> > &NextF,
                              const std::vector<std::vector<int> > &NextE,
                              std::vector<bool> IsBorder,
                              std::vector<std::vector<Pos3> > &StarMeshes)
    {

        StarMeshes.clear();
        StarMeshes.resize(VertPos.size());
        for (size_t i=0;i<Faces.size();i++)
            for (size_t j=0;j<3;j++)
            {

                //then get the star of vertices
                Pos3 TestPos(i,j,Faces[i][j]);
                int CenterV=TestPos.V;

                //check if already initialized or not
                if (StarMeshes[CenterV].size()>0)continue;
                //border do nothing
                if (IsBorder[CenterV])continue;

                std::vector<Pos3> VertFan;
                bool HasBorder;
                GetSortedFacesFaceFun(TestPos,Faces,NextF,NextE,VertFan,HasBorder);
                assert(!HasBorder);

                if (VertFan.size()==0)continue;
                StarMeshes[CenterV]=VertFan;
            }
    }

public:

    static void Smooth(std::vector<Geo::Point3<ScalarType > > &VertPos,
                       std::vector<std::vector<int> > &Faces,
                       UVSmoothParam UVP=UVSmoothParam())
    {

        //initialize constraints
        std::set<std::pair<int,int> > AllFeatures;
        for (size_t i=0;i<UVP.FeatureLines.size();i++)
        {
            int IndexV0=UVP.FeatureLines[i].first;
            int IndexV1=UVP.FeatureLines[i].second;
            std::pair<int,int> entry(std::min(IndexV0,IndexV1),
                                     std::max(IndexV0,IndexV1));
            AllFeatures.insert(entry);
        }

        //then add borders
        std::vector<std::vector<int> > NextF,NextE;
        Geo::getFFAdjacency(Faces,NextF,NextE);


        for (size_t i=0;i<NextF.size();i++)
        {
            assert(NextF[i].size()==3);//only triangles
            for (size_t j=0;j<NextF[i].size();j++)
            {
                if (NextF[i][j]!=-1)continue;
                int IndexV0=Faces[i][j];
                int IndexV1=Faces[i][(j+1)%NextF[i].size()];
                std::pair<int,int> entry(std::min(IndexV0,IndexV1),
                                         std::max(IndexV0,IndexV1));
                AllFeatures.insert(entry);
            }
        }


        //        std::vector<std::vector<PosType> > ProjBasis;
        //        GetStarMeshes(mesh,StarMeshes,ProjBasis,OrigIndex);
        std::vector<bool> IsBorder;
        Geo::getBorderVert(Faces,NextF,VertPos.size(),IsBorder);
        std::vector<std::vector<int > > BorderVertSeq;
        getBorderSequncesVert(Faces,NextF,NextE,BorderVertSeq);

        //create star meshes
        std::vector<std::vector<Pos3> > StarMeshes;
        GetStarMeshes(VertPos,Faces,NextF,NextE,IsBorder,StarMeshes);

        //and finally smooth
        for (size_t i=0;i<UVP.Steps;i++)
            SmoothStep(VertPos,Faces,NextF,NextE,StarMeshes,IsBorder,BorderVertSeq,UVP);

    }
};

}

#endif // LOCAL_PARAM_SMOOTH
