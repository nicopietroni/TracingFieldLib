#ifndef CROSS_FIELD_FUNCTIONS
#define CROSS_FIELD_FUNCTIONS

#include <cmath>
#include <vector>
#include <stdio.h>
#include "point3.h"
#include "eigen_interface.h"
#include "cross_field.h"
#include "pos3.h"
#include <set>
#include <connectivity.h>
#include "best_rotation.h"
#include "mesh_functions.h"

template < typename ScalarType >
Field::CrossF<ScalarType> InterpolateTangentCrosses(const std::vector<Field::CrossF<ScalarType> > &Crosses,
                                                    const std::vector<ScalarType> &W)
{
    assert(Crosses.size()>0);
    assert(W.size()==Crosses.size());

    Field::CrossF<ScalarType> Cret;

    //get reference direction with highest quality
    Geo::Point3<ScalarType> RefDir=Crosses[0].D[0];
    ScalarType MaxQ=Crosses[0].Q;
    for (size_t i=1;i<Crosses.size();i++)
    {
        if (Crosses[i].Q>MaxQ)
        {
            MaxQ=Crosses[i].Q;
            RefDir=Crosses[i].D[0];
        }
        assert(Crosses[0].IsTangent(Crosses[i]));
    }
    for (size_t i=1;i<Crosses.size();i++)
    {
        assert(Crosses[0].IsTangent(Crosses[i]));
    }

    //get the first one as reference direction
    Geo::Point3<ScalarType> AvgDir(0,0,0);
    Geo::Point3<ScalarType> TargetN=Crosses[0].OrthoDir();
    ScalarType SumW=0;
    for (size_t i=0;i<Crosses.size();i++)
    {
        AvgDir+=Crosses[i].BestDirection(RefDir)*W[i];
        Cret.Pos+=Crosses[i].Pos;
        SumW+=W[i];
    }
    assert(SumW>0);

    AvgDir.Normalize();
    Cret.D[0]=AvgDir;
    Cret.D[1]=Cross(Cret.D[0],TargetN);
    Cret.D[1].Normalize();
    Cret.Pos/=(ScalarType)Crosses.size();

    Cret.K[0]=0;
    Cret.K[1]=0;
    for (size_t i=0;i<Crosses.size();i++)
    {
        int IndexD0=(Cret.BestDirectionI(Crosses[i].D[0])%2);
        int IndexD1=(IndexD0+1)%2;
        Cret.K[IndexD0]+=Crosses[i].K[0];
        Cret.K[IndexD1]+=Crosses[i].K[1];
    }

    Cret.K[0]/=(ScalarType)Crosses.size();
    Cret.K[1]/=(ScalarType)Crosses.size();

    return Cret;
}

template <class ScalarType>
void MakeDirectionNormalCoherent(const std::vector<Geo::Point3<ScalarType> > &TargetNormals,
                                 std::vector<Field::CrossF<ScalarType> > &Crosses)
{
    assert(Crosses.size()==TargetNormals.size());
    for (size_t i=0;i<Crosses.size();i++)
        Crosses[i].MakeDirectionNormalCoherent(TargetNormals[i]);
}

template < typename ScalarType >
void InitCrossFieldQualityAsAnisotropy(std::vector<Field::CrossF<ScalarType> > &Crosses,
                                       ScalarType filter_percentile=0.05,
                                       ScalarType minV=0.00001)
{

    if (Crosses.size()==0)return;

    ScalarType trimUp=1;
    ScalarType trimDown=1;

    std::vector<ScalarType> KVal;
    for (size_t i=0;i<Crosses.size();i++)
    {
        KVal.push_back(Crosses[i].K[0]);
        KVal.push_back(Crosses[i].K[1]);
    }
    assert(KVal.size()>0);

    std::sort(KVal.begin(),KVal.end());
    int percUp=int(floor((ScalarType)KVal.size()*(1-filter_percentile)+0.5));
    percUp=std::min(percUp,(int)KVal.size()-1);
    int percDown=int(floor((ScalarType)KVal.size()*filter_percentile+0.5));
    percDown=std::min(percDown,(int)KVal.size()-1);

    if (filter_percentile>0)
    {
        trimUp=KVal[percUp];
        trimDown=KVal[percDown];
    }
    else
    {
        trimUp=KVal.back();
        trimDown=KVal[0];
    }

    for (size_t i=0;i<Crosses.size();i++)
    {
        ScalarType N0=(Crosses[i].K[0]);
        ScalarType N1=(Crosses[i].K[1]);

        //clamp
        if (filter_percentile>0)
        {
            N0=std::min(N0,trimUp);
            N0=std::max(N0,trimDown);
            N1=std::min(N1,trimUp);
            N1=std::max(N1,trimDown);
        }

        ScalarType NMax=std::max(N0,N1)/trimUp;
        ScalarType NMin=std::min(N0,N1)/trimUp;

        if (NMax>1)NMax=1;
        if (NMax<-1)NMax=-1;

        if (NMin>1)NMin=1;
        if (NMin<-1)NMin=-1;

        ScalarType Ani=(NMax-NMin)/2;
        Ani=std::max(Ani,minV);

        Crosses[i].Q=Ani;
    }
    //    vcg::tri::UpdateQuality<MeshType>::FaceFromVertex(mesh);
}


template <class ScalarType>
void SetFaceCrossVectorFromVert(const std::vector<std::vector<int> > &Faces,
                                const std::vector<Geo::Point3<ScalarType> > &FacesNormals,
                                const std::vector<Field::CrossF<ScalarType> > &VertCurv,
                                std::vector<Field::CrossF<ScalarType> > &FaceCurv,
                                bool weightByQuality=true)
{
    //    MakeDirectionNormalCoherent()
    FaceCurv.resize(Faces.size());
    for (size_t i=0;i<Faces.size();i++)
    {
        std::vector<Field::CrossF<ScalarType> > VertCross;
        std::vector<ScalarType> W(Faces[i].size(),1);
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int IndexV=Faces[i][j];
            assert(IndexV<(int)VertCurv.size());
            VertCross.push_back(VertCurv[IndexV]);
            if (weightByQuality)
                W[j]=VertCurv[IndexV].Q;
        }
        //save old quality
        ScalarType OldQ=FaceCurv[i].Q;
        FaceCurv[i]=InterpolateCross(VertCross,FacesNormals[i],W);
        assert(fabs(FaceCurv[i].D[0].Norm()-1)<0.00001);
        assert(fabs(FaceCurv[i].D[1].Norm()-1)<0.00001);
        //restore old quality
        FaceCurv[i].Q=OldQ;
    }
    //make them normal coherent
    MakeDirectionNormalCoherent<ScalarType>(FacesNormals,FaceCurv);
}

template <class ScalarType>
void SetVertCrossFromFace(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                          const std::vector<std::vector<int> > &Faces,
                          const std::vector<Geo::Point3<ScalarType> > &VertNormals,
                          const std::vector<Field::CrossF<ScalarType> > &FaceCurv,
                          std::vector<Field::CrossF<ScalarType> > &VertCurv,
                          bool weightByQuality=true)
{
    VertCurv.resize(VertPos.size());
    std::vector<std::vector<size_t> > VertFaces(VertPos.size());
    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int IndexV=Faces[i][j];
            VertFaces[IndexV].push_back(i);
        }

    for (size_t i=0;i<VertFaces.size();i++)
    {
        std::vector<Field::CrossF<ScalarType> > FaceCross;
        std::vector<ScalarType> W(VertFaces[i].size(),1);
        for (size_t j=0;j<VertFaces[i].size();j++)
        {
            int IndexF=VertFaces[i][j];
            assert(IndexF<(int)FaceCurv.size());
            FaceCross.push_back(FaceCurv[IndexF]);
            if (weightByQuality)
            {
                W[j]=FaceCurv[IndexF].Q;
            }
        }
        //save old Position and quality
        Geo::Point3<ScalarType> OldPos=VertCurv[i].Pos;
        ScalarType OldQ=VertCurv[i].Q;

        VertCurv[i]=InterpolateCross(FaceCross,VertNormals[i],W);

        //restore
        VertCurv[i].Pos=OldPos;
        VertCurv[i].Q=OldQ;
    }
    //make them normal coherent
    //std::cout<<"Test0"<<std::endl;
    MakeDirectionNormalCoherent<ScalarType>(VertNormals,VertCurv);
    //std::cout<<"Test1"<<std::endl;
}


template <class ScalarType>
int CrossSingularValue(const std::vector<Field::CrossF<ScalarType> > &SortedFaceCrosses)
{
    int CurrD=0;
    int sizeCross=SortedFaceCrosses.size();
    //std::cout<<"size:"<<SortedFaceCrosses.size()<<std::endl;
    for (size_t i=0;i<SortedFaceCrosses.size();i++)
    {
        Field::CrossF<ScalarType> currF=SortedFaceCrosses[i];
        Field::CrossF<ScalarType> NextF=SortedFaceCrosses[(i+1)%sizeCross];
        currF.MakeTangentTo(NextF.OrthoDir());

        Geo::Point3<ScalarType> TargetDir=currF.GetDirection(CurrD);
        CurrD=NextF.BestDirectionI(TargetDir);
    }
    return (CurrD);
}

template <class ScalarType>
void GetFaceSingularities(const std::vector<Field::CrossF<ScalarType> > &FaceCrosses,
                          const std::vector<Geo::Point3<ScalarType> > &VertPos,
                          const std::vector<std::vector<int> > &Faces,
                          const std::vector<std::vector<int> > &NextF,
                          const std::vector<std::vector<int> > &NextE,
                          std::vector<int> &SingVert,
                          std::vector<int> &SingValue)
{
    assert(NextF.size()==Faces.size());
    assert(NextE.size()==Faces.size());

    std::vector<size_t> BorderV;
    //GetBorderVertices(Faces,NextF,NextE,BorderV);
    Geo::getBorderVert(Faces,NextF,BorderV);
    std::set<size_t> BorderVSet(BorderV.begin(),BorderV.end());

    SingVert.clear();
    SingValue.clear();
    std::vector<bool> Visited(VertPos.size(),false);
    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[j].size();j++)
        {
            int IndexV=Faces[i][j];

            if (Visited[IndexV])continue;
            Visited[IndexV]=true;

            //if border no singularity
            if (BorderVSet.count(IndexV)>0)continue;

            assert(IndexV<(int)Visited.size());

            Pos3 currP(i,j,IndexV);

            bool HasBorder;
            std::vector<Pos3 > PosFan;
            GetSortedFacesFaceFun(currP,Faces,NextF,NextE,PosFan,HasBorder);

            std::vector<Field::CrossF<ScalarType> > SortedFan;
            for (size_t k=0;k<PosFan.size();k++)
            {
                int IndexF=PosFan[k].F;
                assert(PosFan[k].V==IndexV);
                SortedFan.push_back(FaceCrosses[IndexF]);
            }
            int SingVal=CrossSingularValue(SortedFan);
            if (SingVal!=0)
            {
                SingVert.push_back(IndexV);
                SingValue.push_back(SingVal);
            }
            assert(!HasBorder);
        }
}

template <class ScalarType>
void AdaptFieldToMesh(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                      const std::vector<std::vector<int> > &Faces,
                      std::vector<Field::CrossF<ScalarType> > &VertCurv,
                      std::vector<Field::CrossF<ScalarType> > &FaceCurv)
{
    assert(VertCurv.size()==VertPos.size());
    assert(FaceCurv.size()==Faces.size());

    std::vector<Geo::Point3<ScalarType> > FaceNormals;
    std::vector<Geo::Point3<ScalarType> > VertNormals;
    Geo::ComputeNormals(VertPos,Faces,FaceNormals,VertNormals);

    for (size_t i=0;i<VertCurv.size();i++)
    {
        //reset center
        VertCurv[i].Pos=VertPos[i];
        //make it tangent
        VertCurv[i].MakeTangentTo(VertNormals[i]);
    }

    std::vector<Geo::Point3<ScalarType> > BaryF;
    ComputeFaceBarycenters(VertPos,Faces,BaryF);

    for (size_t i=0;i<Faces.size();i++)
    {
        //reset center
        FaceCurv[i].Pos=BaryF[i];
        //make it tangent
        FaceCurv[i].MakeTangentTo(FaceNormals[i]);
    }
}

template <class ScalarType>
void RotateCrossesToUV(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                       const std::vector<std::vector<int> > &Faces,
                       const std::vector<Geo::Point2<ScalarType> > &UVPos,
                       const std::vector<Field::CrossF<ScalarType> > &Crosses,
                       std::vector<Field::CrossF<ScalarType> > &UVCrosses)
{
    UVCrosses.clear();
    assert(Faces.size()==Crosses.size());
    for (size_t i=0;i<Faces.size();i++)
    {
        std::vector<Geo::Point3<ScalarType> > P3D;
        std::vector<Geo::Point3<ScalarType> > PUV;
        Geo::Point3<ScalarType> BaryF(0,0,0);
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int IndexV=Faces[i][j];
            P3D.push_back(VertPos[IndexV]);
            PUV.push_back(Geo::Point3<ScalarType>(UVPos[IndexV].X,UVPos[IndexV].Y,0));
            BaryF+=PUV.back();
        }
        BaryF/=Faces[i].size();
        Eigen::Matrix3<ScalarType> Matr;
        Geo::PointMatch<ScalarType>::ComputeBestRotation(PUV,P3D,Matr);
        Geo::Point3<ScalarType> D0Cross=Crosses[i].D[0];
        Geo::Point3<ScalarType> D1Cross=Crosses[i].D[1];
        D0Cross=Matr*D0Cross;
        D1Cross=Matr*D1Cross;
        ScalarType K0=Crosses[i].K[0];
        ScalarType K1=Crosses[i].K[1];
        Field::CrossF<ScalarType> RotCross(D0Cross,D1Cross,K0,K1,BaryF);
        RotCross.Q=Crosses[i].Q;
        UVCrosses.push_back(RotCross);
    }

}

template <class ScalarType>
void RedistributeAnisotropy(const std::vector<Field::CrossF<ScalarType> > &Source,
                            std::vector<Field::CrossF<ScalarType> > &Target)
{
    assert(Source.size()==Target.size());
    for (size_t i=0;i<Source.size();i++)
    {
//        Geo::Point3<ScalarType> Dir0S=Source[i].D[0];
//        Geo::Point3<ScalarType> Dir1S=Source[i].D[1];
//        Dir0S.Normalize();
//        Dir1S.Normalize();
//        ScalarType K0S=Source[i].K[0];
//        ScalarType K1S=Source[i].K[1];

//        Geo::Point3<ScalarType> Dir0T=Target[i].D[0];
//        Geo::Point3<ScalarType> Dir1T=Target[i].D[1];
//        Dir0T.Normalize();
//        Dir1T.Normalize();
//        ScalarType K0T=0;
//        ScalarType K1T=0;

//        ScalarType Dot0=fabs(Dot(Dir0S,Dir0T));
//        ScalarType Dot1=fabs(Dot(Dir0S,Dir1T));
//        ScalarType alpha0=Dot0/(Dot0+Dot1);
//        K0T+=alpha0*K0S;
//        K1T+=(1-alpha0)*K0S;

//        ScalarType Dot2=fabs(Dot(Dir1S,Dir0T));
//        ScalarType Dot3=fabs(Dot(Dir1S,Dir1T));
//        ScalarType alpha1=Dot2/(Dot2+Dot3);
//        K0T+=alpha1*K1S;
//        K1T+=(1-alpha1)*K1S;

//        Target[i].K[0]=K0T;
//        Target[i].K[1]=K1T;

        Geo::Point3<ScalarType> Dir0S=Source[i].D[0];
        Geo::Point3<ScalarType> Dir1S=Source[i].D[1];
        Dir0S.Normalize();
        Dir1S.Normalize();
        ScalarType K0S=Source[i].K[0];
        ScalarType K1S=Source[i].K[1];

        Geo::Point3<ScalarType> Dir0T=Target[i].D[0];
        Geo::Point3<ScalarType> Dir1T=Target[i].D[1];
        Dir0T.Normalize();
        Dir1T.Normalize();

        ScalarType Dot0=fabs(Dot(Dir0S,Dir0T));
        ScalarType Dot1=fabs(Dot(Dir0S,Dir1T));

        if (Dot0>Dot1)
        {
             Target[i].K[0]=K0S;
             Target[i].K[1]=K1S;
        }
        else
        {
            Target[i].K[0]=K1S;
            Target[i].K[1]=K0S;
        }
    }
}

#endif
