#ifndef QUADRANGULATION
#define QUADRANGULATION

#include <cmath>
#include <vector>
#include "box2.h"
#include "point3.h"
#include "Space_Query/triangle_grid_2D.h"
#include <map>
#include <connectivity.h>
#include <box3.h>

namespace Geo
{

template <class ScalarType>
class Quadrangulation
{
public:

    static void QuadrangulateRaster(const std::vector<Geo::Point2<ScalarType> > &vert,
                                    const std::vector<std::vector<int> > &face,
                                    std::vector<Geo::Point3<ScalarType> > &bary,
                                    std::vector<int> &baryFIdx,
                                    std::vector<std::vector<int> > &quadF,
                                    ScalarType edgeSize=-1,bool IncludeExternal=true)
    {
        baryFIdx.clear();
        bary.clear();
        quadF.clear();

        Geo::Box2<ScalarType> BB_2D;
        BB_2D.Init(vert);
        if (edgeSize<=0)
            edgeSize=BB_2D.Diag()*0.05;
        BB_2D.Inflate(edgeSize);

        Geo::TriangleSpatialIndex2D<ScalarType> Grid_2D;
        Grid_2D.Init(vert,face);

        int DivX=ceil(BB_2D.DimX()/edgeSize);
        int DivY=ceil(BB_2D.DimY()/edgeSize);
        DivX=std::max(DivX,1);
        DivY=std::max(DivY,1);

        std::map<Geo::Point2<int>,int> GridVert;
        for (int x_i=0;x_i<DivX;x_i++)
            for (int y_i=0;y_i<DivY;y_i++)
            {
                Geo::Point2<int> GridIDX(x_i,y_i);

                Geo::Point2<ScalarType> CurrP=BB_2D.Min;
                CurrP.X+=x_i*edgeSize;
                CurrP.Y+=y_i*edgeSize;

                std::vector<int> FacesIn;
                std::vector<Geo::Point3<ScalarType> > BaryCoords;
                Grid_2D.TrianglesOverPoint(vert,face,CurrP,FacesIn,BaryCoords);

                //no overlapping, should be bijective
                assert(FacesIn.size()<2);
                if (FacesIn.size()==0)continue;

                bary.push_back(BaryCoords[0]);
                baryFIdx.push_back(FacesIn[0]);
                GridVert[GridIDX]=bary.size()-1;
            }

        //then add the other on border if needed
        if (IncludeExternal)
        {
            std::vector<Geo::Point2<int> > AdditionalGridIdx;
            std::vector<int> AdditionalVertIdx;

            for (int x_i=0;x_i<DivX;x_i++)
                for (int y_i=0;y_i<DivY;y_i++)
                {
                    Geo::Point2<int> GridIDX(x_i,y_i);
                    if (GridVert.count(GridIDX)>0)continue;

                    //get closest neighbour to interpolate
                    Geo::Point2<int> ClosIdx(-1,-1);
                    int MinD=3;
                    for (int i=-1;i<=1;i++)
                        for (int j=-1;j<=1;j++)
                        {
                            Geo::Point2<int> CurrIdx=GridIDX+Geo::Point2<int>(i,j);
                            if (GridVert.count(CurrIdx)==0)continue;
                            int CurrD=fabs(i)+fabs(j);
                            if (CurrD>MinD)continue;
                            MinD=CurrD;
                            ClosIdx=CurrIdx;
                        }

                    //if no neighbour do nothing
                    if (ClosIdx==Geo::Point2<int>(-1,-1))continue;

                    assert(GridVert.count(ClosIdx)>0);

                    //otherwise get the closest barycenter coords
                    int IndexVKernel=GridVert[ClosIdx];
                    int IndexF=baryFIdx[IndexVKernel];
                    assert(IndexF>=0);
                    assert(IndexF<face.size());

                    Geo::Point2<ScalarType> P0=vert[face[IndexF][0]];
                    Geo::Point2<ScalarType> P1=vert[face[IndexF][1]];
                    Geo::Point2<ScalarType> P2=vert[face[IndexF][2]];
                    Traingle2<ScalarType> T(P0,P1,P2);

                    Geo::Point2<ScalarType> CurrP=BB_2D.Min;
                    CurrP.X+=x_i*edgeSize;
                    CurrP.Y+=y_i*edgeSize;


                    Geo::Point3<ScalarType> baryV;
                    BarycentricCoords(T,CurrP,baryV.X,baryV.Y,baryV.Z);

                    baryFIdx.push_back(IndexF);
                    bary.push_back(baryV);

                    AdditionalGridIdx.push_back(GridIDX);
                    AdditionalVertIdx.push_back(bary.size()-1);
                }

            for (size_t i=0;i<AdditionalVertIdx.size();i++)
                GridVert[AdditionalGridIdx[i]]=AdditionalVertIdx[i];

        }

        //then compose the quads
        for (int x_i=0;x_i<DivX;x_i++)
            for (int y_i=0;y_i<DivY;y_i++)
            {
                Geo::Point2<int> IdxV0(x_i,y_i);
                Geo::Point2<int> IdxV1(x_i+1,y_i);
                Geo::Point2<int> IdxV2(x_i+1,y_i+1);
                Geo::Point2<int> IdxV3(x_i,y_i+1);

                if (GridVert.count(IdxV0)==0)continue;
                int IndexV0=GridVert[IdxV0];

                if (GridVert.count(IdxV1)==0)continue;
                int IndexV1=GridVert[IdxV1];

                if (GridVert.count(IdxV2)==0)continue;
                int IndexV2=GridVert[IdxV2];

                if (GridVert.count(IdxV3)==0)continue;
                int IndexV3=GridVert[IdxV3];

                quadF.resize(quadF.size()+1);
                assert(IndexV0!=IndexV1);
                assert(IndexV1!=IndexV2);
                assert(IndexV2!=IndexV3);
                quadF.back().push_back(IndexV0);
                quadF.back().push_back(IndexV1);
                quadF.back().push_back(IndexV2);
                quadF.back().push_back(IndexV3);
            }
    }

    static void QuadrangulateRaster(const std::vector<Geo::Point2<ScalarType> > &vert_2D,
                                    const std::vector<Geo::Point3<ScalarType> > &vert_3D,
                                    const std::vector<std::vector<int> > &face,
                                    std::vector<Geo::Point2<ScalarType> > &quad_vert2D,
                                    std::vector<Geo::Point3<ScalarType> > &quad_vert3D,
                                    std::vector<std::vector<int> > &quadF,ScalarType edgeSize=-1,
                                    bool IncludeExternal=true)
    {
        std::vector<Geo::Point3<ScalarType> > bary;
        std::vector<int> baryFIdx;
        QuadrangulateRaster(vert_2D,face,bary,baryFIdx,quadF,edgeSize,IncludeExternal);
        assert(bary.size()==baryFIdx.size());

        //then interpolate coords
        quad_vert2D.clear();
        quad_vert3D.clear();
        for (size_t i=0;i<bary.size();i++)
        {
            int IndexF=baryFIdx[i];
            assert(IndexF>=0);
            assert(IndexF<face.size());

            int IndexV0=face[IndexF][0];
            int IndexV1=face[IndexF][1];
            int IndexV2=face[IndexF][2];

            Geo::Point3<ScalarType> P_3D[3];
            P_3D[0]=vert_3D[IndexV0];
            P_3D[1]=vert_3D[IndexV1];
            P_3D[2]=vert_3D[IndexV2];

            Geo::Point2<ScalarType> P_2D[3];
            P_2D[0]=vert_2D[IndexV0];
            P_2D[1]=vert_2D[IndexV1];
            P_2D[2]=vert_2D[IndexV2];

            Geo::Point3<ScalarType> baryV=bary[i];
            Geo::Point3<ScalarType> Pos3D=P_3D[0]*baryV.X+
                                          P_3D[1]*baryV.Y+
                                          P_3D[2]*baryV.Z;

            Geo::Point2<ScalarType> Pos2D=P_2D[0]*baryV.X+
                                          P_2D[1]*baryV.Y+
                                          P_2D[2]*baryV.Z;

            quad_vert2D.push_back(Pos2D);
            quad_vert3D.push_back(Pos3D);
        }
    }

    static void MultiComponentQuadrangulateRaster(const std::vector<Geo::Point2<ScalarType> > &vert_2D,
                                                  const std::vector<Geo::Point3<ScalarType> > &vert_3D,
                                                  const std::vector<std::vector<int> > &face,
                                                  std::vector<Geo::Point2<ScalarType> > &quad_vert2D,
                                                  std::vector<Geo::Point3<ScalarType> > &quad_vert3D,
                                                  std::vector<std::vector<int> > &quadF,ScalarType edgeSize=-1,
                                                  bool IncludeExternal=true)
    {
        quad_vert2D.clear();
        quad_vert3D.clear();
        quadF.clear();

        //split into connected component
        std::vector<std::vector<Geo::Point3<ScalarType> > > SubMeshesV;
        std::vector<std::vector<std::vector<int> > > SubMeshesF;
        std::vector<std::vector<int> > OriginalV,OriginalF;
        SplitByFFConnectedComponenent<ScalarType>(vert_3D,face,SubMeshesV,SubMeshesF,OriginalV,OriginalF);

        //check
        assert(SubMeshesF.size()==SubMeshesV.size());
        assert(SubMeshesF.size()==OriginalV.size());

        if (edgeSize==-1)
        {
            Geo::Box3<ScalarType> BB;
            BB.Init(vert_3D);
            edgeSize=BB.Diag()*0.03;
        }

        for (size_t i=0;i<SubMeshesV.size();i++)
        {
            std::vector<Geo::Point2<ScalarType> > vert_2DComp;
            for (size_t j=0;j<OriginalV[i].size();j++)
                vert_2DComp.push_back(vert_2D[OriginalV[i][j]]);

            assert(vert_2DComp.size()==SubMeshesV[i].size());

            std::vector<Geo::Point2<ScalarType> > quad_vert2D_Comp;
            std::vector<Geo::Point3<ScalarType> > quad_vert3D_Comp;
            std::vector<std::vector<int> > quadFComp;

            //quadrangulate components
            QuadrangulateRaster(vert_2DComp,SubMeshesV[i],SubMeshesF[i],quad_vert2D_Comp,quad_vert3D_Comp,quadFComp,edgeSize,IncludeExternal);


            //sum up exhisting vertices
            for (size_t j=0;j<quadFComp.size();j++)
                for (size_t k=0;k<quadFComp[j].size();k++)
                    quadFComp[j][k]+=quad_vert3D.size();

            quadF.insert(quadF.end(),quadFComp.begin(),quadFComp.end());

            //copy back to original
            quad_vert2D.insert(quad_vert2D.end(),quad_vert2D_Comp.begin(),quad_vert2D_Comp.end());
            quad_vert3D.insert(quad_vert3D.end(),quad_vert3D_Comp.begin(),quad_vert3D_Comp.end());

        }
    }
};

}
#endif
