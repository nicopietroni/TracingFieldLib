#ifndef QUADRANGULATION
#define QUADRANGULATION

#include "IO/io_obj.h"
#include "Space_Query/triangle_grid_2D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "UV_param.h"
#include "UV_utils.h"
#include "box2.h"
#include "mesh_functions.h"
#include "point2.h"
#include "point3.h"
#include "quadrangulationLUT.h"
#include "segment2.h"
#include "triangle3.h"
#include "box3.h"
#include <cassert>
#include <cmath>
#include "connectivity.h"
#include <cstddef>
#include <iostream>
#include "mesh_mapping.h"
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType> class Quadrangulation {
public:
  // static void QuadrangulateRaster(const std::vector<Geo::Point2<ScalarType> >
  // &vert,
  //                                 const std::vector<std::vector<int> > &face,
  //                                 std::vector<Geo::Point3<ScalarType> >
  //                                 &bary, std::vector<int> &baryFIdx,
  //                                 std::vector<std::vector<int> > &quadF,
  //                                 ScalarType edgeSize=-1,bool
  //                                 IncludeExternal=true)
  // {
  //     baryFIdx.clear();
  //     bary.clear();
  //     quadF.clear();

  //     Geo::Box2<ScalarType> BB_2D;
  //     BB_2D.Init(vert);
  //     if (edgeSize<=0)
  //         edgeSize=BB_2D.Diag()*0.05;
  //     BB_2D.Inflate(edgeSize/2);

  //     Geo::TriangleSpatialIndex2D<ScalarType> Grid_2D;
  //     Grid_2D.Init(vert,face);

  //     int DivX=ceil(BB_2D.DimX()/edgeSize);
  //     int DivY=ceil(BB_2D.DimY()/edgeSize);
  //     DivX=std::max(DivX,1);
  //     DivY=std::max(DivY,1);

  //     std::map<Geo::Point2<int>,int> GridVert;
  //     for (int x_i=0;x_i<=DivX;x_i++)
  //         for (int y_i=0;y_i<=DivY;y_i++)
  //         {
  //             Geo::Point2<int> GridIDX(x_i,y_i);

  //             Geo::Point2<ScalarType> CurrP=BB_2D.Min;
  //             CurrP.X+=x_i*edgeSize;
  //             CurrP.Y+=y_i*edgeSize;

  //             std::vector<int> FacesIn;
  //             std::vector<Geo::Point3<ScalarType> > BaryCoords;
  //             Grid_2D.TrianglesOverPoint(vert,face,CurrP,FacesIn,BaryCoords);

  //             //no overlapping, should be bijective
  //             assert(FacesIn.size()<2);
  //             if (FacesIn.size()==0)continue;

  //             bary.push_back(BaryCoords[0]);
  //             baryFIdx.push_back(FacesIn[0]);
  //             GridVert[GridIDX]=bary.size()-1;
  //         }

  //     //then add the other on border if needed
  //     if (IncludeExternal)
  //     {
  //         std::vector<Geo::Point2<int> > AdditionalGridIdx;
  //         std::vector<int> AdditionalVertIdx;

  //         for (int x_i=0;x_i<=DivX;x_i++)
  //             for (int y_i=0;y_i<=DivY;y_i++)
  //             {
  //                 Geo::Point2<int> GridIDX(x_i,y_i);
  //                 if (GridVert.count(GridIDX)>0)continue;

  //                 //get closest neighbour to interpolate
  //                 Geo::Point2<int> ClosIdx(-1,-1);
  //                 int MinD=3;
  //                 for (int i=-1;i<=1;i++)
  //                     for (int j=-1;j<=1;j++)
  //                     {
  //                         Geo::Point2<int>
  //                         CurrIdx=GridIDX+Geo::Point2<int>(i,j); if
  //                         (GridVert.count(CurrIdx)==0)continue; int
  //                         CurrD=fabs(i)+fabs(j); if (CurrD>MinD)continue;
  //                         MinD=CurrD;
  //                         ClosIdx=CurrIdx;
  //                     }

  //                 //if no neighbour do nothing
  //                 if (ClosIdx==Geo::Point2<int>(-1,-1))continue;

  //                 assert(GridVert.count(ClosIdx)>0);

  //                 //otherwise get the closest barycenter coords
  //                 int IndexVKernel=GridVert[ClosIdx];
  //                 int IndexF=baryFIdx[IndexVKernel];
  //                 assert(IndexF>=0);
  //                 assert(IndexF<face.size());

  //                 Geo::Point2<ScalarType> P0=vert[face[IndexF][0]];
  //                 Geo::Point2<ScalarType> P1=vert[face[IndexF][1]];
  //                 Geo::Point2<ScalarType> P2=vert[face[IndexF][2]];
  //                 Traingle2<ScalarType> T(P0,P1,P2);

  //                 Geo::Point2<ScalarType> CurrP=BB_2D.Min;
  //                 CurrP.X+=x_i*edgeSize;
  //                 CurrP.Y+=y_i*edgeSize;

  //                 Geo::Point3<ScalarType> baryV;
  //                 BarycentricCoords(T,CurrP,baryV.X,baryV.Y,baryV.Z);

  //                 baryFIdx.push_back(IndexF);
  //                 bary.push_back(baryV);

  //                 AdditionalGridIdx.push_back(GridIDX);
  //                 AdditionalVertIdx.push_back(bary.size()-1);
  //             }

  //         for (size_t i=0;i<AdditionalVertIdx.size();i++)
  //             GridVert[AdditionalGridIdx[i]]=AdditionalVertIdx[i];

  //     }

  //     //then compose the quads
  //     for (int x_i=0;x_i<DivX;x_i++)
  //         for (int y_i=0;y_i<DivY;y_i++)
  //         {
  //             Geo::Point2<int> IdxV0(x_i,y_i);
  //             Geo::Point2<int> IdxV1(x_i+1,y_i);
  //             Geo::Point2<int> IdxV2(x_i+1,y_i+1);
  //             Geo::Point2<int> IdxV3(x_i,y_i+1);

  //             if (GridVert.count(IdxV0)==0)continue;
  //             int IndexV0=GridVert[IdxV0];

  //             if (GridVert.count(IdxV1)==0)continue;
  //             int IndexV1=GridVert[IdxV1];

  //             if (GridVert.count(IdxV2)==0)continue;
  //             int IndexV2=GridVert[IdxV2];

  //             if (GridVert.count(IdxV3)==0)continue;
  //             int IndexV3=GridVert[IdxV3];

  //             quadF.resize(quadF.size()+1);
  //             assert(IndexV0!=IndexV1);
  //             assert(IndexV1!=IndexV2);
  //             assert(IndexV2!=IndexV3);
  //             quadF.back().push_back(IndexV0);
  //             quadF.back().push_back(IndexV1);
  //             quadF.back().push_back(IndexV2);
  //             quadF.back().push_back(IndexV3);
  //         }
  // }

  // static void QuadrangulateRaster(const std::vector<Geo::Point2<ScalarType> >
  // &vert_2D,
  //                                 const std::vector<Geo::Point3<ScalarType> >
  //                                 &vert_3D, const
  //                                 std::vector<std::vector<int> > &face,
  //                                 std::vector<Geo::Point2<ScalarType> >
  //                                 &quad_vert2D,
  //                                 std::vector<Geo::Point3<ScalarType> >
  //                                 &quad_vert3D, std::vector<std::vector<int>
  //                                 > &quadF,ScalarType edgeSize=-1, bool
  //                                 IncludeExternal=true)
  // {
  //     std::vector<Geo::Point3<ScalarType> > bary;
  //     std::vector<int> baryFIdx;
  //     QuadrangulateRaster(vert_2D,face,bary,baryFIdx,quadF,edgeSize,IncludeExternal);
  //     assert(bary.size()==baryFIdx.size());

  //     //then interpolate coords
  //     quad_vert2D.clear();
  //     quad_vert3D.clear();
  //     for (size_t i=0;i<bary.size();i++)
  //     {
  //         int IndexF=baryFIdx[i];
  //         assert(IndexF>=0);
  //         assert(IndexF<face.size());

  //         int IndexV0=face[IndexF][0];
  //         int IndexV1=face[IndexF][1];
  //         int IndexV2=face[IndexF][2];

  //         Geo::Point3<ScalarType> P_3D[3];
  //         P_3D[0]=vert_3D[IndexV0];
  //         P_3D[1]=vert_3D[IndexV1];
  //         P_3D[2]=vert_3D[IndexV2];

  //         Geo::Point2<ScalarType> P_2D[3];
  //         P_2D[0]=vert_2D[IndexV0];
  //         P_2D[1]=vert_2D[IndexV1];
  //         P_2D[2]=vert_2D[IndexV2];

  //         Geo::Point3<ScalarType> baryV=bary[i];
  //         Geo::Point3<ScalarType> Pos3D=P_3D[0]*baryV.X+
  //                                       P_3D[1]*baryV.Y+
  //                                       P_3D[2]*baryV.Z;

  //         Geo::Point2<ScalarType> Pos2D=P_2D[0]*baryV.X+
  //                                       P_2D[1]*baryV.Y+
  //                                       P_2D[2]*baryV.Z;

  //         quad_vert2D.push_back(Pos2D);
  //         quad_vert3D.push_back(Pos3D);
  //     }
  // }

  // static void MultiComponentQuadrangulateRaster(const
  // std::vector<Geo::Point2<ScalarType> > &vert_2D,
  //                                               const
  //                                               std::vector<Geo::Point3<ScalarType>
  //                                               > &vert_3D, const
  //                                               std::vector<std::vector<int>
  //                                               > &face,
  //                                               std::vector<Geo::Point2<ScalarType>
  //                                               > &quad_vert2D,
  //                                               std::vector<Geo::Point3<ScalarType>
  //                                               > &quad_vert3D,
  //                                               std::vector<std::vector<int>
  //                                               > &quadF,ScalarType
  //                                               edgeSize=-1, bool
  //                                               IncludeExternal=true)
  // {
  //     quad_vert2D.clear();
  //     quad_vert3D.clear();
  //     quadF.clear();

  //     //split into connected component
  //     std::vector<std::vector<Geo::Point3<ScalarType> > > SubMeshesV;
  //     std::vector<std::vector<std::vector<int> > > SubMeshesF;
  //     std::vector<std::vector<int> > OriginalV,OriginalF;
  //     SplitByFFConnectedComponenent<ScalarType>(vert_3D,face,SubMeshesV,SubMeshesF,OriginalV,OriginalF);

  //     //check
  //     assert(SubMeshesF.size()==SubMeshesV.size());
  //     assert(SubMeshesF.size()==OriginalV.size());

  //     if (edgeSize==-1)
  //     {
  //         Geo::Box3<ScalarType> BB;
  //         BB.Init(vert_3D);
  //         edgeSize=BB.Diag()*0.03;
  //     }

  //     for (size_t i=0;i<SubMeshesV.size();i++)
  //     {
  //         std::vector<Geo::Point2<ScalarType> > vert_2DComp;
  //         for (size_t j=0;j<OriginalV[i].size();j++)
  //             vert_2DComp.push_back(vert_2D[OriginalV[i][j]]);

  //         assert(vert_2DComp.size()==SubMeshesV[i].size());

  //         std::vector<Geo::Point2<ScalarType> > quad_vert2D_Comp;
  //         std::vector<Geo::Point3<ScalarType> > quad_vert3D_Comp;
  //         std::vector<std::vector<int> > quadFComp;

  //         //quadrangulate components
  //         QuadrangulateRaster(vert_2DComp,SubMeshesV[i],SubMeshesF[i],quad_vert2D_Comp,quad_vert3D_Comp,quadFComp,edgeSize,IncludeExternal);

  //         //sum up exhisting vertices
  //         for (size_t j=0;j<quadFComp.size();j++)
  //             for (size_t k=0;k<quadFComp[j].size();k++)
  //                 quadFComp[j][k]+=quad_vert3D.size();

  //         quadF.insert(quadF.end(),quadFComp.begin(),quadFComp.end());

  //         //copy back to original
  //         quad_vert2D.insert(quad_vert2D.end(),quad_vert2D_Comp.begin(),quad_vert2D_Comp.end());
  //         quad_vert3D.insert(quad_vert3D.end(),quad_vert3D_Comp.begin(),quad_vert3D_Comp.end());

  //     }
  // }

  static void SampleRegular(
      const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
      const std::vector<Geo::Point3<ScalarType>> &PatchVert3D,
      const std::vector<std::vector<int>> &Patch_Faces, const int SizeX,
      const int SizeY,
      const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
      const std::vector<std::vector<Geo::Point2<ScalarType>>> &UVSideSubd,
      std::vector<Geo::Point3<ScalarType>> &QuadVert,
      std::vector<std::vector<int>> &QuadFaces)

  {
    std::vector<Geo::Segment2<ScalarType>> Seg_U, Seg_V;

    for (size_t i = 0; i < SizeX; i++) {
      Geo::Point2<ScalarType> P0 = UVSideSubd[0][i];
      Geo::Point2<ScalarType> P1 = UVSideSubd[2][SizeX - i - 1];
      Geo::Segment2<ScalarType> Seg2D(P0, P1);
      Seg_U.push_back(Seg2D);
    }

    for (size_t i = 0; i < SizeY; i++) {
      Geo::Point2<ScalarType> P0 = UVSideSubd[1][i];
      Geo::Point2<ScalarType> P1 = UVSideSubd[3][SizeY - i - 1];
      Geo::Segment2<ScalarType> Seg2D(P0, P1);
      Seg_V.push_back(Seg2D);
    }

    std::vector<std::vector<Geo::Point3<ScalarType>>> Sampled3D;
    std::vector<std::vector<int>> SampledIdx;
    int NumRow = Seg_V.size();
    int NumCol = Seg_U.size();

    Sampled3D.resize(NumRow);
    SampledIdx.resize(NumRow);
    for (size_t i = 0; i < Sampled3D.size(); i++) {
      Sampled3D[i].resize(NumCol, Geo::Point3<ScalarType>(0, 0, 0));
      SampledIdx[i].resize(NumCol, -1);
    }

    // add the first row
    Sampled3D[0] = SideSubd[0];
    // add the second row
    Sampled3D.back() = SideSubd[2];
    std::reverse(Sampled3D.back().begin(), Sampled3D.back().end());

    // add the first column
    for (size_t u = 0; u < SideSubd[1].size(); u++)
      Sampled3D[u].back() = SideSubd[1][u];

    // and the third
    for (size_t u = 0; u < SideSubd[3].size(); u++)
      Sampled3D[u][0] = SideSubd[3][SideSubd[3].size() - u - 1];

    // add the last column
    // WriteOBJ("test.obj",VertUV,Patch_Faces);
    Geo::TriangleSpatialIndex2D<ScalarType> GridUV;
    GridUV.Init(PatchVertUV, Patch_Faces);

    for (size_t u = 1; u < Seg_U.size() - 1; u++)
      for (size_t v = 1; v < Seg_V.size() - 1; v++) {
        Geo::Segment2<ScalarType> S0 = Seg_U[u];
        Geo::Segment2<ScalarType> S1 = Seg_V[v];
        Geo::Point2<ScalarType> IntP;
        bool do_inters = Intersect(S0, S1, IntP);
        assert(do_inters);

        // then retrieve the 3D point
        Geo::Point3<ScalarType> BaryCoords;
        int IndexF = -1;
        GridUV.FirstTrianglesOverPoint(PatchVertUV, Patch_Faces, IntP, IndexF,
                                       BaryCoords);
        if (IndexF == -1) {
          IndexF = Geo::UVUtils<ScalarType>::InsideFace(
              PatchVertUV, Patch_Faces, IntP, BaryCoords);
        }

        assert(IndexF != -1);

        int IndexV0 = Patch_Faces[IndexF][0];
        int IndexV1 = Patch_Faces[IndexF][1];
        int IndexV2 = Patch_Faces[IndexF][2];

        Geo::Point3<ScalarType> P0 = PatchVert3D[IndexV0];
        Geo::Point3<ScalarType> P1 = PatchVert3D[IndexV1];
        Geo::Point3<ScalarType> P2 = PatchVert3D[IndexV2];

        Geo::Point3<ScalarType> InterpP =
            P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
        Sampled3D[v][u] = InterpP;
      }

    for (size_t i = 0; i < Sampled3D.size(); i++)
      for (size_t j = 0; j < Sampled3D[i].size(); j++) {
        QuadVert.push_back(Sampled3D[i][j]);
        SampledIdx[i][j] = QuadVert.size() - 1;
      }

    // then add faces
    for (size_t i = 0; i < SampledIdx.size() - 1; i++)
      for (size_t j = 0; j < SampledIdx[i].size() - 1; j++) {
        int IndexV0 = SampledIdx[i][j];
        int IndexV1 = SampledIdx[i + 1][j];
        int IndexV2 = SampledIdx[i + 1][j + 1];
        int IndexV3 = SampledIdx[i][j + 1];
        QuadFaces.resize(QuadFaces.size() + 1);
        QuadFaces.back().push_back(IndexV0);
        QuadFaces.back().push_back(IndexV1);
        QuadFaces.back().push_back(IndexV2);
        QuadFaces.back().push_back(IndexV3);
      }
  }

  static bool SampleWithLUT(
      const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
      const std::vector<Geo::Point3<ScalarType>> &PatchVert3D,
      const std::vector<std::vector<int>> &Patch_Faces,
      const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
      const std::vector<std::vector<Geo::Point2<ScalarType>>> &UVSideSubd,
      const QuadrangulationLUT<ScalarType> &QLut,
      std::vector<Geo::Point3<ScalarType>> &QuadVert,
      std::vector<std::vector<int>> &QuadFaces)

  {

    std::vector<Geo::Point2<ScalarType>> QuadVertUV;
    bool found = QLut.GetQuadrangulation(UVSideSubd, QuadVertUV, QuadFaces);

    if (!found)
      return false;

    // std::cout<<"DEDE"<<std::endl;
    //     WriteOBJ("Test0.obj",PatchVertUV, Patch_Faces);
    //     WriteOBJ("Test1.obj",QuadVertUV, QuadFaces);
    //     exit(0);
    std::vector<Geo::Point3<ScalarType>> AllSide3D;
    for (size_t i = 0; i < SideSubd.size(); i++)
      AllSide3D.insert(AllSide3D.end(), SideSubd[i].begin(), SideSubd[i].end());

    AllSide3D.erase(std::unique(AllSide3D.begin(), AllSide3D.end()),
                    AllSide3D.end());
    AllSide3D.pop_back();

    QuadVert.resize(QuadVertUV.size(), Geo::Point3<ScalarType>(0, 0, 0));

    // set boundary
    assert(AllSide3D.size() <= QuadVert.size());

    int BorderNum = AllSide3D.size();
    // std::cout<<"Size0 "<<BorderNum<<std::endl;
    // std::cout<<"Size1 "<<QuadVert.size()<<std::endl;

    for (size_t i = 0; i < BorderNum; i++)
      QuadVert[i] = AllSide3D[i];

    // initialize the grid for sampling
    Geo::TriangleSpatialIndex2D<ScalarType> GridUV;
    GridUV.Init(PatchVertUV, Patch_Faces);
    for (size_t i = BorderNum; i < QuadVert.size(); i++) {
      Geo::Point2<ScalarType> UVPos = QuadVertUV[i];

      Geo::Point3<ScalarType> BaryCoords;
      int IndexF = -1;
      GridUV.FirstTrianglesOverPoint(PatchVertUV, Patch_Faces, QuadVertUV[i],
                                     IndexF, BaryCoords);
      if (IndexF == -1) {
        std::cout << "DEDE" << std::endl;
        WriteOBJ("Test0.obj", PatchVertUV, Patch_Faces);
        WriteOBJ("Test1.obj", QuadVertUV, QuadFaces);
        // QuadVertUV
        exit(0);
      }
      assert(IndexF != -1);

      int IndexV0 = Patch_Faces[IndexF][0];
      int IndexV1 = Patch_Faces[IndexF][1];
      int IndexV2 = Patch_Faces[IndexF][2];

      Geo::Point3<ScalarType> P0 = PatchVert3D[IndexV0];
      Geo::Point3<ScalarType> P1 = PatchVert3D[IndexV1];
      Geo::Point3<ScalarType> P2 = PatchVert3D[IndexV2];

      Geo::Point3<ScalarType> InterpP =
          P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
      QuadVert[i] = InterpP;
    }

    // // then set the rest
    // for (size_t i = 0; i < AllSide3D.size(); i++)
    //   QuadVert[i] = AllSide3D[i];
    return true;
  }

  static bool QuadrangulateQuadPatchWithSubs(
      const std::vector<Geo::Point3<ScalarType>> &PatchVert,
      const std::vector<std::vector<int>> &Patch_Faces,
      const std::vector<int> &PatchCorners,
      const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
      const QuadrangulationLUT<ScalarType> &QLut,
      std::vector<Geo::Point3<ScalarType>> &QuadVert,
      std::vector<std::vector<int>> &QuadFaces) {

    QuadVert.clear();
    QuadFaces.clear();

    assert(PatchCorners.size() == 4);
    // if not a quad return false
    if (SideSubd.size() != 4)
      return false;
    if (PatchCorners.size() != 4)
      return false;
    std::vector<Geo::Point2<ScalarType>> VertUV;

    // paremetrize the patch also considering corners
    Geo::UVParam<ScalarType>::TutteEmbedding(PatchVert, Patch_Faces,
                                             PatchCorners, VertUV);
    Geo::TriangleSpatialIndex<ScalarType> Grid;
    Grid.Init(PatchVert, Patch_Faces);
    Geo::Box3<ScalarType> BB;
    BB.Init(PatchVert);

    std::vector<std::vector<Geo::Point2<ScalarType>>> UVSideSubd;
    for (size_t i = 0; i < SideSubd.size(); i++) {
      UVSideSubd.resize(UVSideSubd.size() + 1);
      for (size_t j = 0; j < SideSubd[i].size(); j++) {
        Geo::Point3<ScalarType> Pos3D = SideSubd[i][j];
        int closestF;
        Geo::Point3<ScalarType> closestPt;
        Grid.GridClosest(PatchVert, Patch_Faces, Pos3D, BB.Diag(), closestF,
                         closestPt);

        Geo::Traingle3<ScalarType> CloseT =
            TriangleOfFace<ScalarType>(PatchVert, Patch_Faces, closestF);
        ScalarType alpha, beta, gamma;
        Geo::BarycentricCoords(CloseT, closestPt, alpha, beta, gamma);

        assert(closestF >= 0);
        assert(closestF < Patch_Faces.size());
        int IndexV0 = Patch_Faces[closestF][0];
        int IndexV1 = Patch_Faces[closestF][1];
        int IndexV2 = Patch_Faces[closestF][2];
        Geo::Point2<ScalarType> PosUV0 = VertUV[IndexV0];
        Geo::Point2<ScalarType> PosUV1 = VertUV[IndexV1];
        Geo::Point2<ScalarType> PosUV2 = VertUV[IndexV2];
        Geo::Point2<ScalarType> UVInterp =
            PosUV0 * alpha + PosUV1 * beta + PosUV2 * gamma;
        UVSideSubd.back().push_back(UVInterp);
      }
    }

    size_t Size0 = SideSubd[0].size();
    size_t Size1 = SideSubd[1].size();
    size_t Size2 = SideSubd[2].size();
    size_t Size3 = SideSubd[3].size();

    bool IsRegular = ((Size0 == Size2) && (Size1 == Size3));
    bool has_sampled = true;
    if (IsRegular) {
      SampleRegular(VertUV, PatchVert, Patch_Faces, Size0, Size1, SideSubd,
                    UVSideSubd, QuadVert, QuadFaces);

      has_sampled = true;
    } else {
      has_sampled = SampleWithLUT(VertUV, PatchVert, Patch_Faces, SideSubd,
                                  UVSideSubd, QLut, QuadVert, QuadFaces);
    }
    return has_sampled;
  }

public:
  static void QuadrangulateQuadPatchsWithSubs(
      const std::vector<Geo::Point3<ScalarType>> &Vert_3D,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<int>> &PatchFacesIdx,
      const std::vector<std::vector<int>> &PatchCornerIdx,
      const std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &SideSamples3D,
      const QuadrangulationLUT<ScalarType> &QLut,
      std::vector<Geo::Point3<ScalarType>> &QuadVert,
      std::vector<std::vector<int>> &QuadFaces, std::vector<int> &QuadPatch,
      bool save_debug = false) {

    QuadVert.clear();
    QuadFaces.clear();
    QuadPatch.clear();

    assert(PatchFacesIdx.size() == PatchCornerIdx.size());

    for (size_t i = 0; i < PatchFacesIdx.size(); i++) {
      // for (size_t i = 0; i <1; i++) {
      assert(PatchCornerIdx[i].size() == SideSamples3D[i].size());
      if (PatchCornerIdx[i].size() != 4)
        continue;

      std::vector<Geo::Point3<ScalarType>> CurrPatchVert;
      std::vector<std::vector<int>> CurrPatchFaces;
      std::vector<int> OriginalVert, OriginalElem;
      CopySubMesh<ScalarType>(Vert_3D, Faces, PatchFacesIdx[i], CurrPatchVert,
                              CurrPatchFaces, OriginalVert, OriginalElem);

      std::vector<int> CurrPatchCorners;
      for (size_t j = 0; j < OriginalVert.size(); j++) {
        std::vector<int>::const_iterator it =
            std::find(PatchCornerIdx[i].begin(), PatchCornerIdx[i].end(),
                      OriginalVert[j]);
        if (it == PatchCornerIdx[i].end())
          continue;
        CurrPatchCorners.push_back(j);
      }
      assert(CurrPatchCorners.size() == PatchCornerIdx[i].size());
      std::vector<Geo::Point3<ScalarType>> Curr_Patch_QuadVert;
      std::vector<std::vector<int>> Curr_Patch_QuadFaces;

      QuadrangulateQuadPatchWithSubs(CurrPatchVert, CurrPatchFaces,
                                     CurrPatchCorners, SideSamples3D[i], QLut,
                                     Curr_Patch_QuadVert, Curr_Patch_QuadFaces);
      if (save_debug) {
        std::string out_file = std::string("quad_mesh_") +
                               std::to_string(i) + 
                               std::string(".obj");
        WriteOBJ(out_file.c_str(),Curr_Patch_QuadVert, Curr_Patch_QuadFaces);
      }
      
      AppendMesh<ScalarType>(QuadVert, QuadFaces, Curr_Patch_QuadVert,
                             Curr_Patch_QuadFaces);

      QuadPatch.resize(QuadFaces.size(), i);
    }
  }
};

} // namespace Geo
#endif
