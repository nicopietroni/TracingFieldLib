#ifndef SMOOTH_WITH_PATCHES
#define SMOOTH_WITH_PATCHES

#include "IO/io_obj.h"
#include "UV_param.h"
#include "UV_utils.h"
#include "connectivity.h"
#include "point3.h"
#include "rectify_features.h"
#include "reproject_mesh.h"
#include "smooth_mesh.h"
#include "solve/least_squares_lin_constraint_solver.h"
#include "solve/linear_constraint.h"
#include "solve/linear_constraints_functions.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cotangent.h>
#include <cstddef>
#include <iostream>
#include <limits>
#include <mesh_functions.h>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class InPatchSmooth {

  static void
  FindVertPatches(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces,
                  const std::vector<int> &FacePatch,
                  std::vector<std::vector<int>> &VertPatches) {
    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VFIndex;
    getVFAdjacency(Faces, VF, VFIndex);
    std::vector<bool> IsBorderV;
    getBorderVert(Faces, IsBorderV);

    std::set<std::vector<int>> VertPatchesSet;

    for (size_t i = 0; i < VF.size(); i++) {
      std::set<int> CurrVPatches;
      for (size_t j = 0; j < VF[i].size(); j++) {
        int IndexF = VF[i][j];
        assert(IndexF >= 0);
        assert(IndexF < FacePatch.size());
        CurrVPatches.insert(FacePatch[IndexF]);
      }
      if (CurrVPatches.size() < 2)
        continue;
      bool IsVert = (CurrVPatches.size() > 2) ||
                    ((CurrVPatches.size() == 2) && IsBorderV[i]);
      if (!IsVert)
        continue;
      std::vector<int> CurrP(CurrVPatches.begin(), CurrVPatches.end());
      std::sort(CurrP.begin(), CurrP.end());
      VertPatchesSet.insert(CurrP);
    }

    VertPatches = std::vector<std::vector<int>>(VertPatchesSet.begin(),
                                                VertPatchesSet.end());
  }

  // static void
  // FindVertFromPatches(const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                     const std::vector<std::vector<int>> &Faces,
  //                     const std::vector<int> &FacePatch,
  //                     std::vector<int> &VertIdx) {
  //   std::vector<std::vector<int>> VF;
  //   std::vector<std::vector<int>> VFIndex;
  //   getVFAdjacency(Faces, VF, VFIndex);
  //   std::vector<bool> IsBorderV;
  //   getBorderVert(Faces, IsBorderV);

  //   for (size_t i = 0; i < VF.size(); i++) {
  //     std::set<int> CurrVPatches;
  //     for (size_t j = 0; j < VF[i].size(); j++) {
  //       int IndexF = VF[i][j];
  //       assert(IndexF >= 0);
  //       assert(IndexF < FacePatch.size());
  //       CurrVPatches.insert(FacePatch[IndexF]);
  //     }
  //     if (CurrVPatches.size() < 2)
  //       continue;
  //     bool IsVert = (CurrVPatches.size() > 2) ||
  //                   ((CurrVPatches.size() == 2) && IsBorderV[i]);
  //     if (!IsVert)
  //       continue;

  //     VertIdx.push_back(i);
  //   }
  // }

  static void SmoothRectify(std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &QuadFaces,
                            const std::vector<std::pair<int, int>> &Features,
                            const std::vector<int> &QuadCorners) {
    // WriteOBJ("test_0.obj",VertPos, QuadFaces);
    //  bool hasF=(Features.size()>0);
    //  if (!hasF)return;
    // get the corners
    std::vector<int> BorderV;
    getBorderVert(QuadFaces, BorderV);

    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VFIndex;
    getVFAdjacency(QuadFaces, VF, VFIndex);

    std::vector<int> RectifyCorners = QuadCorners;
    for (size_t i = 0; i < BorderV.size(); i++) {
      int IndexV = BorderV[i];
      if (VF[IndexV].size() == 2)
        continue;
      RectifyCorners.push_back(IndexV);
    }
    RectifyCorners.insert(RectifyCorners.end(), QuadCorners.begin(),
                          QuadCorners.end());

    std::vector<std::vector<int>> TriFaces;
    SplitQuadIntoTris<ScalarType>(VertPos, QuadFaces, TriFaces);

    // std::vector<std::pair<int,int>> FeaturesB;
    // getBorderEdges(TriFaces,FeaturesB);

    std::vector<std::pair<int, int>> RectifyF; //=FeaturesB;
    RectifyF.insert(RectifyF.end(), Features.begin(), Features.end());

    for (size_t i = 0; i < RectifyF.size(); i++) {
      if (RectifyF[i].second > RectifyF[i].first)
        std::swap(RectifyF[i].first, RectifyF[i].second);
    }
    std::sort(RectifyF.begin(), RectifyF.end());
    RectifyF.erase(std::unique(RectifyF.begin(), RectifyF.end()),
                   RectifyF.end());

    std::vector<Geo::Point3<ScalarType>> OldVertPos = VertPos;

    // RectifyFeatures<ScalarType>(VertPos,TriFaces,RectifyF,
    //                 RectifyCorners);

    // then smooth
    Geo::SmoothParam<ScalarType> SParamQuadSmooth;

    // SParamQuadSmooth.ReprojVert = &VertPos;
    // SParamQuadSmooth.ReprojFaces = &TriFaces;
    SParamQuadSmooth.Damp = 0.5;
    SParamQuadSmooth.post_project = false;
    SParamQuadSmooth.NumIte = 1;
    SParamQuadSmooth.SType = Geo::SMLaplacian;
    SParamQuadSmooth.Features = RectifyF;

    Geo::getBorderVert(QuadFaces, SParamQuadSmooth.FixedVert);
    SParamQuadSmooth.FixedVert.insert(SParamQuadSmooth.FixedVert.end(),
                                      QuadCorners.begin(), QuadCorners.end());

    std::vector<Geo::Point3<ScalarType>> SmoothVertPos = VertPos;
    Geo::SmoothLaplacian<ScalarType>(SmoothVertPos, QuadFaces,
                                     SParamQuadSmooth);

    // project back to the original mesh
    std::vector<int> IndexF;
    std::vector<Geo::Point3<ScalarType>> BaryCoord;
    ReprojectBasis(SmoothVertPos, VertPos, TriFaces, IndexF, BaryCoord);

    // //then deform it back
    // WriteOBJ("test_1.obj",SmoothVertPos,QuadFaces);

    for (size_t i = 0; i < SmoothVertPos.size(); i++) {
      int currF = IndexF[i];
      Geo::Point3<ScalarType> CurrB = BaryCoord[i];
      int IndexV0 = TriFaces[currF][0];
      int IndexV1 = TriFaces[currF][1];
      int IndexV2 = TriFaces[currF][2];

      Geo::Point3<ScalarType> OldP0 = OldVertPos[IndexV0];
      Geo::Point3<ScalarType> OldP1 = OldVertPos[IndexV1];
      Geo::Point3<ScalarType> OldP2 = OldVertPos[IndexV2];

      SmoothVertPos[i] = OldP0 * CurrB.X + OldP1 * CurrB.Y + OldP2 * CurrB.Z;
    }
    VertPos = SmoothVertPos;
    // //then restore border and corners
    // WriteOBJ("test_2.obj",SmoothVertPos, QuadFaces);

    // //then remap back

    // exit(0);
  }

  // static void
  // SmoothQuadPatches(std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                   const std::vector<std::vector<int>> &QuadFaces,
  //                   const std::vector<std::pair<int,int>> &QuadFeatures,
  //                   const std::vector<int> &QuadCorners,
  //                   const std::vector<std::vector<int>> &QuadPatchFaces,
  //                   const std::vector<int> &To_SmoothPatches,
  //                   const int smooth_steps = 1) {

  //   if (To_SmoothPatches.size() == 0)
  //     return;

  //   std::vector<int> CurrPatchFacesIdx;

  //   for (size_t i = 0; i < To_SmoothPatches.size(); i++) {
  //     int IndexP = To_SmoothPatches[i];
  //     assert(IndexP < QuadPatchFaces.size());
  //     CurrPatchFacesIdx.insert(CurrPatchFacesIdx.end(),
  //                              QuadPatchFaces[IndexP].begin(),
  //                              QuadPatchFaces[IndexP].end());
  //   }

  //   // get the mesh for that vertex
  //   std::vector<Geo::Point3<ScalarType>> SubMeshVert;
  //   std::vector<std::vector<int>> SubMeshesFaces;
  //   std::vector<int> OriginalVert,OriginalElem;

  //   CopySubMesh(VertPos, QuadFaces,
  //             CurrPatchFacesIdx,
  //             SubMeshVert,
  //             SubMeshesFaces,
  //             OriginalVert,OriginalElem);

  //   //map the features to the submeshes
  //   std::map<int,int> VertMap;
  //   for (size_t i=0;i<OriginalVert.size();i++)
  //     VertMap[OriginalVert[i]]=i;

  //   std::vector<std::pair<int,int>> SubMeshFeatures;
  //   std::vector<int> SubMeshCorners;

  //   for (size_t i=0;i<QuadFeatures.size();i++)
  //   {
  //     int IndexV0=QuadFeatures[i].first;
  //     int IndexV1=QuadFeatures[i].second;

  //     if (VertMap.count(IndexV0)==0)continue;
  //     if (VertMap.count(IndexV1)==0)continue;

  //     IndexV0=VertMap[IndexV0];
  //     IndexV1=VertMap[IndexV1];
  //     std::pair<int,int> KeyF(std::min(IndexV0,IndexV1),
  //                             std::max(IndexV0,IndexV1));
  //     SubMeshFeatures.push_back(KeyF);
  //   }

  //   for (size_t i=0;i<QuadCorners.size();i++)
  //   {
  //     int IndexV=QuadCorners[i];
  //     if (VertMap.count(IndexV)==0)continue;
  //     IndexV=VertMap[IndexV];
  //     SubMeshCorners.push_back(IndexV);
  //   }

  //   Geo::SmoothParam<ScalarType> SParamQuadSmooth;

  //   std::vector<Geo::Point3<ScalarType>> ReprojV=SubMeshVert;
  //   std::vector<std::vector<int>> ReprojF;
  //  //if (IsAllQuad(SubMeshesFaces))
  //   SplitQuadIntoTris<ScalarType>(ReprojV,SubMeshesFaces,ReprojF);

  //   SParamQuadSmooth.ReprojVert = &ReprojV;
  //   SParamQuadSmooth.ReprojFaces = &ReprojF;
  //   SParamQuadSmooth.Damp = 0.5;
  //   SParamQuadSmooth.post_project = false;
  //   SParamQuadSmooth.NumIte = smooth_steps;
  //   SParamQuadSmooth.SType = Geo::SMLaplacian;
  //   SParamQuadSmooth.Features = SubMeshFeatures;

  //   Geo::getBorderVert(SubMeshesFaces,SParamQuadSmooth.FixedVert);
  //   SParamQuadSmooth.FixedVert.insert(SParamQuadSmooth.FixedVert.end(),
  //                                     SubMeshCorners.begin(),
  //                                     SubMeshCorners.end());
  //   Geo::SmoothLaplacian<ScalarType>(SubMeshVert, SubMeshesFaces,
  //   SParamQuadSmooth);

  //   //copy back values
  //   for (size_t i=0;i<SubMeshVert.size();i++)
  //   {
  //     int OrigIndex=OriginalVert[i];
  //     VertPos[OrigIndex]=SubMeshVert[i];
  //   }
  // }

  static void
  SmoothQuadPatches(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &QuadFaces,
                    const std::vector<std::pair<int, int>> &QuadFeatures,
                    const std::vector<int> &QuadCorners,
                    const std::vector<std::vector<int>> &QuadPatchFaces,
                    const std::vector<int> &To_SmoothPatches,
                    const int smooth_steps = 1) {

    if (To_SmoothPatches.size() == 0)
      return;

    std::vector<int> CurrPatchFacesIdx;

    for (size_t i = 0; i < To_SmoothPatches.size(); i++) {
      int IndexP = To_SmoothPatches[i];
      assert(IndexP < QuadPatchFaces.size());
      CurrPatchFacesIdx.insert(CurrPatchFacesIdx.end(),
                               QuadPatchFaces[IndexP].begin(),
                               QuadPatchFaces[IndexP].end());
    }

    // get the mesh for that vertex
    std::vector<Geo::Point3<ScalarType>> SubMeshVert;
    std::vector<std::vector<int>> SubMeshesFaces;
    std::vector<int> OriginalVert, OriginalElem;

    CopySubMesh(VertPos, QuadFaces, CurrPatchFacesIdx, SubMeshVert,
                SubMeshesFaces, OriginalVert, OriginalElem);

    // map the features to the submeshes
    std::map<int, int> VertMap;
    for (size_t i = 0; i < OriginalVert.size(); i++)
      VertMap[OriginalVert[i]] = i;

    std::vector<std::pair<int, int>> SubMeshFeatures;
    std::vector<int> SubMeshCorners;

    for (size_t i = 0; i < QuadFeatures.size(); i++) {
      int IndexV0 = QuadFeatures[i].first;
      int IndexV1 = QuadFeatures[i].second;

      if (VertMap.count(IndexV0) == 0)
        continue;
      if (VertMap.count(IndexV1) == 0)
        continue;

      IndexV0 = VertMap[IndexV0];
      IndexV1 = VertMap[IndexV1];
      std::pair<int, int> KeyF(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      SubMeshFeatures.push_back(KeyF);
    }

    for (size_t i = 0; i < QuadCorners.size(); i++) {
      int IndexV = QuadCorners[i];
      if (VertMap.count(IndexV) == 0)
        continue;
      IndexV = VertMap[IndexV];
      SubMeshCorners.push_back(IndexV);
    }

    SmoothRectify(SubMeshVert, SubMeshesFaces, SubMeshFeatures, SubMeshCorners);

    //   Geo::SmoothParam<ScalarType> SParamQuadSmooth;

    //   std::vector<Geo::Point3<ScalarType>> ReprojV=SubMeshVert;
    //   std::vector<std::vector<int>> ReprojF;
    //  //if (IsAllQuad(SubMeshesFaces))
    //   SplitQuadIntoTris<ScalarType>(ReprojV,SubMeshesFaces,ReprojF);

    //   SParamQuadSmooth.ReprojVert = &ReprojV;
    //   SParamQuadSmooth.ReprojFaces = &ReprojF;
    //   SParamQuadSmooth.Damp = 0.5;
    //   SParamQuadSmooth.post_project = false;
    //   SParamQuadSmooth.NumIte = smooth_steps;
    //   SParamQuadSmooth.SType = Geo::SMLaplacian;
    //   SParamQuadSmooth.Features = SubMeshFeatures;

    //   Geo::getBorderVert(SubMeshesFaces,SParamQuadSmooth.FixedVert);
    //   SParamQuadSmooth.FixedVert.insert(SParamQuadSmooth.FixedVert.end(),
    //                                     SubMeshCorners.begin(),
    //                                     SubMeshCorners.end());
    //   Geo::SmoothLaplacian<ScalarType>(SubMeshVert, SubMeshesFaces,
    //   SParamQuadSmooth);

    // copy back values
    for (size_t i = 0; i < SubMeshVert.size(); i++) {
      int OrigIndex = OriginalVert[i];
      VertPos[OrigIndex] = SubMeshVert[i];
    }
  }

public:
  static void SmoothPatchesWithQuad(
      std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
      const std::vector<std::vector<int>> &QuadFaces,
      const std::vector<std::pair<int, int>> &QuadFeatures,
      const std::vector<int> &QuadCorners, const std::vector<int> &QuadPatch,
      const std::vector<Geo::Point3<ScalarType>> &ReprojVertPos,
      const std::vector<std::vector<int>> &ReprojTriFaces,
      const std::vector<std::pair<int, int>> &ReprojectTriFeatures,
      const int smooth_steps = 10) {
    // initialize the reprojection structures
    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    Geo::SegmentSpatialIndex<ScalarType> SegGrid;
    Geo::InitProjStructures<ScalarType>(ReprojVertPos, ReprojTriFaces,
                                        ReprojectTriFeatures, TriGrid,
                                        FeatureSegs, SegGrid);

    assert(IsAllQuad(QuadFaces));

    // find per patch face indexes
    std::vector<std::vector<int>> QuadPatchFaces;
    Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
        QuadPatch, QuadPatchFaces);

    // split the quads into tris
    std::vector<std::vector<int>> TriFaces;
    SplitQuadIntoTris<ScalarType>(QuadVertPos, QuadFaces, TriFaces);

    // reupdate the patch per triangle
    std::vector<int> TriPatch;
    for (size_t i = 0; i < QuadPatch.size(); i++) {
      TriPatch.push_back(QuadPatch[i]);
      TriPatch.push_back(QuadPatch[i]);
    }

    // find for each vertex the surrounding patches
    std::vector<std::vector<int>> VertPatches;
    FindVertPatches(QuadVertPos, TriFaces, TriPatch, VertPatches);

    // std::vector<std::pair<int, int>> EdgePatches;
    // FindEdgePatches(TriFaces, TriPatch, EdgePatches);

    // then does different steps
    for (size_t s = 0; s < smooth_steps; s++) {

      // smooth around all vertices
      for (size_t v = 0; v < VertPatches.size(); v++) {
        SmoothQuadPatches(QuadVertPos, QuadFaces, QuadFeatures, QuadCorners,
                          QuadPatchFaces, VertPatches[v]);
      }

      // ///backproject
      Geo::AttractProject<ScalarType>(ReprojVertPos, QuadVertPos, TriFaces);

      // reproject
      Geo::ProjectOn(QuadVertPos, QuadFaces, QuadFeatures, ReprojVertPos,
                     ReprojTriFaces, TriGrid, FeatureSegs, SegGrid);
    }
  }

  // static void SmoothPatchesUV(
  //     std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
  //     const std::vector<std::vector<int>> &QuadFaces,
  //     const std::vector<int> &QuadCorners, const std::vector<int> &QuadPatch,
  //     const std::vector<Geo::Point3<ScalarType>> &ReprojVertPos,
  //     const std::vector<std::vector<int>> &ReprojTriFaces,
  //     const std::vector<std::pair<int, int>> &ReprojectTriFeatures) {

  //   // find for each vertex the surrounding patches
  //   std::vector<std::vector<int>> QuadPatchFaces;
  //   Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
  //       QuadPatch, QuadPatchFaces);

  //   for (size_t i = 0; i < QuadPatchFaces.size(); i++) {
  //     std::vector<Geo::Point3<ScalarType>> SubMeshVert;
  //     std::vector<std::vector<int>> SubMeshesFaces;
  //     std::vector<int> OriginalVert, OriginalElem;

  //     CopySubMesh(QuadVertPos, QuadFaces, QuadPatchFaces[i], SubMeshVert,
  //                 SubMeshesFaces, OriginalVert, OriginalElem);

  //     std::vector<std::vector<int>> VF;
  //     std::vector<std::vector<int>> VFIndex;
  //     getVFAdjacency(SubMeshesFaces, VF, VFIndex);
  //     std::vector<int> SubMCornerIDX;
  //     for (size_t j = 0; j < VF.size(); j++) {
  //       if (VF[j].size() == 1)
  //         SubMCornerIDX.push_back(j);
  //     }

  //     // split into tris
  //     std::vector<std::vector<int>> SubMeshesTriFaces;
  //     SplitQuadIntoTris<ScalarType>(SubMeshVert, SubMeshesFaces,
  //                                   SubMeshesTriFaces);
      
  //     std::vector<std::pair<int,int>> BorderE;
  //     getBorderEdges(SubMeshesFaces, BorderE);
  //     std::vector<Geo::Point3<ScalarType>> SubMeshVertOld=SubMeshVert;
  //     RectifyFeatures<ScalarType>(SubMeshVert, SubMeshesFaces, 
  //                       BorderE,SubMCornerIDX);

  //     std::vector<Geo::Point2<ScalarType>> VertUV;
  //     // Geo::UVParam<ScalarType>::TutteEmbedding(SubMeshVert, SubMeshesTriFaces,
  //     //                                           SubMCornerIDX, VertUV,true);
  //     Geo::UVParam<ScalarType>::Arap(SubMeshVert, SubMeshesTriFaces,VertUV);
  //     // then smooth
  //     // SmoothLaplacianImplicit<ScalarType>(int &VertPos, const int &Faces,
  //     // const int &FixedVert)
  //     std::vector<int> BorderV;
  //     getBorderVert(SubMeshesTriFaces, BorderV);
  //     SmoothParam<ScalarType> SParam;
  //     SParam.FixedVert = SubMCornerIDX; //
  //     SParam.NumIte = 10;
  //     SParam.Damp = 0.5;
  //     SParam.Features=BorderE;
  //     //getBorderEdges(SubMeshesFaces, SParam.Features);
  //     std::vector<Geo::Point2<ScalarType>> SmoothVertUV = VertUV;
  //     Geo::SmoothLaplacian<ScalarType>(SmoothVertUV, SubMeshesFaces, SParam);
  //     //WriteOBJ("testUV0.obj",SmoothVertUV,SubMeshesFaces);
  //     // WriteOBJ("testUV1.obj",SmoothVertUV,SubMeshesFaces);
  //     // return;
  //     std::vector<int> faceIDx;
  //     std::vector<Geo::Point3<ScalarType>> bary;
  //     Geo::UVUtils<ScalarType>::MapBetween(SmoothVertUV, VertUV,
  //                                          SubMeshesTriFaces, faceIDx, bary);

      
  //     assert(faceIDx.size() == SmoothVertUV.size());
  //     assert(faceIDx.size() == bary.size());
  //     assert(faceIDx.size() == OriginalVert.size());

  //     // then reassign the values
  //     for (size_t j = 0; j < OriginalVert.size(); j++) {
  //       int OrigVertIDX = OriginalVert[j];
  //       assert(OrigVertIDX >= 0);
  //       assert(OrigVertIDX < QuadVertPos.size());

  //       int IndexF = faceIDx[j];
  //       assert(IndexF >= 0);
  //       assert(IndexF < SubMeshesTriFaces.size());
  //       int IndexV0 = SubMeshesTriFaces[IndexF][0];
  //       int IndexV1 = SubMeshesTriFaces[IndexF][1];
  //       int IndexV2 = SubMeshesTriFaces[IndexF][2];

  //       assert(IndexV0 >= 0);
  //       assert(IndexV0 < SubMeshVert.size());

  //       assert(IndexV1 >= 0);
  //       assert(IndexV1 < SubMeshVert.size());

  //       assert(IndexV2 >= 0);
  //       assert(IndexV2 < SubMeshVertOld.size());

  //       Geo::Point3<ScalarType> P0 = SubMeshVertOld[IndexV0];
  //       Geo::Point3<ScalarType> P1 = SubMeshVertOld[IndexV1];
  //       Geo::Point3<ScalarType> P2 = SubMeshVertOld[IndexV2];

  //       QuadVertPos[OrigVertIDX] = P0 * bary[j].X + P1 * bary[j].Y + P2 * bary[j].Z;
  //     }
  //     //return;
  //   }
  //    // ///backproject
  //   //   Geo::AttractProject<ScalarType>(ReprojVertPos, QuadVertPos, ReprojTriFaces);

  //   //   // reproject
  //   //    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  //   // std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
  //   // Geo::SegmentSpatialIndex<ScalarType> SegGrid;
  //   // Geo::InitProjStructures<ScalarType>(ReprojVertPos, ReprojTriFaces,
  //   //                                     ReprojectTriFeatures, TriGrid,
  //   //                                     FeatureSegs, SegGrid);

  //   //   Geo::ProjectOn(QuadVertPos, QuadFaces, QuadFeatures, ReprojVertPos,
  //   //                  ReprojTriFaces, TriGrid, FeatureSegs, SegGrid);

  // }



  static void SmoothPatchesUV(
      std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
      const std::vector<std::vector<int>> &QuadFaces,
      const std::vector<int> &QuadCorners, const std::vector<int> &QuadPatch,
      const std::vector<Geo::Point3<ScalarType>> &ReprojVertPos,
      const std::vector<std::vector<int>> &ReprojTriFaces,
      const std::vector<std::pair<int, int>> &ReprojectTriFeatures) {

    // find for each vertex the surrounding patches
    std::vector<std::vector<int>> QuadPatchFaces;
    Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
        QuadPatch, QuadPatchFaces);

    for (size_t i = 0; i < QuadPatchFaces.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> SubMeshVert;
      std::vector<std::vector<int>> SubMeshesFaces;
      std::vector<int> OriginalVert, OriginalElem;

      CopySubMesh(QuadVertPos, QuadFaces, QuadPatchFaces[i], SubMeshVert,
                  SubMeshesFaces, OriginalVert, OriginalElem);

      std::vector<std::vector<int>> VF;
      std::vector<std::vector<int>> VFIndex;
      getVFAdjacency(SubMeshesFaces, VF, VFIndex);
      std::vector<int> SubMCornerIDX;
      for (size_t j = 0; j < VF.size(); j++) {
        if (VF[j].size() == 1)
          SubMCornerIDX.push_back(j);
      }

      // split into tris
      std::vector<std::vector<int>> SubMeshesTriFaces;
      SplitQuadIntoTris<ScalarType>(SubMeshVert, SubMeshesFaces,
                                    SubMeshesTriFaces);
      
      std::vector<std::pair<int,int>> BorderE;
      getBorderEdges(SubMeshesFaces, BorderE);
      std::vector<Geo::Point3<ScalarType>> SubMeshVertOld=SubMeshVert;
      RectifyFeatures<ScalarType>(SubMeshVert, SubMeshesFaces, 
                        BorderE,SubMCornerIDX);

      std::vector<Geo::Point2<ScalarType>> VertUV;
      // Geo::UVParam<ScalarType>::TutteEmbedding(SubMeshVert, SubMeshesTriFaces,
      //                                           SubMCornerIDX, VertUV,true);
      Geo::UVParam<ScalarType>::Arap(SubMeshVert, SubMeshesTriFaces,VertUV);
      // then smooth
      // SmoothLaplacianImplicit<ScalarType>(int &VertPos, const int &Faces,
      // const int &FixedVert)
      std::vector<int> BorderV;
      getBorderVert(SubMeshesTriFaces, BorderV);
      SmoothParam<ScalarType> SParam;
      SParam.FixedVert = SubMCornerIDX; //
      SParam.NumIte = 10;
      SParam.Damp = 0.5;
      SParam.Features=BorderE;
      //getBorderEdges(SubMeshesFaces, SParam.Features);
      std::vector<Geo::Point2<ScalarType>> SmoothVertUV = VertUV;
      Geo::SmoothLaplacian<ScalarType>(SmoothVertUV, SubMeshesFaces, SParam);
      //WriteOBJ("testUV0.obj",SmoothVertUV,SubMeshesFaces);
      // WriteOBJ("testUV1.obj",SmoothVertUV,SubMeshesFaces);
      // return;
      std::vector<int> faceIDx;
      std::vector<Geo::Point3<ScalarType>> bary;
      Geo::UVUtils<ScalarType>::MapBetween(SmoothVertUV, VertUV,
                                           SubMeshesTriFaces, faceIDx, bary);

      
      assert(faceIDx.size() == SmoothVertUV.size());
      assert(faceIDx.size() == bary.size());
      assert(faceIDx.size() == OriginalVert.size());

      // then reassign the values
      for (size_t j = 0; j < OriginalVert.size(); j++) {
        int OrigVertIDX = OriginalVert[j];
        assert(OrigVertIDX >= 0);
        assert(OrigVertIDX < QuadVertPos.size());

        int IndexF = faceIDx[j];
        assert(IndexF >= 0);
        assert(IndexF < SubMeshesTriFaces.size());
        int IndexV0 = SubMeshesTriFaces[IndexF][0];
        int IndexV1 = SubMeshesTriFaces[IndexF][1];
        int IndexV2 = SubMeshesTriFaces[IndexF][2];

        assert(IndexV0 >= 0);
        assert(IndexV0 < SubMeshVert.size());

        assert(IndexV1 >= 0);
        assert(IndexV1 < SubMeshVert.size());

        assert(IndexV2 >= 0);
        assert(IndexV2 < SubMeshVertOld.size());

        Geo::Point3<ScalarType> P0 = SubMeshVertOld[IndexV0];
        Geo::Point3<ScalarType> P1 = SubMeshVertOld[IndexV1];
        Geo::Point3<ScalarType> P2 = SubMeshVertOld[IndexV2];

        QuadVertPos[OrigVertIDX] = P0 * bary[j].X + P1 * bary[j].Y + P2 * bary[j].Z;
      }
      //return;
    }
     // ///backproject
    //   Geo::AttractProject<ScalarType>(ReprojVertPos, QuadVertPos, ReprojTriFaces);

    //   // reproject
    //    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    // std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    // Geo::SegmentSpatialIndex<ScalarType> SegGrid;
    // Geo::InitProjStructures<ScalarType>(ReprojVertPos, ReprojTriFaces,
    //                                     ReprojectTriFeatures, TriGrid,
    //                                     FeatureSegs, SegGrid);

    //   Geo::ProjectOn(QuadVertPos, QuadFaces, QuadFeatures, ReprojVertPos,
    //                  ReprojTriFaces, TriGrid, FeatureSegs, SegGrid);

  }
};

} // namespace Geo

#endif
