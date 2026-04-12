/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

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
#include "quad2.h"
#include "quadrangulationLUT.h"
#include "segment2.h"
#include "triangle3.h"
#include <box3.h>
#include <cassert>
#include <cmath>
#include <connectivity.h>
#include <cstddef>
#include <iostream>
#include <mesh_mapping.h>
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType> class Quadrangulation {
public:
  static bool SampleRegularSubPatch(
      const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
      const Geo::TriangleSpatialIndex2D<ScalarType> &PatchUVGrid,
      const std::vector<Geo::Point3<ScalarType>> &PatchVert3D,
      const std::vector<std::vector<int>> &Patch_Faces, const int SizeX,
      const int SizeY, const Geo::Quad2<ScalarType> &DomainQuadUV,
      std::vector<Geo::Point2<ScalarType>> &QuadVertUV,
      std::vector<Geo::Point3<ScalarType>> &QuadVert3D,
      std::vector<std::vector<int>> &QuadFaces)

  {
    std::vector<std::vector<Geo::Point2<ScalarType>>> SampledPosUV;
    std::vector<std::vector<Geo::Point3<ScalarType>>> Sampled3D;
    std::vector<std::vector<int>> IndexSampled;

    QuadVertUV.clear();
    QuadVert3D.clear();
    QuadFaces.clear();
    if (SizeX == 0 || SizeY == 0)
      return false;

    // then update the indexes
    DomainQuadUV.SampleBilinear(SizeX, SizeY, SampledPosUV);
    IndexSampled.resize(SampledPosUV.size());
    // then set the index of the sampled UVs
    for (size_t x = 0; x < SampledPosUV.size(); x++) {
      IndexSampled[x].resize(SampledPosUV[x].size(), -1);
      for (size_t y = 0; y < SampledPosUV[x].size(); y++) {
        QuadVertUV.push_back(SampledPosUV[x][y]);
        IndexSampled[x][y] = QuadVertUV.size() - 1;
      }
    }

    // then retrieve the 3D position for each sampled UV
    Sampled3D.resize(SampledPosUV.size());
    for (size_t x = 0; x < SampledPosUV.size(); x++) {
      Sampled3D[x].resize(SampledPosUV[x].size());
      for (size_t y = 0; y < SampledPosUV[x].size(); y++) {
        Geo::Point2<ScalarType> P = SampledPosUV[x][y];
        Geo::Point3<ScalarType> BaryCoords;
        int IndexF = -1;
        bool found = PatchUVGrid.FirstTrianglesOverPoint(
            PatchVertUV, Patch_Faces, P, IndexF, BaryCoords);
        if (!found)
          return false;

        assert(IndexF != -1);

        int IndexV0 = Patch_Faces[IndexF][0];
        int IndexV1 = Patch_Faces[IndexF][1];
        int IndexV2 = Patch_Faces[IndexF][2];

        Geo::Point3<ScalarType> P0 = PatchVert3D[IndexV0];
        Geo::Point3<ScalarType> P1 = PatchVert3D[IndexV1];
        Geo::Point3<ScalarType> P2 = PatchVert3D[IndexV2];

        Geo::Point3<ScalarType> InterpP =
            P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
        Sampled3D[x][y] = InterpP;
        QuadVert3D.push_back(InterpP);
      }
    }

    // then retrieve the 3D point
    QuadFaces.clear();
    for (size_t i = 0; i < IndexSampled.size() - 1; i++) {
      for (size_t j = 0; j < IndexSampled[i].size() - 1; j++) {
        int IndexV0 = IndexSampled[i][j];
        int IndexV1 = IndexSampled[i + 1][j];
        int IndexV2 = IndexSampled[i + 1][j + 1];
        int IndexV3 = IndexSampled[i][j + 1];
        QuadFaces.resize(QuadFaces.size() + 1);
        QuadFaces.back().push_back(IndexV0);
        QuadFaces.back().push_back(IndexV1);
        QuadFaces.back().push_back(IndexV2);
        QuadFaces.back().push_back(IndexV3);
      }
    }
    
    return true;
  }

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

  static bool OutofLUTBounds(
      const QuadrangulationLUT<ScalarType> &QLut,
      const std::vector<std::vector<Geo::Point2<ScalarType>>> &UVSideSubd) {
    assert(UVSideSubd.size() == 4);
    int size0 = UVSideSubd[0].size();
    int size1 = UVSideSubd[1].size();
    int size2 = UVSideSubd[2].size();
    int size3 = UVSideSubd[3].size();
    assert(size0 >= 0);
    assert(size1 >= 0);
    assert(size2 >= 0);
    assert(size3 >= 0);
    return ((size0 > QLut.MaxSubdVal) || (size1 > QLut.MaxSubdVal) ||
            (size2 > QLut.MaxSubdVal) || (size3 > QLut.MaxSubdVal));
  }

  // static bool QuadrangulateExtendLUT(const QuadrangulationLUT<ScalarType>
  // &QLut,
  //                                   const
  //                                   std::vector<std::vector<Geo::Point2<ScalarType>>>
  //                                   &UVSideSubd)
  // {
  //   int Diff0_U=UVSideSubd[0].size()-QLut.MaxSubdVal;
  //   int Diff1_U=UVSideSubd[2].size()-QLut.MaxSubdVal;
  // }

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
    bool found = false;

    if (!OutofLUTBounds(QLut, UVSideSubd))
      found = QLut.GetQuadrangulation(UVSideSubd, QuadVertUV, QuadFaces);
    else {
      std::cout << "MaxSizd:" << QLut.MaxSubdVal << std::endl;
      std::cout << "Size 0:" << UVSideSubd[0].size() << std::endl;
      std::cout << "Size 1:" << UVSideSubd[1].size() << std::endl;
      std::cout << "Size 2:" << UVSideSubd[2].size() << std::endl;
      std::cout << "Size 3:" << UVSideSubd[3].size() << std::endl;
      return false;
    }
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

    // the first one will be the one along the border, sorted accordignly
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
        std::cout << "Error Quadrangulation 3D Mapping" << std::endl;
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
                                             PatchCorners, VertUV, true);
    Geo::TriangleSpatialIndex<ScalarType> Grid;
    Grid.Init(PatchVert, Patch_Faces);
    Geo::Box3<ScalarType> BB;
    BB.Init(PatchVert);

    // get the UV of side subdivision
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
  // return the amount of non quadrangulated patches
  static int QuadrangulateQuadPatchsWithSubs(
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

    int problematic_patches = 0;
    for (size_t i = 0; i < PatchFacesIdx.size(); i++) {
      // for (size_t i = 0; i <1; i++) {
      assert(PatchCornerIdx[i].size() == SideSamples3D[i].size());

      // only if the patch corner side is 4
      if (PatchCornerIdx[i].size() != 4) {
        std::cout << "Warning Quadrangulation of Patch with Corner Num !=4"
                  << std::endl;
        problematic_patches++;
        continue;
      }
      // check eveness
      int TotSub = 0;
      for (size_t j = 0; j < SideSamples3D[i].size(); j++)
        TotSub += SideSamples3D[i][j].size();
      if ((TotSub % 2) != 0) {
        problematic_patches++;
        continue;
      }
      // copy the submesh
      std::vector<Geo::Point3<ScalarType>> CurrPatchVert;
      std::vector<std::vector<int>> CurrPatchFaces;
      std::vector<int> OriginalVert, OriginalElem;
      CopySubMesh<ScalarType>(Vert_3D, Faces, PatchFacesIdx[i], CurrPatchVert,
                              CurrPatchFaces, OriginalVert, OriginalElem);

      // then map the corners on the subpatch
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

      // then quadrangulate returning also the side samples
      bool has_quadrangulated = QuadrangulateQuadPatchWithSubs(
          CurrPatchVert, CurrPatchFaces, CurrPatchCorners, SideSamples3D[i],
          QLut, Curr_Patch_QuadVert, Curr_Patch_QuadFaces);
      if (!has_quadrangulated) {
        problematic_patches++;
        std::cout << "WARNING: Problem Quadrangulating patch " << i
                  << std::endl;
      }
      if (save_debug) {
        std::string out_file =
            std::string("quad_mesh_") + std::to_string(i) + std::string(".obj");
        WriteOBJ(out_file.c_str(), Curr_Patch_QuadVert, Curr_Patch_QuadFaces);
      }

      AppendMesh<ScalarType>(QuadVert, QuadFaces, Curr_Patch_QuadVert,
                             Curr_Patch_QuadFaces);

      QuadPatch.resize(QuadFaces.size(), i);
    }
    return problematic_patches;
  }
};

} // namespace Geo
#endif
