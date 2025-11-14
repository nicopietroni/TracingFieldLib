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

#ifndef MESH_SUBDIVIDE
#define MESH_SUBDIVIDE

#include "Space_Query/triangle_grid_3D.h"
#include "mesh_refine.h"
#include "point3.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <subdivide.h>
#include <sys/_types/_size_t.h>
#include <tangent_space_smooth.h>
#include <utility>
#include <vector>

namespace Geo{
template <class ScalarType> class SubdivideMesh {
public:
  static void MidPointStep(std::vector<Geo::Point3<ScalarType>> &Pos,
                           std::vector<std::vector<int>> &Faces,
                           std::vector<int> &OldF) {
    std::vector<std::pair<int, int> > MeshE;
    GetAllMeshEdges<ScalarType>(Faces, MeshE);

    std::vector<ScalarType> Alpha(MeshE.size(), 0.5);

    RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, MeshE, Alpha, OldF);
  }

  static void
  MidPointReproject(const std::vector<Geo::Point3<ScalarType> > &TargetPos,
                    const std::vector<std::vector<int> > &TargetFaces,
                    std::vector<Geo::Point3<ScalarType> > &Pos,
                    std::vector<std::vector<int> > &Faces,
                    std::vector<int> &OldF, 
                    int subd_steps) 
  {
    Geo::TriangleSpatialIndex<ScalarType> TriIdx;
    TriIdx.Init(TargetPos, TargetFaces);
    OldF.clear();
    for (size_t i = 0; i < subd_steps; i++) {
      std::vector<int> OldFStep;
      MidPointStep(Pos, Faces, OldFStep);
      
      // reproject
      ScalarType maxDist = TriIdx.bbox.Diag();

      for (size_t j = 0; j < Pos.size(); j++) {
        Geo::Point3<ScalarType> CurrP = Pos[j];
        int closestF;
        Geo::Point3<ScalarType> closestPt;
        TriIdx.GridClosest(TargetPos, TargetFaces, CurrP, maxDist, closestF,
                           closestPt);
        Pos[j] = closestPt;
      }
      //update old faces
      if (OldF.size() == 0) {
        OldF=OldFStep;
      } else {
        for (size_t j = 0; j < OldFStep.size(); j++) {
          OldFStep[j]=OldF[OldFStep[j]];
        }
        OldF=OldFStep;
      }
    }
  }
};
}

#endif
