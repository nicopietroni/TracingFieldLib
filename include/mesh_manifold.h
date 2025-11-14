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

#ifndef MESH_MANIFOLD
#define MESH_MANIFOLD

// #include <string>
// #include <fstream>
// #include <iostream>
// #include <vector>
// #include <map>
// #include <set>
// #include "point3.h"
// #include "point2.h"
// #include "segment2.h"
// #include "plane3.h"
// #include "connectivity.h"
// #include <pos3.h>
// #include "Space_Query/point_grid_3D.h"
// #include "smooth_mesh.h"

#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>

//#include "eigen_interface.h"
// #include "Space_Query/triangle_grid_2D.h"

template <class ScalarType>
static bool EdgeManifold(std::vector<Geo::Point3<ScalarType> > &Pos3D,
                        std::vector<std::vector<int> > &Faces)
{
    Eigen::MatrixXi F;
    typename EigenInterface<ScalarType>::MatrixXm Vf;

    EigenInterface<ScalarType>::GetTriMeshData(Pos3D, Faces, F, Vf);
    Eigen::MatrixXd V = Vf.template cast<double>();

     bool is_manifold = igl::is_edge_manifold(F);
     return is_manifold;
}

template <class ScalarType>
static bool VertManifold(std::vector<Geo::Point3<ScalarType> > &Pos3D,
                        std::vector<std::vector<int> > &Faces)
{
    Eigen::MatrixXi F;
    typename EigenInterface<ScalarType>::MatrixXm Vf;

    EigenInterface<ScalarType>::GetTriMeshData(Pos3D, Faces, F, Vf);
    Eigen::MatrixXd V = Vf.template cast<double>();

    bool is_manifold = igl::is_vertex_manifold(F);
    return is_manifold;
}


#endif
