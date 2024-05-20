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
