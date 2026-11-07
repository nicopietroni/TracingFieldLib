// //Non-Commercial License Agreement
// //Copyright (c) 2025, Nico Pietroni
// //
// //1. Permission is hereby granted to use, copy, modify, and distribute
// //   this software for academic research, teaching, and personal
// //   non-commercial purposes, provided that this copyright notice is
// //   included in all copies or substantial portions of the software.
// //
// //2. Commercial use of this software, including but not limited to use in
// //   commercial products, services, tools, or consulting, is strictly
// //   prohibited without prior written permission from the copyright holder.
// //
// //3. This software is provided "as is", without warranty of any kind,
// //   express or implied, including but not limited to the warranties of
// //   merchantability, fitness for a particular purpose, and
// //   noninfringement. In no event shall the authors or copyright holders
// //   be liable for any claim, damages, or other liability, arising from,
// //   out of, or in connection with the software or the use or other
// //   dealings in the software.
// //
// //4. For commercial licensing, please contact:
// //   Nico Pietroni
// //   Email: nico.pietroni@gmail.com

// sym_align: finds the reflective symmetry plane of a mesh and, if the mesh
// is symmetric, moves it (rotation + translation) so that the plane is the
// YZ plane. Output: <name>_sym.obj next to the input, one row in
// symmetry.csv in the folder of the executable.
// Exit code: 0 symmetric (aligned mesh saved), 1 not symmetric, 2 error.

#include <IO/file_path_functions.h>
#include <IO/io_obj.h>
#include <IO/io_stl.h>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mesh_functions.h>
#include <point3.h>
#include <string>
#include <symmetry_detection.h>
#include <vector>

typedef double ScalarType;

bool IsSTLFile(const std::string &path) {
  if (path.size() < 4)
    return false;
  std::string ext = path.substr(path.size() - 4);
  for (size_t i = 0; i < ext.size(); i++)
    ext[i] = (char)std::tolower((unsigned char)ext[i]);
  return (ext == std::string(".stl"));
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "usage: sym_align <mesh.obj|mesh.stl> [-tol 0.005] "
                 "[-checktol 0.002] [-checkedge 0] [-normangle 40] [-maxout 0.001] [-dirs 600] "
                 "[-samples 2000] [-verbose]"
              << std::endl;
    return 2;
  }
  std::string PathMesh = std::string(argv[1]);
  Geo::SymmetryParam<ScalarType> Param;
  for (int i = 2; i < argc; i++) {
    std::string arg(argv[i]);
    bool has_val = (i < argc - 1);
    if ((arg == "-tol") && has_val)
      Param.tolerance = std::atof(argv[++i]);
    else if ((arg == "-checktol") && has_val)
      Param.check_tolerance = std::atof(argv[++i]);
    else if ((arg == "-checkedge") && has_val)
      Param.check_edge_ratio = std::atof(argv[++i]);
    else if ((arg == "-normangle") && has_val)
      Param.max_normal_angle = std::atof(argv[++i]);
    else if ((arg == "-maxout") && has_val)
      Param.max_outliers = std::atof(argv[++i]);
    else if ((arg == "-dirs") && has_val)
      Param.num_directions = std::atoi(argv[++i]);
    else if ((arg == "-samples") && has_val)
      Param.fine_samples = std::atoi(argv[++i]);
    else if (arg == "-verbose")
      Param.verbose = true;
    else
      std::cout << "Warning: unknown option " << arg << std::endl;
  }

  // LOAD
  std::vector<Geo::Point3<ScalarType>> VertPos;
  std::vector<std::vector<int>> Faces;
  bool loaded = IsSTLFile(PathMesh) ? LoadSTL(PathMesh, VertPos, Faces, true)
                                    : LoadOBJ(PathMesh, VertPos, Faces, true);
  if ((!loaded) || Faces.empty()) {
    std::cout << "Error: cannot load the mesh " << PathMesh << std::endl;
    return 2;
  }
  if (!IsAllTris(Faces)) {
    std::vector<Geo::Point3<ScalarType>> VertPosTris;
    std::vector<std::vector<int>> FacesTris;
    SplitIntoTris(VertPos, Faces, VertPosTris, FacesTris);
    VertPos = VertPosTris;
    Faces = FacesTris;
  }

  // SEARCH
  Geo::SymmetryPlane<ScalarType> Sym;
  bool symmetric = Geo::SymmetryDetection<ScalarType>::FindSymmetryPlane(
      VertPos, Faces, Param, Sym);
  const Geo::Point3<ScalarType> &N = Sym.Plane.Direction();
  std::cout << "SYMMETRIC " << (symmetric ? 1 : 0) << " outliers "
            << Sym.Outliers << " (by normal " << Sym.NormalOutliers
            << ") max_angle " << Sym.MaxAngle << " max_err " << Sym.MaxErr << " normal " << N.X
            << " " << N.Y << " " << N.Z << " offset " << Sym.Plane.Offset()
            << std::endl;

  // ALIGN AND SAVE
  if (symmetric) {
    Geo::SymmetryDetection<ScalarType>::AlignPlaneToYZ(VertPos, Sym.Plane);
    std::string Folder = GetDirectory(PathMesh);
    if (Folder.empty())
      Folder = std::string(".");
    std::string OutFile = Folder + std::string("/") +
                          GetFilanameNoExtension(PathMesh) +
                          std::string("_sym.obj");
    std::cout << "Saving on:" << OutFile << std::endl;
    WriteOBJ(OutFile, VertPos, Faces);
  }

  // CSV (folder of the executable):
  // mesh,symmetric,outliers,max_err,sample_score,nx,ny,nz,offset
  std::string ExeFolder =
      std::filesystem::absolute(std::filesystem::path(argv[0]))
          .parent_path()
          .string();
  std::string CsvFile = ExeFolder + std::string("/symmetry.csv");
  bool new_file = !std::filesystem::exists(CsvFile);
  std::ofstream Csv(CsvFile.c_str(), std::ios::app);
  if (new_file)
    Csv << "mesh,symmetric,outliers,max_err,sample_score,nx,ny,nz,offset"
        << std::endl;
  Csv << GetFilaname(PathMesh) << "," << (symmetric ? 1 : 0) << ","
      << Sym.Outliers << "," << Sym.MaxErr << "," << Sym.SampleScore << ","
      << N.X << "," << N.Y << "," << N.Z << "," << Sym.Plane.Offset()
      << std::endl;

  return symmetric ? 0 : 1;
}
