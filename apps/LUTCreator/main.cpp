#include <IO/file_path_functions.h>
#include <IO/io_obj.h>
#include <cassert>
#include <chrono>
#include <connectivity.h>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <mesh_functions.h>
#include <point3.h>
#include <quadrangulationLUT.h>
#include <sstream>
#include <string>
#include <vector>

typedef double ScalarType;

int firstIndex(const std::vector<int> &quadruple) {
  std::vector<int> best = quadruple;
  std::vector<int> rotated = quadruple;

  int FirstIdx = 0;
  for (int i = 1; i < 4; ++i) {
    rotate(rotated.begin(), rotated.begin() + 1, rotated.end());
    if (rotated < best) {
      best = rotated;
      FirstIdx = i;
    }
  }
  return FirstIdx;
}

void SortConnectivity(std::vector<Geo::Point3<ScalarType>> &VertPos,
                      std::vector<std::vector<int>> &Connectivity,
                      std::vector<int> &CodeE) {
  ScalarType TolAngle = 10;

  std::vector<std::vector<int>> BorderSeq;
  getBorderSequncesVert(Connectivity, BorderSeq);
  assert(BorderSeq.size() == 1);
  int sizeSeq = BorderSeq[0].size();
  std::reverse(BorderSeq[0].begin(), BorderSeq[0].end());
  std::vector<int> IndexCorners;
  for (size_t i = 0; i < BorderSeq[0].size(); i++) {
    int SeqIdx0 = (i + sizeSeq - 1) % sizeSeq;
    int SeqIdx1 = i;
    int SeqIdx2 = (i + 1) % sizeSeq;

    int IndexV0 = BorderSeq[0][SeqIdx0];
    int IndexV1 = BorderSeq[0][SeqIdx1];
    int IndexV2 = BorderSeq[0][SeqIdx2];

    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
    Geo::Point3<ScalarType> P2 = VertPos[IndexV2];
    Geo::Point3<ScalarType> Dir0 = P0 - P1;
    Geo::Point3<ScalarType> Dir1 = P2 - P1;
    Dir0.Normalize();
    Dir1.Normalize();
    ScalarType AngleDir = fabs(Geo::AngleDeg(Dir0, Dir1));
    if (fabs(AngleDir - 180) > TolAngle)
      IndexCorners.push_back(i);
  }
  // std::cout<<"Borders:"<<BorderSeq.size()<<std::endl;
  // std::cout<<"Border Size:"<<BorderSeq[0].size()<<std::endl;
  // std::cout<<"Num Corners:"<<IndexCorners.size()<<std::endl;

  // get the starting point (one with bigger interval)
  CodeE.clear();
  for (size_t i = 0; i < IndexCorners.size(); i++) {
    int Index0 = IndexCorners[i];
    int Index1 = IndexCorners[(i + 1) % IndexCorners.size()];
    int Dist = (Index1 - Index0 + sizeSeq) % sizeSeq;
    CodeE.push_back(Dist);
  }

  int startingIdxE = firstIndex(CodeE);
  int startingIdx = IndexCorners[startingIdxE];

  std::rotate(CodeE.begin(), CodeE.begin() + startingIdxE, CodeE.end());

  std::rotate(BorderSeq[0].begin(), BorderSeq[0].begin() + startingIdx,
              BorderSeq[0].end());

  // then resort the mesh
  std::map<int, int> VertMap;
  std::vector<Geo::Point3<ScalarType>> NewVertPos;

  // start from the border
  for (size_t i = 0; i < BorderSeq[0].size(); i++) {
    int IndexOld = BorderSeq[0][i];
    int IndexNew = i;
    NewVertPos.push_back(VertPos[IndexOld]);
    VertMap[IndexOld] = IndexNew;
  }

  // add the rest
  for (size_t i = 0; i < VertPos.size(); i++) {
    // already inserted
    if (VertMap.count(i) > 0)
      continue;
    NewVertPos.push_back(VertPos[i]);
    VertMap[i] = NewVertPos.size() - 1;
  }
  VertPos = NewVertPos;

  // reupdate indexes of faces
  for (size_t i = 0; i < Connectivity.size(); i++)
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int IndexVOld = Connectivity[i][j];
      assert(VertMap.count(IndexVOld) > 0);
      int IndexVNew = VertMap[IndexVOld];
      Connectivity[i][j] = IndexVNew;
    }
}

bool LoadMesh(const std::string &path,
              std::vector<typename Geo::Point3<ScalarType>> &VertPos,
              std::vector<std::vector<int>> &Connectivity) {
  bool loadedMesh = LoadOBJ(path, VertPos, Connectivity, true);
  if (loadedMesh) {
    // std::cout << "Loaded " << Connectivity.size() << " faces " <<
    // VertPos.size() << " vertices" << std::endl;
    return true;
  } else {
    std::cout << "Mesh not loaded" << std::endl;
    return false;
  }
}

bool WriteLUT(const std::string &path,
              const std::vector<std::vector<std::vector<int>>> &Connectivity,
              const std::vector<std::vector<int>> &CodeE) {
  std::ofstream fout(path.c_str());
  if (!fout)
    return false;
  assert(Connectivity.size() == CodeE.size());

  fout << Connectivity.size() << "\n";

  for (size_t i = 0; i < Connectivity.size(); i++) {
    for (size_t j = 0; j < CodeE[i].size() - 1; j++) {
      fout << CodeE[i][j] << " ";
    }
    fout << CodeE[i][CodeE[i].size() - 1] << "\n";
    fout << Connectivity[i].size() << "\n";
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      // should be a quad
      assert(Connectivity[i][j].size() == 4);
      int IndexV0 = Connectivity[i][j][0];
      int IndexV1 = Connectivity[i][j][1];
      int IndexV2 = Connectivity[i][j][2];
      int IndexV3 = Connectivity[i][j][3];
      fout << "(" << IndexV0 << "," << IndexV1 << "," << IndexV2 << ","
           << IndexV3 << ")";
      fout << std::endl;
    } 
  }
  // for (size_t i=0;i<Connectivity.size();i++)
  // {
  //     fout<<"f";
  //     for (size_t j=0;j<Connectivity[i].size();j++)
  //         fout<<" "<<Connectivity[i][j]+1;
  //     fout<<"\n";
  // }

  fout.close();
  return true;
}

std::string commandLine;

int main(int argc, char *argv[]) {

  // std::string OutLUTFile("Side4.lut");
  // if (argc < 2) {
  //   std::cout << "You should pass at least a directory" << std::endl;
  //   exit(0);
  // }

  // namespace fs = std::filesystem;

  // std::vector<std::string> obj_files;
  // std::string directory_path = argv[1];
  // try {
  //   for (const auto &entry : fs::directory_iterator(directory_path)) {
  //     if (entry.is_regular_file() && entry.path().extension() == ".obj") {
  //       obj_files.push_back(entry.path().string());
  //     }
  //   }
  // } catch (const fs::filesystem_error &e) {
  //   std::cerr << "Filesystem error: " << e.what() << std::endl;
  // }
  // std::cout << "Processing " << obj_files.size() << " OBJ Meshes" << std::endl;

  // std::vector<std::vector<int>> Codes;
  // std::vector<std::vector<std::vector<int>>> Connectivities;
  // for (size_t i = 0; i < obj_files.size(); i++) {
  //   std::vector<typename Geo::Point3<ScalarType>> VertPos;
  //   std::vector<std::vector<int>> Connectivity;
  //   LoadMesh(obj_files[i], VertPos, Connectivity);
  //   std::vector<int> CodeSub;
  //   SortConnectivity(VertPos, Connectivity, CodeSub);
  //   Codes.push_back(CodeSub);
  //   // std::cout<<std::endl;
  //   // for (size_t j=0;j<CodeSub.size();j++)
  //   //   std::cout<<CodeSub[j]<<" ";
  //   // std::cout<<std::endl;

  //   Connectivities.push_back(Connectivity);
  // }

  // WriteLUT("Lut.txt", Connectivities, Codes);
  Geo::QuadrangulationLUT<ScalarType> Q_Lut;
  Q_Lut.Initialize();
  
  Geo::Point2<ScalarType> CornTest[4];
  ScalarType SideL = 1;
  CornTest[0] = Geo::Point2<ScalarType>(0, 0);
  CornTest[1] = Geo::Point2<ScalarType>(SideL, 0);
  CornTest[2] = Geo::Point2<ScalarType>(SideL, SideL);
  CornTest[3] = Geo::Point2<ScalarType>(0, SideL);
  bool found = false;
  int subS[4];
  srand(time(0));

  int MaxDiv=20;
  do {
    for (size_t i = 0; i < 4; i++) {
      subS[i] = rand() % (MaxDiv+1);
      subS[i] = std::max(subS[i], 1);
    }
    // check if ok
    int Sum = 0;
    for (size_t i = 0; i < 4; i++)
      Sum += subS[i];

    // off subdivision
    if ((Sum % 2) == 1)
      continue;

    // equal opp side
    if ((subS[0] == subS[2]) && (subS[1] == subS[3]))
      continue;

    found = true;
  } while (!found);
  std::cout << "Testing Code: " << subS[0] << "," << subS[1] << "," << subS[2]
            << "," << subS[3] << std::endl;

  std::vector<std::vector<Geo::Point2<ScalarType>>> SidePos;
  SidePos.resize(4);
  for (size_t i = 0; i < 4; i++) {
    Geo::Point2<ScalarType> Corn0 = CornTest[i];
    Geo::Point2<ScalarType> Corn1 = CornTest[(i + 1) % 4];

    SidePos[i].push_back(Corn0);
    int ToAddPoints = subS[i] - 1;

    if (ToAddPoints > 0) {
      ScalarType InterpW = ((ScalarType)1) / ((ScalarType)(ToAddPoints + 1));
      //std::cout<<"Int W"<<InterpW<<std::endl;
      for (int j = 0; j < ToAddPoints; j++) {
        ScalarType CurrW = 1 - (InterpW * (j+1));
        //std::cout<<"Curr W"<<CurrW<<std::endl;
        Geo::Point2<ScalarType> IntPoint = Corn0 * CurrW + Corn1 * (1 - CurrW);
        SidePos[i].push_back(IntPoint);
      }
    }
    SidePos[i].push_back(Corn1);
  }
  std::vector<Geo::Point2<ScalarType>> VertPos;
  std::vector<std::vector<int>> Connectivity;
  
  // Test
  
  Q_Lut.GetQuadrangulation(SidePos, VertPos, Connectivity);

  WriteOBJ("testOUT.obj",VertPos, Connectivity);

  return 0;
}