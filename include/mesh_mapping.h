#ifndef MESH_MAPPING
#define MESH_MAPPING


#include "point3.h"
#include "triangle3.h"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <pos3.h>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
int FindVertCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &FromVert,
    const std::vector<Geo::Point3<ScalarType>> &ToVert,
    std::vector<int> &VerteMap)
{
  VerteMap.clear();
  VerteMap.resize(FromVert.size(), -1);
  std::map<Geo::Point3<ScalarType>, int> ToVertMap;
  for (size_t i=0;i<ToVert.size();i++)
    ToVertMap[ToVert[i]]=i;
  
  int Mapped = 0;
  for (size_t i=0;i<FromVert.size();i++)
  {
    if (ToVertMap.count(FromVert[i])==0)continue;
    int IndexV=ToVertMap[FromVert[i]];
    assert(IndexV>=0);
    assert(IndexV<ToVert.size());
    VerteMap[i]=IndexV;
    Mapped++;
  }
  return Mapped;
}

template <class ScalarType>
int FindFacesCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &FromVert,
    const std::vector<std::vector<int>> &FromFaces,
    const std::vector<Geo::Point3<ScalarType>> &ToVert,
    const std::vector<std::vector<int>> &ToFaces, std::vector<int> &FaceMap) {
  FaceMap.clear();
  FaceMap.resize(FromFaces.size(), -1);

  std::map<Geo::Traingle3<ScalarType>, int> ToFacesMap;
  for (size_t i = 0; i < ToFaces.size(); i++) {
    assert(ToFaces[i].size() == 3);
    std::vector<Geo::Point3<ScalarType>> TrPos(3, Geo::Point3<ScalarType>(0, 0, 0));
    Geo::Point3<ScalarType> P0 = ToVert[ToFaces[i][0]];
    Geo::Point3<ScalarType> P1 = ToVert[ToFaces[i][1]];
    Geo::Point3<ScalarType> P2 = ToVert[ToFaces[i][2]];
    TrPos.push_back(P0);
    TrPos.push_back(P1);
    TrPos.push_back(P2);
    std::sort(TrPos.begin(), TrPos.end());
    Geo::Traingle3<ScalarType> T(P0, P1, P2);
    ToFacesMap.insert(std::pair<Geo::Traingle3<ScalarType>, int>(T, i));
  }

  int Mapped = 0;
  for (size_t i = 0; i < FromFaces.size(); i++) {
    assert(FromFaces[i].size() == 3);
    std::vector<Geo::Point3<ScalarType>> TrPos(
        3, Geo::Point3<ScalarType>(0, 0, 0));
    Geo::Point3<ScalarType> P0 = FromVert[FromFaces[i][0]];
    Geo::Point3<ScalarType> P1 = FromVert[FromFaces[i][1]];
    Geo::Point3<ScalarType> P2 = FromVert[FromFaces[i][2]];
    TrPos.push_back(P0);
    TrPos.push_back(P1);
    TrPos.push_back(P2);
    std::sort(TrPos.begin(), TrPos.end());
    Geo::Traingle3<ScalarType> T(P0, P1, P2);
    if (ToFacesMap.count(T) > 0) {
      int Index = ToFacesMap[T];
      assert(Index >= 0);
      assert(Index < ToFaces.size());
      FaceMap[i] = Index;
      Mapped++;
    }
  }
  return Mapped;
}


template <class ScalarType>
void FindEdgeFacesCorrespondences(
    const std::vector<Geo::Point3<ScalarType> > &Vert,
    const std::vector<std::vector<int> > &Faces,
    std::map<std::pair<int,int>,std::vector<std::pair<int,int> > > &EdgeMap) 
{
  EdgeMap.clear();
  for (size_t i=0;i<Faces.size();i++)
  for (size_t j=0;j<Faces[i].size();j++)
  {
    int IndexV0=Faces[i][j];
    int IndexV1=Faces[i][(j+1)%Faces[i].size()];
    std::pair<int,int> Key(std::min(IndexV0,IndexV1),
                           std::max(IndexV0,IndexV1));
    std::pair<int,int> FaceEdge(i,j);
    EdgeMap[Key].push_back(FaceEdge);
  }
}

template <class ScalarType>
void FindEdgeFacesCorrespondences(
    const std::vector<Geo::Point3<ScalarType> > &Vert,
    const std::vector<std::vector<int> > &Faces,
    const std::vector<std::pair<int,int> > &Edges,
    std::vector<std::vector<std::pair<int,int> > > &EdgeFaces) 
{
  EdgeFaces.clear();
  if (Edges.size()==0)return;
  std::map<std::pair<int,int>,std::vector<std::pair<int,int> > > EdgeMap;
  FindEdgeFacesCorrespondences<ScalarType>(Vert,Faces,EdgeMap);
  for (size_t i=0;i<Edges.size();i++)
  {
    assert(EdgeMap.count(Edges[i])>0);
    EdgeFaces.push_back(EdgeMap[Edges[i]]);
  }
}

};
#endif
