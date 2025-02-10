#ifndef PATCH_SIDE_DISCRETIZATION
#define PATCH_SIDE_DISCRETIZATION

#include "mesh_functions.h"
#include "point3.h"
#include "pos3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class SideDiscretization {
  struct PatchSide {
    std::vector<std::vector<Geo::Point3<ScalarType> > > SubSides;
    // std::vector<std::vector<std::array<int, 3> > > PPAdj;
    // std::vector<std::vector<std::array<int, 3> > > PEAdj;
    std::vector<std::array<int, 3> > PPAdj;
    std::vector<int> Variables;
  };

  struct PatchInfo {
    std::vector<PatchSide> Sides;
  };

  std::vector<PatchInfo> PatchSides;

  std::vector<std::vector<std::array<int, 3> > > VarSide;
  int NumVariables;

  void SortFromSeparator(std::vector<Geo::Point3<ScalarType>> &SideSeq,
                         const std::set<Geo::Point3<ScalarType>> &Separator) {
    // sort the sequence from the first corner
    int FistCorner = -1;
    for (size_t j = 0; j < SideSeq.size(); j++) {
      Geo::Point3<ScalarType> VertPos = SideSeq[j];
      if (Separator.count(VertPos) == 0)
        continue;
      FistCorner = j;
    }
    assert(FistCorner != -1);
    std::rotate(SideSeq.begin(), SideSeq.begin() + FistCorner, SideSeq.end());
  }

  void SplitIntoSubSequences(
      const std::vector<Geo::Point3<ScalarType>> &SideSeq,
      const std::set<Geo::Point3<ScalarType>> &Separator,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SplitSeq,
      bool circular=true) {

    assert(Separator.count(SideSeq[0])> 0);

    if (SideSeq.empty())
      return;

    std::vector<Geo::Point3<ScalarType>> TestSeq = SideSeq;
    // repeat the first to the enc
    TestSeq.push_back(TestSeq[0]);

    //add one sequence
    SplitSeq.resize(1);
    //set the first by defualt
    SplitSeq.back().push_back(SideSeq[0]);
    //start from second
    for (size_t i = 1; i < SideSeq.size(); ++i) {
      Geo::Point3<ScalarType> point = SideSeq[i];
      
      bool AddNewSide=(Separator.count(point)>0);

      if (!circular)
        AddNewSide&=(i != (SideSeq.size()-1));

      if (AddNewSide)
      {
        //repeat the separator at the end
        SplitSeq.back().push_back(point);
        //add new sequence
        SplitSeq.resize(SplitSeq.size()+1);
      }

      //add the point
      SplitSeq.back().push_back(point);
    }
    //re add the first one at the end
    if (circular)
        SplitSeq.back().push_back(SideSeq[0]);
    // // If the first element is a separator, make it circular
    // if (!SplitSeq.empty() && Separator.find(SideSeq.front()) !=
    // Separator.end()) {
    //     SplitSeq.front().insert(SplitSeq.front().begin(), SideSeq.front());
    // }
  }

  void InitPatchAdj() 
  {
    VarSide.clear();
    VarSide.resize(NumVariables);
    //resize with max value
    for (size_t i=0;i<NumPatches();i++)
    for (size_t j=0;j<NumSides(i);j++)
    {
      int CurrSubSides=NumSubSides(i,j);
      PatchSides[i].Sides[j].PPAdj.resize(CurrSubSides,std::array<int,3>{-1,-1,-1});

      for (size_t k=0;k<CurrSubSides;k++)
      {
        int CurrVar=Variable(i,j,k);
        std::array<int, 3> Entry;
        Entry[0]=i;
        Entry[1]=j;
        Entry[2]=k;
        VarSide[CurrVar].push_back(Entry);
      }
    }

    for (size_t i=0;i<VarSide.size();i++)
    {
        std::vector<std::array<int, 3> > CurrSide=VarSide[i];

        //manifold check
        assert((CurrSide.size()==1)||(CurrSide.size()==2));

        //in this case is border
        if (CurrSide.size()==1)continue;

        std::array<int, 3> Pside0=CurrSide[0];
        int IndexPatch0=Pside0[0];
        int IndexSide0=Pside0[1];
        int IndexSubS0=Pside0[2];
        
        std::array<int, 3> Pside1=CurrSide[1];
        int IndexPatch1=Pside1[0];
        int IndexSide1=Pside1[1];
        int IndexSubS1=Pside1[2];
        
        PatchSides[IndexPatch0].Sides[IndexSide0].PPAdj[IndexSubS0]=Pside1;
        PatchSides[IndexPatch1].Sides[IndexSide1].PPAdj[IndexSubS1]=Pside0;
    }

    //write a bit of stats
    int BorderSide=0;
    for (size_t i=0;i<NumPatches();i++)
    for (size_t j=0;j<NumSides(i);j++)
    for (size_t k=0;k<NumSubSides(i,j);k++)
    {
        std::array<int, 3> Pside=PatchSides[i].Sides[j].PPAdj[k];
        int IndexPatch=Pside[0];
        int IndexSide=Pside[1];
        int IndexSubS=Pside[2];
        if (Pside==std::array<int, 3>{-1,-1,-1})BorderSide++;
    }
    std::cout<<"Num Border Patch:"<<BorderSide<<std::endl;
  }

  void InitVariables() 
  {
    for (size_t i=0;i<NumPatches();i++)
    for (size_t j=0;j<NumSides(i);j++)
    {
      int CurrSubNum=NumSubSides(i,j);
      PatchSides[i].Sides[j].Variables=std::vector<int>(CurrSubNum,-1);
    }

    std::map<std::vector<Geo::Point3<ScalarType> >,int > SideVarMap;
    int CurrVar=0;
    for (size_t i=0;i<NumPatches();i++)
    for (size_t j=0;j<NumSides(i);j++)
    {
      std::vector<std::vector<Geo::Point3<ScalarType> > > SubSides;
      GetSidePoints(i,j,SubSides);
      for (size_t k=0;k<SubSides.size();k++)
      {
        std::vector<Geo::Point3<ScalarType> > KeySide=SubSides[k];
        std::sort(KeySide.begin(),KeySide.end());
        if (SideVarMap.count(KeySide)==0)
        {
            SideVarMap[KeySide]=CurrVar;
            Variable(i,j,k)=CurrVar;
            CurrVar++;
        }
        else {
            Variable(i,j,k)=SideVarMap[KeySide];
        }
      }
    }
    NumVariables=CurrVar;
    std::cout<<"There are "<<NumVariables<<" Variables"<<std::endl;
  }

public:

  int NumPatches()const
  {
    return PatchSides.size();
  }

  int NumSides(const int &PatchIdx)const
  {
    assert(PatchIdx>=0);
    assert(PatchIdx<NumPatches());
    return PatchSides[PatchIdx].Sides.size();
  }
  
  int NumSubSides(const int &PatchIdx,
                  const int &SideIdx)const
  {
    assert(PatchIdx>=0);
    assert(PatchIdx<NumPatches());
    assert(SideIdx>=0);
    assert(SideIdx<NumSides(PatchIdx));
    return (PatchSides[PatchIdx].Sides[SideIdx].SubSides.size());
  }
  
  int &Variable(int IndexP,int IndexS,int IndexSubS)
  {
    assert(IndexP>=0);
    assert(IndexS>=0);
    assert(IndexSubS>=0);
    assert(IndexP<NumPatches());
    assert(IndexS<NumSides(IndexP));
    assert(IndexSubS<NumSubSides(IndexP,IndexS));
    return(PatchSides[IndexP].Sides[IndexS].Variables[IndexSubS]);
  }

  void GetSidePoints(const int &PatchIdx,
                     const int &SideIdx,
                     std::vector<std::vector<Geo::Point3<ScalarType> > > &SubSides)const
  {
    assert(PatchIdx>=0);
    assert(PatchIdx<NumPatches());
    assert(SideIdx>=0);
    assert(SideIdx<NumSides(PatchIdx));
    SubSides=PatchSides[PatchIdx].Sides[SideIdx].SubSides;
  }

  void InitFromLabels(const std::vector<Geo::Point3<ScalarType>> &Vert,
                      const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::vector<int>> &PatchFacesIdx,
                      const std::vector<std::vector<int>> &PatchCorners) {

    std::set<Geo::Point3<ScalarType>> AllCornerPos;
    std::vector<std::set<Geo::Point3<ScalarType>>> PatchCornerPos;
    PatchCornerPos.resize(PatchCorners.size());

    for (size_t i = 0; i < PatchCorners.size(); i++)
      for (size_t j = 0; j < PatchCorners[i].size(); j++) {
        int IndexV = PatchCorners[i][j];
        assert(IndexV >= 0);
        assert(IndexV < Vert.size());
        Geo::Point3<ScalarType> CornPos = Vert[IndexV];
        AllCornerPos.insert(CornPos);
        PatchCornerPos[i].insert(CornPos);
      }

    PatchSides.resize(PatchFacesIdx.size());

    // then get all boundaries
    for (size_t i = 0; i < PatchFacesIdx.size(); i++) {
      std::vector<std::vector<int>> SubMeshFaces;
      std::vector<Geo::Point3<ScalarType>> SubMeshVert;
      std::vector<int> OriginalF;
      std::vector<int> OriginalV;

      // make a copy of the mesh
      CopySubMesh<ScalarType>(Vert, Faces, PatchFacesIdx[i], SubMeshVert,
                              SubMeshFaces, OriginalV, OriginalF);

      // then get the boundaries
      std::vector<std::vector<int>> BorderVIdxAll;
      getBorderSequncesVert(SubMeshFaces, BorderVIdxAll);
      if (BorderVIdxAll.size() != 1) {
        std::cout << "MULTIPLE BORDERS" << std::endl;
        continue;
      }
      std::vector<int> BorderVIdx = BorderVIdxAll[0];

      std::vector<Geo::Point3<ScalarType>> BorderVPos;
      for (size_t j = 0; j < BorderVIdx.size(); j++)
        BorderVPos.push_back(SubMeshVert[BorderVIdx[j]]);

      // sort from the first
      SortFromSeparator(BorderVPos, PatchCornerPos[i]);

      // split in sequences
      std::vector<std::vector<Geo::Point3<ScalarType>>> SideSeq;
      SplitIntoSubSequences(BorderVPos, PatchCornerPos[i], SideSeq);

      if (SideSeq.size()!=4)
      {
        std::cout<<"WARNING Side Size: "<<SideSeq.size()<<std::endl;
        std::cout<<"Corner Size: "<<PatchCornerPos[i].size()<<std::endl;
      }

      //then split subsequences
      PatchSides[i].Sides.resize(SideSeq.size());

      for (size_t j=0;j<SideSeq.size();j++)
      {
        std::vector<std::vector<Geo::Point3<ScalarType>>> SubSideSeq;
        SplitIntoSubSequences(SideSeq[j], AllCornerPos, SubSideSeq,false);
        PatchSides[i].Sides[j].SubSides=SubSideSeq;
        if (SubSideSeq.size()!=1)
            std::cout<<"Sub Seq Size: "<<SubSideSeq.size()<<std::endl;
      }
    //   std::cout<<"Seq Size: "<<SideSeq.size()<<std::endl;
    }
    InitVariables();
    InitPatchAdj();
  }

  SideDiscretization()
  {
    NumVariables=0;
  }
};

} // namespace Geo
#endif
