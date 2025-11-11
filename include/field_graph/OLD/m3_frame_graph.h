#ifndef M3_FRAME_GRAPH
#define M3_FRAME_GRAPH

#include "triangle2.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <set>
#include <stdio.h>
//#include <sys/_types/_size_t.h>
#include <unordered_set>
#include <utility>
#include <vector>


namespace Geo {
template <class ScalarType> class UndirectedNodeWeightedGraph {
  std::vector<bool> EnabledN;
  //for each node for each label if exit or nor
  std::vector<std::vector<bool> > ExitLab;
  std::vector<ScalarType> W;
  std::vector<std::vector<int>> L;
  std::vector<std::vector<int>> Label;
  std::vector<std::vector<bool>> EnabledL;
  std::vector<NType> NT;

public:

  const NType &NodeType(const size_t &IndexN) const {
    assert(IndexN >= 0);
    assert(IndexN < NT.size());
    return NT[IndexN];
  }

  NType &NodeType(const size_t &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < NT.size());
    return NT[IndexN];
  }

  ScalarType &Weight(const size_t &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < NT.size());
    return W[IndexN];
  }

  ScalarType Weight(const size_t &IndexN)const {
    assert(IndexN >= 0);
    assert(IndexN < NT.size());
    return W[IndexN];
  }

  size_t NumNodes() const {
    assert(L.size() == EnabledN.size());
    assert(L.size() == Label.size());
    assert(L.size() == EnabledL.size());
    assert(L.size() == NT.size());
    return (L.size());
  }

  int NumDisabledN()
  {
    int NumDis=0;
    for (size_t i=0;i<EnabledN.size();i++)
    if (!EnabledN[i])NumDis++;
    return NumDis;
  }

  ScalarType AverageVal()
  {
    int TotL=0;
    for (size_t i=0;i<L.size();i++)
      TotL+=L[i].size();
    return ((ScalarType)TotL/(ScalarType)L.size());
  }

  ScalarType AverageEnabledVal()
  {
    int TotL=0;
    for (size_t i=0;i<L.size();i++)
    {
      for (size_t j=0;j<L[i].size();j++)
      {
        if (!EnabledL[i][j])continue;
        if (!EnabledN[i])continue;
        TotL++;
      }
    }
    return ((ScalarType)TotL/(ScalarType)L.size());
  }

  void GetNodeOfType(const NType &TypeNode, std::vector<int> &IdxNode) const {
    IdxNode.clear();
    for (size_t i = 0; i < NumNodes(); i++) {
      if (NodeType(i) != TypeNode)
        continue;
      IdxNode.push_back(i);
    }
  }

  bool IsEnabled(const int &IndexN) const {
    assert(IndexN >= 0);
    assert(IndexN < EnabledN.size());
    return (EnabledN[IndexN]);
  }

  const bool IsExitLabel(const int &IndexN,const int & IndexLab) const {
    assert(IndexN >= 0);
    assert(IndexN < ExitLab.size());
    assert(IndexLab >= 0);
    assert(IndexLab < ExitLab[IndexN].size());
    return (ExitLab[IndexN][IndexLab]);
  }

  void SetExitLabel(const int &IndexN,const int &IndexLab){
    assert(IndexN >= 0);
    assert(IndexN < ExitLab.size());
    assert(IndexLab >= 0);
    assert(IndexLab < ExitLab[IndexN].size());
    ExitLab[IndexN][IndexLab]=true;
  }

  bool IsEnabledL(const int &IndexN,int IndexL) const {
    assert(IndexN >= 0);
    assert(IndexN < EnabledL.size());
    assert(IndexL >= 0);
    assert(IndexL < EnabledL[IndexN].size());
    return (EnabledL[IndexN][IndexL]);
  }

  void EnableNode(const int &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < EnabledN.size());
    EnabledN[IndexN] = true;
  }

  void DisableNode(const int &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < EnabledN.size());
    EnabledN[IndexN] = false;
  }

  int NumNeigh(const int &IndexN)const
  {
    assert(IndexN>=0);
    assert(IndexN<L.size());
    return (L[IndexN].size());
  }

  int NeighNode(const int &IndexN,const int &NeighIdx)const
  {
    assert(IndexN>=0);
    assert(IndexN<L.size());
    assert(NeighIdx>=0);
    assert(NeighIdx<L[IndexN].size());
    return (L[IndexN][NeighIdx]);
  }

  int NeighLabel(const int &IndexN,const int &NeighIdx)const
  {
    assert(IndexN>=0);
    assert(IndexN<L.size());
    assert(NeighIdx>=0);
    assert(NeighIdx<Label[IndexN].size());
    return (Label[IndexN][NeighIdx]);
  }

  void AddSymmLink(const int &IndexN0, 
                  const int &IndexN1, 
                  const int &Label0,
                  const int &Label1,
                  bool Enabled=true) {
    assert(IndexN0 >= 0);
    assert(IndexN0 < (int)L.size());
    assert(IndexN0 < (int)EnabledL.size());

    assert(IndexN1 >= 0);
    assert(IndexN1 < (int)L.size());
    assert(IndexN1 < (int)EnabledL.size());

    assert(L[IndexN0].size() == EnabledL[IndexN0].size());
    assert(L[IndexN1].size() == EnabledL[IndexN1].size());

    // add it symmetrically
    L[IndexN0].push_back(IndexN1);
    L[IndexN1].push_back(IndexN0);
    Label[IndexN0].push_back(Label0);
    Label[IndexN1].push_back(Label1);
    EnabledL[IndexN0].push_back(Enabled);
    EnabledL[IndexN1].push_back(Enabled);
  }

  void MakeUniqueConnection() {

    struct LinkData
    {
      int Label;
      bool Enabled;

      LinkData(int _Label,bool _Enabled)
      {
        Label=_Label;
        Enabled=_Enabled;
      }

      LinkData()
      {
        Label=-1;
        Enabled=false;
      }

      bool operator ==(const LinkData &LD)const
      {
        if (Label!=LD.Label)return false;
        if (Enabled!=LD.Enabled)return false;
        return true;
      }
    };

    std::map<std::pair<int, int>, LinkData> GraphLinks;
    assert(L.size() == Label.size());
    assert(L.size() == EnabledL.size());

    for (size_t i = 0; i < L.size(); i++) {
      assert(L[i].size() == Label[i].size());
      assert(L[i].size() == EnabledL[i].size());
      
      for (size_t j = 0; j < L[i].size(); j++) 
      {
        int N0=i;
        int N1=L[i][j];
        int Lab=Label[i][j];
        bool IsEnabled=EnabledL[i][j];
        LinkData LData(Lab,IsEnabled);

        std::pair<int, int> Key(N0,N1);
        if (GraphLinks.count(Key)==0)
           GraphLinks[Key]=LData;

        assert(GraphLinks[Key]==LData);
      }
    }

    for (size_t i = 0; i < L.size(); i++)
    {
      L[i].clear();
      Label[i].clear();
      EnabledL[i].clear();
    }

    typename std::map<std::pair<int, int>, LinkData>::const_iterator IteMap;
    for (IteMap=GraphLinks.begin();IteMap!=GraphLinks.end();IteMap++)
    {
      int N0=(*IteMap).first.first;
      int N1=(*IteMap).first.second;
      int Lab=(*IteMap).second.Label;
      bool IsEn=(*IteMap).second.Enabled;

      L[N0].push_back(N1);
      Label[N0].push_back(Lab);
      EnabledL[N0].push_back(IsEn);
    }
  }

  void AllocateNodes(int NumNodes) {
    L.clear();
    L.resize(NumNodes);
    Label.clear();
    Label.resize(NumNodes);
    ExitLab.clear();
    W.clear();
    W.resize(NumNodes, -1);

    NT.clear();
    NT.resize(NumNodes, Internal);

    EnabledN.clear();
    EnabledN.resize(NumNodes, true);
    EnabledL.clear();
    EnabledL.resize(NumNodes);
    ExitLab.resize(NumNodes,std::vector<bool>(3,false));
  }

  void WriteStats()
  {
    std::vector<int> TangentN,InternalN;

    std::cout<<"Number of Nodes: "<<NumNodes()<<std::endl;
    GetNodeOfType(Geo::TangentDomain,TangentN);

    GetNodeOfType(Geo::TangentDomain,TangentN);
    std::cout<<"Number of Tangent Nodes: "<<TangentN.size()<<std::endl;

    GetNodeOfType(Geo::Internal,InternalN);
    std::cout<<"Number of Internal Nodes: "<<InternalN.size()<<std::endl;

    int disNodes=NumDisabledN();
    std::cout<<"Number of Disabled Nodes: "<<disNodes<<std::endl;

    std::cout<<"Average Valence: "<<AverageVal()<<std::endl;
    std::cout<<"Average Enabled Valence: "<<AverageEnabledVal()<<std::endl;
  }

};

} // namespace Geo
#endif
