#ifndef MESH_FUNCTIONS
#define MESH_FUNCTIONS

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "point3.h"
#include "point2.h"
#include "connectivity.h"
#include <pos3.h>

template <class ScalarType>
bool HasUnfererencedVert(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                         const std::vector<std::vector<int> > &Faces)
{
    std::vector<bool> Referred(VertPos.size(),false);
    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[i].size();j++)
            Referred[Faces[i][j]]=true;

    for (size_t i=0;i<Referred.size();i++)
        if (!Referred[i])return true;

    return false;
}

template <class ScalarType>
void RemoveUnfererencedVert(std::vector<Geo::Point3<ScalarType> > &VertPos,
                            std::vector<std::vector<int> > &Faces)
{
    std::vector<bool> Referred(VertPos.size(),false);
    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[i].size();j++)
            Referred[Faces[i][j]]=true;

    std::vector<int> NewIdx(VertPos.size(),-1);
    std::vector<Geo::Point3<ScalarType> > NewVertPos;

    int currIdx=0;
    for (size_t i=0;i<Referred.size();i++)
        if (Referred[i])
        {
            NewIdx[i]=currIdx;
            currIdx++;
            NewVertPos.push_back(VertPos[i]);
        }

    VertPos=NewVertPos;
    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int currIdx=Faces[i][j];
            assert(currIdx<NewIdx.size());
            int newIdx=NewIdx[currIdx];
            assert(newIdx>=0);
            assert(newIdx<VertPos.size());
            Faces[i][j]=newIdx;
        }
}


void NumFaceValences(const std::vector<std::vector<int> > &Faces,
                     std::map<int,int> &ValMap)
{
    ValMap.clear();
    for (size_t i=0;i<Faces.size();i++)
    {
        int Valence=Faces[i].size();
        if (ValMap.find(Valence) == ValMap.end())
            ValMap[Valence]=0;

        ValMap[Valence]++;
    }
}

bool IsAllQuad(const std::vector<std::vector<int> > &Faces)
{
    std::map<int,int> ValMap;
    NumFaceValences(Faces,ValMap);
    return(ValMap[4]==(int)Faces.size());
}

bool IsAllTris(const std::vector<std::vector<int> > &Faces)
{
    std::map<int,int> ValMap;
    NumFaceValences(Faces,ValMap);
    return(ValMap[3]==(int)Faces.size());
}

void VertValences(const std::vector<std::vector<int> > &Faces,
                  const int NumV,
                  std::vector<int> &VertValence)
{
    VertValence.clear();
    VertValence.resize(NumV,0);
    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int IndexV=Faces[i][j];
            assert(IndexV>=0);
            assert(IndexV<(int)VertValence.size());
            VertValence[IndexV]++;
        }
}

void SingularVertices(const std::vector<std::vector<int> > &Faces,
                      const std::vector<bool> &BorderV,
                      std::vector<int>  &SingularVert)
{
    std::vector<int> VertValence;
    VertValences(Faces,BorderV.size(),VertValence);
    for (size_t i=0;i<VertValence.size();i++)
    {
        if (BorderV[i]&&(VertValence[i]>2))SingularVert.push_back(i);
        if (!BorderV[i]&&(VertValence[i]!=4))SingularVert.push_back(i);
    }
}

enum MType{MT_Triangle,MT_Tetra,MT_Other};

MType MeshType(const std::vector<std::vector<int> > &Elems)
{
    bool AllTris=true;
    bool AllTetra=true;
    for (size_t i=0;i<Elems.size();i++)
    {
        AllTris&=(Elems[i].size()==3);
        AllTetra&=(Elems[i].size()==4);
        if ((!AllTris)&&(!AllTetra))break;
    }
    assert(!(AllTris && AllTetra));
    if (AllTris)return MT_Triangle;
    if (AllTetra)return MT_Tetra;
    return MT_Other;
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                const std::vector<std::vector<int> > &Faces,
                const int &IndexF)
{
    assert(Faces[IndexF].size()==3);
    Geo::Point3<ScalarType> P0=VertPos[Faces[IndexF][0]];
    Geo::Point3<ScalarType> P1=VertPos[Faces[IndexF][1]];
    Geo::Point3<ScalarType> P2=VertPos[Faces[IndexF][2]];
    ScalarType FaceArea=(Geo::Cross((P1-P0),(P2-P0))).Norm()/2;
    return FaceArea;
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point2<ScalarType> > &VertPos,
                const std::vector<std::vector<int> > &Faces,
                const int &IndexF,bool unsignedArea=true)
{
    assert(Faces[IndexF].size()==3);
    Geo::Point2<ScalarType> P0=VertPos[Faces[IndexF][0]];
    Geo::Point2<ScalarType> P1=VertPos[Faces[IndexF][1]];
    Geo::Point2<ScalarType> P2=VertPos[Faces[IndexF][2]];
    ScalarType FaceArea=(Geo::Cross((P1-P0),(P2-P0))/2);

    if (unsignedArea)
        return fabs(FaceArea);
    else
        return (FaceArea);
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                const std::vector<std::vector<int> > &Faces)
{
    ScalarType AreaVal=0;
    for (size_t i=0; i<Faces.size(); i++)
        AreaVal+=Area(VertPos,Faces,i);
    return AreaVal;
}

template <class ScalarType>
ScalarType Area(const std::vector<Geo::Point2<ScalarType> > &VertPos,
                const std::vector<std::vector<int> > &Faces,
                bool unsignedArea=true)
{
    ScalarType AreaVal=0;
    for (size_t i=0; i<Faces.size(); i++)
        AreaVal+=Area(VertPos,Faces,i,unsignedArea);
    return AreaVal;
}

template <class ScalarType>
int ClosestVert(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                const Geo::Point3<ScalarType> &TargetPos)
{
    ScalarType DistMin=std::numeric_limits<ScalarType>::max();
    int IndexV=-1;
    for (size_t i=0;i<VertPos.size();i++)
    {
        ScalarType TestDist=(VertPos[i]-TargetPos).Norm();
        if (TestDist>DistMin)continue;
        DistMin=TestDist;
        IndexV=i;
    }
    return IndexV;
}

template <class ScalarType>
ScalarType EdgeLenght(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                      const std::vector<std::vector<int> > &Faces,
                      const int &IndexF,const int &IndexE)
{
    assert((IndexF>=0)&&(IndexF<(int)Faces.size()));
    assert((IndexE>=0)&&(IndexE<(int)Faces[IndexF].size()));

    int IndexV0=Faces[IndexF][IndexE];
    int IndexV1=Faces[IndexF][(IndexE+1)%Faces[IndexF].size()];
    Geo::Point3<ScalarType> P0=VertPos[IndexV0];
    Geo::Point3<ScalarType> P1=VertPos[IndexV1];
    return (P1-P0).Norm();
}

template <class ScalarType>
ScalarType EdgeLenght(const std::vector<Geo::Point2<ScalarType> > &VertPos,
                      const std::vector<std::vector<int> > &Faces,
                      const int &IndexF,const int &IndexE)
{
    assert((IndexF>=0)&&(IndexF<(int)Faces.size()));
    assert((IndexE>=0)&&(IndexE<(int)Faces[IndexF].size()));

    int IndexV0=Faces[IndexF][IndexE];
    int IndexV1=Faces[IndexF][(IndexE+1)%Faces[IndexF].size()];
    Geo::Point2<ScalarType> P0=VertPos[IndexV0];
    Geo::Point2<ScalarType> P1=VertPos[IndexV1];
    return (P1-P0).Norm();
}


template <class ScalarType>
void LongestEdge(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                 const std::vector<std::vector<int> > &Faces,
                 int &IndexF,int &IndexE)
{
    ScalarType LargestELenght=0;
    IndexF=-1;
    IndexE=-1;
    for (size_t i=0; i<Faces.size(); i++)
        for (size_t j=0; j<Faces[i].size(); j++)
        {
            //            int IndexV0=Faces[i][j];
            //            int IndexV1=Faces[i][(j+1)%Faces[i].size()];
            //            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
            //            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
            //            ScalarType CurrELenght=(P0-P1).Norm();

            ScalarType CurrELenght=EdgeLenght(VertPos,Faces,i,j);

            if (LargestELenght>CurrELenght)continue;

            LargestELenght=CurrELenght;
            IndexF=i;
            IndexE=j;
        }
}

template <class ScalarType>
ScalarType AvgEdgeLen(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                      const std::vector<std::vector<int> > &Faces)
{
    ScalarType EdgeSum=0;
    int Num=0;
    for (size_t i=0; i<Faces.size(); i++)
        for (size_t j=0; j<Faces[i].size(); j++)
        {
            ScalarType CurrELenght=EdgeLenght(VertPos,Faces,i,j);
            //should consider also border but become slow
            EdgeSum+=CurrELenght;
            Num++;
        }
    return EdgeSum/Num;
}

template <class ScalarType>
ScalarType FaceVertAngleDeg(const std::vector<Geo::Point2<ScalarType> > &VertPos,
                            const std::vector<std::vector<int> > &Faces,
                            const int &IndexF,const int &IndexV)
{
    assert((IndexF>=0)&&(IndexF<(int)Faces.size()));
    assert((IndexV>=0)&&(IndexV<(int)Faces[IndexF].size()));

    size_t NumV=Faces[IndexF].size();

    int IndexV0=Faces[IndexF][(IndexV+NumV-1)%NumV];
    int IndexV1=Faces[IndexF][IndexV];
    int IndexV2=Faces[IndexF][(IndexV+1)%NumV];

    Geo::Point2<ScalarType> p0=VertPos[IndexV0];
    Geo::Point2<ScalarType> p1=VertPos[IndexV1];
    Geo::Point2<ScalarType> p2=VertPos[IndexV2];

    Geo::Point2<ScalarType>  dir0=p0-p1;
    Geo::Point2<ScalarType>  dir1=p2-p1;
    dir0.Normalize();
    dir1.Normalize();
    ScalarType angle=Geo::AngleDeg(dir0,dir1);
    return angle;
}

template <class ScalarType>
ScalarType FaceVertAngleDeg(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                            const std::vector<std::vector<int> > &Faces,
                            const int &IndexF,const int &IndexV)
{
    assert((IndexF>=0)&&(IndexF<(int)Faces.size()));
    assert((IndexV>=0)&&(IndexV<(int)Faces[IndexF].size()));

    size_t NumV=Faces[IndexF].size();

    int IndexV0=Faces[IndexF][(IndexV+NumV-1)%NumV];
    int IndexV1=Faces[IndexF][IndexV];
    int IndexV2=Faces[IndexF][(IndexV+1)%NumV];

    Geo::Point3<ScalarType> p0=VertPos[IndexV0];
    Geo::Point3<ScalarType> p1=VertPos[IndexV1];
    Geo::Point3<ScalarType> p2=VertPos[IndexV2];

    Geo::Point3<ScalarType>  dir0=p0-p1;
    Geo::Point3<ScalarType>  dir1=p2-p1;
    dir0.Normalize();
    dir1.Normalize();
    ScalarType angle=Geo::AngleDeg(dir0,dir1);
    return angle;
}

template <class ScalarType>
void InterpolatePerFace(const std::vector<Geo::Point3<ScalarType> > &VertVal,
                        const std::vector<std::vector<int> > &Faces,
                        std::vector<Geo::Point3<ScalarType> > &FaceVal)
{
    FaceVal.clear();
    FaceVal.resize(Faces.size(),Geo::Point3<ScalarType>(0,0,0));
    for (size_t i=0; i<Faces.size(); i++)
    {
        for (size_t j=0; j<Faces[i].size(); j++)
        {
            int IndexV=Faces[i][j];
            assert(IndexV<(int)VertVal.size());
            FaceVal[i]+=VertVal[IndexV];
        }
        FaceVal[i]/=Faces[i].size();
    }
}

template <class ScalarType>
void SplitQuadIntoTris(const std::vector<Geo::Point3<ScalarType> > &Pos,
                       const std::vector<std::vector<int> > &Quads,
                       std::vector<std::vector<int> > &Tris)
{
    Tris.clear();
    //2 triangles for each quad
    Tris.resize((int)Quads.size()*2,std::vector<int> (3,-1));
    for (size_t i=0;i<Quads.size();i++)
    {
        assert(Quads[i].size()==4);
        int IndexV0=Quads[i][0];
        int IndexV1=Quads[i][1];
        int IndexV2=Quads[i][2];
        int IndexV3=Quads[i][3];

        Geo::Point3<ScalarType> PosV0=Pos[IndexV0];
        Geo::Point3<ScalarType> PosV1=Pos[IndexV1];
        Geo::Point3<ScalarType> PosV2=Pos[IndexV2];
        Geo::Point3<ScalarType> PosV3=Pos[IndexV3];

        //split along shortest distance
        ScalarType D0=(PosV0-PosV2).Norm();
        ScalarType D1=(PosV1-PosV3).Norm();

        int IndexTriF0=i*2;
        int IndexTriF1=(i*2)+1;
        if (D0<D1)
        {
            Tris[IndexTriF0][0]=IndexV0;
            Tris[IndexTriF0][1]=IndexV1;
            Tris[IndexTriF0][2]=IndexV2;
            Tris[IndexTriF1][0]=IndexV0;
            Tris[IndexTriF1][1]=IndexV2;
            Tris[IndexTriF1][2]=IndexV3;
        }
        else
        {
            Tris[IndexTriF0][0]=IndexV0;
            Tris[IndexTriF0][1]=IndexV1;
            Tris[IndexTriF0][2]=IndexV3;
            Tris[IndexTriF1][0]=IndexV1;
            Tris[IndexTriF1][1]=IndexV2;
            Tris[IndexTriF1][2]=IndexV3;
        }
    }
}

int WhichVert(const std::vector<int> &Face,
              const int &VertIndex)
{
    for (size_t i=0;i<Face.size();i++)
        if (Face[i]==VertIndex)return i;
    return -1;
}

void FFConnectedComponenent(const std::vector<std::vector<int> > &Elements,
                            const std::vector<std::vector<int> > &NextP,
                            const std::vector<std::vector<int> > &NextE,
                            std::vector<std::vector<int> > &ComponentsF)
{
    ComponentsF.clear();

    //find partitions
    std::vector<bool> Has_Explored(Elements.size(),false);

    for (size_t i=0;i<Elements.size();i++)
    {
        if (Has_Explored[i])continue;

        ComponentsF.resize(ComponentsF.size()+1);

        std::vector<size_t> StackP;
        StackP.push_back(i);
        do
        {
            int CurrP=StackP.back();
            StackP.pop_back();

            if (Has_Explored[CurrP])continue;

            Has_Explored[CurrP]=true;
            ComponentsF.back().push_back(CurrP);


            for (size_t n=0;n<NextP[CurrP].size();n++)
            {
                //get next elem
                int AdjP=NextP[CurrP][n];
                assert(AdjP!=CurrP);

                //get edge as seen from the other face
                int AdjE=NextE[CurrP][n];

                if (AdjP==-1)continue;
                assert(AdjE!=-1);

                bool has_to_diffuse=(!Has_Explored[AdjP]);
                if (has_to_diffuse)
                {
                    StackP.push_back(AdjP);
                }
            }
        }while (StackP.size()>0);
    }

}

void FFConnectedComponenent(const std::vector<std::vector<int> > &Elements,
                            const std::vector<std::vector<int> > &NextP,
                            const std::vector<std::vector<int> > &NextE,
                            std::vector<std::vector<int> > &ComponentsF,
                            std::vector<std::vector<int> > &ComponentsV)
{
    FFConnectedComponenent(Elements,NextP,NextE,ComponentsF);
    ComponentsV.resize(ComponentsF.size());

    for (size_t c=0;c<ComponentsF.size();c++)
    {
        for (size_t i=0;i<ComponentsF[c].size();i++)
        {
            int IndexF=ComponentsF[c][i];
            assert(IndexF>=0);
            assert(IndexF<(int)Elements.size());
            for (size_t j=0;j<Elements[IndexF].size();j++)
                ComponentsV[c].push_back(Elements[IndexF][j]);
        }

        std::sort(ComponentsV[c].begin(),ComponentsV[c].end());
        std::vector<int>::iterator last = std::unique(ComponentsV[c].begin(),ComponentsV[c].end());
        ComponentsV[c].erase(last, ComponentsV[c].end());
    }
}

void FFConnectedComponenent(const std::vector<std::vector<int> > &Elements,
                            std::vector<std::vector<int> > &ComponentsF,
                            std::vector<std::vector<int> > &ComponentsV)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Elements, NextF,NextE);
    FFConnectedComponenent(Elements,NextF,NextE,ComponentsF,ComponentsV);
}

template <class ScalarType>
void CopySubMesh(const std::vector<Geo::Point3<ScalarType> > &Vert,
                 const std::vector<std::vector<int> > &Elements,
                 const std::vector<int> &ToCopyEl,
                 std::vector<Geo::Point3<ScalarType> > &CopiedMeshesVert,
                 std::vector<std::vector<int> > &CopiedMeshesElem,
                 std::vector<int> &OriginalVert,
                 std::vector<int> &OriginalElem)
{
    OriginalVert.clear();
    OriginalElem.clear();
    CopiedMeshesVert.clear();
    CopiedMeshesElem.clear();

    //copy faces
    for (size_t i=0;i<ToCopyEl.size();i++)
    {
        int IndexEl=ToCopyEl[i];
        assert(IndexEl>=0);
        assert(IndexEl<Elements.size());

        CopiedMeshesElem.push_back(Elements[IndexEl]);
        OriginalElem.push_back(IndexEl);
    }

    //then copy/reassign vertices
    std::map<int,int> VertToVert;
    for (size_t i=0;i<CopiedMeshesElem.size();i++)
        for (size_t j=0;j<CopiedMeshesElem[i].size();j++)
        {
            int IndexV=CopiedMeshesElem[i][j];
            assert(IndexV>=0);
            assert(IndexV<Vert.size());
            if (VertToVert.count(IndexV)==0)
            {
                CopiedMeshesVert.push_back(Vert[IndexV]);
                OriginalVert.push_back(IndexV);
                VertToVert[IndexV]=CopiedMeshesVert.size()-1;
            }
            int NewIndexV=VertToVert[IndexV];
            CopiedMeshesElem[i][j]=NewIndexV;
    }
}

template <class ScalarType>
void SplitByFFConnectedComponenent(const std::vector<Geo::Point3<ScalarType> > &Vert,
                            const std::vector<std::vector<int> > &Elements,
                            const std::vector<std::vector<int> > &NextP,
                            const std::vector<std::vector<int> > &NextE,
                            std::vector<std::vector<Geo::Point3<ScalarType> > > &SubMeshesV,
                            std::vector<std::vector<std::vector<int> > > &SubMeshesF,
                            std::vector<std::vector<int> > &OriginalV,
                            std::vector<std::vector<int> > &OriginalF)
{
    std::vector<std::vector<int> > Components;
    FFConnectedComponenent(Elements,NextP,NextE,Components);

    SubMeshesV.clear();
    SubMeshesV.resize(Components.size());
    SubMeshesF.clear();
    SubMeshesF.resize(Components.size());


    OriginalV.clear();
    OriginalV.resize(Components.size());
    OriginalF.clear();
    OriginalF.resize(Components.size());
    for (size_t i=0;i<Components.size();i++)
        CopySubMesh(Vert,Elements,Components[i],
                    SubMeshesV[i],SubMeshesF[i],
                    OriginalV[i],OriginalF[i]);
}


template <class ScalarType>
void SplitByFFConnectedComponenent(const std::vector<Geo::Point3<ScalarType> > &Vert,
                            const std::vector<std::vector<int> > &Elements,
                            std::vector<std::vector<Geo::Point3<ScalarType> > > &SubMeshesV,
                            std::vector<std::vector<std::vector<int> > > &SubMeshesF,
                            std::vector<std::vector<int> > &OriginalV,
                            std::vector<std::vector<int> > &OriginalF)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Elements, NextF,NextE);
    SplitByFFConnectedComponenent(Vert,Elements,NextF,NextE,SubMeshesV,
                                  SubMeshesF,OriginalV,OriginalF);
}

template <class ScalarType>
void SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType> > &Pos,
                        std::vector<std::vector<int> > &Faces,
                        const  std::vector<std::vector<int> > &NextF,
                        const  std::vector<std::vector<int> > &NextE,
                        const std::vector<std::pair<int,int> > &SplitEdges)
{
    std::set<std::pair<int,int> > SplitEdgeSet(SplitEdges.begin(),SplitEdges.end());
    //    std::vector<std::vector<int> > NextF,NextE;
    //    Geo::getFFAdjacency(Faces, NextF,NextE);

    std::vector<Geo::Point3<ScalarType> > NewPos=Pos;
    std::vector<std::vector<int> > NewFaces=Faces;

    std::vector<bool> CheckedV(Pos.size(),false);

    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[i].size();j++)
        {
            //if border do nothing
            if (NextF[i][j]==-1)continue;

            int IndexV0=Faces[i][j];
            //already checked

            int IndexV1=Faces[i][(j+1)%Faces[i].size()];
            std::pair<int,int> Entry0(IndexV0,IndexV1);
            std::pair<int,int> Entry1(IndexV1,IndexV0);
            bool HasSplit= ((SplitEdgeSet.count(Entry0)>0) ||
                            (SplitEdgeSet.count(Entry1)>0));

            //no split
            if (!HasSplit)continue;

            if (CheckedV[IndexV0])continue;
            CheckedV[IndexV0]=true;

            Pos3 StartP(i,j,IndexV0);
            std::vector<Pos3> PosFan;
            StartP.FlipE(Faces);
            bool HasBorder;

            GetSortedFacesFaceFun(StartP,Faces,NextF,NextE,PosFan,HasBorder);

            std::vector<std::vector<Pos3> > PosPartitions;
            PosPartitions.resize(1);
            for (size_t k=0;k<PosFan.size();k++)
            {
                PosPartitions.back().push_back(PosFan[k]);

                //do not check the last one because moight be back on first marked edge
                if ((!HasBorder)&&(k==PosFan.size()-1))continue;

                //check if has a border on other side
                Pos3 TestPos=PosFan[k];
                //TestPos.FlipE(Faces);
                int IndexVTest0=TestPos.V;
                TestPos.FlipV(Faces);
                int IndexVTest1=TestPos.V;
                std::pair<int,int> Entry0(IndexVTest0,IndexVTest1);
                std::pair<int,int> Entry1(IndexVTest1,IndexVTest0);
                bool HasSplit= ((SplitEdgeSet.count(Entry0)>0) ||
                                (SplitEdgeSet.count(Entry1)>0));

                if (HasSplit)
                    PosPartitions.resize(PosPartitions.size()+1);
            }

            //no need to split the vertex
            if (PosPartitions.size()==1)continue;

            Geo::Point3<ScalarType> VertPos=Pos[StartP.V];
            //first partition remain the same then from the secon to other change vertex
            for(size_t k=1;k<PosPartitions.size();k++)
            {
                NewPos.push_back(VertPos);
                int newVIndex=NewPos.size()-1;
                for(size_t w=0;w<PosPartitions[k].size();w++)
                {
                    int IndexF=PosPartitions[k][w].F;
                    int IndexV=PosPartitions[k][w].V;
                    int LocV=WhichVert(Faces[IndexF],IndexV);
                    assert(LocV>=0);
                    assert(LocV<=(int)Faces[IndexF].size());
                    NewFaces[IndexF][LocV]=newVIndex;
                }
            }

        }
    Pos=NewPos;
    Faces=NewFaces;
}


#endif
