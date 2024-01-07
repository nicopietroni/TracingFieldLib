#ifndef PATCH_MANAGING
#define PATCH_MANAGING

#include <vector>
#include "../point3.h"
#include "../mesh_functions.h"
#include "connectivity.h"
#include "graph_cross_field.h"
//#include "../math_utils.h"
#include <unordered_map>
#include <pos3.h>
#include "graph.h"
#include "color_functions.h"
#include "../IO/io_obj.h"
#include "mesh_refine.h"
#include <UV_param.h>
#include "../tangent_space_smooth.h"

namespace Geo{

template <class ScalarType>
class PatchManaging
{

public:
    enum CornerType{CTypeFlat,CTypeConvex,CTypeSpike,CTypeConcave};

    //per patch per loop corners
    std::vector<std::vector<std::vector<int> > > PatchCorners;
    std::vector<std::vector<std::vector<CornerType> > > TypeOfCorner;

    //for each patch, for each loop, for each edge the sequence
    std::vector<std::vector<std::vector<std::vector<std::pair<int,int> > > > > EdgePatchSides;
    //per patch per loop per side
    std::vector<std::vector<std::vector<ScalarType> > > SideLenghts;

    bool UpdateVisual;
    //ONLY FOR VISUALIZATION
    std::vector<std::vector<Geo::Point3<ScalarType> > > VisualCorners;
    std::vector<std::vector<CornerType> > VisualCornerType;

    //for each patch, for each side , the vector of points
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType> > > > VisualPaths;

    //SUBMESHES
    std::vector<std::vector<Geo::Point3<ScalarType> > > SubPatchVertPos;
    std::vector<std::vector<std::vector<int> > > SubPatchFaces;

    //their mapping to original one
    std::vector<std::vector<int > > SubPatchVertToOriginal;
    std::vector<std::vector<int > > SubPatchFacesToOriginal;
    //END SUBMESHES

    //the original mesh
    std::vector<Geo::Point3<ScalarType> > &VertPos;
    std::vector<std::vector<int> > &Faces;
    const Geo::TracingGraph<ScalarType> &TGraph;

    //mapping patch to faces and viceversa
    std::vector<std::vector<int> > PatchOriginalFaces;
    std::vector<int> OriginalFaceToPatch;

private:

    //adjacency on the mesh
    std::vector<std::vector<int> > NextF;
    std::vector<std::vector<int> > NextE;

    //the vactor with nodes composing the paths
    const std::vector<std::vector<int> > &PathNodes;

    //the old nodes used to check if need to update
    std::vector<std::vector<int> > OLDPathNodes;

    //for each edge the nodes that are passed with respect to the first node
    // std::unordered_map<std::pair<int,int>,int[2],Geo::hashPairIntFunction > PathsEdges;
    std::map<std::pair<int,int>,std::vector<int> > PathsEdges;

    //this one store for each edge the faces that are on it
    std::map<std::pair<int,int>,std::vector<int> > EdgesToFace;


    void InitVisualCorners(int IndexP)
    {
        VisualCorners[IndexP].clear();
        VisualCornerType[IndexP].clear();

        std::map<int,std::pair<Geo::Point3<ScalarType>,int> > CornersAvgPos;
        std::map<int,CornerType> PerCornerType;
        //initialize to cumulate per corner
        std::vector<int> CornerVert;
        std::vector<CornerType> CurrType;
        GetAllCorners(IndexP,CornerVert,CurrType);
        assert(CornerVert.size()==CurrType.size());
        for (size_t k=0;k<CornerVert.size();k++)
        {
            int IndexCornerV=CornerVert[k];
            CornersAvgPos[IndexCornerV]=std::pair<Geo::Point3<ScalarType>,int>(Geo::Point3<ScalarType>(0,0,0),0);
            PerCornerType[IndexCornerV]=CurrType[k];
        }

        std::vector<Geo::Point3<ScalarType> > BaryF;
        ComputeFaceBarycenters(SubPatchVertPos[IndexP],SubPatchFaces[IndexP],BaryF);

        //go over all faces
        for (size_t i=0;i<SubPatchFaces[IndexP].size();i++)
        {
            size_t sizeF=SubPatchFaces[IndexP][i].size();
            assert(sizeF==3);

            //sum up if there is the vertex
            for (size_t j=0;j<sizeF;j++)
            {
                int IndexV=SubPatchFaces[IndexP][i][j];
                if (CornersAvgPos.count(IndexV)==0)continue;
                CornersAvgPos[IndexV].first+=BaryF[i];
                CornersAvgPos[IndexV].second++;
            }
        }

        typename std::map<int,std::pair<Geo::Point3<ScalarType>,int> >::iterator IteMapCorn;
        for (IteMapCorn=CornersAvgPos.begin();IteMapCorn!=CornersAvgPos.end();IteMapCorn++)
        {
            Geo::Point3<ScalarType> AvgP=(*IteMapCorn).second.first;
            int SumP=(*IteMapCorn).second.second;
            Geo::Point3<ScalarType> PosCorn=AvgP/(ScalarType)SumP;
            VisualCorners[IndexP].push_back(PosCorn);
            int IndexV=(*IteMapCorn).first;
            assert(PerCornerType.count(IndexV)>0);
            VisualCornerType[IndexP].push_back(PerCornerType[IndexV]);
        }

    }

    void GetAllSideVert(int IndexP,std::vector<int> &SideVert)
    {
        SideVert.clear();
        for (size_t i=0;i<EdgePatchSides[IndexP].size();i++)
            for (size_t j=0;j<EdgePatchSides[IndexP][i].size();j++)
                for (size_t k=0;k<EdgePatchSides[IndexP][i][j].size();k++)
                {
                    int IndexV0=EdgePatchSides[IndexP][i][j][k].first;
                    int IndexV1=EdgePatchSides[IndexP][i][j][k].second;
                    SideVert.push_back(IndexV0);
                    SideVert.push_back(IndexV1);
                }

        // remove consecutive (adjacent) duplicates
        std::vector<int>::iterator last = std::unique(SideVert.begin(), SideVert.end());
        SideVert.erase(last, SideVert.end());
    }

    void InitVisualSides(int IndexP)
    {
        VisualPaths[IndexP].clear();

        std::map<int,std::pair<Geo::Point3<ScalarType>,int> > SidePatchPos;
        std::vector<int> SideVert;
        GetAllSideVert(IndexP,SideVert);
        for (size_t i=0;i<SideVert.size();i++)
        {
            int IndexV=SideVert[i];
            SidePatchPos[IndexV]=std::pair<Geo::Point3<ScalarType>,int>(Geo::Point3<ScalarType>(0,0,0),0);
        }

        std::vector<Geo::Point3<ScalarType> > BaryF;
        ComputeFaceBarycenters(SubPatchVertPos[IndexP],SubPatchFaces[IndexP],BaryF);

        //for each subface
        for (size_t i=0;i<SubPatchFaces[IndexP].size();i++)
        {
            size_t sizeF=SubPatchFaces[IndexP][i].size();
            assert(sizeF==3);

            for (size_t j=0;j<sizeF;j++)
            {
                int IndexV=SubPatchFaces[IndexP][i][j];
                if (SidePatchPos.count(IndexV)==0)continue;
                SidePatchPos[IndexV].first+=BaryF[i];
                SidePatchPos[IndexV].second++;
            }
        }

        //for each loop
        for (size_t i=0;i<EdgePatchSides[IndexP].size();i++)
            //for each side
            for (size_t j=0;j<EdgePatchSides[IndexP][i].size();j++)
            {
                //for each side
                std::vector<Geo::Point3<ScalarType> > SidePos;
                for (size_t k=0;k<EdgePatchSides[IndexP][i][j].size();k++)
                {
                    //get the vertex on the side
                    int IndexV0=EdgePatchSides[IndexP][i][j][k].first;
                    assert(SidePatchPos.count(IndexV0)>0);
                    //average the cumulated posistion
                    Geo::Point3<ScalarType> AvgPos=SidePatchPos[IndexV0].first;
                    int NumPos=SidePatchPos[IndexV0].second;
                    SidePos.push_back(AvgPos/(ScalarType)NumPos);

                    //add the last one
                    if (k==EdgePatchSides[IndexP][i][j].size()-1)
                    {
                        int IndexV1=EdgePatchSides[IndexP][i][j][k].second;
                        assert(SidePatchPos.count(IndexV1)>0);
                        Geo::Point3<ScalarType> AvgPos=SidePatchPos[IndexV1].first;
                        int NumPos=SidePatchPos[IndexV1].second;
                        SidePos.push_back(AvgPos/(ScalarType)NumPos);
                    }
                }
                VisualPaths[IndexP].push_back(SidePos);
            }
    }

    ScalarType ComputeSideLenght(int IndexP,int IndexL,int IndexS,int smooth_steps=3)
    {
        assert(IndexP<EdgePatchSides.size());
        assert(IndexL<EdgePatchSides[IndexP].size());
        assert(IndexS<EdgePatchSides[IndexP][IndexL].size());
        if (EdgePatchSides[IndexP][IndexL][IndexS].size()==0)return 0;

        std::vector<Geo::Point3<ScalarType> > EdgePos;
        for (size_t i=0;i<EdgePatchSides[IndexP][IndexL][IndexS].size();i++)
        {
            int IndexV= EdgePatchSides[IndexP][IndexL][IndexS][i].first;
            assert(IndexV>=0);
            assert(IndexV<SubPatchVertPos[IndexP].size());
            EdgePos.push_back(SubPatchVertPos[IndexP][IndexV]);
        }

        //then add the last one
        int IndexVLast= EdgePatchSides[IndexP][IndexL][IndexS].back().second;
        assert(IndexVLast>=0);
        assert(IndexVLast<SubPatchVertPos[IndexP].size());
        EdgePos.push_back(SubPatchVertPos[IndexP][IndexVLast]);

        //then smooth a bit
        for (size_t s=0;s<smooth_steps;s++)
        {
            std::vector<Geo::Point3<ScalarType> > EdgePosNew=EdgePos;
            for (int i=1;i<(int)(EdgePos.size()-1);i++)
            {
                Geo::Point3<ScalarType> P0=EdgePos[i-1];
                Geo::Point3<ScalarType> P1=EdgePos[i];
                Geo::Point3<ScalarType> P2=EdgePos[i+1];
                Geo::Point3<ScalarType> Avg=(P0+P1+P2)/3;
                EdgePosNew[i]=Avg;
            }
            EdgePos=EdgePosNew;
        }

        //finally compute lenghts
        ScalarType LenghtP=0;
        for (int i=0;i<(int)(EdgePos.size()-1);i++)
        {
            Geo::Point3<ScalarType> P0=EdgePos[i];
            Geo::Point3<ScalarType> P1=EdgePos[i+1];
            LenghtP+=(P1-P0).Norm();
        }
        return LenghtP;
    }

    void UpdateSideLenghs(int IndexP)
    {
        assert(IndexP<SideLenghts.size());
        SideLenghts[IndexP].resize(EdgePatchSides[IndexP].size());

        for (size_t i=0;i<EdgePatchSides[IndexP].size();i++)
        {
            for (size_t j=0;j<EdgePatchSides[IndexP][i].size();j++)
            {
                ScalarType SideL=ComputeSideLenght(IndexP,i,j);
                SideLenghts[IndexP][i].push_back(SideL);
            }
        }
    }

    void UpdateSideLenghs(const std::vector<int> &PatchesIdx)
    {
        for (size_t i=0;i<PatchesIdx.size();i++)
        {
            UpdateSideLenghs(PatchesIdx[i]);
        }
    }

    void InitVisualCornersEdges(const std::vector<int> &PatchesIdx)
    {
        for (size_t i=0;i<PatchesIdx.size();i++)
        {
            InitVisualCorners(PatchesIdx[i]);
            InitVisualSides(PatchesIdx[i]);
        }
    }

    void GetTangentNodesOn(const int &IndexV,
                           std::vector<int> &TangNodeIdx)
    {
        TangNodeIdx.clear();
        std::vector<int> NodeIdx;
        Field::CrossFieldGraph::NodesOfElem(IndexV,NodeIdx);
        for (size_t i=0;i<NodeIdx.size();i++)
            if (TGraph.NodeType(NodeIdx[i])==TangentDomain)
                TangNodeIdx.push_back(NodeIdx[i]);
    }

    void GetPairTangentNodes(const int &IndexV0,
                             const int &IndexV1,
                             int &IndexN0,
                             int &IndexN1)
    {

        std::vector<int> TangNodeIdx0;
        GetTangentNodesOn(IndexV0,TangNodeIdx0);

        std::vector<int> TangNodeIdx1;
        GetTangentNodesOn(IndexV1,TangNodeIdx1);

        assert(TangNodeIdx0.size()==2);
        assert(TangNodeIdx1.size()==2);

        //check with one only border
        IndexN0=-1;
        IndexN1=-1;
        ScalarType MinW=std::numeric_limits<ScalarType>::max();

        for (size_t i=0;i<TangNodeIdx0.size();i++)
        {
            int Test_IndexN0=TangNodeIdx0[i];
            for (size_t j=0;j<TangNodeIdx1.size();j++)
            {
                int Test_IndexN1=TangNodeIdx1[j];
                ScalarType ConnW=TGraph.WeightConnectection(Test_IndexN0,Test_IndexN1,false);
                if (ConnW<0)continue;
                if (ConnW>MinW)continue;
                MinW=ConnW;
                IndexN0=Test_IndexN0;
                IndexN1=Test_IndexN1;
            }
        }
        assert(IndexN0>=0);
        assert(IndexN1>=0);
    }

    void GetNodeOfType(const int IndexV,const Geo::NType TypeNode,
                       std::vector<size_t> &IndexN)
    {
        IndexN.clear();
        std::vector<int> Elems;
        Field::CrossFieldGraph::NodesOfElem(IndexV,Elems);

        for (size_t k=0;k<Elems.size();k++)
        {
            if (TGraph.NodeType(Elems[k])!=TypeNode)continue;
            IndexN.push_back(Elems[k]);
        }
    }

    void AddBorderPaths()
    {
        //get all border pos
        std::vector<std::vector<Pos3 > > BorderPos;
        getBorderSequnces(Faces,NextF,NextE,BorderPos);

        for (size_t i=0;i<BorderPos.size();i++)
        {
            for (size_t j=0;j<BorderPos[i].size();j++)
            {
                Pos3 CurrPos=BorderPos[i][j];
                int IndexV0=CurrPos.V;
                int IndexV1=CurrPos.VFlip(Faces);

                std::vector<size_t> TangentV0;
                GetNodeOfType(IndexV0,Geo::TangentDomain,TangentV0);

                std::vector<size_t> TangentV1;
                GetNodeOfType(IndexV1,Geo::TangentDomain,TangentV1);

                std::pair<int,int> key0(IndexV0,IndexV1);
                std::pair<int,int> key1(IndexV1,IndexV0);

                for (size_t k=0;k<TangentV0.size();k++)
                    PathsEdges[key0].push_back(TangentV0[k]);

                for (size_t k=0;k<TangentV1.size();k++)
                    PathsEdges[key1].push_back(TangentV1[k]);
            }
        }

    }

    void UpdateFaceToPatch()
    {
        //then update the correspondence
        OriginalFaceToPatch.clear();
        OriginalFaceToPatch.resize(Faces.size(),-1);

        for (size_t i=0;i<PatchOriginalFaces.size();i++)
            for (size_t j=0;j<PatchOriginalFaces[i].size();j++)
            {
                int IndexF=PatchOriginalFaces[i][j];
                assert(IndexF>=0);
                assert(IndexF<OriginalFaceToPatch.size());
                OriginalFaceToPatch[IndexF]=i;
            }
    }


    struct VertDirs
    {
        int IndexV;
        int Dir[2];

        VertDirs()
        {
            IndexV=-1;
            Dir[0]=-1;
            Dir[1]=-1;
        }
    };


    CornerType GetCornerType(int CountDir[4],const ScalarType &angleDeg)
    {
        std::set<int> DirM2;
        for (size_t k=0;k<4;k++)
        {
            if (CountDir[k]>0)
                DirM2.insert(k%2);
        }

        if (DirM2.size()==2)
        {
            if (angleDeg<180)
                return CTypeConvex;
            else
                return CTypeConcave;
        }

        if ((CountDir[0]==2)&&(CountDir[2]==0))
            return CTypeSpike;
        if ((CountDir[2]==2)&&(CountDir[0]==0))
            return CTypeSpike;
        if ((CountDir[1]==2)&&(CountDir[3]==0))
            return CTypeSpike;
        if ((CountDir[3]==2)&&(CountDir[1]==0))
            return CTypeSpike;
        //        if ((CountDir[1]==2)&&(CountDir[3]==2))
        //             return CTypeFlat;


        //        return CTypeSpike;
        return CTypeFlat;
    }

    void FindPatchCorners(const int &IndexP)
    {

        if (SubPatchFaces[IndexP].size()==0)return;
        //get the Adjacency of the submesh
        std::vector<std::vector<int> > NextFPatch;
        std::vector<std::vector<int> > NextEPatch;
        getFFAdjacency(SubPatchFaces[IndexP],NextFPatch,NextEPatch);

        //compute the angle on the vertices to cech convex/concave
        std::vector<ScalarType> VertAngle;
        ComputePerVertAngle(SubPatchVertPos[IndexP],SubPatchFaces[IndexP],VertAngle);

        //get the border sequences
        std::vector<std::vector<int > > BorderVIdx;
        getBorderSequncesVert(SubPatchFaces[IndexP],NextFPatch,NextEPatch,BorderVIdx);

        std::set<int> CornerSet;
        assert(IndexP<PatchCorners.size());
        assert(IndexP<TypeOfCorner.size());
        PatchCorners[IndexP].clear();
        TypeOfCorner[IndexP].clear();
        PatchCorners[IndexP].resize(BorderVIdx.size());
        TypeOfCorner[IndexP].resize(BorderVIdx.size());
        //for each loop
        for (size_t i=0;i<BorderVIdx.size();i++)
        {

            int sizeL=BorderVIdx[i].size();

            for (size_t j=0;j<BorderVIdx[i].size();j++)
            {
                int IndexV_Prev=BorderVIdx[i][(j+sizeL-1)%sizeL];
                int IndexV_Curr=BorderVIdx[i][j];
                int IndexV_Next=BorderVIdx[i][(j+1)%sizeL];

                //then map to original
                int IndexV_Prev_Or=MapToOriginal(IndexP,IndexV_Prev);
                int IndexV_Curr_Or=MapToOriginal(IndexP,IndexV_Curr);
                int IndexV_Next_Or=MapToOriginal(IndexP,IndexV_Next);

                //get the two edges as emanating from the current node
                std::pair<int,int> key_prev(IndexV_Curr_Or,IndexV_Prev_Or);
                std::pair<int,int> key_next(IndexV_Curr_Or,IndexV_Next_Or);

                //check they shoudl exhist
                assert(PathsEdges.count(key_prev)>0);
                assert(PathsEdges.count(key_next)>0);

                //get the direction from nodes
                std::vector<int> NodesV_prev=PathsEdges[key_prev];
                std::vector<int> NodesV_next=PathsEdges[key_next];

                //add direction per vertex
                std::vector<int> DirV_prev;
                Field::CrossFieldGraph::DirOfNodes(NodesV_prev,DirV_prev);
                std::vector<int> DirV_next;
                Field::CrossFieldGraph::DirOfNodes(NodesV_next,DirV_next);
                std::vector<int> AllDirsV=DirV_prev;
                AllDirsV.insert(AllDirsV.end(),DirV_next.begin(),DirV_next.end());
                int CountDir[4]={0,0,0,0};
                for (size_t k=0;k<AllDirsV.size();k++)
                {
                    int IndexD=AllDirsV[k];
                    CountDir[IndexD]++;
                }

                //then check the kind of vertex
                ScalarType currAVert=VertAngle[IndexV_Curr];
                CornerType CType=GetCornerType(CountDir,currAVert);
                if (CType!=CTypeFlat)
                {
                    PatchCorners[IndexP][i].push_back(IndexV_Curr);
                    TypeOfCorner[IndexP][i].push_back(CType);
                    CornerSet.insert(IndexV_Curr);
                }
            }
        }

        //then sort with respect to
        //std::cout<<"CornerSet size:"<<CornerSet.size()<<std::endl;
        assert(IndexP<EdgePatchSides.size());
        EdgePatchSides[IndexP].resize(BorderVIdx.size());
        for (size_t i=0;i<BorderVIdx.size();i++)
        {
            int IndexCorner=-1;
            for (int j=0;j<(int)BorderVIdx[i].size();j++)
            {
                if (CornerSet.count(BorderVIdx[i][j])>0)
                {
                    IndexCorner=j;
                    assert(CornerSet.count(BorderVIdx[i][IndexCorner])>0);
                    break;
                }
            }
            if (IndexCorner!=-1)
            {
                std::rotate(BorderVIdx[i].begin(), BorderVIdx[i].begin() + IndexCorner, BorderVIdx[i].end());
                assert(CornerSet.count(BorderVIdx[i][0])>0);
            }

            //retrieve the side vertices
            std::vector<std::vector<std::pair<int,int> > > SideEdges;
            SideEdges.resize(1);

            //makes no sense a border with lenght 1
            //retrieve edge sequences
            assert(BorderVIdx[i].size()>1);
            int BorderL=(int)BorderVIdx[i].size();
            for (int j=0;j<BorderL;j++)
            {
                int IndexV0=BorderVIdx[i][j];
                int IndexV1=BorderVIdx[i][(j+1)%BorderL];
                std::pair<int,int> EdgeV(IndexV0,IndexV1);
                SideEdges.back().push_back(EdgeV);
                if ((CornerSet.count(IndexV1)>0)&&(j!=(BorderL-1)))
                    SideEdges.resize(SideEdges.size()+1);
            }
            EdgePatchSides[IndexP][i]=SideEdges;
        }

    }

//    void FindAddRemovePaths(std::vector<int> &ToAdd,
//                            std::vector<int> &ToRemove)
//    {
//        ToAdd.clear();
//        ToRemove.clear();
//        for (size_t i=0;i<PathNodes.size();i++)
//        {
//            bool found=false;
//            for (size_t j=0;j<OLDPathNodes.size();j++)
//            {
//                if (PathNodes[i]!=OLDPathNodes[j])continue;
//                found=true;
//                break;
//            }
//            if (!found)
//                ToAdd.push_back(i);
//        }
//        for (size_t i=0;i<OLDPathNodes.size();i++)
//        {
//            bool found=false;
//            for (size_t j=0;j<PathNodes.size();j++)
//            {
//                if (OLDPathNodes[i]!=(PathNodes[j]))continue;
//                found=true;
//                break;
//            }
//            if (!found)
//                ToRemove.push_back(i);
//        }
//    }

    void InitEdgeToFace()
    {
        EdgesToFace.clear();
        for (size_t i=0;i<Faces.size();i++)
            for (size_t j=0;j<Faces[i].size();j++)
            {
                int IndexV0=Faces[i][j];
                int IndexV1=Faces[i][(j+1)%3];
                std::pair<int,int> key(std::min(IndexV0,IndexV1),
                                       std::max(IndexV0,IndexV1));
                EdgesToFace[key].push_back(i);
            }
    }

    void RemovePath(const std::vector<int> &ToRem)
    {
        if (ToRem.size()==0)return;
        //loop repeat the last two times so no need to check
        for (int i=0;i<ToRem.size()-1;i++)
        {
            int IndexN0=ToRem[i];
            int IndexN1=ToRem[i+1];
            int IndexV0=Field::CrossFieldGraph::ElemOfNode(IndexN0);
            int IndexV1=Field::CrossFieldGraph::ElemOfNode(IndexN1);
            std::pair<int,int> Key0(IndexV0,IndexV1);
            std::pair<int,int> Key1(IndexV1,IndexV0);
            assert(PathsEdges.count(Key0)>0);
            assert(PathsEdges.count(Key1)>0);
            PathsEdges.erase(Key0);
            PathsEdges.erase(Key1);
        }
    }

    void AddPath(const std::vector<int> &ToAdd)
    {
        //loop repeat the last two times so no need to check
        if (ToAdd.size()==0)return;

        for (int i=0;i<ToAdd.size()-1;i++)
        {
            int IndexN0=ToAdd[i];
            int IndexN1=ToAdd[i+1];
            //int OppN0=Field::CrossFieldGraph::OppositeNodeCrossField(IndexN0);
            int OppN1=Field::CrossFieldGraph::OppositeNodeCrossField(IndexN1);

            //safety check
            int IndexV0=Field::CrossFieldGraph::ElemOfNode(IndexN0);
            int IndexV1=Field::CrossFieldGraph::ElemOfNode(IndexN1);
            std::pair<int,int> Key0(IndexV0,IndexV1);
            std::pair<int,int> Key1(IndexV1,IndexV0);

            assert(PathsEdges.count(Key0)==0);
            assert(PathsEdges.count(Key1)==0);

            PathsEdges[Key0].push_back(IndexN0);
            PathsEdges[Key1].push_back(OppN1);

        }
    }

    void FFConnectedComponenentElements(const std::vector<int> &IndexF,
                                        std::vector<std::vector<int> > &ComponentsF)
    {
        ComponentsF.clear();

        //find partitions
        std::set<int> Explored;

        for (size_t i=0;i<IndexF.size();i++)
        {
            if (Explored.count(IndexF[i])>0)continue;

            ComponentsF.resize(ComponentsF.size()+1);

            std::vector<size_t> StackF;
            StackF.push_back(IndexF[i]);
            do
            {
                int CurrF=StackF.back();
                StackF.pop_back();

                if (Explored.count(CurrF)>0)continue;

                Explored.insert(CurrF);

                ComponentsF.back().push_back(CurrF);


                for (size_t n=0;n<NextF[CurrF].size();n++)
                {
                    //get next elem
                    int AdjF=NextF[CurrF][n];
                    assert(AdjF!=CurrF);

                    //get edge as seen from the other face
                    int AdjE=NextE[CurrF][n];

                    if (AdjF==-1)continue;
                    assert(AdjE!=-1);

                    //check if is an edge
                    int IndexV0=Faces[CurrF][n];
                    int IndexV1=Faces[CurrF][(n+1)%3];

                    std::pair<int,int> Key(IndexV0,IndexV1);

                    //there is a path
                    if (PathsEdges.count(Key)>0)continue;

                    bool has_to_diffuse=(Explored.count(AdjF)==0);
                    if (has_to_diffuse)
                    {
                        StackF.push_back(AdjF);
                    }
                }
            }while (StackF.size()>0);
        }
    }

    void UpdatePatches(const std::vector<int> &ToUpdPatches)
    {
        std::vector<int> FacesToCheck;
        for (size_t i=0;i<ToUpdPatches.size();i++)
        {
            int IndexP=ToUpdPatches[i];
            assert(IndexP>=0);
            assert(IndexP<PatchOriginalFaces.size());
            FacesToCheck.insert(FacesToCheck.end(),
                                PatchOriginalFaces[IndexP].begin(),
                                PatchOriginalFaces[IndexP].end());
        }
        std::vector<std::vector<int> > ComponentsF;
        FFConnectedComponenentElements(FacesToCheck,ComponentsF);

        //clear the old one
        for (size_t i=0;i<ToUpdPatches.size();i++)
        {
            int IndexP=ToUpdPatches[i];
            ClearDataForPatch(IndexP);

        }

        //then add the new ones
        for (size_t i=0;i<ComponentsF.size();i++)
        {
            AllocateNewPatch();
            for (size_t j=0;j<ComponentsF[i].size();j++)
            {
                int IndexF=ComponentsF[i][j];
                OriginalFaceToPatch[IndexF]=PatchOriginalFaces.size()-1;
                PatchOriginalFaces.back().push_back(ComponentsF[i][j]);
            }
        }
    }

    void SplitSubMeshWithInternalEdges(int IndexPatch)
    {
        std::vector<std::vector<int> > NextFPatch,NextEPatch;
        getFFAdjacency(SubPatchFaces[IndexPatch],NextFPatch,NextEPatch);
        std::set<std::pair<int,int> > SplitEdgeSet;
        for (size_t i=0;i<SubPatchFaces[IndexPatch].size();i++)
        {
            int sizeF=SubPatchFaces[IndexPatch][i].size();
            assert(sizeF==3);
            for (size_t j=0;j<3;j++)
            {
                if (NextFPatch[i][j]==-1)continue;

                //get the two vertices and maps to the original
                int IndexV0=SubPatchFaces[IndexPatch][i][j];
                assert(IndexV0<SubPatchVertToOriginal[IndexPatch].size());
                int IndexV0Or=SubPatchVertToOriginal[IndexPatch][IndexV0];

                int IndexV1=SubPatchFaces[IndexPatch][i][(j+1)%3];
                assert(IndexV1<SubPatchVertToOriginal[IndexPatch].size());
                int IndexV1Or=SubPatchVertToOriginal[IndexPatch][IndexV1];

                std::pair<int,int> EdgeKOr(std::min(IndexV0Or,IndexV1Or),std::max(IndexV0Or,IndexV1Or));
                if (PathsEdges.count(EdgeKOr)==0)continue;

                std::pair<int,int> EdgeK(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                SplitEdgeSet.insert(EdgeK);
            }
        }
        std::vector<std::pair<int,int> > SplitEdges(SplitEdgeSet.begin(),SplitEdgeSet.end());

        if (SplitEdges.size()>0)
        {
            //save the original vertex mapping
            std::map<Geo::Point3<ScalarType>, int> PosToIDX;
            for (size_t i=0;i<SubPatchVertPos[IndexPatch].size();i++)
            {
                Geo::Point3<ScalarType> Pos=SubPatchVertPos[IndexPatch][i];
                PosToIDX[Pos]=SubPatchVertToOriginal[IndexPatch][i];
            }

            //split the mesh
            SplitMeshFromEdges(SubPatchVertPos[IndexPatch],SubPatchFaces[IndexPatch],NextFPatch,NextEPatch,SplitEdges);

            //restore the right mapping
            SubPatchVertToOriginal[IndexPatch].clear();
            for (size_t i=0;i<SubPatchVertPos[IndexPatch].size();i++)
            {
                Geo::Point3<ScalarType> Pos=SubPatchVertPos[IndexPatch][i];
                assert(PosToIDX.count(Pos)>0);
                SubPatchVertToOriginal[IndexPatch].push_back(PosToIDX[Pos]);
            }
        }
    }

    void UpdateSubMesh(int IndexPatch)
    {
        assert(IndexPatch>=0);
        std::map<int,int> VertMap;

        assert(IndexPatch<(int)SubPatchVertPos.size());
        assert(IndexPatch<(int)SubPatchFaces.size());
        assert(IndexPatch<(int)SubPatchVertToOriginal.size());
        assert(IndexPatch<(int)SubPatchFacesToOriginal.size());

        SubPatchVertPos[IndexPatch].clear();
        SubPatchFaces[IndexPatch].clear();
        SubPatchVertToOriginal[IndexPatch].clear();
        SubPatchFacesToOriginal[IndexPatch].clear();

        //for each face
        for (size_t i=0;i<PatchOriginalFaces[IndexPatch].size();i++)
        {
            int IndexOrFace=PatchOriginalFaces[IndexPatch][i];

            assert(Faces[IndexOrFace].size()==3);
            SubPatchFaces[IndexPatch].push_back(std::vector<int>(3,-1));
            SubPatchFacesToOriginal[IndexPatch].push_back(IndexOrFace);

            //for each vertex
            for (size_t j=0;j<3;j++)
            {
                int IndexOrV=Faces[IndexOrFace][j];
                int IndexNewV=-1;

                if (VertMap.count(IndexOrV)>0)
                    IndexNewV=VertMap[IndexOrV];
                else
                {
                    SubPatchVertPos[IndexPatch].push_back(VertPos[IndexOrV]);
                    SubPatchVertToOriginal[IndexPatch].push_back(IndexOrV);
                    IndexNewV=SubPatchVertPos[IndexPatch].size()-1;

                    VertMap[IndexOrV]=IndexNewV;
                }
                SubPatchFaces[IndexPatch].back()[j]=IndexNewV;
            }
        }

        //then split internally if needed
        SplitSubMeshWithInternalEdges(IndexPatch);
    }

    void UpdateSubMeshes(const std::vector<int> &IndexPatches)
    {
        for (size_t i=0;i<IndexPatches.size();i++)
            UpdateSubMesh(IndexPatches[i]);
    }

    void FindPatchCorners(const std::vector<int> &IndexPatches)
    {
        for (size_t i=0;i<IndexPatches.size();i++)
            FindPatchCorners(IndexPatches[i]);
    }

    void ClearDataForPatch(const int &IndexPatch)
    {
        assert(IndexPatch>=0);
        assert(IndexPatch<(int)SideLenghts.size());
        assert(IndexPatch<(int)PatchOriginalFaces.size());
        assert(IndexPatch<(int)SubPatchVertPos.size());
        assert(IndexPatch<(int)SubPatchFaces.size());
        assert(IndexPatch<(int)SubPatchVertToOriginal.size());
        assert(IndexPatch<(int)SubPatchFacesToOriginal.size());
        assert(IndexPatch<(int)EdgePatchSides.size());

        SideLenghts[IndexPatch].clear();
        PatchOriginalFaces[IndexPatch].clear();
        SubPatchVertPos[IndexPatch].clear();
        SubPatchFaces[IndexPatch].clear();
        SubPatchVertToOriginal[IndexPatch].clear();
        SubPatchFacesToOriginal[IndexPatch].clear();
        PatchCorners[IndexPatch].clear();
        TypeOfCorner[IndexPatch].clear();
        EdgePatchSides[IndexPatch].clear();
        VisualCorners[IndexPatch].clear();
        VisualCornerType[IndexPatch].clear();
        VisualPaths[IndexPatch].clear();
    }

    void AllocateNewPatch()
    {
        SideLenghts.resize(SideLenghts.size()+1);
        PatchOriginalFaces.resize(PatchOriginalFaces.size()+1);
        SubPatchVertPos.resize(SubPatchVertPos.size()+1);
        SubPatchFaces.resize(SubPatchFaces.size()+1);
        SubPatchVertToOriginal.resize(SubPatchVertToOriginal.size()+1);
        SubPatchFacesToOriginal.resize(SubPatchFacesToOriginal.size()+1);
        PatchCorners.resize(PatchCorners.size()+1);
        TypeOfCorner.resize(TypeOfCorner.size()+1);
        EdgePatchSides.resize(EdgePatchSides.size()+1);
        VisualCorners.resize(VisualCorners.size()+1);
        VisualCornerType.resize(VisualCornerType.size()+1);
        VisualPaths.resize(VisualPaths.size()+1);
    }

    void UpdateSubPatchPos(int IndexPatch)
    {
        assert(IndexPatch>=0);
        assert(IndexPatch<SubPatchVertToOriginal.size());
        for (size_t i=0;i<SubPatchVertToOriginal[IndexPatch].size();i++)
        {
            int IndexVOr=SubPatchVertToOriginal[IndexPatch][i];
            SubPatchVertPos[IndexPatch][i]=VertPos[IndexVOr];
        }
    }

    //redo with saved data
    void GetOriginalBoundaryVert(const int IndexPatch,
                                 std::vector<std::vector<int> > &BoundVIdx)
    {
        assert(IndexPatch<NumPatches());
        std::vector<std::vector<int> > NextF,NextE;
        getFFAdjacency(SubPatchFaces[IndexPatch],NextF,NextE);
        getBorderSequncesVert(SubPatchFaces[IndexPatch],NextF,NextE,BoundVIdx);
        MapToOriginal(IndexPatch,BoundVIdx);
    }

    //reorganise the sequence from the first corner found
    bool SortFromCorner(std::vector<std::vector<int> > &BoundVIdx,
                        const std::set<int> &CornerSet)
    {

        //sort with respect to one vertex
        for (size_t i=0;i<BoundVIdx.size();i++)
        {
            bool found=false;
            for (size_t j=0;j<BoundVIdx[i].size();j++)
            {
                int CurrV=BoundVIdx[i][j];
                if (CornerSet.count(CurrV)>0)
                {
                    std::rotate(BoundVIdx[i].begin(), BoundVIdx[i].begin() + j, BoundVIdx[i].end());
                    found=true;
                    break;
                }
            }
            if (!found)return false;
        }
        return true;
    }

    //return sides considering global indexes
    bool GetGlobalSubsides(const int IndexPatch,const std::set<int> &CornerSet,
                           std::vector<std::vector<std::vector<int> > > &SidesVert)
    {
        //get boundary
        std::vector<std::vector<int> > BoundVIdx;
        GetOriginalBoundaryVert(IndexPatch,BoundVIdx);

        //sort with respect to corner
        bool has_sorted=SortFromCorner(BoundVIdx,CornerSet);
        if (!has_sorted)return false;

        SidesVert.clear();
        SidesVert.resize(BoundVIdx.size());
        for (size_t i=0;i<BoundVIdx.size();i++)
        {
            if (BoundVIdx[i].size()==0)return false;

            int IndexV0=BoundVIdx[i][0];
            assert(CornerSet.count(IndexV0)>0);

            //allocate first side
            SidesVert[i].resize(1);
            //add first vertex
            SidesVert[i].back().push_back(IndexV0);
            //start from second on
            for (size_t j=1;j<BoundVIdx[i].size();j++)
            {
                //add it
                int IndexV=BoundVIdx[i][j];
                SidesVert[i].back().push_back(IndexV);

                //then check if need to allocate a new edge
                if (CornerSet.count(IndexV)>0)
                {
                    //then add a new side
                    SidesVert[i].resize(SidesVert[i].size()+1);
                    //add this as the first of the new edge
                    SidesVert[i].back().push_back(IndexV);
                }
            }
            //add the first one
            SidesVert[i].back().push_back(IndexV0);
        }
        return true;
    }

    //END FUNCTION FOR SUBDIVIDE
    std::map<int,int> SingValence;

public:

    void FindAddRemovePaths(std::vector<int> &ToAdd,
                            std::vector<int> &ToRemove)
    {
        ToAdd.clear();
        ToRemove.clear();
        for (size_t i=0;i<PathNodes.size();i++)
        {
            bool found=false;
            for (size_t j=0;j<OLDPathNodes.size();j++)
            {
                if (PathNodes[i]!=OLDPathNodes[j])continue;
                found=true;
                break;
            }
            if (!found)
                ToAdd.push_back(i);
        }
        for (size_t i=0;i<OLDPathNodes.size();i++)
        {
            bool found=false;
            for (size_t j=0;j<PathNodes.size();j++)
            {
                if (OLDPathNodes[i]!=(PathNodes[j]))continue;
                found=true;
                break;
            }
            if (!found)
                ToRemove.push_back(i);
        }
    }

    int MapFaceToOriginal(const int &IndexPatch,const int &FIndexLocal)const
    {
        assert(IndexPatch<SubPatchFacesToOriginal.size());
        assert(FIndexLocal<SubPatchFacesToOriginal[IndexPatch].size());
        return (SubPatchFacesToOriginal[IndexPatch][FIndexLocal]);
    }


    int MapToOriginal(const int &IndexPatch,const int &VIndexLocal)const
    {
        assert(IndexPatch<SubPatchVertToOriginal.size());
        assert(VIndexLocal<SubPatchVertToOriginal[IndexPatch].size());
        return (SubPatchVertToOriginal[IndexPatch][VIndexLocal]);
    }

    void MapToOriginal(const int &IndexPatch,std::vector<std::vector<int > > &VIndexes)const
    {
        assert(IndexPatch<SubPatchVertToOriginal.size());
        for (size_t i=0;i<VIndexes.size();i++)
            for (size_t j=0;j<VIndexes[i].size();j++)
            {
                int IndexLocalV=VIndexes[i][j];
                VIndexes[i][j]=MapToOriginal(IndexPatch,IndexLocalV);
            }
    }

    void MapToOriginal(const int &IndexPatch,std::vector<int > &VIndexes)const
    {
        assert(IndexPatch<SubPatchVertToOriginal.size());
        for (size_t i=0;i<VIndexes.size();i++)
        {
            int IndexLocalV=VIndexes[i];
            VIndexes[i]=MapToOriginal(IndexPatch,IndexLocalV);
        }
    }

    //return the midpoint for each side
    bool GetSideMidpoints(const int IndexPatch,
                          const std::set<int> &CornerSet,
                          std::vector<std::pair<int,int> > &CornersExtremes,
                          std::vector<std::pair<int,int> > &EdgeInt,
                          std::vector<ScalarType> &Alpha)
    {
        CornersExtremes.clear();
        EdgeInt.clear();
        Alpha.clear();
        //MidPoint.clear();

        //get subsides
        std::vector<std::vector<std::vector<int> > > SidesVert;
        bool isOK=GetGlobalSubsides(IndexPatch,CornerSet,SidesVert);
        if (!isOK)return false;
        for (size_t i=0;i<SidesVert.size();i++)
        {
            for (size_t j=0;j<SidesVert[i].size();j++)
            {
                assert(SidesVert[i][j].size()>=2);
                int indexAvg0=floor(SidesVert[i][j].size()/2)-1;
                indexAvg0=std::max(indexAvg0,0);
                int indexAvg1=indexAvg0+1;
                assert(indexAvg0>=0);
                assert(indexAvg0<SidesVert[i][j].size());
                assert(indexAvg1>=0);
                assert(indexAvg1<SidesVert[i][j].size());
                int IndexV0=(SidesVert[i][j][indexAvg0]);
                int IndexV1=(SidesVert[i][j][indexAvg1]);
                if (IndexV0>IndexV1)
                    std::swap(IndexV0,IndexV1);

                EdgeInt.push_back(std::pair<int,int>(IndexV0,IndexV1));
                Alpha.push_back(0.5);
                //                EdgeInt EInt(IndexV0,IndexV1,0.5);
                //                MidPoint.push_back(EInt);
                std::pair<int,int> SideE(SidesVert[i][j][0],SidesVert[i][j].back());
                CornersExtremes.push_back(SideE);
            }
        }
        return isOK;
    }

    void GetAllCorners(int IndexP,std::vector<int> &CornerVert,std::vector<CornerType> &CurrType)const
    {
        CornerVert.clear();
        CurrType.clear();
        std::map<int,std::pair<Geo::Point3<ScalarType>,int> > CornersPos;
        //for each loop
        for (size_t i=0;i<PatchCorners[IndexP].size();i++)
            //for each side
            for (size_t j=0;j<PatchCorners[IndexP].size();j++)
                //for each corner
                for (size_t k=0;k<PatchCorners[IndexP][j].size();k++)
                {
                    int IndexCornerV=PatchCorners[IndexP][j][k];
                    CornerVert.push_back(IndexCornerV);
                    CurrType.push_back(TypeOfCorner[IndexP][j][k]);
                }
    }


    void GetAllCorners(std::vector<int> &CornerVert)const
    {
        CornerVert.clear();
        for (size_t i=0;i<NumPatches();i++)
        {
            std::vector<int> PCornerVert;
            std::vector<CornerType> PCurrType;
            GetAllCorners(i,PCornerVert,PCurrType);
            MapToOriginal(i,PCornerVert);
            CornerVert.insert(CornerVert.end(),
                              PCornerVert.begin(),
                              PCornerVert.end());
        }
        std::sort(CornerVert.begin(),
                  CornerVert.end());
        std::vector<int>::iterator last = std::unique(CornerVert.begin(),
                                                      CornerVert.end());
        CornerVert.erase(last, CornerVert.end());
    }

    void GetAllCornerPos(std::vector<Geo::Point3<ScalarType >> &CornerPos)const
    {
        CornerPos.clear();
        std::vector<int> CornerVert;
        GetAllCorners(CornerVert);
        for (size_t i=0;i<CornerVert.size();i++)
        {
            int IndexV=CornerVert[i];
            CornerPos.push_back(VertPos[IndexV]);
        }
    }

    //    std::vector<Geo::Point3<ScalarType> > SidesMidP;
    //    bool ExtractDual()
    //    {
    //        CompactEmptyPatches();

    //        SidesMidP.clear();

    //        //first get all corners
    //        std::set<int> CornerSet;
    //        for (int i=0;i<NumPatches();i++)
    //        {
    //            std::vector<int> CornerVert;
    //            std::vector<CornerType> CurrType;
    //            GetAllCorners(i,CornerVert,CurrType);
    //            MapToOriginal(i,CornerVert);
    //            CornerSet.insert(CornerVert.begin(),CornerVert.end());
    //        }

    //        //compute mid point of edge
    //        std::set<std::pair<int,int> > AddedMidP;

    //        //std::set<EdgeInt> InterPEdge;

    //        std::map<std::pair<int,int> , ScalarType> EdgeAlpha;

    //        assert(SubPatchFaces.size()==SubPatchVertPos.size());
    //        for (size_t i=0;i<SubPatchFaces.size();i++)
    //        {
    //            //avoid empty patches
    //            if (SubPatchFaces[i].size()==0)continue;

    //            std::vector<std::pair<int,int> > CornersExtremes;
    //            //std::vector<EdgeInt> CurrMidP;

    //            std::vector<std::pair<int,int> > EdgeInt;
    //            std::vector<ScalarType> Alpha;
    //            bool foundMidP=GetSideMidpoints(i,CornerSet,CornersExtremes,EdgeInt,Alpha);
    //            if (!foundMidP)return false;

    //            assert(CornersExtremes.size()==EdgeInt.size());
    //            assert(CornersExtremes.size()==Alpha.size());

    //            for (size_t j=0;j<CornersExtremes.size();j++)
    //            {
    //                int IndexC0=CornersExtremes[j].first;
    //                int IndexC1=CornersExtremes[j].second;
    //                std::pair<int,int> Entry(std::min(IndexC0,IndexC1),std::max(IndexC0,IndexC1));
    //                if (AddedMidP.count(Entry)>0)continue;
    //                AddedMidP.insert(Entry);

    //                int IndexV0=EdgeInt[j].first;
    //                int IndexV1=EdgeInt[j].second;
    //                assert(IndexV0<IndexV1);
    //                ScalarType CurrAlpha=Alpha[j];
    //                assert(IndexV0>=0);
    //                assert(IndexV1>=0);
    //                assert(IndexV0<VertPos.size());
    //                assert(IndexV1<VertPos.size());
    //                EdgeAlpha[EdgeInt[j]]=CurrAlpha;

    //                Geo::Point3<ScalarType> P0=VertPos[IndexV0];
    //                Geo::Point3<ScalarType> P1=VertPos[IndexV1];
    //                Geo::Point3<ScalarType> NewPos= (P0*CurrAlpha+P1*(1-CurrAlpha));

    //                SidesMidP.push_back(NewPos);
    //            }
    //        }


    //        std::vector<Geo::Point3<ScalarType> > SaveMPos;
    //        std::vector<std::vector<int> > SaveMFaces;

    //        //compute tutte embedding and split faces
    //        for (size_t i=0;i<SubPatchFaces.size();i++)
    //        {
    //            //avoid empty patches
    //            if (SubPatchFaces[i].size()==0)continue;

    //            //compute tutte's embedding
    //            std::vector<Geo::Point2<ScalarType> > VertUV;
    //            bool done_param=Geo::UVParam<ScalarType>::TutteEmbedding(SubPatchVertPos[i],SubPatchFaces[i],VertUV);
    //            assert(done_param);

    //            //get the split segments
    //            //std::vector<Geo::Segment2<ScalarType> > SplitSeg;
    //            std::vector<std::vector<Geo::Point2<ScalarType> > > SplitPoly;
    //            for (size_t j=0;j<SubPatchFaces[i].size();j++)
    //            {
    //                for (size_t k=0;k<3;k++)
    //                {
    //                    int IndexV0=SubPatchFaces[i][j][k];
    //                    int IndexV1=SubPatchFaces[i][j][(k+1)%3];

    //                    Geo::Point2<ScalarType> UV0=VertUV[IndexV0];
    //                    Geo::Point2<ScalarType> UV1=VertUV[IndexV1];

    //                    IndexV0=MapToOriginal(i,IndexV0);
    //                    IndexV1=MapToOriginal(i,IndexV1);

    //                    std::pair<int,int> Key(std::min(IndexV0,IndexV1),
    //                                           std::max(IndexV0,IndexV1));

    //                    if (EdgeAlpha.count(Key)==0)continue;
    //                    ScalarType alpha=EdgeAlpha[Key];

    //                    Geo::Point2<ScalarType> InterpPUV=UV0*alpha+UV1*(1-alpha);
    //                    //inflate a bit
    //                    InterpPUV=InterpPUV*1.01;
    //                    SplitPoly.resize(SplitPoly.size()+1);
    //                    SplitPoly.back().push_back(Geo::Point2<ScalarType>(0,0));
    //                    SplitPoly.back().push_back(InterpPUV);
    ////                    Geo::Segment2<ScalarType> CurrSeg(Geo::Point2<ScalarType>(0,0),InterpPUV);
    ////                    SplitSeg.push_back(CurrSeg);
    //                }
    //            }
    //            //std::cout<<"Size split:"<<SplitP.size()<<std::endl;

    ////            std::vector<Geo::Point3<ScalarType> > TestPos;
    ////            for (size_t i=0;i<VertUV.size();i++)
    ////            {
    ////                Geo::Point2<ScalarType> uvP=VertUV[i];
    ////                TestPos.push_back(Geo::Point3<ScalarType>(uvP.X,uvP.Y,0));
    ////            }
    //            std::vector<std::pair<int,int> > NewEdges;
    //            //SplitMeshWithUVPolyline(TestPos,VertUV,SubPatchFaces[i],SplitSeg,NewEdges);
    //            SplitMeshWithUVPolylines(SubPatchVertPos[i],VertUV,SubPatchFaces[i],SplitPoly,NewEdges);

    //            std::vector<std::vector<int> > NextFP;
    //            std::vector<std::vector<int> > NextEP;
    //            getFFAdjacency(SubPatchFaces[i],NextFP,NextEP);

    //            SplitMeshFromEdges(SubPatchVertPos[i],SubPatchFaces[i],NextFP,NextEP,NewEdges);

    //            AppendMesh(SaveMPos,SaveMFaces,SubPatchVertPos[i],SubPatchFaces[i]);

    ////            std::string NameF= "Test_split" + std::to_string(i) + ".obj";
    ////            WriteOBJ(NameF,SubPatchVertPos[i],SubPatchFaces[i]);


    ////            SplitMeshFromEdges(TestPos,SubPatchFaces[i],NextFP,NextEP,NewEdges);
    ////            std::string NameF= "Test_split" + std::to_string(i) + ".obj";
    ////            WriteOBJ(NameF,TestPos,SubPatchFaces[i]);


    ////            std::vector<std::pair<int,int> > NewEdges;
    ////            SplitMeshWithUVPolyline(SubPatchVertPos[i],VertUV,SubPatchFaces[i],
    ////                                    SplitSeg,NewEdges);

    ////            std::vector<std::vector<int> > NextFP;
    ////            std::vector<std::vector<int> > NextEP;
    ////            getFFAdjacency(SubPatchFaces[i],NextFP,NextEP);

    ////            SplitMeshFromEdges(SubPatchVertPos[i],SubPatchFaces[i],NextFP,NextEP,NewEdges);
    ////            std::string NameF= "Test_split" + std::to_string(i) + ".obj";
    ////            WriteOBJ(NameF,SubPatchVertPos[i],SubPatchFaces[i]);
    //        }
    //        WriteOBJ("./dualized.obj",SaveMPos,SaveMFaces);
    //        return true;
    //    }


    void CornersOfType(const int &IndexPatch,
                       const CornerType &CType,
                       std::vector<int> &SelectCorn)const
    {
        SelectCorn.clear();
        assert(IndexPatch<(int)PatchCorners.size());
        assert(IndexPatch<(int)TypeOfCorner.size());

        for (size_t i=0;i<TypeOfCorner[IndexPatch].size();i++)
            for (size_t j=0;j<TypeOfCorner[IndexPatch][i].size();j++)
                if (TypeOfCorner[IndexPatch][i][j]==CType)
                {
                    int CornerIdx=PatchCorners[IndexPatch][i][j];
                    CornerIdx=MapToOriginal(IndexPatch,CornerIdx);
                    SelectCorn.push_back(CornerIdx);
                }
    }

    void getMeshOfPatches(const std::vector<int> &To_CopyFrom,
                          std::vector<Geo::Point3<ScalarType> > &CopiedMeshesVert,
                          std::vector<std::vector<int> > &CopiedMeshesElem,
                          std::vector<int> &OriginalVert,std::vector<int> &OriginalElem)
    {
        std::vector<int> ToCopyFace;
        for (int i=0;i<To_CopyFrom.size();i++)
        {
            int IndexPatch=To_CopyFrom[i];
            ToCopyFace.insert(ToCopyFace.end(),
                              SubPatchFacesToOriginal[IndexPatch].begin(),
                              SubPatchFacesToOriginal[IndexPatch].end());
        }
        CopiedMeshesVert.clear();
        CopiedMeshesElem.clear();
        OriginalVert.clear();
        OriginalElem.clear();
        CopySubMesh<ScalarType>(VertPos,Faces,ToCopyFace,
                                CopiedMeshesVert,CopiedMeshesElem,
                                OriginalVert,OriginalElem);
    }

    void ComposeMeshFromPatches(std::vector<Geo::Point3<ScalarType> > &VertPos,
                                std::vector<std::vector<int> >  &Faces)const
    {
        VertPos.clear();
        Faces.clear();

        //append all faces
        std::vector<std::pair<int,int> >  PatchF;
        for (size_t i=0;i<NumPatches();i++)
        {
            if (isEmpty(i))continue;

            AppendMesh(VertPos,Faces,SubPatchVertPos[i],SubPatchFaces[i]);
            for (size_t j=0;j<SubPatchFaces[i].size();j++)
                PatchF.push_back(std::pair<int,int>(i,j));
        }

        assert(PatchF.size()==Faces.size());

        //then reorder
        std::vector<std::vector<int> >  NewFaces=Faces;
        for (size_t i=0;i<PatchF.size();i++)
        {
            int IndexP=PatchF[i].first;
            int IndexF=PatchF[i].second;
            assert(IndexP<NumPatches());
            assert(IndexP>=0);
            assert(IndexF<SubPatchFacesToOriginal[IndexP].size());
            int NewFIdx=SubPatchFacesToOriginal[IndexP][IndexF];
            NewFaces[NewFIdx]=Faces[i];
        }
        Faces=NewFaces;
    }

    void AddTouchedPatchesToSet(const std::vector<int> &PathN,
                                std::set<int> &PatchIdx)
    {
        if (PathN.size()==0)return;
        //loop repeat the last two times so no need to check
        for (int i=0;i<PathN.size()-1;i++)
        {
            int IndexN0=PathN[i];
            int IndexN1=PathN[i+1];
            int IndexV0=Field::CrossFieldGraph::ElemOfNode(IndexN0);
            int IndexV1=Field::CrossFieldGraph::ElemOfNode(IndexN1);
            std::pair<int,int> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            assert(EdgesToFace.count(key)>0);
            for (size_t j=0;j<EdgesToFace[key].size();j++)
            {
                int IndexF=EdgesToFace[key][j];
                assert(IndexF>=0);
                assert(IndexF<(int)Faces.size());
                assert(IndexF<(int)OriginalFaceToPatch.size());
                PatchIdx.insert(OriginalFaceToPatch[IndexF]);
            }
        }
    }

    void GetTouchedPatchesFromPathIndex(const std::vector<int> &PathIDX,
                                         std::set<int> &PatchIdx)
    {
        PatchIdx.clear();
        for (size_t i=0;i<PathIDX.size();i++)
        {
            int IndexPath=PathIDX[i];
            assert(IndexPath>=0);
            assert(IndexPath<PathNodes.size());
            AddTouchedPatchesToSet(PathNodes[IndexPath],PatchIdx);
        }
    }

    void GetFaceOnPath(const int IndexPath,std::vector<int> &Faces)
    {
        Faces.clear();
        assert(IndexPath>=0);
        assert(IndexPath<PathNodes.size());
        for (int i=0;i<PathNodes[IndexPath].size()-1;i++)
        {
            int IndexN0=PathNodes[IndexPath][i];
            int IndexN1=PathNodes[IndexPath][i+1];
            int IndexV0=Field::CrossFieldGraph::ElemOfNode(IndexN0);
            int IndexV1=Field::CrossFieldGraph::ElemOfNode(IndexN1);
            std::pair<int,int> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            assert(EdgesToFace.count(key)>0);
            for (size_t j=0;j<EdgesToFace[key].size();j++)
            {
                int IndexF=EdgesToFace[key][j];
                assert(IndexF>=0);
                assert(IndexF<(int)Faces.size());
                assert(IndexF<(int)OriginalFaceToPatch.size());
                Faces.push_back(OriginalFaceToPatch[IndexF]);
            }
        }
    }

    void GetFaceOnPaths(const std::vector<int> &IndexPaths,std::vector<int> &Faces)
    {
       Faces.clear();
       for (size_t i=0;i<IndexPaths.size();i++)
       {
           std::vector<int> PathF;
           GetFaceOnPath(IndexPaths[i],PathF);
           Faces.insert(PathF.end(),PathF.begin(),PathF.end());
       }
    }

    void CompactEmptyPatches(std::map<int,int> &PatchIdxRemap)
    {
        std::vector<std::vector<std::vector<ScalarType> > > SideLenghts_Swap;
        std::vector<std::vector<int> > PatchOriginalFaces_Swap;
        std::vector<std::vector<Geo::Point3<ScalarType> > > SubPatchVertPos_Swap;
        std::vector<std::vector<std::vector<int> > > SubPatchFaces_Swap;
        std::vector<std::vector<int > > SubPatchVertToOriginal_Swap;
        std::vector<std::vector<int > > SubPatchFacesToOriginal_Swap;
        std::vector<std::vector<std::vector<int> > > PatchCorners_Swap;
        std::vector<std::vector<std::vector<CornerType> > > TypeOfCorner_Swap;
        std::vector<std::vector<std::vector<std::vector<std::pair<int,int> > > > > EdgePatchSides_Swap;
        std::vector<std::vector<Geo::Point3<ScalarType> > > VisualCorners_Swap;
        std::vector<std::vector<CornerType> > VisualCornerType_Swap;
        std::vector<std::vector<std::vector<Geo::Point3<ScalarType> > > > VisualPaths_Swap;


        for (size_t i=0;i<PatchOriginalFaces.size();i++)
        {
            if (PatchOriginalFaces[i].size()==0)continue;

            SideLenghts_Swap.push_back(SideLenghts[i]);
            PatchOriginalFaces_Swap.push_back(PatchOriginalFaces[i]);
            SubPatchVertPos_Swap.push_back(SubPatchVertPos[i]);
            SubPatchFaces_Swap.push_back(SubPatchFaces[i]);
            SubPatchVertToOriginal_Swap.push_back(SubPatchVertToOriginal[i]);
            SubPatchFacesToOriginal_Swap.push_back( SubPatchFacesToOriginal[i]);
            PatchCorners_Swap.push_back(PatchCorners[i]);
            TypeOfCorner_Swap.push_back(TypeOfCorner[i]);
            EdgePatchSides_Swap.push_back(EdgePatchSides[i]);

            PatchIdxRemap[i]=EdgePatchSides_Swap.size()-1;

            if (UpdateVisual)
            {
                VisualCorners_Swap.push_back(VisualCorners[i]);
                VisualCornerType_Swap.push_back(VisualCornerType[i]);
                VisualPaths_Swap.push_back(VisualPaths[i]);
            }
        }

        SideLenghts=SideLenghts_Swap;
        PatchOriginalFaces=PatchOriginalFaces_Swap;
        SubPatchVertPos=SubPatchVertPos_Swap;
        SubPatchFaces=SubPatchFaces_Swap;
        SubPatchVertToOriginal=SubPatchVertToOriginal_Swap;
        SubPatchFacesToOriginal=SubPatchFacesToOriginal_Swap;
        PatchCorners=PatchCorners_Swap;
        TypeOfCorner=TypeOfCorner_Swap;
        EdgePatchSides=EdgePatchSides_Swap;
        if (UpdateVisual)
        {
            VisualCorners=VisualCorners_Swap;
            VisualCornerType=VisualCornerType_Swap;
            VisualPaths=VisualPaths_Swap;
        }

        UpdateFaceToPatch();
    }

    bool HasMultipleSingularities(int IndexPatch)const
    {
        int numSing=0;
        for (size_t i=0;i<SubPatchVertToOriginal[IndexPatch].size();i++)
        {
            int IndexVOr=SubPatchVertToOriginal[IndexPatch][i];
            if (SingValence.count(IndexVOr)==0)continue;
            numSing++;
            if (numSing>1)return true;
        }
        return false;
    }

    bool MatchSingularity(int IndexPatch)const
    {
        assert(IndexPatch>=0);
        assert(IndexPatch<SubPatchVertToOriginal.size());

        int TargetValence=-1;
        for (size_t i=0;i<SubPatchVertToOriginal[IndexPatch].size();i++)
        {
            int IndexVOr=SubPatchVertToOriginal[IndexPatch][i];
            if (SingValence.count(IndexVOr)==0)continue;
            if (TargetValence!=-1)return false;//multiple singularities
            TargetValence=SingValence.at(IndexVOr);
        }
        //std::cout<<"Target Val 0 "<<TargetValence<<std::endl;

        //CHECK , NEED A BETTER VALENCE CALCULATION
        switch(TargetValence) {
        case -1:
            TargetValence=4;
            break;
        case 1:
            TargetValence=5;
            break;
        case 2:
            TargetValence=6;
            break;
        case 3:
            TargetValence=3;
            break;
        default:
            TargetValence=4;
        }

        return (NumSides(IndexPatch)==TargetValence);
    }



    void SetSingularities(const std::vector<int> &SingIndex,
                          const std::vector<int> &SingValue)
    {
        SingValence.clear();
        assert(SingIndex.size()==SingValue.size());
        for (size_t i=0;i<SingIndex.size();i++)
            SingValence[SingIndex[i]]=SingValue[i];
    }

    int NumSides(const int &IndexPatch)const
    {
        int sideNum=0;
        for (size_t i=0;i<EdgePatchSides[IndexPatch].size();i++)
            sideNum+=EdgePatchSides[IndexPatch][i].size();
        return sideNum;
    }

    int NumCorners(const int &IndexPatch)const
    {
        int cornNum=0;
        for (size_t i=0;i<PatchCorners[IndexPatch].size();i++)
            cornNum+=PatchCorners[IndexPatch][i].size();
        return cornNum;
    }

    int NumBorders(const int &IndexPatch)const
    {
        assert(IndexPatch<NumPatches());
        return EdgePatchSides[IndexPatch].size();
    }

    int Genus(const int &IndexPatch)const
    {
        assert(IndexPatch<NumPatches());
        int NHoles=NumBorders(IndexPatch);
        int NEdges=NumEdges(SubPatchFaces[IndexPatch]);
        int NVert=SubPatchVertPos[IndexPatch].size();
        int NFaces=SubPatchFaces[IndexPatch].size();

        int Gen=PolyMeshGenus(NVert,NEdges,NFaces,NHoles,1);
        return Gen;
    }

    bool isDiskLike(const int &IndexPatch)const
    {
        //if (NumBorders(IndexPatch)!=1)return false;
        if (Genus(IndexPatch)!=0)return false;
        return true;
    }

    bool isEmpty(const int &IndexPatch)const
    {
        assert(IndexPatch<NumPatches());
        return(PatchOriginalFaces[IndexPatch].empty());
    }

    void InitVisualCornersEdges()
    {
        for (size_t i=0;i<NumPatches();i++)
        {
            InitVisualCorners(i);
            InitVisualSides(i);
        }
    }

    void UpdateSubPatchPos()
    {
        for (size_t i=0;i<NumPatches();i++)
            UpdateSubPatchPos(i);
    }

    //    bool SavePatchesLayout(const std::string meshFile)
    //    {
    //        CompactEmptyPatches();
    //        std::ofstream fout(meshFile.c_str());
    //        if(!fout) return false;

    //        //first set the number of patches
    //        fout<<"#p"<<EdgePatchSides.size()<<"\n";

    //        //for each patch
    //        for (size_t pi=0;pi<EdgePatchSides.size();pi++)
    //        {
    //            //write patch index
    //            fout<<"p "<<pi<<"\n";

    //            //for each loops
    //            fout<<"#l "<<EdgePatchSides[pi].size()<<"\n";
    //            assert(EdgePatchSides[pi].size()==PatchCorners[pi].size());

    //            //for each loop
    //            for (size_t li=0;li<PatchCorners[pi].size();li++)
    //            {
    //                fout<<"l "<<li<<"\n";
    //                //number of sides equal nummber of corners
    //                assert(PatchCorners[pi][li].size()==EdgePatchSides[pi][li].size());

    //                //write corners
    //                fout<<"#c "<<PatchCorners[pi][li].size()<<"\n";

    //                for (size_t ci=0;ci<PatchCorners[pi][li].size();ci++)
    //                {
    //                    int IndexC=PatchCorners[pi][li][ci];
    //                    assert(IndexC>=0);
    //                    assert(IndexC<VertPos.size());
    //                    fout<<"c "<<IndexC<<"\n";
    //                }

    //                //write number of sides
    //                fout<<"#s "<<EdgePatchSides[pi][li].size()<<"\n";
    //                for (size_t si=0;si<EdgePatchSides[pi][li].size();si++)
    //                {
    //                    fout<<"s "<<si<<"\n";
    //                    int sizeL=EdgePatchSides[pi][li][si].size();
    //                    fout<<"#e "<<sizeL<<"\n";

    //                    //for each edge
    //                    for (size_t ei=0;ei<EdgePatchSides[pi][li][si].size();ei++)
    //                    {
    //                        int IndexV0=EdgePatchSides[pi][li][si][ei].first;
    //                        int IndexV1=EdgePatchSides[pi][li][si][ei].second;
    //                        assert(IndexV0>=0);
    //                        assert(IndexV0<VertPos.size());
    //                        assert(IndexV1>=0);
    //                        assert(IndexV1<VertPos.size());
    //                        fout<<"e "<<IndexV0<<" "<<IndexV1<<"\n";
    //                    }
    //                }
    //            }
    //        }
    //        fout.close();
    //        return true;
    //    }


    //    void SaveSubMeshes(const std::string PathSave)
    //    {
    //        UpdateSubPatchPos();

    //        std::vector<Geo::Point3<ScalarType> > PosVect;
    //        std::vector<std::vector<int> >  FaceVect;
    //        assert(SubPatchVertPos.size()==SubPatchFaces.size());

    //        for (size_t i=0;i<SubPatchVertPos.size();i++)
    //        {
    //            int baseIdx=PosVect.size();
    //            PosVect.insert(PosVect.end(),SubPatchVertPos[i].begin(),SubPatchVertPos[i].end());

    //            for (size_t j=0;j<SubPatchFaces[i].size();j++)
    //            {
    //                FaceVect.push_back(std::vector<int>(3,-1));
    //                assert(SubPatchFaces[i][j].size()==3);
    //                for (size_t k=0;k<3;k++)
    //                {
    //                    int IndexV=SubPatchFaces[i][j][k];
    //                    FaceVect.back()[k]=IndexV+baseIdx;
    //                }
    //            }
    //        }
    //        WriteOBJ(PathSave.c_str(),PosVect,FaceVect);
    //    }

    void FindPatchCorners()
    {
        for (size_t i=0;i<NumPatches();i++)
            FindPatchCorners(i);
    }


    void GetColorFaceByPatch(std::vector<Geo::Point3<ScalarType> > &FaceCol)
    {
        if (NumPatches()==0)return;
        FaceCol.clear();
        FaceCol.resize(Faces.size(),Geo::Point3<ScalarType>(1,1,1));
        for (size_t i=0;i<Faces.size();i++)
            FaceCol[i]=Geo::ScatterColor<ScalarType>(NumPatches(),OriginalFaceToPatch[i]);
    }

    void GetColorFaceIfDiskLike(std::vector<Geo::Point3<ScalarType> > &FaceCol)
    {
        if (NumPatches()==0)return;
        FaceCol.clear();
        FaceCol.resize(Faces.size(),Geo::Point3<ScalarType>(1,1,1));
        for (size_t i=0;i<PatchOriginalFaces.size();i++)
        {
            Geo::Point3<ScalarType> Col;
            if (isDiskLike(i))
                Col=Geo::Point3<ScalarType>(0,1,0);
            else
                Col=Geo::Point3<ScalarType>(1,0,0);
            for (size_t j=0;j<PatchOriginalFaces[i].size();j++)
            {
                int IndexF=PatchOriginalFaces[i][j];
                //std::cout<<IndexF<<std::endl;

                FaceCol[IndexF]=Col;
            }
        }
    }

    void GetColorByValence(std::vector<Geo::Point3<ScalarType> > &FaceCol)
    {
        if (NumPatches()==0)return;
        FaceCol.clear();
        FaceCol.resize(Faces.size(),Geo::Point3<ScalarType>(1,1,1));
        assert(PatchOriginalFaces.size()==PatchCorners.size());
        assert(PatchOriginalFaces.size()==TypeOfCorner.size());
        for (size_t i=0;i<PatchOriginalFaces.size();i++)
        {
            int sideN=NumSides(i);
            int cornN=NumCorners(i);
            assert(sideN==cornN);
            Geo::Point3<ScalarType> Col=Geo::Point3<ScalarType>(1,0,0);
            if (cornN==3)
                Col=Geo::Point3<ScalarType>(1,1,0);
            if (cornN==4)
                Col=Geo::Point3<ScalarType>(0,1,0);
            if (cornN==5)
                Col=Geo::Point3<ScalarType>(0,1,1);
            if (cornN==6)
                Col=Geo::Point3<ScalarType>(1,1,1);

            for (size_t j=0;j<PatchOriginalFaces[i].size();j++)
            {
                int IndexF=PatchOriginalFaces[i][j];
                assert(IndexF>=0);
                assert(IndexF<FaceCol.size());
                //std::cout<<"Face Idx:"<<IndexF<<std::endl;
                FaceCol[IndexF]=Col;
            }
        }
    }

    void GetColorByMatchSingularity(std::vector<Geo::Point3<ScalarType> > &FaceCol)
    {
        if (NumPatches()==0)return;
        FaceCol.clear();
        FaceCol.resize(Faces.size(),Geo::Point3<ScalarType>(1,1,1));
        assert(PatchOriginalFaces.size()==PatchCorners.size());
        assert(PatchOriginalFaces.size()==TypeOfCorner.size());

        for (size_t i=0;i<PatchOriginalFaces.size();i++)
        {
            bool MatchS=MatchSingularity(i);
            Geo::Point3<ScalarType> Col;
            if (MatchS)
                Col=Geo::Point3<ScalarType>(0,1,0);
            else
                Col=Geo::Point3<ScalarType>(1,0,0);

            for (size_t j=0;j<PatchOriginalFaces[i].size();j++)
            {
                int IndexF=PatchOriginalFaces[i][j];
                assert(IndexF>=0);
                assert(IndexF<FaceCol.size());
                FaceCol[IndexF]=Col;
            }
        }
    }

    void Update(const std::vector<int> &AddedPaths,
                const std::vector<int> &RemovedPaths,
                std::vector<int> &UpdatedPatches)
    {
        UpdatedPatches.clear();

//        std::vector<int> ToAdd,ToRemove;
        std::set<int> ToUpdPatches;
//        FindAddRemovePaths(ToAdd,ToRemove);

        //        std::cout<<"Adding "<<ToAdd.size()<<" Paths"<<std::endl;
        //        std::cout<<"Removing "<<ToRemove.size()<<" Paths"<<std::endl;

        for (size_t i=0;i<RemovedPaths.size();i++)
        {
            int IndexPath=RemovedPaths[i];
            AddTouchedPatchesToSet(OLDPathNodes[IndexPath],ToUpdPatches);
        }

        for (size_t i=0;i<AddedPaths.size();i++)
        {
            int IndexPath=AddedPaths[i];
            AddTouchedPatchesToSet(PathNodes[IndexPath],ToUpdPatches);
        }

        //        std::cout<<"To Update "<<ToUpdPatches.size()<<" Patches"<<std::endl;

        for (size_t i=0;i<RemovedPaths.size();i++)
        {
            int IndexPath=RemovedPaths[i];
            RemovePath(OLDPathNodes[IndexPath]);
        }

        for (size_t i=0;i<AddedPaths.size();i++)
        {
            int IndexPath=AddedPaths[i];
            AddPath(PathNodes[IndexPath]);
        }

        std::vector<int> ToUpdPatchesV(ToUpdPatches.begin(),ToUpdPatches.end());

        int NumPatch0=NumPatches();
        UpdatePatches(ToUpdPatchesV);
        int NumPatch1=NumPatches();
        assert(NumPatch1<=NumPatch1);

        std::vector<int> ToUpdatePatches;

        for (int i=NumPatch0;i<NumPatch1;i++)
            ToUpdatePatches.push_back(i);


        UpdateSubMeshes(ToUpdatePatches);

        FindPatchCorners(ToUpdatePatches);

        UpdateSideLenghs(ToUpdatePatches);

        if (UpdateVisual)
            InitVisualCornersEdges(ToUpdatePatches);

        UpdatedPatches=ToUpdatePatches;

        OLDPathNodes=PathNodes;

        assert(SideLenghts.size()==NumPatches());
        assert(PatchOriginalFaces.size()==NumPatches());
        assert(SubPatchVertPos.size()==NumPatches());
        assert(SubPatchFaces.size()==NumPatches());
        assert(SubPatchVertToOriginal.size()==NumPatches());
        assert(SubPatchFacesToOriginal.size()==NumPatches());
        assert(PatchCorners.size()==NumPatches());
        assert(TypeOfCorner.size()==NumPatches());
        assert(EdgePatchSides.size()==NumPatches());
    }


    void Update(std::vector<int> &UpdatedPatches)
    {
        std::vector<int> ToAdd,ToRemove;
        FindAddRemovePaths(ToAdd,ToRemove);
        Update(ToAdd,ToRemove,UpdatedPatches);
    }

    int NumPatches()const
    {
        return PatchOriginalFaces.size();
    }

    int NumFaces()const
    {
        return Faces.size();
    }

    void WriteStats()
    {
        std::cout<<"There are :"<<NumPatches()<<" Patches"<<std::endl;
        //std::cout<<"With 0 Corners:"<<PatchOriginalFaces.size()<<std::endl;
    }

    void Reset()
    {
        PatchCorners.clear();
        TypeOfCorner.clear();
        SideLenghts.clear();
        PatchOriginalFaces.clear();
        OriginalFaceToPatch.clear();

        SubPatchVertPos.clear();
        SubPatchFaces.clear();

        SubPatchVertToOriginal.clear();
        SubPatchVertToOriginal.clear();

        VisualCorners.clear();
        VisualCornerType.clear();
        VisualPaths.clear();

        std::vector<std::vector<int> > OLDPathNodes=PathNodes;
    }

    void Init()
    {
        Reset();

        getFFAdjacency(Faces,NextF,NextE);

        std::cout<<"- Updating Patches"<<std::endl;

        //make one patch for each connected components
        FFConnectedComponenent(Faces,NextF,NextE,PatchOriginalFaces);
        UpdateFaceToPatch();

        //get the submeshes if needed
        std::cout<<"- Updating SubMeshes"<<std::endl;
        SplitByFFConnectedComponenent(VertPos,Faces,SubPatchVertPos,
                                      SubPatchFaces,SubPatchVertToOriginal,
                                      SubPatchFacesToOriginal);
        //safety checks
        assert(PatchOriginalFaces.size()==SubPatchVertPos.size());
        assert(PatchOriginalFaces.size()==SubPatchFaces.size());
        assert(PatchOriginalFaces.size()==SubPatchVertToOriginal.size());
        assert(PatchOriginalFaces.size()==SubPatchFacesToOriginal.size());
        assert(PatchOriginalFaces.size()==SubPatchFacesToOriginal.size());

        SideLenghts.resize(PatchOriginalFaces.size());
        PatchCorners.resize(PatchOriginalFaces.size());
        TypeOfCorner.resize(PatchOriginalFaces.size());
        EdgePatchSides.resize(PatchOriginalFaces.size());
        VisualCorners.resize(PatchOriginalFaces.size());
        VisualCornerType.resize(PatchOriginalFaces.size());
        VisualPaths.resize(PatchOriginalFaces.size());

        //then set the corners as the bundary convex
        PathsEdges.clear();
        std::cout<<"- Add Border Paths"<<std::endl;
        AddBorderPaths();

        std::cout<<"- Find Patch Corners"<<std::endl;
        FindPatchCorners();

        InitEdgeToFace();

        if (UpdateVisual)
        {
            std::cout<<"- Find Visual Data"<<std::endl;
            InitVisualCornersEdges();
        }
        std::cout<<"DONE Update"<<std::endl;
        WriteStats();
    }

    PatchManaging(std::vector<Geo::Point3<ScalarType> > &_VertPos,
                  std::vector<std::vector<int> > &_Faces,
                  const Geo::TracingGraph<ScalarType> &_TGraph,
                  const std::vector<std::vector<int> > &_PathNodes):
        VertPos(_VertPos),Faces(_Faces),
        TGraph(_TGraph),PathNodes(_PathNodes)
    {
        UpdateVisual=true;
        //        UpdateCorners=true;
        //        UpdateSubMeshes=true;
    }
};

}
#endif
