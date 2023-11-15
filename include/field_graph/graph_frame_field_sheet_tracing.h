#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING
#define FRAME_FIELD_GRAPH_SHEET_TRACING

//#include "../GeoLibrary/frame_field.h"
#include "graph.h"
#include "graph_tracing.h"
#include "./frame_field.h"
#include "graph_frame_field.h"
#include <connectivity.h>
//#define TRACE_DEBUG

template <class ScalarType>
class SheetTracer
{
public:

    struct Sheet
    {
        std::vector<int> TetraIdxSide0;
        std::vector<int> TetraFaceSide0;
        std::vector<int> TetraIdxSide1;
        std::vector<int> TetraFaceSide1;
        std::vector<int> NodeSide0;
        std::vector<int> NodeSide1;
        std::vector<int> SourceNode0;
        std::vector<int> SourceNode1;


        std::vector<ScalarType> NodeSide0Dist,NodeSide1Dist;
        std::vector<std::vector<int> > PathNodeSide0;
        std::vector<std::vector<int> > PathNodeSide1;

        void Clear()
        {
            TetraIdxSide0.clear();
            TetraFaceSide0.clear();
            TetraIdxSide1.clear();
            TetraFaceSide1.clear();
            NodeSide0.clear();
            NodeSide1.clear();
            SourceNode0.clear();
            SourceNode1.clear();
            NodeSide0Dist.clear();
            NodeSide1Dist.clear();
            PathNodeSide0.clear();
            PathNodeSide1.clear();
        }

        bool ProperSources()
        {
            int IndexPoly0=PolyOfNode(SourceNode0[0]);
            int IndexPoly1=PolyOfNode(SourceNode1[0]);

            for (size_t i=1;i<SourceNode0.size();i++)
                if(IndexPoly0!=PolyOfNode(SourceNode0[i]))return false;

            for (size_t i=1;i<SourceNode0.size();i++)
                if(IndexPoly1!=PolyOfNode(SourceNode1[i]))return false;

            return true;
        }

        bool IsProper()const
        {
            for (size_t i=0;i<SourceNode0.size();i++)
                for (size_t j=0;j<SourceNode1.size();j++)
                    if(SourceNode0[i]==SourceNode1[j])return false;

            if(NodeSide0.size()!=NodeSide1.size())return false;
            if(NodeSide0.size()!=TetraFaceSide0.size())return false;
            if(NodeSide0.size()!=TetraFaceSide1.size())return false;
            if(NodeSide0.size()!=TetraIdxSide0.size())return false;
            if(NodeSide0.size()!=TetraIdxSide1.size())return false;
            if(NodeSide0.size()!=NodeSide0Dist.size())return false;
            if(NodeSide0.size()!=NodeSide1Dist.size())return false;
            if(NodeSide0.size()!=PathNodeSide0.size())return false;
            if(NodeSide0.size()!=PathNodeSide1.size())return false;
            return true;
        }

        int WhichSource(int IndexN)const
        {
            std::vector<int>::const_iterator it;
            it = std::find(SourceNode0.begin(), SourceNode0.end(), IndexN);
            if (it != SourceNode0.end()) return 0;
            it = std::find(SourceNode1.begin(), SourceNode1.end(), IndexN);
            if (it != SourceNode1.end()) return 1;
            return -1;
        }

        void CopyEntry(const Sheet &S1,const int IndexE)
        {
            assert(IndexE>=0);

            assert(IndexE<S1.NodeSide0.size());
            assert(IndexE<S1.NodeSide1.size());
            assert(IndexE<S1.TetraFaceSide0.size());
            assert(IndexE<S1.TetraFaceSide1.size());
            assert(IndexE<S1.TetraIdxSide0.size());
            assert(IndexE<S1.TetraIdxSide1.size());
            assert(IndexE<S1.NodeSide1Dist.size());
            assert(IndexE<S1.NodeSide1Dist.size());
            assert(IndexE<S1.PathNodeSide0.size());
            assert(IndexE<S1.PathNodeSide0.size());


            NodeSide0.push_back(S1.NodeSide0[IndexE]);
            NodeSide1.push_back(S1.NodeSide1[IndexE]);
            TetraFaceSide0.push_back(S1.TetraFaceSide0[IndexE]);
            TetraFaceSide1.push_back(S1.TetraFaceSide1[IndexE]);
            TetraIdxSide0.push_back(S1.TetraIdxSide0[IndexE]);
            TetraIdxSide1.push_back(S1.TetraIdxSide1[IndexE]);

            NodeSide0Dist.push_back(S1.NodeSide0Dist[IndexE]);
            NodeSide1Dist.push_back(S1.NodeSide1Dist[IndexE]);
            PathNodeSide0.push_back(S1.PathNodeSide0[IndexE]);
            PathNodeSide1.push_back(S1.PathNodeSide1[IndexE]);
        }
    };

private:

    //    static void SideSheetSeparetrix(const Geo::TracingGraph<ScalarType> &TGraph,
    //                                    const int &T0,const Field::FrameField<ScalarType> &FT0,
    //                                    const int &T1,const Field::FrameField<ScalarType> &FT1,
    //                                    const Sheet &S,
    //                                    int &SetT0,int &SetT1,
    //                                    int &Node0,int &Node1)
    //    {
    //        std::vector<int> TetraN0;
    //        NodesOfPoly(T0,TetraN0);
    //        std::vector<int> TetraN1;
    //        NodesOfPoly(T1,TetraN1);
    //        assert(TetraN0.size()==6);
    //        assert(TetraN0.size()==TetraN1.size());

    //        SetT0=-1;SetT1=-1;
    //        Node0=-1;Node1=-1;
    //        for (size_t i=0;i<TetraN0.size();i++)
    //        {
    //            int IndexN0=TetraN0[i];
    //            int SourceN0=-1;

    //            if (TGraph.IsMarked(IndexN0))
    //                SourceN0=TGraph.Source(IndexN0);

    //            if(SourceN0==-1)continue;

    //            int SourceSetN0=S.WhichSource(SourceN0);
    //            assert(SourceSetN0!=-1);
    //            assert((SourceSetN0==0)||(SourceSetN0==1));
    //            int IndexFDir0=FrameOfNode(IndexN0);
    //            for (size_t j=0;j<TetraN1.size();j++)
    //            {
    //                int IndexN1=TetraN1[j];
    //                assert(IndexN0!=IndexN1);

    //                int SourceN1=-1;
    //                if (TGraph.IsMarked(IndexN1))
    //                    SourceN1=TGraph.Source(IndexN1);

    //                if(SourceN1==-1)continue;

    //                int SourceSetN1=S.WhichSource(SourceN1);
    //                if (SourceSetN1==-1)
    //                {
    //                    std::cout<<"Node Test"<<IndexN1<<std::endl;
    //                    std::cout<<"Source Test 0"<<std::endl;
    //                    for (size_t k=0;k<S.SourceNode0.size();k++)
    //                        std::cout<<S.SourceNode0[k]<<std::endl;
    //                    std::cout<<"Node Test"<<IndexN1<<std::endl;
    //                    std::cout<<"Source Test 1"<<std::endl;
    //                    for (size_t k=0;k<S.SourceNode1.size();k++)
    //                        std::cout<<S.SourceNode1[k]<<std::endl;
    //                    assert(0);
    //                }
    //                assert((SourceSetN1==0)||(SourceSetN1==1));

    //                //if same source set nothing to do
    //                if (SourceSetN0==SourceSetN1)continue;

    //                //otherwise check directions
    //                int IndexFDir1=FrameOfNode(IndexN1);

    //                //if they are orthogonal does not care
    //                Field::RelationType Rel0=Field::FrameRelationDir(FT0,IndexFDir0,FT1,IndexFDir1);
    //                Field::RelationType Rel1=Field::FrameRelationDir(FT1,IndexFDir1,FT0,IndexFDir0);
    //                //                if (Rel0!=Rel1)
    //                //                {
    //                //                   std::cout<<"R0:"<<Rel0<<" R1:"<<Rel1<<std::endl;
    //                //                   assert(0);
    //                //                }
    //                assert(Rel0==Rel1);
    //                //                if (Rel0==Field::NodeOrthogonal)continue;
    //                //                if (Rel0==Field::NodeOpposite)continue;
    //                //if (Rel0!=Field::NodeOpposite)continue;
    //                if (Rel0!=Field::NodeCollinear)continue;
    //                //return;
    //                SetT0=SourceSetN0;
    //                SetT1=SourceSetN1;
    //                Node0=IndexN0;
    //                Node1=IndexN1;
    //                return;
    //            }
    //        }
    //    }

    static ScalarType PathDistances(const ScalarType &D0,const ScalarType &D1)
    { return (D0*D0+D1*D1);}

    static void SideSheetSeparetrix(const Geo::TracingGraph<ScalarType> &TGraph,
                                    const int &T0,const Field::FrameField<ScalarType> &FT0,
                                    const int &T1,const Field::FrameField<ScalarType> &FT1,
                                    const Sheet &S,
                                    int &SetT0,int &SetT1,
                                    int &Node0,int &Node1)
    {
        std::vector<int> TetraN0;
        NodesOfPoly(T0,TetraN0);
        std::vector<int> TetraN1;
        NodesOfPoly(T1,TetraN1);
        assert(TetraN0.size()==6);
        assert(TetraN0.size()==TetraN1.size());

        SetT0=-1;SetT1=-1;
        Node0=-1;Node1=-1;
        ScalarType MinD=std::numeric_limits<ScalarType>::max();

        for (size_t i=0;i<TetraN0.size();i++)
        {
            int IndexN0=TetraN0[i];

            //should both be selected
            if (!TGraph.IsSelected(IndexN0,0))continue;
            if (!TGraph.IsSelected(IndexN0,1))continue;

            int SourceN0=-1;

            if (TGraph.IsMarked(IndexN0))
                SourceN0=TGraph.Source(IndexN0);

            if(SourceN0==-1)continue;

            int SourceSetN0=S.WhichSource(SourceN0);
            assert(SourceSetN0!=-1);
            assert((SourceSetN0==0)||(SourceSetN0==1));
            int IndexFDir0=FrameOfNode(IndexN0);
            for (size_t j=0;j<TetraN1.size();j++)
            {
                int IndexN1=TetraN1[j];

                //should both be selected
                if (!TGraph.IsSelected(IndexN1,0))continue;
                if (!TGraph.IsSelected(IndexN1,1))continue;

                assert(IndexN0!=IndexN1);

                int SourceN1=-1;
                if (TGraph.IsMarked(IndexN1))
                    SourceN1=TGraph.Source(IndexN1);

                if(SourceN1==-1)continue;

                int SourceSetN1=S.WhichSource(SourceN1);
                if (SourceSetN1==-1)
                {
                    assert(0);
                }
                assert((SourceSetN1==0)||(SourceSetN1==1));

                //if same source set nothing to do
                if (SourceSetN0==SourceSetN1)continue;

                //otherwise check directions
                int IndexFDir1=FrameOfNode(IndexN1);

                //if they are orthogonal does not care
                Field::RelationType Rel0=Field::FrameRelationDir(FT0,IndexFDir0,FT1,IndexFDir1);
                Field::RelationType Rel1=Field::FrameRelationDir(FT1,IndexFDir1,FT0,IndexFDir0);
                assert(Rel0==Rel1);

                if (Rel0!=Field::NodeCollinear)continue;

                ScalarType D0=TGraph.Distance(IndexN0);
                ScalarType D1=TGraph.Distance(IndexN1);
                ScalarType PathD=PathDistances(D0,D1);
                if (PathD<MinD)
                {
                    SetT0=SourceSetN0;
                    SetT1=SourceSetN1;
                    Node0=IndexN0;
                    Node1=IndexN1;
                    MinD=PathD;
                }
                //return;
            }
        }
    }

public:

    //    static void ExtractSheetFrom(Geo::TracingGraph<ScalarType> &TGraph,
    //                                 const std::vector<std::vector<int> > &NextPoly,
    //                                 const std::vector<std::vector<int> > &NextPolyF,
    //                                 Sheet &S,bool WriteDbg=true)
    //    {

    //        S.NodeSide0.clear();
    //        S.NodeSide1.clear();
    //        S.TetraFaceSide0.clear();
    //        S.TetraFaceSide1.clear();
    //        S.TetraIdxSide0.clear();
    //        S.TetraIdxSide1.clear();

    //        //#ifdef TRACE_DEBUG
    //        assert(S.NodeSide0.size()>0);
    //        assert(S.NodeSide1.size()>0);
    //        assert(S.NodeSide0[0]!=S.NodeSide1[0]);

    //        int IndexPoly0=PolyOfNode(S.SourceNode0[0]);
    //        int IndexPoly1=PolyOfNode(S.SourceNode1[0]);

    //        for (size_t i=1;i<S.SourceNode0.size();i++)
    //        {
    //            assert(IndexPoly0==PolyOfNode(S.SourceNode0[0]));
    //        }

    //        for (size_t i=1;i<S.SourceNode1.size();i++)
    //        {
    //            assert(IndexPoly1==PolyOfNode(S.SourceNode1[0]));
    //        }
    //        //#endif
    //        //set to make further checks

    //        //trace from the two sets contemporanily
    //        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    //        SParam.stop_at_S0=false;
    //        SParam.stop_at_S1=false;
    //        //SParam.AvoidIntersections=true;

    //        std::vector<int> SourceNodes(S.SourceNode0);
    //        SourceNodes.insert(SourceNodes.end(),S.SourceNode1.begin(),S.SourceNode1.end());

    //        std::set<int> SourceNodeSet0(S.SourceNode0.begin(),S.SourceNode0.end());
    //        std::set<int> SourceNodeSet1(S.SourceNode1.begin(),S.SourceNode1.end());
    //        std::set<int> SourceNodeSetAll(SourceNodes.begin(),SourceNodes.end());

    //        //compute distance to nodes
    //        int LastN;
    //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,SourceNodes,SParam,LastN);

    //        //then check all tetra
    //        assert(NextPoly.size()==NextPolyF.size());

    //        int AllPairs=0;
    //        int NoSource=0;
    //        int SameSource=0;
    //        // for each tetra
    //        for (int i=0;i<(int)NextPoly.size();i++)
    //            for (int j=0;j<(int)NextPoly[i].size();j++)
    //            {
    //                int T0=i;
    //                int T0_F=j;
    //                int T1=NextPoly[T0][T0_F];
    //                int T1_F=NextPolyF[T0][T0_F];

    //                assert(T0!=T1);

    //                //if border do nothing
    //                if (T1==-1)continue;

    //                //check only once for each pair
    //                if (T0>T1)continue;

    //                AllPairs++;

    //                //get the nodes
    //                std::vector<int> TetraN;
    //                NodesOfPoly(T0,TetraN);

    //                //get the source of the tetra
    //                int SourceT0=GetSourceTetra(TGraph,T0,SourceNodeSetAll);

    //                //not reached during the diffusion
    //                if (SourceT0==-1){NoSource++;continue;}


    //                //get the source of the other tetra
    //                int SourceT1=GetSourceTetra(TGraph,T1,SourceNodeSetAll);

    //                //not reached during the diffusion
    //                if (SourceT1==-1){NoSource++;continue;}

    //                //see if they belong to different sets
    //                bool IsT0_Set[2];
    //                bool IsT1_Set[2];
    //                IsT0_Set[0]=(SourceNodeSet0.count(SourceT0)>0);
    //                IsT0_Set[1]=(SourceNodeSet1.count(SourceT0)>0);
    //                IsT1_Set[0]=(SourceNodeSet0.count(SourceT1)>0);
    //                IsT1_Set[1]=(SourceNodeSet1.count(SourceT1)>0);

    //                //should belong to at least one set
    //                assert(IsT0_Set[0]||IsT0_Set[1]);
    //                assert(IsT1_Set[0]||IsT1_Set[1]);

    //                //but not at both
    //                assert(!(IsT0_Set[0]&&IsT0_Set[1]));
    //                assert(!(IsT1_Set[0]&&IsT1_Set[1]));

    //                //same set
    //                if (IsT0_Set[0] && IsT1_Set[0]){SameSource++;continue;}
    //                if (IsT0_Set[1] && IsT1_Set[1]){SameSource++;continue;}

    //                //check if they belong to different sets
    //                if (IsT0_Set[0]&&IsT1_Set[1])
    //                {
    //                    S.TetraIdxSide0.push_back(T0);
    //                    S.TetraFaceSide0.push_back(T0_F);
    //                    S.TetraIdxSide1.push_back(T1);
    //                    S.TetraFaceSide1.push_back(T1_F);
    //                }

    //                if (IsT0_Set[1]&&IsT1_Set[0])
    //                {
    //                    S.TetraIdxSide0.push_back(T1);
    //                    S.TetraFaceSide0.push_back(T1_F);
    //                    S.TetraIdxSide1.push_back(T0);
    //                    S.TetraFaceSide1.push_back(T0_F);
    //                }
    //            }
    //        if (WriteDbg)
    //        {
    //            std::cout<<"All Tested Pairs: "<<AllPairs<<std::endl;
    //            std::cout<<"NoSource: "<<NoSource<<std::endl;
    //            std::cout<<"SameSource: "<<SameSource<<std::endl;
    //        }
    //    }

    //    static void CompleteSheet(Geo::TracingGraph<ScalarType> &TGraph,
    //                              const std::vector<Field::FrameField<ScalarType> > &TetraField,
    //                              const std::vector<std::vector<int> > &NextPoly,
    //                              const std::vector<std::vector<int> > &NextPolyF,
    //                              Sheet &S,bool WriteDbg=true)
    //    {

    //        S.TetraIdxSide0.clear();
    //        S.TetraFaceSide0.clear();
    //        S.TetraIdxSide1.clear();
    //        S.TetraFaceSide1.clear();
    //        S.NodeSide0.clear();
    //        S.NodeSide1.clear();

    //        //check
    //        assert(S.SourceNode0.size()>0);
    //        assert(S.SourceNode1.size()>0);
    //        assert(S.IsProper());
    //        assert(S.FromSameTetraSource());


    //        //trace from the two sets contemporanily
    //        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    //        SParam.stop_at_S0=false;
    //        SParam.stop_at_S1=false;
    //        //SParam.AvoidIntersections=true;

    //        std::vector<int> SourceNodes(S.SourceNode0);
    //        SourceNodes.insert(SourceNodes.end(),S.SourceNode1.begin(),S.SourceNode1.end());

    //        //compute distance to nodes
    //        int LastN;
    //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,SourceNodes,SParam,LastN);

    //        //then check all tetra
    //        assert(NextPoly.size()==NextPolyF.size());

    //        int AllPairs=0;
    //        int NoSource=0;
    //        int SameSource=0;
    //        // for each tetra
    //        for (int i=0;i<(int)NextPoly.size();i++)
    //            for (int j=0;j<(int)NextPoly[i].size();j++)
    //            {
    //                int T0=i;
    //                int T0_F=j;
    //                int T1=NextPoly[T0][T0_F];
    //                int T1_F=NextPolyF[T0][T0_F];

    //                assert(T0!=T1);

    //                //if border do nothing
    //                if (T1==-1)continue;

    //                //check only once for each pair
    //                if (T0>T1)continue;

    //                AllPairs++;

    //                int T0Set=-1;
    //                int T1Set=-1;
    //                int IndexN0=-1;
    //                int IndexN1=-1;
    //                SideSheetSeparetrix(TGraph,T0,TetraField[T0],T1,TetraField[T1],S,T0Set,T1Set,IndexN0,IndexN1);

    //                if (T0Set==-1)
    //                {
    //                    assert(T1Set==-1);
    //                    continue;
    //                }

    //                assert((T0Set==0)||(T0Set==1));
    //                assert((T1Set==0)||(T1Set==1));
    //                assert(T0Set!=T1Set);
    //                if (T0Set==0)
    //                {
    //                    S.TetraIdxSide0.push_back(T0);
    //                    S.TetraFaceSide0.push_back(T0_F);
    //                    S.NodeSide0.push_back(IndexN0);
    //                    S.TetraIdxSide1.push_back(T1);
    //                    S.TetraFaceSide1.push_back(T1_F);
    //                    S.NodeSide1.push_back(IndexN1);
    //                }
    //                else
    //                {
    //                    assert(T1Set==0);
    //                    S.TetraIdxSide0.push_back(T1);
    //                    S.TetraFaceSide0.push_back(T1_F);
    //                    S.NodeSide0.push_back(IndexN1);
    //                    S.TetraIdxSide1.push_back(T0);
    //                    S.TetraFaceSide1.push_back(T0_F);
    //                    S.NodeSide1.push_back(IndexN0);
    //                }
    //            }
    //        if (WriteDbg)
    //        {
    //            std::cout<<"All Tested Pairs: "<<AllPairs<<std::endl;
    //            std::cout<<"NoSource: "<<NoSource<<std::endl;
    //            std::cout<<"SameSource: "<<SameSource<<std::endl;
    //        }
    //    }

    static void SelectReachableNodes(Geo::TracingGraph<ScalarType> &TGraph,Sheet &S)
    {
        TGraph.DeSelectAll(0);
        TGraph.DeSelectAll(1);

        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.stop_at_S0=false;
        SParam.stop_at_S1=false;

        //propagate the first one
        int LastN;
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,S.SourceNode0,SParam,LastN);
        for (size_t i=0;i<TGraph.NumNodes();i++)
            if (TGraph.IsMarked(i))
            {
                assert(TGraph.Distance(i)>=0);
                TGraph.Select(i,0);
            }

        //propagate the second one
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,S.SourceNode1,SParam,LastN);
        for (size_t i=0;i<TGraph.NumNodes();i++)
            if (TGraph.IsMarked(i))
            {
                assert(TGraph.Distance(i)>=0);
                TGraph.Select(i,1);
            }
    }


    static void AddFacesAdjTetras(const std::vector<Field::FrameField<ScalarType> > &TetraField,
                                  Geo::TracingGraph<ScalarType> &TGraph,
                                  const std::vector<std::vector<int> > &NextPoly,
                                  const std::vector<std::vector<int> > &NextPolyF,
                                  Sheet &S)
    {
        int AllPairs=0;
        //        int NoSource=0;
        //        int SameSource=0;
        // for each tetra
        for (int i=0;i<(int)NextPoly.size();i++)
            for (int j=0;j<(int)NextPoly[i].size();j++)
            {
                int T0=i;
                int T0_F=j;
                int T1=NextPoly[T0][T0_F];
                int T1_F=NextPolyF[T0][T0_F];

                assert(T0!=T1);

                //if border do nothing
                if (T1==-1)continue;

                //check only once for each pair
                if (T0>T1)continue;

                AllPairs++;

                int T0Set=-1;
                int T1Set=-1;
                int IndexN0=-1;
                int IndexN1=-1;
                SideSheetSeparetrix(TGraph,T0,TetraField[T0],T1,TetraField[T1],S,T0Set,T1Set,IndexN0,IndexN1);

                if (T0Set==-1)
                {
                    assert(T1Set==-1);
                    continue;
                }

                assert((T0Set==0)||(T0Set==1));
                assert((T1Set==0)||(T1Set==1));
                assert(T0Set!=T1Set);
                if (T0Set==0)
                {
                    S.TetraIdxSide0.push_back(T0);
                    S.TetraFaceSide0.push_back(T0_F);
                    S.NodeSide0.push_back(IndexN0);
                    S.TetraIdxSide1.push_back(T1);
                    S.TetraFaceSide1.push_back(T1_F);
                    S.NodeSide1.push_back(IndexN1);
                }
                else
                {
                    assert(T1Set==0);
                    S.TetraIdxSide0.push_back(T1);
                    S.TetraFaceSide0.push_back(T1_F);
                    S.NodeSide0.push_back(IndexN1);
                    S.TetraIdxSide1.push_back(T0);
                    S.TetraFaceSide1.push_back(T0_F);
                    S.NodeSide1.push_back(IndexN0);
                }
            }
    }

    static void FindSeparatricesNodes(const std::vector<Field::FrameField<ScalarType> > &TetraField,
                                      Geo::TracingGraph<ScalarType> &TGraph,
                                      const std::vector<std::vector<int> > &NextPoly,
                                      const std::vector<std::vector<int> > &NextPolyF,
                                      Sheet &S)
    {
        SelectReachableNodes(TGraph,S);

        //trace from the two sets contemporanily
        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.stop_at_S0=false;
        SParam.stop_at_S1=false;
        //SParam.AvoidIntersections=true;

        std::vector<int> SourceNodes(S.SourceNode0);
        SourceNodes.insert(SourceNodes.end(),S.SourceNode1.begin(),S.SourceNode1.end());

        //compute distance to nodes
        int LastN;
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,SourceNodes,SParam,LastN);

        //then check all tetra
        assert(NextPoly.size()==NextPolyF.size());

        //then add the nodes from adj tetras
        AddFacesAdjTetras(TetraField,TGraph,NextPoly,NextPolyF,S);


    }


    //    static void CompleteSheet(const std::vector<Field::FrameField<ScalarType> > &TetraField,
    //                              Geo::TracingGraph<ScalarType> &TGraph,
    //                              const std::vector<std::vector<int> > &NextPoly,
    //                              const std::vector<std::vector<int> > &NextPolyF,
    //                              Sheet &S,bool WriteDbg=true)
    //    {

    //        S.TetraIdxSide0.clear();
    //        S.TetraFaceSide0.clear();
    //        S.TetraIdxSide1.clear();
    //        S.TetraFaceSide1.clear();
    //        S.NodeSide0.clear();
    //        S.NodeSide1.clear();
    //        S.NodeSide0Dist.clear();
    //        S.NodeSide1Dist.clear();
    //        S.PathNodeSide0.clear();
    //        S.PathNodeSide1.clear();

    //        //check
    //        assert(S.SourceNode0.size()>0);
    //        assert(S.SourceNode1.size()>0);
    //        assert(S.IsProper());
    //        assert(S.ProperSources());



    //        //trace from the two sets contemporanily
    //        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    //        SParam.stop_at_S0=false;
    //        SParam.stop_at_S1=false;
    //        //SParam.AvoidIntersections=true;

    //        std::vector<int> SourceNodes(S.SourceNode0);
    //        SourceNodes.insert(SourceNodes.end(),S.SourceNode1.begin(),S.SourceNode1.end());

    //        //compute distance to nodes
    //        int LastN;
    //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,SourceNodes,SParam,LastN);

    //        //then check all tetra
    //        assert(NextPoly.size()==NextPolyF.size());

    //        int AllPairs=0;
    //        int NoSource=0;
    //        int SameSource=0;
    //        // for each tetra
    //        for (int i=0;i<(int)NextPoly.size();i++)
    //            for (int j=0;j<(int)NextPoly[i].size();j++)
    //            {
    //                int T0=i;
    //                int T0_F=j;
    //                int T1=NextPoly[T0][T0_F];
    //                int T1_F=NextPolyF[T0][T0_F];

    //                assert(T0!=T1);

    //                //if border do nothing
    //                if (T1==-1)continue;

    //                //check only once for each pair
    //                if (T0>T1)continue;

    //                AllPairs++;

    //                int T0Set=-1;
    //                int T1Set=-1;
    //                int IndexN0=-1;
    //                int IndexN1=-1;
    //                SideSheetSeparetrix(TGraph,T0,TetraField[T0],T1,TetraField[T1],S,T0Set,T1Set,IndexN0,IndexN1);

    //                if (T0Set==-1)
    //                {
    //                    assert(T1Set==-1);
    //                    continue;
    //                }

    //                assert((T0Set==0)||(T0Set==1));
    //                assert((T1Set==0)||(T1Set==1));
    //                assert(T0Set!=T1Set);
    //                if (T0Set==0)
    //                {
    //                    S.TetraIdxSide0.push_back(T0);
    //                    S.TetraFaceSide0.push_back(T0_F);
    //                    S.NodeSide0.push_back(IndexN0);
    //                    S.TetraIdxSide1.push_back(T1);
    //                    S.TetraFaceSide1.push_back(T1_F);
    //                    S.NodeSide1.push_back(IndexN1);
    //                }
    //                else
    //                {
    //                    assert(T1Set==0);
    //                    S.TetraIdxSide0.push_back(T1);
    //                    S.TetraFaceSide0.push_back(T1_F);
    //                    S.NodeSide0.push_back(IndexN1);
    //                    S.TetraIdxSide1.push_back(T0);
    //                    S.TetraFaceSide1.push_back(T0_F);
    //                    S.NodeSide1.push_back(IndexN0);
    //                }
    //            }
    //        if (WriteDbg)
    //        {
    //            std::cout<<"All Tested Pairs: "<<AllPairs<<std::endl;
    //            std::cout<<"NoSource: "<<NoSource<<std::endl;
    //            std::cout<<"SameSource: "<<SameSource<<std::endl;
    //        }
    //    }


    //    static bool ExtractSheetFrom(Geo::TracingGraph<ScalarType> &TGraph,
    //                                 const std::vector<std::vector<int> > &NextPoly,
    //                                 const std::vector<std::vector<int> > &NextPolyF,
    //                                 const int &IndexT,const int &IndexDir,Sheet &S,
    //                                 bool writeDBG=true)
    //    {
    //        assert(IndexT>=0);
    //        assert(IndexDir>=0);
    //        assert(IndexDir<6);

    //        S.Clear();

    //        int IndexN0=NodeOfPolyFrame(IndexT,IndexDir);
    //        if(!TGraph.EnabledN[IndexN0])return false;

    //        //get oposite node
    //        int IndexN1=-1;
    //        int IndexOppN0=OppositeNodeFrameField(IndexN0);
    //        //then get the closest node
    //        IndexN1=TGraph.FirstActiveNeigh(IndexOppN0);
    //        if (IndexN1==-1)return false;
    //        assert(TGraph.EnabledN[IndexN1]);


    ////        //then get the ortho nodes
    //        OrthoNodesFrameField(IndexN0,S.SourceNode0);
    //        OrthoNodesFrameField(IndexN1,S.SourceNode1);

    //        //filter to have only internals
    //        S.SourceNode0.push_back(IndexN0);
    //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode0);

    //        S.SourceNode1.push_back(IndexN1);
    ////        NodeSet1.push_back(IndexOppN0);
    //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode1);

    //        if (S.SourceNode0.size()==0)return false;
    //        if (S.SourceNode1.size()==0)return false;

    //        if (writeDBG)
    //        {
    //            std::cout<<"Size of Node Set 0:"<<S.SourceNode0.size()<<std::endl;
    //            std::cout<<"Size of Node Set 1:"<<S.SourceNode1.size()<<std::endl;
    //        }
    //        ExtractSheetFrom(TGraph,NextPoly,NextPolyF,S);

    //        return true;
    //    }

    //    static bool ExtractSheetFrom(const std::vector<Field::FrameField<ScalarType> > &TetraField,
    //                                 Geo::TracingGraph<ScalarType> &TGraph,
    //                                 const std::vector<std::vector<int> > &NextPoly,
    //                                 const std::vector<std::vector<int> > &NextPolyF,
    //                                 const int &IndexT,const int &IndexDir,Sheet &S,
    //                                 bool useOrtho=true,bool UseOpposite=false,
    //                                 bool writeDBG=true)
    //    {
    //        assert(IndexT>=0);
    //        assert(IndexDir>=0);
    //        assert(IndexDir<6);
    //        assert(useOrtho||UseOpposite);

    //        if (writeDBG)
    //            std::cout<<"Tracing From Tet "<<IndexT<<" Dir "<<IndexDir<<std::endl;

    //        S.Clear();

    //        //get the current node of the tetra
    //        int IndexN0=NodeOfPolyFrame(IndexT,IndexDir);

    //        if (writeDBG)
    //            std::cout<<"Node0 "<<IndexN0<<std::endl;

    //        //if not enabled return true (maybe to change in assert?)
    //        if(!TGraph.EnabledN[IndexN0])return false;

    //        //get oposite node
    //        int IndexN1=-1;
    //        int IndexOppN0=OppositeNodeFrameField(IndexN0);

    //        //then get the closest node
    //        IndexN1=TGraph.FirstActiveNeigh(IndexOppN0);

    //        //if there is not another opposite return false
    //        if (IndexN1==-1)return false;

    //        if (writeDBG)
    //            std::cout<<"Node1 "<<IndexN1<<std::endl;

    //        assert(TGraph.EnabledN[IndexN1]);

    //        //then get the ortho nodes
    //        if (useOrtho)
    //        {
    //            OrthoNodesFrameField(IndexN0,S.SourceNode0);
    //            OrthoNodesFrameField(IndexN1,S.SourceNode1);
    //        }

    //        //add the opposite
    //        if (UseOpposite)
    //        {
    //            S.SourceNode0.push_back(IndexN0);
    //            S.SourceNode1.push_back(IndexN1);
    //        }

    //        std::cout<<"Size 0:"<<S.SourceNode0.size()<<std::endl;
    //        std::cout<<"Size 1:"<<S.SourceNode1.size()<<std::endl;
    //        //filter to have only internals
    //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode0);
    //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode1);

    //        std::cout<<"Size 0 _ new:"<<S.SourceNode0.size()<<std::endl;
    //        std::cout<<"Size 1 _ new:"<<S.SourceNode1.size()<<std::endl;

    //        //if no node remains no possible extraction
    //        if (S.SourceNode0.size()==0)return false;
    //        if (S.SourceNode1.size()==0)return false;

    //        if (writeDBG)
    //        {
    //            std::cout<<"Size of Node Set 0:"<<S.SourceNode0.size()<<std::endl;
    //            std::cout<<"Size of Node Set 1:"<<S.SourceNode1.size()<<std::endl;
    //        }

    //        CompleteSheet(TGraph,TetraField,NextPoly,NextPolyF,S);
    //        assert(S.NodeSide0.size()==S.NodeSide1.size());

    //        return true;
    //    }

    static bool ExtractSheetFrom(const std::vector<Field::FrameField<ScalarType> > &TetraField,
                                 Geo::TracingGraph<ScalarType> &TGraph,
                                 const std::vector<std::vector<int> > &NextPoly,
                                 const std::vector<std::vector<int> > &NextPolyF,
                                 const int &IndexT,const int &IndexDir,Sheet &S,
                                 bool useOrtho=true,bool UseOpposite=false,
                                 bool writeDBG=true)
    {
        assert(IndexT>=0);
        assert(IndexDir>=0);
        assert(IndexDir<6);
        assert(useOrtho||UseOpposite);

        if (writeDBG)
            std::cout<<"Tracing From Tet "<<IndexT<<" Dir "<<IndexDir<<std::endl;

        S.Clear();

        //get the current node of the tetra
        int IndexN0=NodeOfPolyFrame(IndexT,IndexDir);

        //        if (writeDBG)
        //            std::cout<<"Node0 "<<IndexN0<<std::endl;

        //if not enabled return true (maybe to change in assert?)
        if(!TGraph.EnabledN[IndexN0])return false;

        //get oposite node
        int IndexN1=-1;
        int IndexOppN0=OppositeNodeFrameField(IndexN0);

        //then get the closest node
        IndexN1=TGraph.FirstActiveNeigh(IndexOppN0);

        //if there is not another opposite return false
        if (IndexN1==-1)return false;

        //        if (writeDBG)
        //            std::cout<<"Node1 "<<IndexN1<<std::endl;

        assert(TGraph.EnabledN[IndexN1]);

        //then get the ortho nodes
        if (useOrtho)
        {
            OrthoNodesFrameField(IndexN0,S.SourceNode0);
            OrthoNodesFrameField(IndexN1,S.SourceNode1);
        }

        //add the opposite
        if (UseOpposite)
        {
            S.SourceNode0.push_back(IndexN0);
            S.SourceNode1.push_back(IndexN1);
        }

        //        std::cout<<"Size 0:"<<S.SourceNode0.size()<<std::endl;
        //        std::cout<<"Size 1:"<<S.SourceNode1.size()<<std::endl;
        //filter to have only internals
        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode0);
        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode1);


        //if no node remains no possible extraction
        if (S.SourceNode0.size()==0)return false;
        if (S.SourceNode1.size()==0)return false;

        if (writeDBG)
        {
            std::cout<<"Size of Node Set 0:"<<S.SourceNode0.size()<<std::endl;
            std::cout<<"Size of Node Set 1:"<<S.SourceNode1.size()<<std::endl;
        }

        //CompleteSheet(TGraph,TetraField,NextPoly,NextPolyF,S);
        FindSeparatricesNodes(TetraField,TGraph,NextPoly,NextPolyF,S);
        assert(S.NodeSide0.size()==S.NodeSide1.size());

        return true;
    }

    static void GetDistances(const Geo::TracingGraph<ScalarType> &TGraph,
                             const Sheet&S,std::vector<ScalarType> &DistN)
    {
        DistN.clear();
        for (size_t i=0;i<S.NodeSide0.size();i++)
        {
            int IndexN0=S.NodeSide0[i];
            int IndexN1=S.NodeSide1[i];
            assert(IndexN0>=0);
            assert(IndexN1>=0);
            assert(IndexN0<(int)TGraph.NumNodes());
            assert(IndexN1<(int)TGraph.NumNodes());
            assert(IndexN0!=IndexN1);
            ScalarType D0=TGraph.Distance(IndexN0);
            ScalarType D1=TGraph.Distance(IndexN1);
            ScalarType PathD=PathDistances(D0,D1);
            DistN.push_back(PathD);
        }
    }

    static void RetrievePaths(const Sheet&S,
                              const Geo::TracingGraph<ScalarType> &TGraph,
                              std::vector<std::vector<int> > &PathSide0,
                              std::vector<std::vector<int> > &PathSide1)
    {
        //std::set<int> Ed
        std::set<int> startN0(S.SourceNode0.begin(),S.SourceNode0.end());
        std::set<int> startN1(S.SourceNode1.begin(),S.SourceNode1.end());

        //        std::set<int> startN(S.SourceNode0.begin(),S.SourceNode0.end());;
        //        startN.insert(S.SourceNode1.begin(),S.SourceNode1.end());

        PathSide0.clear();
        PathSide1.clear();

        assert(S.NodeSide0.size()==S.NodeSide1.size());

        PathSide0.resize(S.NodeSide0.size());
        PathSide1.resize(S.NodeSide1.size());

        for (size_t i=0;i<S.NodeSide0.size();i++)
        {
            bool found=Geo::GraphQuery<ScalarType>::RetrievePath(TGraph,S.NodeSide0[i],startN0,PathSide0[i]);
            assert(found);
        }

        for (size_t i=0;i<S.NodeSide1.size();i++)
        {
            bool found=Geo::GraphQuery<ScalarType>::RetrievePath(TGraph,S.NodeSide1[i],startN1,PathSide1[i]);
            assert(found);
        }
    }


    static void CleanSheet(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                           const std::vector<std::vector<int> > &TetraConnectivity,
                           const std::vector<Field::FrameField<ScalarType> > &FrameField,
                           const Geo::TracingGraph<ScalarType> &TGraph,
                           Sheet &S)
    {
        assert(S.IsProper());
        //S.CheckProper();
        //extract the surface
        std::vector<Geo::Point3<ScalarType> > SheetMeshVertPos;
        std::vector<std::vector<int> > SheetMeshConnectivity;
        Tetra::TetraFaceMesh(VertPos,TetraConnectivity,S.TetraIdxSide0,
                             S.TetraFaceSide0,SheetMeshVertPos,
                             SheetMeshConnectivity);

        assert(SheetMeshConnectivity.size()==S.NodeSide0.size());
        assert(SheetMeshConnectivity.size()==S.NodeSide1.size());

        //update the connectivity
        std::vector<std::vector<std::vector<int> > > NextF;
        std::vector<std::vector<std::vector<int> > > NextE;
        Geo::getFFAdjacencyNonManifold(SheetMeshConnectivity,NextF,NextE);

        //compute the distance of each face
        std::vector<ScalarType> FaceDist;
        ScalarType MinD=std::numeric_limits<ScalarType>::max();
        int MinI=-1;
        for (size_t i=0;i<S.NodeSide0.size();i++)
        {
            int IndexN0=S.NodeSide0[i];
            int IndexN1=S.NodeSide1[i];
            assert(IndexN0>=0);
            assert(IndexN1>=0);
            assert(IndexN0<(int)TGraph.NumNodes());
            assert(IndexN1<(int)TGraph.NumNodes());
            assert(IndexN0!=IndexN1);
            ScalarType D0=TGraph.Distance(IndexN0);
            ScalarType D1=TGraph.Distance(IndexN1);
            ScalarType PathD=PathDistances(D0,D1);
            FaceDist.push_back(PathD);
            if (PathD<MinD)
            {
                MinD=PathD;
                MinI=i;
            }
        }

        std::vector<bool> SelectedF(S.TetraFaceSide0.size(),false);

        //create the initial stack
        std::vector<std::pair<ScalarType,int> > Heap(1,std::pair<ScalarType,int>(MinD,MinI));
        do
        {
            std::pop_heap(Heap.begin(),Heap.end());
            int CurrIdx=(Heap.back()).second;
            assert(CurrIdx>=0);
            assert(CurrIdx<(int)SelectedF.size());

            SelectedF[CurrIdx]=true;
            Heap.pop_back();

            int CurrN0=S.NodeSide0[CurrIdx];
            int CurrN1=S.NodeSide1[CurrIdx];

            //get the neigbours
            for (size_t e=0;e<NextF[CurrIdx].size();e++)
            {
                int NextNEdge=-1;
                ScalarType MinD=std::numeric_limits<ScalarType>::max();
                for (size_t i=0;i<NextF[CurrIdx][e].size();i++)
                {
                    int NeighIdx=NextF[CurrIdx][e][i];
                    if (NeighIdx==-1)continue;//this is a border
                    assert(NeighIdx!=CurrIdx);

                    //get the two other nodes
                    int IndexN0=S.NodeSide0[NeighIdx];
                    int IndexN1=S.NodeSide1[NeighIdx];

                    bool AreCoherent0=AreDirectionCoherent(FrameField,CurrN0,IndexN0);
                    bool AreCoherent1=AreDirectionCoherent(FrameField,CurrN1,IndexN1);

                    if (!AreCoherent0)continue;
                    if (!AreCoherent1)continue;

                    //otherwise get the minimum along that direction
                    assert(NeighIdx>=0);
                    assert(NeighIdx<(int)FaceDist.size());
                    ScalarType TestD=FaceDist[NeighIdx];

                    if (SelectedF[NeighIdx])
                    {
                        NextNEdge=-1;
                        break;
                    }
                    if (TestD<MinD)
                    {
                        MinD=TestD;
                        NextNEdge=NeighIdx;
                    }
                }
                if (NextNEdge!=-1)
                {
                    //std::cout<<"DE BOUA"<<std::endl;
                    //in such case already been explored
                    if (SelectedF[NextNEdge])continue;
                    assert(NextNEdge>=0);
                    assert(NextNEdge<(int)FaceDist.size());
                    std::pair<ScalarType,int> Entry(FaceDist[NextNEdge],NextNEdge);
                    Heap.push_back(Entry);
                    std::push_heap(Heap.begin(),Heap.end());
                }
            }
        }while (!Heap.empty());

        //remove the marked faces
        Sheet S_Swap;
        S_Swap.SourceNode0=S.SourceNode0;
        S_Swap.SourceNode1=S.SourceNode1;
        for (size_t i=0;i<SelectedF.size();i++)
        {
            if (!SelectedF[i])continue;
            S_Swap.CopyEntry(S,i);
        }
        S=S_Swap;
    }
};
#endif
