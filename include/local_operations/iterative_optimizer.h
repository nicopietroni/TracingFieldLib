#ifndef TRI_ITERATIVE_OPTIMIZER
#define TRI_ITERATIVE_OPTIMIZER

//#include "../adjacency/adjacency.h"
//#include "../adjacency/vert_pair.h"
//#include "../local_operations/local_operation.h"
//#include "../local_operations/edge_collapse.h"
//#include "../local_operations/vert_smooth.h"
#include "mesh_functions.h"
#include "basic_local_operations.h"
#include "edge_flip.h"
#include "vert_smooth.h"
#include "pos3.h"
#include "mesh_aspect_ratio.h"
#include <iostream>
//#include "../local_operations/crease_finder.h"

namespace Geo{
///*! \brief enumeration to list all possible operation type
//*/
enum OperationType{OpSmooth,OpRefine,OpFlip,OpEdgeCollapse};

/*! \brief class used to implement priority operations
*/
template <class ScalarType>
class HeapEntry
{
    typedef Geo::Point3<ScalarType> CoordType;
    //    typedef HeapEntry<ScalarType> HeapEntryType;

    std::vector<int> involved_faces;
    std::vector<int> *FaceMark;
    int TMark;

    //    /*! \brief return the quality of a set of faces
    //     * \param faces the vector of faces on which the quality have to be computed
    //    */
    //    ScalarType FaceQuality(const std::vector<int> &faces)
    //    {
    //        ScalarType currQ=0;
    //        int num=0;

    //        for (size_t i=0;i<faces.size();i++)
    //        {
    //            int indexF=faces[i];
    //            if ((*LOP->DeletedF)[indexF])continue;
    //            std::vector<CoordType> FaceP;
    //            num++;
    //            for (size_t j=0;j<(*LOP->Faces)[indexF].size();j++)
    //            {
    //                int indexV=(*LOP->Faces)[indexF][j];
    //                FaceP.push_back((*LOP->Vert)[indexV]);
    //            }
    //            currQ+=QualityRadii(FaceP);
    //        }
    //        currQ/=(ScalarType)num;
    //        return currQ;
    //    }

    //    /*! \brief return the worst quality of a set of faces
    //     * \param faces the vector of faces on which the quality have to be computed
    //    */
    //    ScalarType WorstQuality(const std::vector<int> &faces)
    //    {
    //        ScalarType WorstQ=1;
    //        for (size_t i=0;i<faces.size();i++)
    //        {
    //            int indexF=faces[i];
    //            if ((*LOP->DeletedF)[indexF])continue;
    //            std::vector<CoordType> FaceP;
    //            for (size_t j=0;j<(*LOP->Faces)[indexF].size();j++)
    //            {
    //                int indexV=(*LOP->Faces)[indexF][j];
    //                FaceP.push_back((*LOP->Vert)[indexV]);
    //            }
    //            ScalarType TempQ=QualityRadii(FaceP);
    //            if (TempQ<WorstQ)WorstQ=TempQ;
    //        }
    //        return WorstQ;
    //    }

    public:

    OperationType OpType;
    LocalOperation<ScalarType> *LOP;

    ScalarType Priority;

    //    /*! \brief return the worst quality of the involved faces in the operation
    //    */
    //    ScalarType WorstQuality()
    //    {
    //        std::vector<int> faces;
    //        LOP->GetInvolvedFaces(faces);
    //        return WorstQuality(faces);
    //    }

    //    /*! \brief execute the current operation
    //    */
    //    void Execute()
    //    {
    //        ///execute operation
    //        LOP->Do();
    //    }

        /*! \brief return true if the current operation is updated and can be performed
        */
        bool IsUpToDate()
        {
//            int V0=this->LOP->VPair.V0;
//            int V1=this->LOP->VPair.V1;

//            assert(V0>=0);
//            assert(V1>=0);
//            assert(V0<LOP->DeletedV.size());
//            assert(V1<LOP->DeletedV.size());

//            if ((LOP->DeletedV)[V0]||
//                (LOP->DeletedV)[V1])return false;

//            std::vector<int> involved_old_faces;
//            //std::cout<<"0"<<std::endl;
//            LOP->GetInvolvedOLDFaces(involved_old_faces);
            //std::cout<<"1"<<std::endl;
            assert(involved_faces.size()>0);
            for (size_t i=0;i<involved_faces.size();i++)
            {
                int indexF=involved_faces[i];
                assert(indexF>=0);
                assert(indexF<LOP->DeletedF.size());
                assert(indexF<(*FaceMark).size());
                if ((LOP->DeletedF)[indexF])return false;
                if ((*FaceMark)[indexF]>TMark)return false;
            }
            return true;
        }

        /*! \brief ordering to perform higher priority operations first
        */
        inline bool operator <(const HeapEntry &H1)const
        {return (Priority<H1.Priority);}

        /*! \brief update the priority after an operation has been performed
         * \param force_collapse the threshold to force collapse if the quality of the face is below
        */
//        void UpdatePriority(ScalarType force_collapse)
//        {
////            //see if need to force collapse
////            ScalarType AddGain=0;
////            if ((OpType==OpEdgeCollapse)&&(force_collapse>0))
////            {
////                int V0=this->LOP->VPair.V0;
////                int V1=this->LOP->VPair.V1;
////                std::vector<int> edge_faces;
////                LOP->Adj->GetSharedEdgeDataFromV(V0,V1,edge_faces);
////                ScalarType QTest=WorstQuality(edge_faces);
////                if (QTest<force_collapse)
////                {
////                    AddGain=1;
////                }
////            }

////            std::vector<int> faces;
////            LOP->GetInvolvedFaces(faces);

////            ///Esteem Aspect Ratio of Faces
////            ScalarType Q0=FaceQuality(faces);
////            //ScalarType Q0=WorstQuality(faces);

////            LOP->SaveOriginalStatus();

////            LOP->Do(false);

////            LOP->SaveNewStatus();

////            ///gets new faces
////            LOP->GetInvolvedFaces(faces);

////            ///Esteem Aspect Ratio of Faces
////            ScalarType Q1=FaceQuality(faces);
////            //ScalarType Q1=WorstQuality(faces);

////            ///restore state
////            LOP->Undo();

////            //Priority=Q1-Q0;
////            if (Q0>0)
////                Priority=AddGain+(Q1-Q0)/Q0;
////            else
////                Priority=AddGain;
//        }

    void UpdatePriority(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                        const std::vector<std::vector<int> > &Faces)
    {
            //see if need to force collapse
//           ScalarType AddGain=0;
//            if ((OpType==OpEdgeCollapse)&&(force_collapse>0))
//            {
//                int V0=this->LOP->VPair.V0;
//                int V1=this->LOP->VPair.V1;
//                std::vector<int> edge_faces;
//                LOP->Adj->GetSharedEdgeDataFromV(V0,V1,edge_faces);
//                ScalarType QTest=WorstQuality(edge_faces);
//                if (QTest<force_collapse)
//                {
//                    AddGain=1;
//                }
//            }


            //std::vector<int> involved_old_faces;
            LOP->GetInvolvedOLDFaces(involved_faces);

//            ///Esteem Aspect Ratio of Faces
//            ScalarType Q0=FaceQuality(faces);

            //Esteem Aspect Ratio of Faces
            ScalarType Q0=Geo::TriMeshAvgQuality(VertPos,Faces,involved_faces);
//            //ScalarType Q0=WorstQuality(faces);
            
            //std::cout<<"Q0:"<<Q0<<std::endl;
            LOP->SaveOriginalStatus(involved_faces);

//            std::cout<<"test flip A 0"<<std::endl;
//            LOP->pos.FlipV(Faces);
//            std::cout<<"test flip A 1"<<std::endl;

//            std::cout<<"POS before"<<std::endl;
//            std::cout<<"pos V "<<LOP->pos.V<<std::endl;
//            std::cout<<"pos E "<<LOP->pos.E<<std::endl;
//            std::cout<<"pos F "<<LOP->pos.F<<std::endl;

//            std::cout<<"V0 "<<Faces[Lop->pos.F][0]<<std::endl;
//            std::cout<<"V1 "<<Faces[Lop->pos.F][1]<<std::endl;
//            std::cout<<"V2 "<<Faces[Lop->pos.F][2]<<std::endl;

            LOP->Do(false);

//            std::cout<<"POS after"<<std::endl;
//            std::cout<<"pos V "<<LOP->pos.V<<std::endl;
//            std::cout<<"pos E "<<LOP->pos.E<<std::endl;
//            std::cout<<"pos F "<<LOP->pos.F<<std::endl;

//            LOP->SaveNewStatus();

            //gets new faces
            std::vector<int> involved_new_faces;
            LOP->GetInvolvedNEWFaces(involved_new_faces);
            ScalarType Q1=Geo::TriMeshAvgQuality(VertPos,Faces,involved_new_faces);
            //std::cout<<"Q1:"<<Q1<<std::endl;
//            ///Esteem Aspect Ratio of Faces
//            ScalarType Q1=FaceQuality(faces);
//            //ScalarType Q1=WorstQuality(faces);

//            ///restore state
            LOP->RestoreOriginalStatus();

//            std::cout<<"test flip B 0"<<std::endl;
//            LOP->pos.FlipV(Faces);
//            std::cout<<"test flip B 1"<<std::endl;

//            //Priority=Q1-Q0;
//            if (Q0>0)
//                Priority=AddGain+(Q1-Q0)/Q0;
//            else
//                Priority=AddGain;
            Priority=(Q1-Q0)/Q0;
    }


        HeapEntry(std::vector<int> *_FaceMark,
                  int _TMark,
                  LocalOperation<ScalarType> *_LOP,
                  OperationType _OpType)//:FaceMark(_FaceMark)
        {
            FaceMark=_FaceMark;
            TMark=_TMark;
            LOP=_LOP;
            OpType=_OpType;
            Priority=0;
        }

    //    /*! \brief copy constructor
    //     * \param _right object to be copied from
    //    */
        HeapEntry(const HeapEntry &_right)//:FaceMark(_right.FaceMark)
        {
            FaceMark=_right.FaceMark;
            TMark=_right.TMark;
            LOP=_right.LOP;
            Priority=_right.Priority;
            OpType=_right.OpType;
            involved_faces=_right.involved_faces;
        }
};

template <class ScalarType>
class IterativeOptimizer
{
    typedef  Geo::Point3<ScalarType> CoordType;
    typedef  HeapEntry<ScalarType> HEntryType;

    std::vector<HEntryType> HeapOperations;

    bool addSmooth;
    bool addFlip;
    //bool addEdgeColl;

    std::vector<std::vector<int> > &Faces;
    std::vector<CoordType > &Vert;
    std::vector<bool> DeletedV;
    std::vector<bool> DeletedF;
    std::vector<std::vector<int> > VF;
    std::vector<std::vector<int> > VFIndex;
    std::vector<std::vector<int> > NextF;
    std::vector<std::vector<int> > NextE;
    std::vector<CoordType> VertNorm;
    std::vector<CoordType> FaceNorm;

    std::vector<int>  FaceMark;

    //    CreaseFinder<ScalarType> *CreaseF;

    /*! \brief the global temporary mark
    */
    int TMark;

    /*! \brief the number of performed steps
    */
    int performed_steps;

    //    ScalarType last_priority;

    ScalarType min_gain;

    //    ScalarType min_quality;

    //    ScalarType force_collapse;

    //    bool EndConditions()
    //    {
    //        if (performed_steps==0)return true;
    //        return false;
    //    }

    //    void AddOperation(VertPair Vp,OperationType OpType)
    //    {
    //        LocalOperation<ScalarType> *Lop=NULL;
    //        if (OpType==OpSmooth)
    //            Lop= new VertSmoothOperation<ScalarType>(VFAdj,CreaseF,&Faces,&Vert,&DeletedF,&DeletedV,Vp);
    //        if (OpType==OpFlip)
    //            Lop= new EdgeFlipOperation<ScalarType>(VFAdj,CreaseF,&Faces,&Vert,&DeletedF,&DeletedV,Vp);
    //        if (OpType==OpEdgeCollapse)
    //            Lop= new EdgeCollapseOperation<ScalarType>(VFAdj,CreaseF,&Faces,&Vert,&DeletedF,&DeletedV,Vp);

    //        assert(Lop!=NULL);


    //        HEntryType HEntry(&FaceMark,TMark,Lop,OpType);
    //        if (HEntry.WorstQuality()>min_quality)
    //        {
    //            delete(Lop);
    //            return;
    //        }
    //        if (!Lop->Feseable())
    //        {
    //            delete(Lop);
    //            return;
    //        }
    //        //
    //        HEntry.UpdatePriority(force_collapse);

    //        //add only if profitable operation
    //        if (HEntry.Priority<min_gain)
    //        {
    //            delete(Lop);
    //            return;
    //        }

    //        HeapOperations.push_back(HEntry);

    //        //make it as an heap
    //        std::push_heap(HeapOperations.begin(),HeapOperations.end());
    //    }

    void AddOperation(const Pos3 &pos,OperationType OpType)
    {
        Geo::LocalOperation<ScalarType> *Lop=NULL;
        if (OpType==OpSmooth)
            //assert(0);
            Lop= new Geo::VertSmoothOperation<ScalarType>(pos,Faces,Vert,
                                                          DeletedV,DeletedF,
                                                          VF,VFIndex,NextF,NextE,
                                                          VertNorm,FaceNorm);
        if (OpType==OpFlip)
            Lop = new Geo::TriangleEdgeFlipOperation<ScalarType>(pos,Faces,Vert,
                                                                 DeletedV,DeletedF,
                                                                 VF,VFIndex,NextF,NextE,
                                                                 VertNorm,FaceNorm);
        //assert(0);
        //Lop= new EdgeFlipOperation<ScalarType>(VFAdj,CreaseF,&Faces,&Vert,&DeletedF,&DeletedV,Vp);
        if (OpType==OpEdgeCollapse)
            assert(0);
        //Lop= new EdgeCollapseOperation<ScalarType>(VFAdj,CreaseF,&Faces,&Vert,&DeletedF,&DeletedV,Vp);

        assert(Lop!=NULL);


        //HEntryType HEntry(FaceMark,TMark,Lop,OpType);
        HEntryType HEntry(&FaceMark,TMark,Lop,OpType);
        //        if (HEntry.WorstQuality()>min_quality)
        //        {
        //            delete(Lop);
        //            return;
        //        }

        if (!Lop->Feasible())
        {
            delete(Lop);
            return;
        }

//
//        std::cout<<"test flip A 0"<<std::endl;
//        HEntry.LOP->pos.FlipV(Faces);
//        std::cout<<"test flip A 1"<<std::endl;

        HEntry.UpdatePriority(Vert,Faces);//force_collapse);

//        std::cout<<"test flip B 0"<<std::endl;
//        HEntry.LOP->pos.FlipV(Faces);
//        std::cout<<"test flip B 1"<<std::endl;

        //add only if profitable operation
        if (HEntry.Priority<min_gain)
        {
            delete(Lop);
            return;
        }


        HeapOperations.push_back(HEntry);

//        //make it as an heap
        std::push_heap(HeapOperations.begin(),HeapOperations.end());
    }

    void AddOperation(const Pos3 &pos)
    {

        if (addSmooth)   AddOperation(pos,OpSmooth);
        if (addFlip)     AddOperation(pos,OpFlip);
        //if (addEdgeColl) AddOperation(Vp,OpEdgeCollapse);
    }

    void UpdateMarks(const std::vector<int> &to_update_mark)
    {
        for (size_t i=0;i<to_update_mark.size();i++)
        {
            //update temporary mark
            int currF=to_update_mark[i];
            assert(currF>=0);
            assert((size_t)currF<FaceMark.size());
            FaceMark[currF]=TMark;
        }
    }

//    void UpdateHeap(std::vector<Pos3> &Operations)
//    {
//        for (size_t i=0;i<Operations.size();i++)
//            AddOperation(Operations[i]);
//    }

    void UpdateHeap(const std::vector<int> &involved_faces)
    {
        for (size_t i=0;i<involved_faces.size();i++)
        {
            int IndexF=involved_faces[i];
            assert(IndexF>=0);
            assert(IndexF<Faces.size());
            assert(IndexF<DeletedF.size());
            assert(!DeletedF[IndexF]);
            assert(IndexF<NextF.size());
            assert(IndexF<NextE.size());
            for (size_t j=0;j<Faces[IndexF].size();j++)
            {
                int IndexE=j;
                int IndexV=Faces[IndexF][j];
                assert(IndexV>=0);
                assert(IndexV<Vert.size());
                assert(IndexV<DeletedV.size());
                assert(!DeletedV[IndexV]);
                Pos3 NexOp(IndexF,IndexE,IndexV);
                AddOperation(NexOp);
            }
        }
    }

    void InitSupportStructures()
    {
        Geo::getFFAdjacency(Faces,NextF,NextE);
        Geo::getVFAdjacency(Faces,VF,VFIndex);

        //        if (VFAdj!=NULL)
        //            delete(VFAdj);


        //        VFAdj= new VFAdjacencyStructure();
        //        VFAdj->Init(Vert.size(),Faces);

        //        if (CreaseF!=NULL)
        //            delete(CreaseF);

        //        CreaseF=new CreaseFinder<ScalarType>(Vert,Faces,Components,*VFAdj);
        //        CreaseF->Init();

    }


    //    /*! \brief initialize the whole heap of operations
    //    */
    //    void InitHeap()
    //    {
    //        ///initialize the heap
    //        HeapOperations.clear();

    //        DeletedV.clear();
    //        DeletedF.clear();
    //        FaceMark.clear();

    //        ///initialize the rest of structures
    //        DeletedV.resize(Vert.size(),false);
    //        DeletedF.resize(Faces.size(),false);
    //        FaceMark.resize(Faces.size(),0);

    //        TMark=0;

    //        std::set<VertPair> inserted;

    //        ///with all possible operations
    //        for (size_t i=0;i<Faces.size();i++)
    //        {
    //            for (size_t j=0;j<Faces[i].size();j++)
    //            {
    //                int size=Faces[i].size();
    //                int IndexV0=Faces[i][j];
    //                int IndexV1=Faces[i][(j+1) % size];
    //                VertPair Vp(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
    //                if (inserted.count(Vp)!=0)continue;
    //                AddOperation(Vp);
    //                inserted.insert(Vp);
    //            }
    //        }
    //        //then check
    //        for (size_t i=0;i<HeapOperations.size();i++)
    //        {
    //            assert(HeapOperations[i].Priority>=min_gain);
    //        }
    //        printf("Initializated heap of size %d \n",(int)HeapOperations.size());
    //        fflush(stdout);

    //    }

    void InitHeap()
    {
        ///initialize the heap
        HeapOperations.clear();

        DeletedV.clear();
        DeletedF.clear();
        FaceMark.clear();

        ///initialize the rest of structures
        DeletedV.resize(Vert.size(),false);
        DeletedF.resize(Faces.size(),false);
        FaceMark.resize(Faces.size(),0);

        TMark=0;

        std::set<std::pair<int,int> > insertedE;
        //std::set<int> insertedV;

        //with all possible operations
        for (size_t i=0;i<Faces.size();i++)
        {
            for (size_t j=0;j<Faces[i].size();j++)
            {
                int size=Faces[i].size();
                int IndexV0=Faces[i][j];
                // int IndexV1=Faces[i][(j+1) % size];
                // std::pair<int,int> Key(std::min(IndexV0,IndexV1),
                //                        std::max(IndexV0,IndexV1));
                // //already inserted in this case
                // if (insertedE.count(Key)>00)continue;
                Pos3 Pos(i,j,IndexV0);
                AddOperation(Pos);
                //insertedE.insert(Key);
            }
        }
        //then check
        for (size_t i=0;i<HeapOperations.size();i++)
        {
            assert(HeapOperations[i].Priority>=min_gain);
        }
        printf("Initializated heap of size %d \n",(int)HeapOperations.size());
        fflush(stdout);

    }

public:

    //    /*! \brief initialize the optimizer
    //     * \param _force_collapse the threshold of quality when a face is forced to collapse independently of priority queue
    //     * \param _min_gain the minimum gain when an operation is considered profitable, the lowest the better the result but become slower
    //     * \param _min_quality the minimum quality of a face if is over this value the face is considered to be ok, the highest the better the result but may become really slow
    //     * \param _addSmooth use smooth operations within the framework
    //     * \param _addFlip use flip operations within the framework
    //     * \param _addEdgeColl use edge collapse operations within the framework
    //    */
    //    void Init(ScalarType _force_collapse=0.05,
    //              ScalarType _min_gain=0.001,
    //              ScalarType _min_quality=0.6,
    //              bool _addSmooth=false,
    //              bool _addFlip=true,
    //              bool _addEdgeColl=true)
    //    {

    //        InitSupportStructures();

    //        force_collapse=_force_collapse;
    //        min_gain=_min_gain;
    //        min_quality=_min_quality;
    //        addSmooth=_addSmooth;
    //        addFlip=_addFlip;
    //        addEdgeColl=_addEdgeColl;

    //        CreaseF->InitCorners();
    //        CreaseF->Smooth();
    //        CreaseF->ClearCorners();

    //        //        InitHeap();
    //    }


    void Init(ScalarType _min_gain=0.001,
              bool _addSmooth=true,
              bool _addFlip=true,
              bool _addEdgeColl=true)
    {

        InitSupportStructures();

        min_gain=_min_gain;
        addSmooth=_addSmooth;
        addFlip=_addFlip;
        //addEdgeColl=_addEdgeColl;

        //        CreaseF->InitCorners();
        //        CreaseF->Smooth();
        //        CreaseF->ClearCorners();
        std::vector<std::vector<int> > Faces0=Faces;
        std::vector<std::vector<int> > VF0=VF;
        std::vector<std::vector<int> > VFIndex0=VFIndex;
        std::vector<std::vector<int> > NextF0=NextF;
        std::vector<std::vector<int> > NextE0=NextE;
        InitHeap();
        assert(Faces0==Faces);
        assert(VF0==VF);
        assert(VFIndex0==VFIndex);
        assert(NextF0==NextF);
        assert(NextE0==NextE);
    }

//    bool OptimizeStep()
//    {

//        performed_steps=0;
//        int print_steps=1000;
//        int not_updated=0;
//        //        int bad_priority=0;
//        //        bad_priority=1;

//        //       int smoothDone=0;
//        int flipDone=0;
//        //       int collapseDone=0;
//        //int tot_ite=0;

//        ScalarType last_priority=0;

//        if (HeapOperations.size()==0)return false;
//        do{
//            std::pop_heap(HeapOperations.begin(),HeapOperations.end());

//            Geo::LocalOperation<ScalarType> *Lop=HeapOperations.back().LOP;

//            bool isUptoDate=HeapOperations.back().IsUpToDate();
//            assert(HeapOperations.back().Priority>=min_gain);
//            //bool notPrio=HeapOperations.back().Priority<0;
//            if (!isUptoDate) {not_updated++;
//            //if (notPrio)bad_priority++;

//            bool modified=false;

//            if ((isUptoDate)&&
//                    (!notPrio)/*&&
//                                (Lop->Feseable())*/)
//            {
//                Lop->Do();

//                performed_steps++;
//                TMark++;
//                modified=true;

//                last_priority=HeapOperations.back().Priority;
//                //if (HeapOperations.back().OpType==OpSmooth)smoothDone++;
//                if (HeapOperations.back().OpType==OpFlip)flipDone++;
//                //if (HeapOperations.back().OpType==OpEdgeCollapse)collapseDone++;
//            }

//            HeapOperations.pop_back();

//            if (modified)
//            {
//                std::vector<int> to_upate_faces;
//                std::vector<Pos3> to_upate_operations;
//                Lop->GetToUpdate(to_upate_faces,to_upate_operations);
//                UpdateMarks(to_upate_faces);
//                UpdateHeap(to_upate_operations);
//            }

//            //erase current operation
//            delete(Lop);

//            if ((modified)&&((performed_steps%print_steps)==0))
//                printf("Heap Size %d \n",(int)HeapOperations.size());

//        }while ((!HeapOperations.empty())&&(last_priority>=min_gain));

//        assert(HeapOperations.empty());
//        printf("Performed %d Optimization Steps \n",performed_steps);
//        //printf("Smooth Steps %d \n",smoothDone);
//        printf("Flip Steps  %d \n",flipDone);
//        //printf("CollapseSteps %d \n",collapseDone);

//        //printf("Bad Priority %d \n",bad_priority);
//        printf("Not Updated  %d \n",not_updated);
//        fflush(stdout);

//        //Compress(Vert,Faces,Components,DeletedF,DeletedV);
//        RemoveDeleted(Vert,Faces,DeletedF,DeletedV);

//        InitSupportStructures();
//        return (performed_steps>0);
//    }

    bool OptimizeStep()
    {

        performed_steps=0;
        int print_steps=1000;
        int not_updated=0;
        //        int bad_priority=0;
        //        bad_priority=1;
        //       int smoothDone=0;
        int flipDone=0;
        //       int collapseDone=0;
        //int tot_ite=0;

        //ScalarType last_priority=0;

        if (HeapOperations.size()==0)return false;
        do{
            std::pop_heap(HeapOperations.begin(),HeapOperations.end());

            //std::cout<<"performed steps"<<performed_steps<<std::endl;

            Geo::LocalOperation<ScalarType> *Lop=HeapOperations.back().LOP;
//            std::cout<<"check is update"<<std::endl;
//            std::cout<<"pos V "<<Lop->pos.V<<std::endl;
//            std::cout<<"pos E "<<Lop->pos.E<<std::endl;
//            std::cout<<"pos F "<<Lop->pos.F<<std::endl;

//            std::cout<<"V0 "<<Faces[Lop->pos.F][0]<<std::endl;
//            std::cout<<"V1 "<<Faces[Lop->pos.F][1]<<std::endl;
//            std::cout<<"V2 "<<Faces[Lop->pos.F][2]<<std::endl;

            bool isUptoDate=HeapOperations.back().IsUpToDate();

            assert(HeapOperations.back().Priority>=min_gain);
            //bool notPrio=HeapOperations.back().Priority<0;
            //if (!isUptoDate) {not_updated++;continue;}
            //if (notPrio)bad_priority++;

            bool has_performed=false;
            if ((isUptoDate)&&(Lop->Feasible()))
            {
                //std::cout<<"check feasible"<<std::endl;
                //assert(Lop->Feasible());
                //std::cout<<"doing"<<std::endl;
                Lop->Do();

                performed_steps++;
                TMark++;
                has_performed=true;

                //last_priority=HeapOperations.back().Priority;
                //if (HeapOperations.back().OpType==OpSmooth)smoothDone++;
                if (HeapOperations.back().OpType==OpFlip)flipDone++;
                //if (HeapOperations.back().OpType==OpEdgeCollapse)collapseDone++;
            }

            HeapOperations.pop_back();

            if (has_performed)
            {
                //std::vector<int> to_upate_faces;
                //std::vector<Pos3> to_upate_operations;
                //Lop->GetToUpdate(to_upate_faces,to_upate_operations);

                std::vector<int> involved_new_faces;
                Lop->GetInvolvedNEWFaces(involved_new_faces);
                UpdateMarks(involved_new_faces);
                UpdateHeap(involved_new_faces);
            }

            //erase current operation
            delete(Lop);

            if ((has_performed)&&((performed_steps%print_steps)==0))
            {
                printf("Heap Size %d \n",(int)HeapOperations.size());
            }

        }while (!HeapOperations.empty());//&&(last_priority>=min_gain));

        assert(HeapOperations.empty());
        printf("Performed %d Optimization Steps \n",performed_steps);
        //printf("Smooth Steps %d \n",smoothDone);
        printf("Flip Steps  %d \n",flipDone);
        //printf("CollapseSteps %d \n",collapseDone);

        //printf("Bad Priority %d \n",bad_priority);
        printf("Not Updated  %d \n",not_updated);
        fflush(stdout);

        //Compress(Vert,Faces,Components,DeletedF,DeletedV);
//        std::cout<<"Size 0 "<<DeletedV.size()<<std::endl;
//        std::cout<<"Size 1 "<<Vert.size()<<std::endl;

        RemoveDeleted(Vert,Faces,DeletedV,DeletedF);

        InitSupportStructures();
        return (performed_steps>0);
    }

    /*! \brief perform the whole optimization
     */
    void Optimize()
    {
        //InitHeap();
        while (OptimizeStep())
        {
            //InitHeap();
        }
    }

    IterativeOptimizer(std::vector<CoordType> &_Vert,
                       std::vector<std::vector<int> > &_Faces):Vert(_Vert),Faces(_Faces)
    {
    }

    ~IterativeOptimizer()
    {
        for (size_t i=0;i<HeapOperations.size();i++)
            delete(HeapOperations[i].LOP);
    }
};
}
#endif 
