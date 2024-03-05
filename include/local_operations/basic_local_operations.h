#ifndef BASIC_LOCAL_OPERATION
#define BASIC_LOCAL_OPERATION

#include <point3.h>
#include <point2.h>
#include <vector>
#include <pos3.h>
#include <connectivity.h>

namespace Geo{

//template <class ScalarType>
//struct MeshData
//{
//    std::vector<Geo::Point3<ScalarType > > &VertPos;
//    std::vector<std::vector<int> > &Faces;
//    std::vector<Geo::Point3<ScalarType > > &VertNorm;
//    std::vector<Geo::Point3<ScalarType > > &FaceNorm;
//    std::vector<bool> DeletedF,DeletedV;

//    std::vector<std::vector<int> > FFAdj,FEAdj;
//    std::vector<std::vector<int> > VF,VFIndex;

//    void Init()
//    {
//        DeletedF.clear();
//        DeletedV.clear();
//        DeletedV.resize(VertPos.size(),false);
//        DeletedF.resize(DeletedF.size(),false);

//        getVFAdjacency(Faces,VF,VFIndex);
//        getFFAdjacency(Faces,FFAdj,FEAdj);
//    }

//    MeshData(std::vector<Geo::Point3<ScalarType > > &_VertPos,
//             std::vector<std::vector<int> > &_Faces,
//             std::vector<Geo::Point3<ScalarType > > &_VertNorm,
//             std::vector<Geo::Point3<ScalarType > > &_FaceNorm):
//             VertPos(_VertPos),Faces(_Faces),
//             VertNorm(_VertNorm),FaceNorm(_FaceNorm)
//    {}
//};

template <class ScalarType>
class LocalOperation
{
    typedef  Geo::Point3<ScalarType> CoordType;

    std::vector<int> Save_FaceIdx;
    std::vector<std::vector<int> > Save_Face;
    std::vector<int > Save_VertIdx;
    std::vector<CoordType > Save_VertPos;
    int SaveFaceNum;
    int SaveVertNum;

public:

    std::vector<std::vector<int> > &Faces;
    std::vector<CoordType > &Vert;
    std::vector<bool> &DeletedF;
    std::vector<bool> &DeletedV;
    std::vector<std::vector<int> > &VF;
    std::vector<std::vector<int> > &VFIndex;
    std::vector<std::vector<int> > &NextF;
    std::vector<std::vector<int> > &NextE;
    std::vector<CoordType> &VertNorm;
    std::vector<CoordType> &FaceNorm;
    Pos3 pos;

    std::vector<int> ModifiedF;


    virtual bool Feasible()const
    {return true;}

    //    virtual ScalarType UpdatePriority(const Pos3 &P,MeshData<ScalarType> &MData)
    //    {
    //        (void)P;
    //        (void)MData;
    //        return 0;
    //    }

    virtual void Do(bool UpdateAdj=true)
    {(void)UpdateAdj;}

    virtual void GetToUpdate(std::vector<int> &to_upate_faces,
                             std::vector<Pos3> &to_upate_operations)
    {
        (void)to_upate_faces;
        (void)to_upate_operations;
    }

    virtual void GetInvolvedOLDFaces(std::vector<int> &involved_f)const
    {
        (void)involved_f;
    }

    void GetInvolvedNEWFaces(std::vector<int> &involved_f)const
    {
       involved_f=ModifiedF;
    }

    void SaveOriginalStatus(const std::vector<int> &involved_f)
    {
        Save_FaceIdx.clear();
        Save_Face.clear();
        Save_VertPos.clear();
        Save_VertIdx.clear();

        std::set<int> Save_VertSet;
        for (size_t i=0;i<involved_f.size();i++)
        {
            int IndexF=involved_f[i];
            assert(IndexF>=0);
            assert(IndexF<DeletedF.size());
            assert(!DeletedF[IndexF]);
            assert(IndexF<Faces.size());
            assert(IndexF<NextF.size());
            assert(IndexF<NextE.size());

            Save_FaceIdx.push_back(IndexF);
            Save_Face.push_back(Faces[IndexF]);
            for (size_t j=0;j<Faces[IndexF].size();j++)
                Save_VertSet.insert(Faces[IndexF][j]);
        }
        Save_VertIdx=std::vector<int>(Save_VertSet.begin(),Save_VertSet.end());
        for (size_t i=0;i<Save_VertIdx.size();i++)
        {
            int IndexV=Save_VertIdx[i];
            assert(IndexV>=0);
            assert(IndexV<DeletedV.size());
            assert(!DeletedV[IndexV]);
            assert(IndexV<Vert.size());
            assert(IndexV<VF.size());
            assert(IndexV<VFIndex.size());
            Save_VertPos.push_back(Vert[IndexV]);
        }

        SaveVertNum=Vert.size();
        SaveFaceNum=Faces.size();
    }

    void RestoreOriginalStatus()
    {

        Vert.resize(SaveVertNum);
        Faces.resize(SaveFaceNum);

        assert(Save_FaceIdx.size()==Save_Face.size());
        assert(Save_VertIdx.size()==Save_VertPos.size());
        for (size_t i=0;i<Save_VertIdx.size();i++)
        {
            int IndexV=Save_VertIdx[i];
            assert(IndexV>=0);
            assert(IndexV<Vert.size());
            assert(IndexV<DeletedV.size());
            assert(IndexV<VF.size());
            assert(IndexV<VFIndex.size());
            DeletedV[IndexV]=false;
            Vert[IndexV]=Save_VertPos[i];
        }

        for (size_t i=0;i<Save_FaceIdx.size();i++)
        {
            int IndexF=Save_FaceIdx[i];
            assert(IndexF>=0);
            assert(IndexF<DeletedF.size());
            assert(!DeletedF[IndexF]);
            assert(IndexF<Faces.size());
            assert(IndexF<NextF.size());
            assert(IndexF<NextE.size());
            DeletedF[IndexF]=false;
            Faces[IndexF]=Save_Face[i];
        }
    }


    LocalOperation(Pos3 _pos,
                   std::vector<std::vector<int> > &_Faces,
                   std::vector<CoordType > &_Vert,
                   std::vector<bool> &_DeletedV,
                   std::vector<bool> &_DeletedF,
                   std::vector<std::vector<int> > &_VF,
                   std::vector<std::vector<int> > &_VFIndex,
                   std::vector<std::vector<int> > &_NextF,
                   std::vector<std::vector<int> > &_NextE,
                   std::vector<CoordType> &_VertNorm,
                   std::vector<CoordType> &_FaceNorm):
        Faces(_Faces),Vert(_Vert),
        DeletedV(_DeletedV),DeletedF(_DeletedF),
        VF(_VF),VFIndex(_VFIndex),
        NextF(_NextF),NextE(_NextE),
        VertNorm(_VertNorm),FaceNorm(_FaceNorm)
    {
        pos=_pos;
    }

};

}
#endif
