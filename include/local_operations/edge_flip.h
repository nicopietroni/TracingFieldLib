#ifndef TRI_EDGE_FLIP
#define TRI_EDGE_FLIP

#include "basic_local_operations.h"
#include <connectivity.h>
#include <mesh_functions.h>
#include <triangle3.h>

namespace  Geo{

template <class ScalarType>
class TriangleEdgeFlip
{
    typedef Geo::Point3<ScalarType> CoordType;

    /*! \brief Return the entities involved in the operation
    *
    *     V1 ___________ V3
    *       |\     F1   |
    *       |  \        |
    *       |    \      |
    *       |      \    |
    *       |  F0    \  |
    *       |__________\|
    *     V2            V0
    * \param F0 the face index as in shown in the ascii art
    * \param F1 the face index as in shown in the ascii art
    * \param V0 the vert index as in shown in the ascii art
    * \param V1 the vert index as in shown in the ascii art
    * \param V2 the vert index as in shown in the ascii art
    * \param V3 the vert index as in shown in the ascii art
    */
    static void FindSets(const Pos3 &FlipP,
                         const std::vector<std::vector<int> > &Faces,
                         const std::vector<std::vector<int> > &NextF,
                         const std::vector<std::vector<int> > &NextEd,
                         int &F0,int &F1,
                         int &V0,int &V1,
                         int &V2,int &V3)
    {
        V0=FlipP.V;
        V1=FlipP.VFlip(Faces);
        assert(V0!=V1);

        Pos3 TestP=FlipP;
        assert(!TestP.IsBorder(NextF));
        F0=TestP.F;

        Pos3 OppP=TestP;
        OppP.FlipF(NextF,NextEd);
        F1=OppP.F;

        TestP.FlipE(Faces);
        TestP.FlipV(Faces);
        V2=TestP.V;

        OppP.FlipE(Faces);
        OppP.FlipV(Faces);
        V3=OppP.V;

        assert((V0!=V1)&&(V1!=V2)&&(V2!=V3)&&
               (V1!=V2)&&(V1!=V3)&&(V2!=V3));
    }

    //    /*! \brief dot product for max normal deviation to avoid invert normal */
    //    ScalarType maxNormDev;

public:

    /*! \brief check if the operation is feseable */
    static bool Feasible(const Pos3 &FlipP,
                         const std::vector<Geo::Point3<ScalarType> > &Pos,
                         const std::vector<std::vector<int> > &Faces,
                         const std::vector<std::vector<int> > &NextF,
                         const std::vector<std::vector<int> > &NextE,
                         std::vector<std::vector<int> > &VF,
                         ScalarType maxAngleDev=45)
    {
        /*! \brief this function perform the actual edge flip
        *
        *     V1 ___________ V3         V1 ___________ V3
        *       |\     F1   |             |          / |
        *       |  \        |             |   F1   /   |
        *       |    \      |   ----->    |      /     |
        *       |      \    |             |    /  F0   |
        *       |  F0    \  |             |  /         |
        *       |__________\|             |/___________|
        *     V2            V0          V2            V0
        * \param UpdateFeatures not used in this context
        */

        assert(FlipP.E>=0);
        assert(FlipP.E<3);
        assert(FlipP.V>=0);
        assert(FlipP.V<Pos.size());
        assert(FlipP.F>=0);
        assert(FlipP.F<Faces.size());
        assert(NextF.size()==Faces.size());
        assert(NextE.size()==Faces.size());
        assert(VF.size()==Pos.size());

        if (FlipP.IsBorder(NextF))return false;

        //check if the edge appearing after flip already exists or not
        int F0,F1,V0,V1,V2,V3;
        FindSets(FlipP,Faces,NextF,
                 NextE,F0,F1,
                 V0,V1,V2,V3);

        int NumShar=NumSharedFaces(V2,V3,VF);
        if (NumShar!=0)return false;//already exist the edge

        assert(V0>=0);
        assert(V1>=0);
        assert(V2>=0);
        assert(V3>=0);
        assert(V0<Pos.size());
        assert(V1<Pos.size());
        assert(V2<Pos.size());
        assert(V3<Pos.size());

        CoordType P0=Pos[V0];
        CoordType P1=Pos[V1];
        CoordType P2=Pos[V2];
        CoordType P3=Pos[V3];

        CoordType NF0=Geo::Normal(P0,P1,P2);
        CoordType NF1=Geo::Normal(P0,P3,P1);
        ScalarType AngleN=(Geo::AngleDeg(NF0,NF1));
        if (AngleN>maxAngleDev)return false;

        //check collinear faces
        ScalarType delta=0.00001;
        ScalarType Angle0=(Geo::AngleRad(P2 - P0, P1 - P0) + Geo::AngleRad(P3 - P0, P1 - P0));
        if (Angle0>= (M_PI-delta))return false;
        ScalarType Angle1=(Geo::AngleRad(P2 - P1, P0 - P1) + Geo::AngleRad(P3 - P1, P0 - P1));
        if (Angle1>= (M_PI-delta))return false;

        CoordType NF0_new=Geo::Normal(P0,P3,P2);
        CoordType NF1_new=Geo::Normal(P2,P3,P1);
        if (Geo::Dot(NF0,NF0_new)<0)return false;
        if (Geo::Dot(NF1,NF0_new)<0)return false;
        if (Geo::Dot(NF0,NF1_new)<0)return false;
        if (Geo::Dot(NF1,NF1_new)<0)return false;

        return true;
    }


    static void Kernel(const Pos3 &FlipP,
                       const std::vector<std::vector<int> > &NextF,
                       const std::vector<std::vector<int> > &NextE,
                       std::vector<int> &KernelF)
    {
        //check if the edge appearing after flip already exists or not
        int F0,F1;
        F0=FlipP.F;
        F1=FlipP.FFlip(NextF,NextE);
        KernelF.clear();
        KernelF.push_back(F0);
        KernelF.push_back(F1);

    }

    static void GetInvolvedFaces(const Pos3 &FlipP,
                                 const std::vector<Geo::Point3<ScalarType> > &Pos,
                                 const std::vector<std::vector<int> > &Faces,
                                 const std::vector<std::vector<int> > &NextF,
                                 const std::vector<std::vector<int> > &NextEd,
                                 std::vector<int > &InvolvedF)
    {
        int F0,F1,V0,V1,V2,V3;
        FindSets(FlipP,Faces,NextF,NextEd,F0,F1,V0,V1,V2,V3);
        InvolvedF.clear();
        InvolvedF.push_back(F0);
        InvolvedF.push_back(F1);
    }

    //    static void FacesBeforeAfter(const Pos3 &FlipP,
    //                                 const std::vector<Geo::Point3<ScalarType> > &Pos,
    //                                const std::vector<std::vector<int> > &Faces,
    //                                const std::vector<std::vector<int> > &NextF,
    //                                const std::vector<std::vector<int> > &NextEd,
    //                                std::vector<std::vector<int> > &VF,
    //                                std::vector<std::vector<Geo::Traingle3<ScalarType> > > &T_before,
    //                                std::vector<std::vector<Geo::Traingle3<ScalarType> > > &T_after)
    //    {
    //        /*! \brief this function perform the actual edge flip
    //        *
    //        *     V1 ___________ V3         V1 ___________ V3
    //        *       |\     F1   |             |          / |
    //        *       |  \        |             |   F1   /   |
    //        *       |    \      |   ----->    |      /     |
    //        *       |      \    |             |    /  F0   |
    //        *       |  F0    \  |             |  /         |
    //        *       |__________\|             |/___________|
    //        *     V2            V0          V2            V0
    //        * \param UpdateFeatures not used in this context
    //        */
    //        assert(Feaseable(FlipP,Pos,Faces,NextF,NextEd,VF,180));

    //        //check if the edge appearing after flip already exists or not
    //        int F0,F1,V0,V1,V2,V3;
    //        FindSets(FlipP,Faces,NextF,NextEd,F0,F1,V0,V1,V2,V3);

    //        CoordType P0=Pos[V0];
    //        CoordType P1=Pos[V1];
    //        CoordType P2=Pos[V2];
    //        CoordType P3=Pos[V3];

    //        Geo::Traingle3<ScalarType> T0_before(P0,P1,P2);
    //        Geo::Traingle3<ScalarType> T1_before(P0,P3,P1);

    //        Geo::Traingle3<ScalarType> T0_after(P0,P3,P2);
    //        Geo::Traingle3<ScalarType> T1_after(P2,P3,P1);

    //        T_before.clear();
    //        T_after.clear();

    //        T_before.push_back(T0_before);
    //        T_before.push_back(T1_before);

    //        T_after.push_back(T0_after);
    //        T_after.push_back(T1_after);
    //    }

    /*! \brief this function perform the actual edge flip
    *
    *     V1 ___________ V3         V1 ___________ V3
    *       |\     F1   |             |          / |
    *       |  \        |             |   F1   /   |
    *       |    \      |   ----->    |      /     |
    *       |      \    |             |    /  F0   |
    *       |  F0    \  |             |  /         |
    *       |__________\|             |/___________|
    *     V2            V0          V2            V0
    * \param UpdateFeatures not used in this context
    */
    static void Do(const Pos3 &FlipP,
                   std::vector<std::vector<int> > &Faces,
                   std::vector<std::vector<int> > &NextF,
                   std::vector<std::vector<int> > &NextE,
                   std::vector<std::vector<int> > &VF,
                   std::vector<std::vector<int> > &VFIndex,
                   bool UpdateAdj=true)
    {
        int F0,F1,V0,V1,V2,V3;

        FindSets(FlipP,Faces,NextF,NextE,
                 F0,F1,V0,V1,V2,V3);

        Faces[F0][0]=V0;
        Faces[F0][1]=V3;
        Faces[F0][2]=V2;

        Faces[F1][0]=V2;
        Faces[F1][1]=V3;
        Faces[F1][2]=V1;

        if (UpdateAdj)
        {
            std::vector<int> To_UpdF;
            To_UpdF.push_back(F0);
            To_UpdF.push_back(F1);
            UpdateVF(Faces,To_UpdF,VF,VFIndex);
            UpdateFF(Faces,To_UpdF,NextF,NextE);
        }

//        FlipP.E=1;
//        FlipP.V=V3;
    }

};


template <class ScalarType>
class TriangleEdgeFlipOperation:public Geo::LocalOperation<ScalarType>
{
    typedef  Geo::Point3<ScalarType> CoordType;

public:

    bool Feasible()const override
    {
        return (TriangleEdgeFlip<ScalarType>::Feasible(this->pos,this->Vert,
                                                       this->Faces,this->NextF,
                                                       this->NextE,this->VF));
    }

    virtual void Do(bool UpdateAdj=true) override
    {
        this->ModifiedF.clear();
        this->ModifiedF.push_back(this->pos.F);
        this->ModifiedF.push_back(this->pos.FFlip(this->NextF,this->NextE));

        TriangleEdgeFlip<ScalarType>::Do(this->pos,this->Faces,this->NextF,
                                         this->NextE,this->VF,this->VFIndex,UpdateAdj);
    }

    void GetInvolvedOLDFaces(std::vector<int> &involved_f)const override
    {
        TriangleEdgeFlip<ScalarType>::GetInvolvedFaces(this->pos,this->Vert,
                                                       this->Faces,this->NextF,
                                                       this->NextE,
                                                       involved_f);
    }

//    void GetInvolvedNEWFaces(std::vector<int> &involved_f)const override
//    {
//        GetInvolvedNEWFaces(std::vector<int> &involved_f)const
////        TriangleEdgeFlip<ScalarType>::GetInvolvedFaces(this->pos,this->Vert,
////                                                       this->Faces,this->NextF,
////                                                       this->NextE,
////                                                       involved_f);
//    }

    TriangleEdgeFlipOperation(Pos3 _pos,
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
        Geo::LocalOperation<ScalarType>(_pos,_Faces,_Vert,
                                        _DeletedV,_DeletedF,
                                        _VF,_VFIndex,
                                        _NextF,_NextE,
                                        _VertNorm,_FaceNorm)
    {}
};
}

#endif 
