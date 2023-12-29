#ifndef UV_PARAM
#define UV_PARAM

#include <point2.h>
#include <point3.h>
#include <igl/lscm.h>
#include <igl/arap.h>
#include <igl/harmonic.h>
#include <eigen_interface.h>
#include <box3.h>
#include <box2.h>
#include <mesh_functions.h>
#include <connectivity.h>
#include <UV_utils.h>

namespace Geo{

template <class ScalarType >
class UVParam{

    static int SelectBestConstraintOrientedAlong(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                                 const std::vector<int> &Candidates,
                                                 const int &Vert0)
    {
        assert(Vert0>=0);
        assert(Vert0<VertPos.size());
        //        Geo::Point3<ScalarType> DirTarget=Geo::Point3<ScalarType>(1,0,0),
        //        Geo::Point3<ScalarType> DirOrtho=Geo::Point3<ScalarType>(0,1,0)
        Geo::Box3<ScalarType> MeshBB;
        MeshBB.Init(VertPos);

        //get two directions, better to do with PCA in the future
        std::vector<std::pair<ScalarType,Geo::Point3<ScalarType> > > DimDir;
        std::pair<ScalarType,Geo::Point3<ScalarType> > P0(MeshBB.DimX(),Geo::Point3<ScalarType>(1,0,0));
        std::pair<ScalarType,Geo::Point3<ScalarType> > P1(MeshBB.DimY(),Geo::Point3<ScalarType>(0,1,0));
        std::pair<ScalarType,Geo::Point3<ScalarType> > P2(MeshBB.DimZ(),Geo::Point3<ScalarType>(0,0,1));
        DimDir.push_back(P0);
        DimDir.push_back(P1);
        DimDir.push_back(P2);
        std::sort(DimDir.begin(),DimDir.end());
        Geo::Point3<ScalarType> DirTarget=DimDir[2].second;
        Geo::Point3<ScalarType> DirOrtho=DimDir[1].second;

        ScalarType BestDot=std::numeric_limits<double>::lowest();;
        int Vert1=-1;

        for (int i=0;i<(int)Candidates.size();i++)
        {
            int IndexVTest=Candidates[i];
            if (IndexVTest==Vert0)continue;

            Geo::Point3<ScalarType> DirTest(VertPos[IndexVTest]-VertPos[Vert0]);
            ScalarType AlignDot=Dot(DirTest,DirTarget);
            ScalarType OrthoDot=fabs(Geo::Dot(DirTest,DirOrtho));
            ScalarType TestDot=AlignDot-OrthoDot;

            if (TestDot<BestDot)continue;

            BestDot=TestDot;
            Vert1=IndexVTest;
        }
        assert(Vert1>=0);
        return Vert1;
    }

    //    static int ClosestVert(const std::vector<Geo::Point3<ScalarType> > &VertPos,
    //                    const Geo::Point3<ScalarType> &TargetPos)
    //    {
    //        ScalarType DistMin=std::numeric_limits<ScalarType>::max();
    //        int IndexV=-1;
    //        for (size_t i=0;i<VertPos.size();i++)
    //        {
    //            ScalarType TestDist=(VertPos[i]-TargetPos).Norm();
    //            if (TestDist>DistMin)continue;
    //            DistMin=TestDist;
    //            IndexV=i;
    //        }
    //        return IndexV;
    //    }

    static void GetDefaultPairVertexToFix_AlongX(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                                 std::vector<size_t> &FixedV,
                                                 std::vector<Geo::Point2<ScalarType> > &FixedUV)
    {

        Geo::Box3<ScalarType> BB;
        BB.Init(VertPos);

        int Vert0=ClosestVertFromVect<ScalarType>(VertPos,BB.Center());
        std::vector<int> Candidates;
        for (int i=0;i<(int)VertPos.size();i++)
            Candidates.push_back(i);

        int Vert1=SelectBestConstraintOrientedAlong(VertPos,Candidates,Vert0);

        ScalarType Dist=(VertPos[Vert0]-VertPos[Vert1]).Norm();
        //the conformal is folded so check this then
        FixedV.push_back(Vert0);
        FixedUV.push_back(Geo::Point2<ScalarType>(0,0));
        FixedV.push_back(Vert1);
        FixedUV.push_back(Geo::Point2<ScalarType>(Dist,0));
    }

    static void GetDefaultEdgeToFix(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                    const std::vector<std::vector<int> > &Faces,
                                    std::vector<size_t> &FixedV,
                                    std::vector<Geo::Point2<ScalarType> > &FixedUV)
    {

        Geo::Box3<ScalarType> BB;
        BB.Init(VertPos);

        int Vert0=ClosestVert(VertPos,BB.Center());

        //then get the edge
        std::vector<std::vector<int> > Kernel;
        Geo::FindVertexStar(Faces,VertPos.size(),Kernel);

        std::vector<int> Candidates=Kernel[Vert0];
        int Vert1=SelectBestConstraintOrientedAlong(VertPos,Candidates,Vert0);

        ScalarType Dist=(VertPos[Vert0]-VertPos[Vert1]).Norm();
        //the conformal is folded so check this then
        FixedV.push_back(Vert0);
        FixedUV.push_back(Geo::Point2<ScalarType>(0,0));
        FixedV.push_back(Vert1);
        FixedUV.push_back(Geo::Point2<ScalarType>(Dist,0));
    }

public:

    static bool IsFoldedParam(const std::vector<std::vector<int> > &Faces,
                              const std::vector<Geo::Point2<ScalarType> > &VertUV)
    {
        ScalarType meshArea2D = Area(VertUV,Faces,false);
        return (meshArea2D<0);
    }

    static void FoldParam(std::vector<Geo::Point2<ScalarType> > &VertUV)
    {
        Geo::Box2<ScalarType> BBoxUV;
        BBoxUV.Init(VertUV);
        Geo::Point2<ScalarType> CenterUV=BBoxUV.Center();
        for(size_t i=0;i<VertUV.size();i++)
            VertUV[i]=CenterUV+(CenterUV-VertUV[i]);
    }

    static void ScaleToMatchAreas(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                  const std::vector<std::vector<int> > &Faces,
                                  std::vector<Geo::Point2<ScalarType> > &VertUV)
    {
        //Rescale the parametrization to match the 3D area
        ScalarType meshArea2D = Area(VertUV,Faces);
        ScalarType meshArea3D = Area(VertPos,Faces);
        /*
        std::cout<<"Area 2D:"<<meshArea2D<<std::endl;
        std::cout<<"Area 3D:"<<meshArea3D<<std::endl;*/
        ScalarType scaleFact = std::sqrt(meshArea3D / meshArea2D);

        for (size_t i=0; i<VertUV.size(); i++)
            VertUV[i]= VertUV[i]*scaleFact;
    }

    static bool LeastSquaresConformal(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                      const std::vector<std::vector<int> > &Faces,
                                      std::vector<size_t> &FixedV,
                                      std::vector<Geo::Point2<ScalarType> > &FixedUV,
                                      std::vector<Geo::Point2<ScalarType> > &VertUV)
    {
        assert(FixedV.size()==FixedUV.size());

        Eigen::MatrixXi F;
        typename EigenInterface<ScalarType>::MatrixXm Vf;

        EigenInterface<ScalarType>::GetTriMeshData(VertPos,Faces,F,Vf);
        Eigen::MatrixXd V = Vf.template cast<double>();

        Eigen::VectorXi b;
        Eigen::MatrixXd bc;
        Eigen::MatrixXd V_uv;

        if (FixedUV.size()==0)
            GetDefaultPairVertexToFix_AlongX(VertPos,FixedV,FixedUV);
        //GetDefaultEdgeToFix(VertPos,Faces,FixedV,FixedUV);

        b.resize(FixedV.size());
        bc.resize(FixedV.size(),2);
        for (int i=0; i<(int)FixedV.size(); i++)
        {
            b(i) = FixedV[i];
            bc(i,0) = FixedUV[i].X;
            bc(i,1) = FixedUV[i].Y;
        }
        // apply Least Square Conformal Maps
        bool solved=igl::lscm(V,F,b,bc,V_uv);
        if (!solved)return false;
        VertUV.resize(VertPos.size());
        // copy results back to mesh
        for (int i=0; i<(int)VertUV.size(); i++) {
            VertUV[i].X = V_uv(i,0);
            VertUV[i].Y = V_uv(i,1);
        }

        return true;
    }

    static bool LeastSquaresConformal(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                      const std::vector<std::vector<int> > &Faces,
                                      std::vector<Geo::Point2<ScalarType> > &VertUV)
    {
        std::vector<size_t> FixedV;
        std::vector<Geo::Point2<ScalarType> > FixedUV;
        return LeastSquaresConformal(VertPos,Faces,FixedV,FixedUV,VertUV);
    }

    static bool Arap(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                     const std::vector<std::vector<int> > &Faces,
                     std::vector<size_t> &FixedV,
                     std::vector<Geo::Point2<ScalarType> > &FixedUV,
                     std::vector<Geo::Point2<ScalarType> > &VertUV,
                     bool InitializeConformal=true, int MaxIte=100)
    {

        bool solved=true;
        bool Has_No_Constr=(FixedV.size()==0);
        if (InitializeConformal)
        {
            solved=LeastSquaresConformal(VertPos,Faces,FixedV,FixedUV,VertUV);
            if (!solved)return false;
            if (Has_No_Constr)
            {
                FixedV.clear();
                FixedUV.clear();
            }
        }
        else
        {
            assert(VertUV.size()==VertPos.size());
        }

        Eigen::MatrixXd V_uv;
        V_uv.resize(VertPos.size(),2);
        for (unsigned i=0;i<V_uv.rows();++i)
        {
            V_uv(i,0)=VertUV[i].X;
            V_uv(i,1)=VertUV[i].Y;
        }

        Eigen::MatrixXi F;
        typename EigenInterface<ScalarType>::MatrixXm Vf;

        //Eigen::MatrixXd PD1,PD2,PV1,PV2;
        EigenInterface<ScalarType>::GetTriMeshData(VertPos,Faces,F,Vf);
        Eigen::MatrixXd V = Vf.template cast<double>();

        Eigen::VectorXi b;
        Eigen::MatrixXd bc;

        if (Has_No_Constr)
        {
            FixedV.push_back(0);
            FixedUV.push_back(VertUV[0]);
        }
        //            FixedV.push_back(1);
        //            FixedUV.push_back(VertUV[1]);
        ////            GetDefaultPairVertexToFix_AlongX(VertPos,FixedV,FixedUV);
        //        }

        b.resize(FixedV.size());
        bc.resize(FixedV.size(),2);
        for (int i=0; i<(int)FixedV.size(); i++)
        {
            b(i) = FixedV[i];
            bc(i,0) = FixedUV[i].X;
            bc(i,1) = FixedUV[i].Y;
        }
        // Add dynamic regularization to avoid to specify boundary conditions
        igl::ARAPData arap_data;
        arap_data.max_iter = MaxIte;
        arap_data.with_dynamics = false;

        // compute ARAP parametrization
        solved=igl::arap_precomputation(V, F, 2, b, arap_data);
        if (!solved)return false;

        igl::arap_solve(bc, arap_data, V_uv);
        VertUV.resize(VertPos.size());
        // copy results back to mesh
        for (int i=0; i<(int)VertUV.size(); i++) {
            VertUV[i].X = V_uv(i,0);
            VertUV[i].Y = V_uv(i,1);
        }

        return true;
    }

    static bool Arap(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                     const std::vector<std::vector<int> > &Faces,
                     std::vector<Geo::Point2<ScalarType> > &VertUV,
                     bool InitializeConformal=true, int MaxIte=100)
    {
        std::vector<size_t> FixedV;
        std::vector<Geo::Point2<ScalarType> > FixedUV;
        return Arap(VertPos,Faces,FixedV,FixedUV,VertUV,InitializeConformal,MaxIte);
    }

    static bool TutteEmbedding(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                               const std::vector<std::vector<int> > &Faces,
                               std::vector<Geo::Point2<ScalarType> > &VertUV)
    {
        VertUV.clear();

        std::vector<std::vector<int> > NextF,NextE;
        getFFAdjacency(Faces,NextF,NextE);
        std::vector<std::vector<int> > BorderVIdx;
        getBorderSequncesVert(Faces,NextF,NextE,BorderVIdx);
        if (BorderVIdx.size()==0)return false;//multiple boundaries

        ScalarType AreaM=Area(VertPos,Faces);
        ScalarType Radius=sqrt(AreaM/M_PI);
        std::vector<Geo::Point2<ScalarType> > vert_UV;
        Geo::UVUtils<ScalarType>::MapToACircle(VertPos,BorderVIdx[0],Radius,vert_UV);

        //check if there is at least one unknown or not
        if (BorderVIdx[0].size()==VertPos.size())
        {
            VertUV.resize(VertPos.size());
            for (size_t i=0;i<BorderVIdx[0].size();i++)
            {
                int CurrIndexV=BorderVIdx[0][i];
                assert(CurrIndexV>=0);
                assert(CurrIndexV<VertUV.size());
                VertUV[CurrIndexV]=vert_UV[i];
            }
            return true;
        }

        Eigen::MatrixXi F;
        typename EigenInterface<ScalarType>::MatrixXm Vf;
        EigenInterface<ScalarType>::GetTriMeshData(VertPos,Faces,F,Vf);

        Eigen::VectorXi b;
        b.resize(BorderVIdx[0].size());
        Eigen::MatrixXd bc;
        bc.resize(BorderVIdx[0].size(),2);
        for (int i=0; i<(int)BorderVIdx[0].size(); i++)
        {
            b(i) = BorderVIdx[0][i];
            bc(i,0) = vert_UV[i].X;
            bc(i,1) = vert_UV[i].Y;
        }


        Eigen::MatrixXd uv_harm;
        bool solved=igl::harmonic(F, b, bc ,1, uv_harm);
        if (!solved)return false;

        VertUV.resize(VertPos.size());
        // copy results back to mesh
        for (int i=0; i<(int)VertUV.size(); i++) {
            VertUV[i].X = uv_harm(i,0);
            VertUV[i].Y = uv_harm(i,1);
        }
        return true;

    }
};

}

#endif
