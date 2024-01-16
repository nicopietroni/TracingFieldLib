#ifndef BEST_ROTATION
#define BEST_ROTATION

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "point3.h"
#include "eigen_interface.h"

namespace  Geo {

template <class ScalarType>
class PointMatch
{
    /*! \brief Computes the best fitting rigid transformations to align two sets of corresponding points
 *
 * Ref:
 * Olga Sorkine-Hornung and Michael Rabinovich
 * Least-Squares Rigid Motion Using SVD
 */
    static void ComputeLeastSquaresRigidMotionEigen(const std::vector<Eigen::Vector3<ScalarType> >  &pFix,
                                                    const std::vector<Eigen::Vector3<ScalarType> >  &pMov,
                                                    Eigen::Matrix3<ScalarType> &R,
                                                    Eigen::Vector3<ScalarType> &t)
    {
        if (pFix.size() != pMov.size())
        {
            assert(0);
        }

        Eigen::Matrix3X<ScalarType> p(3, pMov.size()); // moving
        Eigen::MatrixX3<ScalarType> q(pFix.size(), 3); // fixed

        for (size_t i=0; i<pMov.size(); i++)
        {
            assert(pMov[i].size()==3);
            Eigen::Vector3<ScalarType> v=pMov[i];
            p.col(i) = v;
        }
        Eigen::Vector3<ScalarType> avgP = p.rowwise().mean();
        p.colwise() -= avgP;

        for (size_t i=0; i<pFix.size(); i++)
        {
            assert(pFix[i].size()==3);
            Eigen::Vector3<ScalarType> v=pFix[i];
            q.row(i) = v;
        }

        Eigen::Vector3<ScalarType> avgQ = q.colwise().mean();
        q.rowwise() -= avgQ.transpose();

        Eigen::Matrix3<ScalarType> cov = p * q;
        Eigen::JacobiSVD<Eigen::Matrix3<ScalarType> > svd;
        svd.compute(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::Matrix3<ScalarType> d = Eigen::Matrix3<ScalarType>::Identity();
        d(2,2) = (svd.matrixV() * svd.matrixU().transpose()).determinant() > 0 ? 1 : -1;

        R = (svd.matrixV() * d * svd.matrixU().transpose());
        t = avgQ - R * avgP;
    }

public:


    static void ComputeLeastSquaresRigidMotion(const std::vector<Point3<ScalarType> >  &pFix,
                                               const std::vector<Point3<ScalarType> >  &pMov,
                                               Eigen::Matrix3<ScalarType> &R,
                                               Eigen::Vector3<ScalarType> &T)
    {
        assert(pFix.size()==pMov.size());
        std::vector<Eigen::Vector3<ScalarType> >  pFixE,pMovE;

        pFixE=EigenInterface<ScalarType>::ToEigenVect(pFix);
        pMovE=EigenInterface<ScalarType>::ToEigenVect(pMov);

        ComputeLeastSquaresRigidMotionEigen(pFixE,pMovE,R,T);
    }

    //make the rotation matrix from V0 to V1
    static Eigen::Matrix3<ScalarType> BestRotationMatrix(const Geo::Point3<ScalarType> &V0,
                                                         const Geo::Point3<ScalarType> &V1)
    {
        Eigen::Matrix3<ScalarType> RotM;

        const ScalarType eps=(ScalarType)0.000000001;
        ScalarType dot=Dot(V0,V1);
        ///control if there is no rotation
        if (dot>(ScalarType(1)-eps))
        {
            RotM.setIdentity();
            return RotM;
        }

        //find the axis of rotation
        Geo::Point3<ScalarType> axis;

        //if dot = -1 rotating to opposite vertex
        //the problem is underdefined, so choose axis such that division is more stable
        //alternative solution at http://cs.brown.edu/research/pubs/pdfs/1999/Moller-1999-EBA.pdf
        if (dot < ScalarType(-1) + eps)
        {
            ScalarType maxV;
            size_t maxInd;
            V0.GetMaxAbsCoord(maxV,maxInd);
            axis.V(maxInd) = - (V0.cV((maxInd+2) % 3) / V0.cV(maxInd));
            axis.V((maxInd+1) % 3) = 0;
            axis.V((maxInd+2) % 3) = 1;

            dot = ScalarType(-1);
        }
        else
        {
            axis=Cross(V0,V1);
        }

        axis.Normalize();

        ///construct rotation matrix
        ScalarType u=axis.X;
        ScalarType v=axis.Y;
        ScalarType w=axis.Z;
        ScalarType phi=acos(dot);
        ScalarType rcos = cos(phi);
        ScalarType rsin = sin(phi);

        RotM(0,0) =      rcos + u*u*(1-rcos);
        RotM(1,0) =  w * rsin + v*u*(1-rcos);
        RotM(2,0) = -v * rsin + w*u*(1-rcos);
        RotM(0,1) = -w * rsin + u*v*(1-rcos);
        RotM(1,1) =      rcos + v*v*(1-rcos);
        RotM(2,1) =  u * rsin + w*v*(1-rcos);
        RotM(0,2) =  v * rsin + u*w*(1-rcos);
        RotM(1,2) = -u * rsin + v*w*(1-rcos);
        RotM(2,2) =      rcos + w*w*(1-rcos);

        return RotM;
    }

    static void ApplyTrasformRigidMatch(std::vector<Point3<ScalarType> >  &p,
                                        const Eigen::Matrix3d &R,
                                        const Eigen::Vector3d &t)
    {
        for (size_t i=0;i<p.size();i++)
        {
            p[i]=R*p[i];
            p[i].X=p[i].X+t(0);
            p[i].Y=p[i].Y+t(1);
            p[i].Z=p[i].Z+t(2);
        }
    }

    static void TrasformRigidMatch(const std::vector<Point3<ScalarType> >  &pFix,
                                   std::vector<Point3<ScalarType> >  &pMov)
    {
        assert(pFix.size()==pMov.size());


        Eigen::Matrix3d R;
        Eigen::Vector3d t;
        ComputeLeastSquaresRigidMotion(pFix,pMov,R,t);

        ApplyTrasformRigidMatch(pMov,R,t);
    }


//    static void ComputeBestRotation(const Point3<ScalarType>  &PCentFix,
//                                    const Point3<ScalarType>  &PCentMov,
//                                    const std::vector<Point3<ScalarType> >  &pFix,
//                                    const std::vector<Point3<ScalarType> >  &pMov,
//                                    Eigen::Matrix3<ScalarType> &R)
    static void ComputeBestRotation(const std::vector<Point3<ScalarType> >  &pFix,
                                    const std::vector<Point3<ScalarType> >  &pMov,
                                    Eigen::Matrix3<ScalarType> &R)
    {
        Point3<ScalarType> BaryFix(0,0,0);
        for (size_t i=0;i<pFix.size();i++)
            BaryFix+=pFix[i];
        BaryFix/=pFix.size();

        Point3<ScalarType> BaryMov(0,0,0);
        for (size_t i=0;i<pMov.size();i++)
            BaryMov+=pMov[i];
        BaryMov/=pMov.size();

        std::vector<Point3<ScalarType> > pFix0=pFix;
        std::vector<Point3<ScalarType> > pMov0=pMov;

        for (size_t i=0;i<pFix0.size();i++)
            pFix0[i]=pFix0[i]-BaryFix;//PCentFix;

        for (size_t i=0;i<pMov0.size();i++)
            pMov0[i]=pMov0[i]-BaryMov;//PCentMov;

        Eigen::Vector3<ScalarType> T;
        ComputeLeastSquaresRigidMotion(pFix0,pMov0,R,T);
        //std::cout<<T(0)<<","<<T(1)<<","<<T(2)<<std::endl;
    }

    static void ComputeCovarianceMatrix(const std::vector<Geo::Point3<ScalarType> > &pointVec,
                                 Geo::Point3<ScalarType> &barycenter, Eigen::Matrix<ScalarType,3,3> &m)
    {
      // first cycle: compute the barycenter
      barycenter=Geo::Point3<ScalarType>(0,0,0);
      typename  std::vector< Point3<ScalarType> >::const_iterator pit;
      for( pit = pointVec.begin(); pit != pointVec.end(); ++pit) barycenter+= (*pit);
      barycenter/=pointVec.size();

      // second cycle: compute the covariance matrix
      m.setZero();
      Eigen::Matrix<ScalarType,3,1> p;
      for(pit = pointVec.begin(); pit != pointVec.end(); ++pit) {
          Geo::Point3<ScalarType> Diff=((*pit)-barycenter);
          p=EigenInterface<ScalarType>::ToEigenVect(Diff);
        //((*pit)-barycenter).ToEigenVector(p);
        m+= p*p.transpose(); // outer product
      }
    }

    static void FitPlaneToPointSet(const std::vector< Geo::Point3<ScalarType> > & pointVec,
                            Geo::Plane3<ScalarType> & plane)
    {
      Eigen::Matrix<ScalarType,3,3> covMat = Eigen::Matrix<ScalarType,3,3>::Zero();
      Geo::Point3<ScalarType> b;
      ComputeCovarianceMatrix(pointVec,b,covMat);

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ScalarType,3,3> > eig(covMat);
      Eigen::Matrix<ScalarType,3,1> eval = eig.eigenvalues();
      Eigen::Matrix<ScalarType,3,3> evec = eig.eigenvectors();
      eval = eval.cwiseAbs();
      int minInd;
      eval.minCoeff(&minInd);
      Point3<ScalarType> d;
      d.X = evec(0,minInd);
      d.Y = evec(1,minInd);
      d.Z = evec(2,minInd);

      plane.Init(b,d);
    }
};

}

#endif
