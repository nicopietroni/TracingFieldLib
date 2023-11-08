#ifndef LEAST_SQUARES_SOLVE
#define LEAST_SQUARES_SOLVE

#include <Eigen/Sparse>
#include "linear_constraint.h"

namespace Geo
{
template <class ScalarType>
void InitSparse(const std::vector<std::pair<int,int> > &Index,
                const std::vector<ScalarType> &Values,
                const int m,
                const int n,
                Eigen::SparseMatrix<ScalarType>& A)
{
    A.setZero();

    assert(Index.size()==Values.size());

    std::vector<Eigen::Triplet<ScalarType> > IJV;
    IJV.reserve(Index.size());

    for(size_t i= 0;i<Index.size();i++)
    {
        int row=Index[i].first;
        int col=Index[i].second;
        ScalarType val=Values[i];

        assert(row<m);
        assert(col<n);

        IJV.push_back(Eigen::Triplet<ScalarType>(row,col,val));
    }
    A.resize(m,n);
    A.setFromTriplets(IJV.begin(),IJV.end());
}

template <class ScalarType>
void InitFromConstraints(const std::vector<LinearConstraint<ScalarType> > &LinConstr,
                         Eigen::SparseMatrix<ScalarType>& A0,
                         Eigen::SparseMatrix<ScalarType>& A)
{
    std::vector<std::pair<int,int> > index;
    std::vector<ScalarType> entry;
    int MaxVarCoord=-1;
    for (size_t i=0;i<LinConstr.size();i++)
    {
        int CurrRow=i;
        for (size_t j=0;j<LinConstr[i].LHSide.size();j++)
        {
            size_t IndexVar=LinConstr[i].LHSide[j].IndexVar;
            size_t IndexCoord=LinConstr[i].LHSide[j].IndexCoord;
            size_t VarCoord=IndexVar*3+IndexCoord;
            MaxVarCoord=std::max(MaxVarCoord,(int)VarCoord);
            ScalarType W=LinConstr[i].LHSide[j].W;
            if (W==0)continue;
            index.push_back(std::pair<int,int>(CurrRow,VarCoord));
            entry.push_back(W*LinConstr[i].Weight);
        }
    }

    //int sizeM=Kernels.size()*3+FixedVert.size()*3+LinConstr.size();
    int sizeRow=LinConstr.size();
    int sizeCol=MaxVarCoord+1;
    InitSparse(index,entry,sizeRow,sizeCol,A0);
    //then

    A=A0.transpose()*A0;
}

template <class ScalarType>
class LinearConstraintLeastSQSolver
{
    Eigen::SparseMatrix<ScalarType> A0,A;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType > > *solver;
    std::vector<LinearConstraint<ScalarType> > LinConstr;
    typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;
    MatrixXm V;

public:

    void Init(const std::vector<LinearConstraint<ScalarType> > &_LinConstr)
    {
        LinConstr=_LinConstr;
        InitFromConstraints<ScalarType>(LinConstr,A0,A);
        solver= new Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType > > (A);
        assert(solver->info() == Eigen::Success);
    }

    void Solve(std::vector<ScalarType> &Solution,const std::vector<ScalarType> &RHS=std::vector<ScalarType>())
    {
        bool UseOriginalRHS=(RHS.size()!=LinConstr.size());
        int NumRows=LinConstr.size();

        V=MatrixXm(NumRows,1);
        V.setZero();
        for (int i=0;i<NumRows;i++)
        {
            int CurrRow=i;

            if (UseOriginalRHS)
                V(CurrRow)=LinConstr[i].RHSide*LinConstr[i].Weight;
            else
                V(CurrRow)=RHS[i]*LinConstr[i].Weight;
        }

        V=A0.transpose()*V;
        //solve the system
        V = solver->solve(V).eval();

        //update the position
        Solution.resize(V.rows());

        for (int i=0;i<V.rows();i++)
            Solution[i]=V(i,0);

    }
};
}

#endif
