#ifndef INTEGER_LEAST_SQUARES_SOLVE
#define INTEGER_LEAST_SQUARES_SOLVE

#include "least_squares_lin_constraint_solver.h"
#include "linear_constraint.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <eigen_interface.h>
#include <iostream>
#include <utility>
#include <vector>

#define WeightFIX 100000

namespace Geo {
template <class ScalarType> class IntegerLinearConstraintLeastSQSolver {
  LinearConstraintLeastSQSolver<ScalarType> LQSolver;
  std::vector<bool> FixedVar;
  std::vector<int> FixedInt;
  // std::vector<ScalarType> FixedW;

  enum SolveType { ST_IterativeFix, ST_Increment };
  enum FixOrderType { ST_ClosestInt, ST_LowError, StQRPivoting };

  SolveType SType;
  FixOrderType FOrdType;

  int MinValue;
  bool PrioritizeMinVal;
  int InitialSize;

  Eigen::VectorXi StQRorder;

  ScalarType computeEnergyChange(const std::vector<ScalarType> &Residual,
                                 int indexVar, ScalarType x_i_new,
                                 ScalarType x_i_old) {

    Eigen::VectorX<ScalarType> r = EigenInterface<ScalarType>::toEigenVector(Residual);
    Eigen::VectorX<ScalarType> delta_r = LQSolver.OriginalA().col(indexVar) * (x_i_new - x_i_old);
    Eigen::SparseMatrix<ScalarType> A0=LQSolver.OriginalA();
    // std::cout<<"rows A "<<A0.rows()<<std::endl;
    // std::cout<<"rows r "<<r.rows()<<std::endl;
    // std::cout<<"rows delta_r "<<delta_r.rows()<<std::endl;
    // std::cout<<"Size 1 "<<r.rows()<<std::endl;
    ScalarType deltaE = 2 * r.dot(delta_r) + delta_r.squaredNorm();
    return deltaE;
  }

  void GetLessThenMinVal(const std::vector<ScalarType> &Solution,
                         std::vector<int> &LessIdx) {
    LessIdx.clear();
    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i])
        continue;
      if (floor(Solution[i] + 0.5) < MinValue) {
        LessIdx.push_back(i);
        // std::cout<<" Value: "<<Solution[i]<<std::endl;
      }
    }
    // if (LessIdx.size()>0)
    //   std::cout<<" Less than Minimal: "<<LessIdx.size()<<" out of
    //   "<<Solution.size()<<std::endl;
  }

  ScalarType RoundedValue(const ScalarType &CurrVal) {
    return std::max((ScalarType)MinValue, (floor(CurrVal + 0.5)));
  }

  ScalarType FindErrApproximate(const std::vector<ScalarType> &Solution,
                                const int &IndexVar) {
    return (fabs(Solution[IndexVar] - RoundedValue(Solution[IndexVar])));
  }

  ScalarType DiffEnApproximate(const std::vector<ScalarType> &Solution,
                               const std::vector<ScalarType> &Residual,
                               const int &IndexVar) {
    ScalarType OldVal = Solution[IndexVar];
    ScalarType NewVal = RoundedValue(OldVal);
    return (computeEnergyChange(Residual,IndexVar,NewVal,OldVal));
    //return (fabs(Solution[IndexVar] - RoundedValue(Solution[IndexVar])));
  }

  void SortByRound(const std::vector<ScalarType> &Solution,
                   std::vector<std::pair<int, int>> &RoundRank) {
    RoundRank.clear();
    std::vector<std::pair<ScalarType, int>> RoundPrio;

    std::vector<ScalarType> Residual;
    if (FOrdType == ST_LowError) {
      LQSolver.Residual(Solution, Residual);
    }

    // if (FOrdType == StQRPivoting) {
    //   Eigen::ColPivHouseholderQR<Eigen::MatrixX<ScalarType>> qr(LQSolver.OriginalA());
    //   Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = qr.colsPermutation();
    //   order = P.indices();
    //   assert(order.size()==Solution.size());
    // }

    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i])
        continue;

      ScalarType Diff = 0;
      if (FOrdType == ST_ClosestInt)
        Diff = FindErrApproximate(Solution, i); // fabs(Solution[i] - floor(Solution[i] + 0.5));//
      if (FOrdType == ST_LowError)
        Diff = DiffEnApproximate(Solution,Residual, i);
      if (FOrdType == StQRPivoting)
      {
        assert(i<StQRorder.size());
        Diff= StQRorder[i];
      }
      RoundPrio.push_back(std::pair<ScalarType, int>(Diff,i));
      // if (Solution[i]<MinValue)
      // std::cout<<"dede"<<std::endl;
    }
    std::sort(RoundPrio.begin(), RoundPrio.end());
    // std::cout<<"Diff:"<<RoundPrio[0].first<<std::endl;
    for (size_t i = 0; i < RoundPrio.size(); i++) {
      int VarIdx = RoundPrio[i].second;
      int IntValue = floor(Solution[VarIdx] + 0.5);
      RoundRank.push_back(std::pair<int, int>(VarIdx, IntValue));
    }
  }

  std::vector<LinearConstraint<ScalarType>> LinConstr;

public:

  void Init(const std::vector<LinearConstraint<ScalarType>> &_LinConstr,
            const int _MinValue) {
    LinConstr = _LinConstr;
    InitialSize = LinConstr.size();

    MinValue = _MinValue;
    int MaxVar = 0;
    int MaxComp = 0;
    for (size_t i = 0; i < LinConstr.size(); i++)
      for (size_t j = 0; j < LinConstr[i].LHSide.size(); j++) {
        int IndexVar = LinConstr[i].LHSide[j].IndexVar;
        int IndexComp = LinConstr[i].LHSide[j].IndexCoord;
        MaxVar = std::max(MaxVar, IndexVar);
        MaxComp = std::max(MaxComp, IndexComp);
      }
    int NumVar = (1 + MaxVar) * (MaxComp + 1);
    // std::cout<<"Var Max "<<NumVar<<std::endl;
    FixedVar.resize(NumVar, false);
    FixedInt.resize(NumVar, -1);

  }

  bool SolveIterationIterativeFix(std::vector<ScalarType> &Solution) {
    LQSolver.Init(LinConstr);

    if ((FOrdType == StQRPivoting)&&(StQRorder.size()==0)) 
    {
      Eigen::ColPivHouseholderQR<Eigen::MatrixX<ScalarType>> qr(LQSolver.OriginalA());
      Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = qr.colsPermutation();
      StQRorder = P.indices();
      //assert(StQRorder.size()==Solution.size());
    }

    Solution.clear();
    LQSolver.Solve(Solution);

    std::vector<int> LessIdx;
    if (PrioritizeMinVal)
      GetLessThenMinVal(Solution, LessIdx);

    std::vector<int> FixVar;
    std::vector<int> FixValue;
    if (LessIdx.size() > 0) {
      FixVar = LessIdx;
      FixValue.resize(FixVar.size(), MinValue);
    } else {
      std::vector<std::pair<int, int>> RoundRank;
      SortByRound(Solution, RoundRank);
      if (RoundRank.size() > 0) {
        // fix the one on top
        FixVar.push_back(RoundRank[0].first);
        FixValue.push_back(RoundRank[0].second);
      }
    }
    if (FixVar.size() == 0)
      return true;

    for (size_t i = 0; i < FixVar.size(); i++) {
      int IndexVar = FixVar[i];
      int ValueVar = FixValue[i];
      ValueVar = std::max(MinValue, ValueVar);
      FixedVar[IndexVar] = true;
      FixedInt[IndexVar] = ValueVar;

      LinearConstraint<ScalarType> LConstr;
      LConstr.Weight = WeightFIX;
      LConstr.LHSide.push_back(VarTerm<ScalarType>(IndexVar, 0, 1));
      LConstr.RHSide = ValueVar;
      LinConstr.push_back(LConstr);
    }
    return false;
  }

  void SolveIterationIncrement(std::vector<ScalarType> &Solution,
                               const ScalarType &CurrW) {
    // remove the added ones
    LinConstr.resize(InitialSize);

    for (size_t i = 0; i < Solution.size(); i++) {
      LinearConstraint<ScalarType> LConstr;
      int ValueVar = floor(Solution[i] + 0.5);
      if (ValueVar >= MinValue) {
        LConstr.Weight = CurrW;
        LConstr.LHSide.push_back(VarTerm<ScalarType>(i, 0, 1));

        LConstr.RHSide = ValueVar;
        LinConstr.push_back(LConstr);
      } else {
        LConstr.Weight = WeightFIX;
        LConstr.LHSide.push_back(VarTerm<ScalarType>(i, 0, 1));
        LConstr.RHSide = MinValue;
        LinConstr.push_back(LConstr);
      }
    }

    LQSolver.Init(LinConstr);
    Solution.clear();
    LQSolver.Solve(Solution);
  }

  void Solve(std::vector<int> &Solution) {
    bool completed = false;
    std::vector<ScalarType> Float_Solution;
    if (SType == ST_IterativeFix) {
      do {
        completed = SolveIterationIterativeFix(Float_Solution);
        // std::cout<<"Iteration Completed "<<std::endl;
      } while (!completed);
    } else {
      size_t max_steps = 10;
      LQSolver.Init(LinConstr);
      Solution.clear();
      LQSolver.Solve(Float_Solution);
      ScalarType CurrW = 1;
      for (size_t s = 0; s < max_steps; s++) {
        SolveIterationIncrement(Float_Solution, CurrW);
        CurrW *= 2;
      }
    }

    // then round
    Solution.clear();
    for (size_t i = 0; i < Float_Solution.size(); i++) {
      int IntValue = floor(Float_Solution[i] + 0.5);
      IntValue = std::max(MinValue, IntValue);

      Solution.push_back(IntValue);
    }
  }

  IntegerLinearConstraintLeastSQSolver() {
    SType = ST_IterativeFix;
    FOrdType = ST_ClosestInt;//StQRPivoting;;//;//ST_LowError;//
    PrioritizeMinVal = true;
  }
};

} // namespace Geo

#endif
