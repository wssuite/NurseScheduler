/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#ifndef SRC_SOLVERS_MP_CONSTRAINTS_ASSIGNMENTCONSTRAINTS_H_
#define SRC_SOLVERS_MP_CONSTRAINTS_ASSIGNMENTCONSTRAINTS_H_

#include "ConstraintsMP.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "solvers/mp/sp/rcspp/RCGraph.h"


class RosterAssignmentConstraint : public ConstraintMP {
 public:
  explicit RosterAssignmentConstraint(MasterProblem *pMaster);

  // update the dual values of the constraints based on the current solution
  void updateDuals() override;

  // update the dual values of the constraints randomly
  void randomUpdateDuals(bool useInputData, int nPerturbations) override;

  // return the dual cost of a stretch based on its consumption of
  // the constraints
  double getDualCost(int nurseNum,
                     const Stretch &st,
                     const PAbstractShift &prevS) const override;

// add a given constraint to the column
  void addConsToCol(std::vector<MyCons *> *cons,
                    std::vector<double> *coeffs,
                    const Column &col) const override;

  std::string toString() const override;
  std::string toString(int nurseNum, const Stretch &st) const override;

  const vector<MyVar*>& getVariables() const {
    return feasibilityVars_;
  }

  const vector<MyCons*>& getConstraints() const {
    return assignmentCons_;
  }

  double getTotalCost() const override {
    return pModel()->getTotalCost(feasibilityVars_);
  }

 protected:
  // Ensure that is nurse has a roster assigned to her
  std::vector<MyCons *> assignmentCons_;
  // feasibility variables for the assignment constraints
  std::vector<MyVar *> feasibilityVars_;
  // dual values per nurse
  vector<double> dualValues_;

  // build the constraints
  void build();
};


class RCSolution;

class RotationGraphConstraint : public ConstraintMP {
 public:
  RotationGraphConstraint(
      MasterProblem *pMaster,
      vector2D<PBoundedResource> masterRotationGraphResources);

  // update the dual values of the constraints based on the current solution
  void updateDuals() override;

  // update the dual values of the constraints randomly
  void randomUpdateDuals(bool useInputData, int nPerturbations) override;

  // return the dual cost of a stretch based on its consumption of
  // the constraints
  double getDualCost(int nurseNum,
                     const Stretch &st,
                     const PAbstractShift &prevS) const override;

  // add a given constraint to the column
  void addConsToCol(std::vector<MyCons *> *cons,
                    std::vector<double> *coeffs,
                    const Column &col) const override;

  std::string toString(int nurseNum, const Stretch &st) const override;

  // return the cost of the initial state by cost
  std::map<CostType, double> getInitialStateCost() const;

  // return the cost of the variables by cost
  std::map<CostType, double> getCosts() const;

  // get a reference to the restsPerDay_ for a Nurse
  const std::vector<MyVar *> &getVariables(PLiveNurse pNurse, int day) const {
    return restsPerDay_[pNurse->num_][day];
  }

  const vector2D<MyVar*>& getVariables() const {
    return restVars_;
  }

  const vector2D<MyCons*>& getConstraints() const {
    return restCons_;
  }

  double getTotalCost() const override {
    return pModel()->getTotalCost(restVars_);
  }

  bool printInSolutionCosts() const override {
    return false;
  }

 protected:
  /* Build each set of constraints
   * Add also the coefficient of a column for each set */
  void build();
  void createRotationNodes(const PRCGraph &pG, const PLiveNurse &pN);
  void createRotationArcs(const PRCGraph &pG, const PLiveNurse &pN);
  void addRotationRestCost(const PRCArc &pArc, const PLiveNurse &pN);
  void createRotationArcsVars(const PRCGraph &pG, const PLiveNurse &pN);
  void createRotationNodesCons(const PRCGraph &pG, const PLiveNurse &pN);

  RCSolution computeCost(const PRCArc &pArc, const PLiveNurse &pN) const;

  // compute the minimum and maximum consecutive rests based on the resources
  // LBs and UBs. Return the value of the LB and the linear cost associated.
  std::pair<int, double> minConsRest(const PLiveNurse &pN);
  std::pair<int, double> maxConsRest(const PLiveNurse &pN);

  // PResources of the rotation graph
  vector2D<PBoundedResource> masterRotationGraphResources_;

  /*
   * Variables
   */
  // rotation graph for each nurse
  vector<PRCGraph> pRotationGraphs_;
  vector2D<PRCNode> firstRestNodes_, maxRestNodes_;
  // stores all the arcs that are resting on a day for each nurse
  vector3D<MyVar *> restsPerDay_;
  // binary variables for the resting arcs in the rotation network
  vector2D<MyVar *> restVars_;
  // contains all the arc variables that are associated to
  // a RC solution of cost != 0
  std::vector<std::map<PRCArc, RCSolution>> rcSolByArc_;

  // stores all the initial solution finishing on the first day
  std::vector<RCSolution> initialStateRCSolutions_;

  /*
   * Constraints
   */
  // transmission of the flow on rotation graph nodes
  // initialization of the flow constraint at the first position of
  // each restFlowCons_[i] (i=nurse)
  vector2D<MyCons *> restCons_;

  // dual values associated to the start and end of a work arc for each nurse
  vector2D<double> startWorkDualValues_, endWorkDualValues_;
};

#endif  // SRC_SOLVERS_MP_CONSTRAINTS_ASSIGNMENTCONSTRAINTS_H_
