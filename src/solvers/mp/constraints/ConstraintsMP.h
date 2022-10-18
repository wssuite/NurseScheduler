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

#ifndef SRC_SOLVERS_MP_CONSTRAINTS_CONSTRAINTSMP_H_
#define SRC_SOLVERS_MP_CONSTRAINTS_CONSTRAINTSMP_H_

#include <memory>
#include <string>
#include <vector>

#include "data/Shift.h"
#include "solvers/mp/modeler/Modeler.h"

class MasterProblem;
class Column;
class Modeler;


class ConstraintMP {
 public:
  explicit ConstraintMP(MasterProblem *pMaster,
                        std::string name,
                        bool impactColumns = false,
                        bool addConstraint = true);

  virtual ~ConstraintMP() = default;

  // update the dual values of the constraints based on the current solution
  virtual void updateDuals() {}

  // update the values of the variables and constraints based
  // on the current model values
  virtual void update() {}

  // update the dual values of the constraints randomly
  virtual void randomUpdateDuals(bool useInputData, int nPerturbations) {}

  // return the dual cost of a stretch based on its consumption of
  // the constraints
  virtual double getDualCost(int nurseNum,
                             const Stretch &st,
                             const PAbstractShift &prevS) const {
    return 0;
  }

// add a given constraint to the column
  virtual void addConsToCol(std::vector<MyCons *> *cons,
                            std::vector<double> *coeffs,
                            const Column &col) const {}

  virtual std::string toString() const { return ""; }
  virtual std::string toString(int nurseNum, const Stretch &st) const {
    return "";
  }

  virtual double getTotalCost() const = 0;

  virtual bool printInSolutionCosts() const {
    return true;
  }

  Modeler * pModel() const;

  const std::string name;

 protected:
  MasterProblem *pMaster_;
  PScenario pScenario_;
};

/*
 * Constraints that count the number of nurses working on a given position for
 * each day and shift based on the selected columns.
 */
class NursePositionCountConstraint : public ConstraintMP {
 public:
  explicit NursePositionCountConstraint(MasterProblem *pMaster);

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

  const vector3D<MyVar*>& getVariables() const {
    return numberOfNursesByPositionVars_;
  }

  const vector3D<MyCons*>& getConstraints() const {
    return numberOfNursesByPositionCons_;
  }

  double getTotalCost() const override {
    return 0;
  }

 protected:
  // count the number of nurses by position on each day, shift
  vector3D<MyVar *> numberOfNursesByPositionVars_;
  // ensure there are enough nurses for numberOfNursesByPositionVars_
  vector3D<MyCons *> numberOfNursesByPositionCons_;
  // dual values per nurse, day and shift
  vector3D<double> dualValues_;

  // build the constraints
  void build();
};


class AllocationConstraint : public ConstraintMP {
 public:
  explicit AllocationConstraint(MasterProblem *pMaster);

  const vector4D<MyVar*>& getVariables() const {
    return skillsAllocVars_;
  }

  const vector3D<MyCons*>& getConstraints() const {
    return feasibleSkillsAllocCons_;
  }

  double getTotalCost() const override {
    return 0;
  }

 protected:
  // makes the allocation of the skills
  vector4D<MyVar *> skillsAllocVars_;
  // ensures that each nurse works with the good skill
  vector3D<MyCons *> feasibleSkillsAllocCons_;

  // build the constraints
  void build();
};

class DemandConstraint : public ConstraintMP {
 public:
  DemandConstraint(MasterProblem *pMaster,
                   bool minDemand,
                   bool soft = false,
                   double weight = 0);

  // update the values of the variables and constraints based
  // on the current model values
  void update() override;

  const vector3D<MyVar*>& getVariables() const {
    return slackVars_;
  }

  const vector3D<MyCons*>& getConstraints() const {
    return demandCons_;
  }

  double getTotalCost() const override {
    return pModel()->getTotalCost(slackVars_) +
        pModel()->getTotalCost(feasNegVars_);
  }

 protected:
  bool minDemand_;
  std::string prefix_;
  bool soft_;
  double weight_;
  // demand constraints per day, shift, skills
  vector3D<MyCons *> demandCons_;
  // slack variables for each constraint
  // if soft has a weight, if hard it's for feasibilty (it has big weight)
  vector3D<MyVar *> slackVars_;
  vector3D<MyVar *> feasNegVars_;  // useful for INRC

  // build the constraints
  void build();

  const vector3D<int>& demand() const;
};

#endif  // SRC_SOLVERS_MP_CONSTRAINTS_CONSTRAINTSMP_H_
