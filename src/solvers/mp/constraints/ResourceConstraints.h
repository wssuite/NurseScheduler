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

#ifndef SRC_SOLVERS_MP_CONSTRAINTS_RESOURCECONSTRAINTS_H_
#define SRC_SOLVERS_MP_CONSTRAINTS_RESOURCECONSTRAINTS_H_

#include "ConstraintsMP.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "solvers/mp/sp/rcspp/resources/ConsWeekendShiftResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalShiftDurationResource.h"
#include "solvers/mp/sp/rcspp/resources/TotalWeekendsResource.h"

template <class T>
class BoundedResourceConstraint : public ConstraintMP {
 public:
  explicit BoundedResourceConstraint(MasterProblem *pMaster);

  // update the dual values of the constraints based on the current solution
  void updateDuals() override;

  // return the dual cost of a stretch based on its consumption of
  // the constraints
  double getDualCost(int nurseNum,
                     const Stretch &st,
                     const PAbstractShift &prevS) const override;

  // add a given constraint to the column
  void addConsToCol(std::vector<MyCons *> *cons,
                    std::vector<double> *coeffs,
                    const Pattern &col) const override;

  std::string toString(int nurseNum, const Stretch &st) const override;

  const std::vector<std::map<T*, vector<MyVar*>>>& getVariables() const {
    return resourceVars_;
  }

  const std::vector<std::map<T*, vector<MyCons*>>>& getConstraints() const {
    return resourceCons_;
  }

 protected:
  // variables associated to the resources
  std::vector<std::map<T*, vector<MyVar*>>> resourceVars_;
  // constraints associated to resources
  std::vector<std::map<T*, vector<MyCons*>>> resourceCons_;
  // dual values per nurse and resources
  std::vector<std::map<T*, double>> dualValues_;

  // update the dual values of the constraints randomly
  void createRandomUpdateDuals(double weight);

  template <class H, class S>
  void build(H *pHR, S *pSR, const LiveNurse &pN) {
    char name[255];
    vector<double> coeffs;
    vector<MyVar*> vars;
    vector<MyCons*> cons;
    // compute bounds and slacks
    int lb = 0, ub = INT_MAX, lb_slack = 0, ub_slack = 0;
    bool lb_soft = false, ub_soft = false;
    // update bounds
    if (pHR) {
      lb = pHR->getLb();
      ub = pHR->getUb();
    }
    if (pSR) {
      if (pSR->getLb() > lb) {
        // update slack
        lb_soft = true;
        lb_slack = pSR->getLb() - lb;
        lb = pSR->getLb();
      }
      if (pSR->getUb() < ub) {
        // update slack
        ub_soft = true;
        // if no hard constraint on the UB, the slack is not bounded
        ub_slack = ub < INT_MAX ? ub - pSR->getUb() : INT_MAX;
        ub = pSR->getUb();
      }
    }

    /* create constraints and slack variables if soft constraint */
    const char* cName = pSR ? pSR->name.c_str() : pHR->name.c_str();
    MyCons *pC;
    // LB slack and constraint
    if (lb_soft) {
      vector<MyVar*> lbVars;
      vector<double>  lbCoeffs;
      if (pSR) {
        MyVar *vLB;
        snprintf(name, sizeof(name), "%s_lb_N%d_slack", cName, pN.num_);
        pModel()->createPositiveVar(
            &vLB, name, pSR->getLbCost(), {}, 0, lb_slack);
        lbVars = {vLB}; lbCoeffs = {1};
        vars.push_back(vLB);
      }
      snprintf(name, sizeof(name), "%s_lb_N%d", cName, pN.num_);
      pModel()->createGEConsLinear(&pC, name, lb, lbVars, lbCoeffs);
      cons.push_back(pC);
    }
    // UB slack and constraint
    if (ub_soft) {
      vector<MyVar*> ubVars;
      vector<double>  ubCoeffs;
      if (pSR) {
        MyVar *vUB;
        snprintf(name, sizeof(name), "%s_ub_N%d_slack", cName, pN.num_);
        pModel()->createPositiveVar(
            &vUB, name, pSR->getUbCost(), {}, 0, ub_slack);
        ubVars = {vUB}; ubCoeffs = {-1};
        vars.push_back(vUB);
      }
      snprintf(name, sizeof(name), "%s_ub_N%d", cName, pN.num_);
      pModel()->createLEConsLinear(&pC, name, ub, ubVars, ubCoeffs);
      cons.push_back(pC);
    }

    // store the variables and constraints
    T *pR = pHR;
    if (!pR) pR = pSR;
    resourceVars_[pN.num_][pR] = vars;
    resourceCons_[pN.num_][pR] = cons;
  }
};

class TotalShiftDurationConstraint :
    public BoundedResourceConstraint<TotalShiftDuration> {
 public:
  explicit TotalShiftDurationConstraint(MasterProblem *pMaster);


  void addConstraintFor(const shared_ptr<HardTotalShiftDurationResource> &pHR,
                        const shared_ptr<SoftTotalShiftDurationResource> &pSR,
                        const LiveNurse &pN) {
    build(pHR.get(), pSR.get(), pN);
  }

  // update the dual values of the constraints randomly
  void randomUpdateDuals(bool useInputData, int nPerturbations) override;
};

class TotalWeekendConstraint :
    public BoundedResourceConstraint<TotalWeekend> {
 public:
  explicit TotalWeekendConstraint(MasterProblem *pMaster);

  void addConstraintFor(const shared_ptr<HardTotalWeekendsResource> &pHR,
                        const shared_ptr<SoftTotalWeekendsResource> &pSR,
                        const LiveNurse &pN) {
    build(pHR.get(), pSR.get(), pN);
  }

  // update the dual values of the constraints randomly
  void randomUpdateDuals(bool useInputData, int nPerturbations) override;
};
#endif  // SRC_SOLVERS_MP_CONSTRAINTS_RESOURCECONSTRAINTS_H_
