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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_MYRCLABEL_H_
#define SRC_SOLVERS_MP_SP_RCSPP_MYRCLABEL_H_


#include <algorithm>
#include <map>
#include <memory>
#include <typeinfo>
#include <utility>
#include <vector>

#include "data/Scenario.h"
#include "data/Shift.h"
#include "tools/MyTools.h"
#include "solvers/mp/sp/rcspp/RCGraph.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;

/** Class describing one label used in the label setting algorithm of the
 * subproblem
*
*/
class Resource;
typedef shared_ptr<Resource> PResource;

struct RCNode;
struct RCArc;


struct ResourceValues {
  int consumption = 0;  // consumption of the resource
  bool readyToConsume = true;  // can be set to false if consumption must
  // be prevented on next arc (e.g. for worked week-end counts)
  double worstLbCost = .0;  // worst costs due to soft lower bound
  double worstUbCost = .0;  // worst costs due to soft upper bounds

  void clear() {
    consumption = 0;
    readyToConsume = true;
    worstLbCost = .0;
    worstUbCost = .0;
  }

  void print() const {
    std::cout << "Consumption=" << consumption
              << " , Ready=" << readyToConsume
              << ", Lb cost=" << worstLbCost
              << ", Ub cost=" << worstUbCost << std::endl;
  }
};

class RCLabelFactory;

class RCLabel {
 private:
  void setNumLabel(int idLabel) {
    num_ = idLabel;
  }

  friend class RCLabelFactory;

 public:
  RCLabel();
  explicit RCLabel(int nResLabels);
  RCLabel(const vector<shared_ptr<Resource>> &resources,
          const State &initialState);
  RCLabel(int idLabel, const RCLabel &l);
  RCLabel(const RCLabel &l);

  void copy(const RCLabel& l);

  int getNum() { return num_; }

  double cost() const { return cost_; }
  double baseCost() const { return baseCost_; }
  double dualCost() const { return dualCost_; }
  double consShiftCost() const { return consShiftCost_; }
  double totalShiftCost() const { return totalShiftCost_; }
  double totalWeekendCost() const { return totalWeekendCost_; }

  void setCost(double c) { cost_ = c; }

  void addCost(double c) { cost_ += c; }

  void addBaseCost(double c) {
    baseCost_ += c;
    cost_ += c;
  }

  void addDualCost(double c) {
    dualCost_ += c;
    cost_ += c;
  }


  void addConsShiftCost(double c) {
    consShiftCost_ += c;
    cost_ += c;
  }


  void addTotalShiftCost(double c) {
    totalShiftCost_ += c,
        cost_ += c;
  }


  void addTotalWeekendCost(double c) {
    totalWeekendCost_ += c;
    cost_ += c;
  }

  int getConsumption(int r) const {
    return resourceValues_[r].consumption;
  }
  int getReadyToConsume(int r) const {
    return resourceValues_[r].readyToConsume;
  }
  double getWorstLbCost(int r) const {
    return resourceValues_[r].worstLbCost;
  }
  double getWorstUbCost(int r) const {
    return resourceValues_[r].worstUbCost;
  }
  ResourceValues& getResourceValues(int r) {
    return resourceValues_[r];
  }

  void setResourceValues(vector<ResourceValues> resourceValues) {
    resourceValues_ = move(resourceValues);
  }

  int nResLabels() const { return nResLabels_; }

  shared_ptr<RCNode> getNode() const { return pNode_; }
  void setNode(shared_ptr<RCNode> pN) {  pNode_ = pN; }

  shared_ptr<RCArc> getInArc() const { return pInArc_; }
  void setInArc(shared_ptr<RCArc> pArc);

  shared_ptr<RCLabel> getPreviousLabel() const {return pParentLabel_;}
  void setParentLabel(shared_ptr<RCLabel> pL) { pParentLabel_ = std::move(pL); }


  bool operator()(const shared_ptr<RCLabel> &pr1, const shared_ptr<RCLabel>
  &pr2) {
    return pr1->cost() < pr2->cost();
  }

 private:
  int num_;  // Id of the label
  int nResLabels_;  // Number of resources represented in the label
  shared_ptr<RCNode> pNode_;  // current residing node of the label
  shared_ptr<RCArc> pInArc_;   // index of the arc where this label
  // was last expanded
  shared_ptr<RCLabel> pParentLabel_;  // Parent label which the current one
  // comes from
  double cost_;   // current cumulated cost of the label
  // containing all the resource values
  vector<ResourceValues> resourceValues_;
  double baseCost_;  // current part of the cost due to base cost
  double dualCost_;  // current part of the cost due to dual cost
  double consShiftCost_;  // current cumulated cost due to consecutive shifts
  double totalShiftCost_;  // current cumulated cost due to total shifts
  double totalWeekendCost_;  // current cumulated cost due to worked weekends
};

typedef shared_ptr<RCLabel> PRCLabel;

/**
 * Factory to create RCLabel
 */

class RCLabelFactory {
 public:
  RCLabelFactory(): numLabels_(0) {}

  template<typename ...Args> PRCLabel makePRCLabel(Args... args) {
    PRCLabel pL(std::make_shared<RCLabel>(args...));
    updateNumlabel(pL);
    return pL;
  }

  void updateNumlabel(PRCLabel pL) {
    pL->setNumLabel(numLabels_++);
  }

 private:
  int numLabels_;  // Total number of labels initialized by the factory
};

typedef shared_ptr<RCLabelFactory> PRCLabelFactory;

/**
* Class used to compare two labels costs in order to sort by increasing cost
*/
class LabelCostIncreasing{
 public:
  bool operator()(const PRCLabel &pr1, const PRCLabel
  &pr2) {
    return pr1->cost() < pr2->cost();
  }
};


/**
 * Class used to compare two labels costs in order to sort by decreasing cost
 */
class LabelCostDecreasing{
 public:
  bool operator()(const PRCLabel &pr1, const PRCLabel
  &pr2) {
    return pr1->cost() > pr2->cost();
  }
};

/**
 * Structure storing the information that is necessary to expand a label
 * through an arc
 */
struct Expander {
  explicit Expander(int rId): resourceId(rId) {}
  virtual ~Expander() = default;

  virtual bool expand(const ResourceValues &vParent,
                      const PRCLabel &pLChild,
                      ResourceValues *vChild) = 0;
  const int resourceId;
};
typedef shared_ptr<Expander> PExpander;

/**
 * Default expander where no information is precomputed before the expansion
 * of labels. Here, the expander only contains the stretch of the arc and the
 * update of the resource consumption is entirely computed during the expansion
 */
struct PlainExpander : public Expander {
  explicit PlainExpander(int rId, const Stretch &stretch) :
  Expander(rId), stretch(stretch) {}

  const Stretch &getStretch() const { return stretch; }

  const Stretch stretch;
};

/**  Base class that describes a resource of the RCGraph
*
*/
class Resource {
 public:
  Resource() : id_(-1) {}

  explicit Resource(int id) : id_(id) {}

  virtual ~Resource() = default;


  void setId(int id) { id_ = id; }

  int id() const { return id_; }

  void initialize(const Shift &prevShift,
                  const Stretch &stretch,
                  const shared_ptr<RCArc> &pArc);

  virtual bool dominates(int conso1,
                         int conso2) {
    return conso1 <= conso2;
  }

  virtual bool isActive(int dayId, const Shift &shift) const { return true; }

  // return true, if to be considered as a hard constraint
  virtual bool isHard() const = 0;

  virtual int getConsumption(const State &initialState) const = 0;

  virtual double getWorstLbCost(int consumption) const = 0;

  virtual double getWorstUbCost(int consumption, int nLeft = 0) const = 0;

 protected:
  int id_;  // id of the resource

  virtual PExpander init(const Shift &prevShift,
                         const Stretch &stretch,
                         const RCArc &arc) = 0;
};

/**
 * Class describing a resource which has both an upper and a lower bound,
 * each being possibly soft with an associated cost
 */
class BoundedResource : public Resource {
 public:
  // Constructor
  //
  BoundedResource(bool isLbSoft, bool isUbSoft,
                  int lb, int ub, double cL = 0.0, double cU = 0.0) :
      isLbSoft_(isLbSoft), isUbSoft_(isUbSoft), lb_(lb), ub_(ub),
      lbCost_(cL), ubCost_(cU) {}

  // true if rl1 dominates rl2, false otherwise
  bool dominates(int conso1, int conso2) override;

  // return true, if at least one bound is hard
  bool isHard() const override {
    return !isLbSoft_ || !isUbSoft_;
  }

  bool isLbSoft() const {
    return isLbSoft_;
  }
  bool isUbSoft() const {
    return isUbSoft_;
  }
  int getLb() const {
    return lb_;
  }
  int getUb() const {
    return ub_;
  }
  double getLbCost() const {
    return lbCost_;
  }
  double getUbCost() const {
    return ubCost_;
  }

  double getLbCost(int consumption) const {
    if (consumption >= lb_) return  0;
    return lbCost_ * (lb_ - consumption);
  }

  double getUbCost(int consumption) const {
    if (consumption <= ub_) return  0;
    return ubCost_ * (consumption - ub_);
  }

 protected:
  bool isLbSoft_;  // true if there is soft lower bound
  bool isUbSoft_;  // true if there is a soft upper bound
  int lb_ = 0;
  int ub_ = 0;
  double lbCost_ = 0.0;
  double ubCost_ = 0.0;
};

/**
 * Bounded resource where both lower and upper bounds are soft
 */
class SoftBoundedResource : public BoundedResource {
 public:
  // Constructor
  SoftBoundedResource(int lb, int ub, double lbCost, double ubCost) :
      BoundedResource(true, true, lb, ub, lbCost, ubCost) {}

  double getWorstLbCost(int consumption) const override {
    return getLbCost(consumption);
  }

  double getWorstUbCost(int consumption, int nLeft = 0) const override {
    return getUbCost(consumption + nLeft);
  }
};

/**
 * Bounded resource where both lower and upper bounds are hard
 */
class HardBoundedResource : public BoundedResource {
 public:
  // Constructor
  HardBoundedResource(int lb, int ub) :
      BoundedResource(true, false, lb, ub) {}

  double getWorstLbCost(int consumption) const override { return .0; }

  double getWorstUbCost(int consumption, int nLeft = 0) const override {
    return .0;
  }
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_MYRCLABEL_H_
