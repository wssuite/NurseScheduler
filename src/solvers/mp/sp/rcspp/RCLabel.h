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

#ifndef SRC_SOLVERS_MP_SP_RCSPP_RCLABEL_H_
#define SRC_SOLVERS_MP_SP_RCSPP_RCLABEL_H_


#include <algorithm>
#include <map>
#include <memory>
#include <typeinfo>
#include <string>
#include <utility>
#include <vector>

#include "data/Scenario.h"
#include "data/Shift.h"
#include "tools/Tools.h"

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
  explicit RCLabel(int nResources);
  RCLabel(const vector<shared_ptr<Resource>> &resources,
          const State &initialState);
  RCLabel(const RCLabel &l);

  void copy(const RCLabel& l);

  double cost() const { return cost_; }

  void addCost(double c) { cost_ += c; }

  void setAsNext(const shared_ptr<RCLabel> &pLPrevious,
                 const shared_ptr<RCArc> &pArc);
  void setAsPrevious(const shared_ptr<RCLabel> &pLNext,
                     const shared_ptr<RCArc> &pArc);

#ifdef DBG
  double baseCost() const { return baseCost_; }
  double dualCost() const { return dualCost_; }
  double consShiftCost() const { return consShiftCost_; }
  double totalShiftCost() const { return totalShiftCost_; }
  double totalWeekendCost() const { return totalWeekendCost_; }
  void addConsShiftCost(double c) {
    consShiftCost_ += c;
  }
  void addTotalShiftCost(double c) {
    totalShiftCost_ += c;
  }
  void addTotalWeekendCost(double c) {
    totalWeekendCost_ += c;
  }
  void addBaseCost(double c) {
    baseCost_ += c;
  }
  void addDualCost(double c) {
    dualCost_ += c;
  }
#endif

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

  vector<ResourceValues> allResourceValues() {
    return resourceValues_;
  }

  int nResources() const { return resourceValues_.size(); }

  shared_ptr<RCNode> getNode() const { return pNode_; }
  void setNode(shared_ptr<RCNode> pN) {  pNode_ = pN; }

  shared_ptr<RCArc> getInArc() const { return pInArc_; }
  shared_ptr<RCArc> getOutArc() const { return pOutArc_; }

  shared_ptr<RCLabel> getPreviousLabel() const {return pPreviousLabel_;}
  void setPreviousLabel(shared_ptr<RCLabel> pL) {
    pPreviousLabel_ = std::move(pL);
  }
  shared_ptr<RCLabel> getNextLabel() const { return pNextLabel_; }
  void setNextLabel(shared_ptr<RCLabel> pL) {pNextLabel_ = std::move(pL);}

  bool operator()(const shared_ptr<RCLabel> &pr1, const shared_ptr<RCLabel>
  &pr2) {
    return pr1->cost() < pr2->cost();
  }

 private:
  int num_;  // Id of the label
  shared_ptr<RCNode> pNode_;  // current residing node of the label
  shared_ptr<RCArc> pInArc_;   // pointer to the arc trough which this label
  // was forward-expanded
  shared_ptr<RCArc> pOutArc_;   // pointer to the arc trough which this label
  // was backward-expanded
  shared_ptr<RCLabel> pPreviousLabel_;  // Label from which the current one
  // was expanded in forward propagation
  shared_ptr<RCLabel> pNextLabel_;  // Label from which the current one
  // was backward expanded
  double cost_;   // current cumulated cost of the label
  vector<ResourceValues> resourceValues_;   // vector containing all the
  // resource values

#ifdef DBG
  double baseCost_;  // current part of the cost due to base cost
  double dualCost_;  // current part of the cost due to dual cost
  double consShiftCost_;  // current cumulated cost due to consecutive shifts
  double totalShiftCost_;  // current cumulated cost due to total shifts
  double totalWeekendCost_;  // current cumulated cost due to worked weekends
#endif
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
* Class used to compare two labels with a view to expand only the K first
 * labels, where K is set by the option param_.NbExpanded_
*/
class SortForFewExpansions {
 public:
  bool operator()(const PRCLabel &pr1, const PRCLabel
  &pr2) {
    return pr1->cost() < pr2->cost();
  }
};

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

  // TODO(AL) : can we remove vChild as accessible from pLChild ?
  virtual bool expand(const PRCLabel &pLChild, ResourceValues *vChild) = 0;
  virtual bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) = 0;
  virtual bool merge(const ResourceValues &vForward,
                     const ResourceValues &vBack,
                     ResourceValues *vMerged,
                     const PRCLabel &pLMerged) = 0;

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
  explicit Resource(std::string _name = "", int id = -1) :
      name(std::move(_name)), id_(-1) {}

  virtual ~Resource() = default;


  void setId(int id) { id_ = id; }

  int id() const { return id_; }

  void initialize(const AbstractShift &prevAShift,
                  const Stretch &stretch,
                  const shared_ptr<RCArc> &pArc);

  virtual bool dominates(const PRCLabel &pL1,
                         const PRCLabel &pL2,
                         double *cost = nullptr);

  virtual void useDefaultDomination() { useDefaultDomination_ = true; }

  virtual void useAltenativeDomination() { useDefaultDomination_ = false; }

  virtual bool isActive(int dayId, const AbstractShift &aShift) const {
    return true;
  }

  // return true, if to be considered as a hard constraint
  virtual bool isHard() const = 0;

  virtual bool isAnyWorkShiftResource() const { return false; }

  virtual int getConsumption(const State &initialState) const = 0;

  const std::string name;   // name of the resource

 protected:
  int id_;  // id of the resource
  // allow to switch between two different domination functions is necessary
  bool useDefaultDomination_ = true;

  virtual PExpander init(const AbstractShift &prevAShift,
                         const Stretch &stretch,
                         const std::shared_ptr<RCArc> &pArc) = 0;
};

/**
 * Class describing a resource which has both an upper and a lower bound
 */
class BoundedResource : public Resource {
 public:
  // Constructor
  //
  BoundedResource(std::string _name, int lb, int ub) :
      Resource(std::move(_name)), lb_(lb), ub_(ub) {}

  // true if rl1 dominates rl2, false otherwise
  bool dominates(const PRCLabel &pL1,
                 const PRCLabel &pL2,
                 double *cost) override;

  int getLb() const {
    return lb_;
  }

  int getUb() const {
    return ub_;
  }

 protected:
  int lb_ = 0;
  int ub_ = 0;
};

/**
 * Bounded resource where both lower and upper bounds are soft
 */
class SoftBoundedResource : public BoundedResource {
 public:
  // Constructor
  SoftBoundedResource(
      std::string _name, int lb, int ub, double lbCost, double ubCost) :
      BoundedResource(std::move(_name), lb, ub),
      lbCost_(lbCost),
      ubCost_(ubCost) {}

  // true if rl1 dominates rl2, false otherwise
  // use worst case for the bounds to determinate if domination
  bool dominates(const PRCLabel &pL1,
                 const PRCLabel &pL2,
                 double *cost) override;

  bool isHard() const override {
    return false;
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

  virtual double getWorstLbCost(int consumption) const {
    return getLbCost(consumption);
  }

  virtual double getWorstUbCost(int consumption, int nLeft = 0) const {
    return getUbCost(consumption + nLeft);
  }

 protected:
  double lbCost_ = 0.0;
  double ubCost_ = 0.0;
};

/**
 * Bounded resource where both lower and upper bounds are hard
 */
class HardBoundedResource : public BoundedResource {
 public:
  // Constructor
  HardBoundedResource(std::string _name, int lb, int ub) :
      BoundedResource(std::move(_name), lb, ub) {}

  // return true. If soft, should be overridden
  bool isHard() const override {
    return true;
  }
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RCLABEL_H_
