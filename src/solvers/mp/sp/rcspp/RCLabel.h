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

#include "data/Nurse.h"
#include "data/Scenario.h"
#include "data/Shift.h"
#include "tools/Tools.h"

using std::shared_ptr;
using std::unique_ptr;
using std::vector;


// Allow to break down the total cost into smaller pieces.
// It is mainly used for:
// - Classifying the cost type of a resource or a variable
// - Retrieving the part of the total cost related to this type
enum CostType {
  // parts of the cost depending on a rotation, i.e. the sum of:
  // CONS_SHIFTS_COST, CONS_WORK_COST, IDENT_WEEKEND_COST, PREFERENCE_COST
  NO_COST,  // for hard constraint if necessary
  ANY_COST,
  CONS_SHIFTS_COST,  // consecutive shifts constraint
  CONS_WEEKEND_SHIFTS_COST,  // consecutive weekend shifts constraint
  CONS_WORK_COST,  // consecutive worked shifts constraint
  IDENT_WEEKEND_COST,  // cost of not working a complete weekend
  FORBIDDEN_PATTERN_COST,
  PREFERENCE_COST,  // cost of not respecting nurses' preferences
  ALTERNATIVE_COST,  // cost of using an alternative skill or shift
  REST_AFTER_SHIFT_COST,
  CONS_REST_COST,  // consecutive rest shifts constraint
  TOTAL_WORK_COST,  // total shifts' duration
  TOTAL_WEEKEND_COST,  // total worked weekends
  END_INDEX_COST  // mark the last index of the enum: NOT TO BE USED
};

static const std::map<std::string, CostType> costTypesByName = {
    {"CONS_SHIFTS_COST", CONS_SHIFTS_COST},
    {"CONS_WORK_COST", CONS_WORK_COST},
    {"CONS_REST_COST", CONS_REST_COST},
    {"CONS_WEEKEND_SHIFTS_COST", CONS_WEEKEND_SHIFTS_COST},
    {"IDENT_WEEKEND_COST", IDENT_WEEKEND_COST},
    {"FORBIDDEN_PATTERN_COST", FORBIDDEN_PATTERN_COST},
    {"PREFERENCE_COST", PREFERENCE_COST},
    {"ALTERNATIVE_COST", ALTERNATIVE_COST},
    {"REST_AFTER_SHIFT_COST", REST_AFTER_SHIFT_COST},
    {"TOTAL_WORK_COST", TOTAL_WORK_COST},
    {"TOTAL_WEEKEND_COST", TOTAL_WEEKEND_COST},
    {"ANY_COST", ANY_COST},
    {"NO_COST", NO_COST}
};

static const std::map<CostType, std::string> namesByCostType =
    Tools::buildNamesByType(costTypesByName);

static const std::map<CostType, std::string> prettyNamesByCostType =
    Tools::buildPrettyNamesByType(costTypesByName);

std::map<CostType, double> initCostPerType();
std::map<int, double> initCostPerIntType();

// Enum to return the domination status
enum DominationStatus {
  NOT_DOMINATED,  // not dominated in any case
  DOMINATED,  // dominated in any case
  UB_DOMINATED  // dominated only if looking at UB
};

/** Class describing one label used in the label setting algorithm of the
 * subproblem
*
*/
class Resource;
typedef shared_ptr<Resource> PResource;

class RCGraph;
typedef shared_ptr<RCGraph> PRCGraph;

struct RCNode;
typedef shared_ptr<RCNode> PRCNode;

struct RCArc;
typedef shared_ptr<RCArc> PRCArc;

struct ResourceValues {
  int consumption = 0;  // consumption of the resource
  int cyclicConsumption = -1;  // consumption at the beginning of the cycle
  bool readyToConsume = true;  // can be set to false if consumption must
  // be prevented on next arc (e.g. for worked weekend counts)
  double worstLbCost = .0;  // worst costs due to soft lower bound
  double worstUbCost = .0;  // worst costs due to soft upper bounds
  vector<int> states = {};  // for resources such as Forbidden patterns,
  // there can be multiple states for the same label, e.g., if a day-shift
  // appears several time in the pattern

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
  void num(int idLabel) {
    num_ = idLabel;
  }

  friend class RCLabelFactory;

 public:
  RCLabel();
  explicit RCLabel(int nResources);
  explicit RCLabel(const vector<PResource> &resources);
  RCLabel(const vector<PResource> &resources,
          const State &initialState);
  RCLabel(const RCLabel &l);

  ~RCLabel() = default;
  void copy(const RCLabel& l);

  // only copy the costs and the resource values
  void copyValues(const RCLabel& l);

  double cost() const { return cost_; }

  void addCost(double c) { cost_ += c; }

  void setAsNext(shared_ptr<RCLabel> pLPrevious,
                 const shared_ptr<RCArc> &pArc);
  void setAsPrevious(shared_ptr<RCLabel> pLNext,
                     const shared_ptr<RCArc> &pArc);
  void setAsMerged(const shared_ptr<RCLabel> &pLForward,
                   const shared_ptr<RCLabel> &pLBackward);

  double baseCost() const { return baseCost_; }

  void addBaseCost(double c) {
    addCost(c);
    baseCost_ += c;
  }

#ifdef NS_DEBUG
  double dualCost() const { return dualCost_; }

  void addDualCost(double c) {
    dualCost_ += c;
  }

  void addConsShiftCost(double c) {
    consShiftCost_ += c;
  }
  void addConsWeekendShiftCost(double c) {
    consWeekendShiftCost_ += c;
  }
  void addTotalShiftCost(double c) {
    totalShiftCost_ += c;
  }
  void addTotalWeekendCost(double c) {
    totalWeekendCost_ += c;
  }

  void addPreferencesCost(double c) {
    preferencesCost_ += c;
  }
#endif

  int getConsumption(int r) const {
    return resourceValues_[r].consumption;
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

  const vector<ResourceValues> &allResourceValues() const {
    return resourceValues_;
  }

  size_t nResources() const { return resourceValues_.size(); }

  RCNode* getNode() const { return pNode_; }
  void setNode(const shared_ptr<RCNode> &pN) {  pNode_ = pN.get(); }

  const shared_ptr<RCArc>& getInArc() const { return pInArc_; }
  const shared_ptr<RCArc>& getOutArc() const { return pOutArc_; }
  const shared_ptr<RCLabel>& getPreviousLabel() const {return pPreviousLabel_;}
  const shared_ptr<RCLabel>& getNextLabel() const { return pNextLabel_; }

  bool operator()(const shared_ptr<RCLabel> &pr1, const shared_ptr<RCLabel>
  &pr2) {
    return pr1->cost() < pr2->cost();
  }

  // put in a string the label representation
  std::string toString(const vector<PResource> &pResources = {}) const;

  // put in a string the whole chain of labels' representation
  std::string toStringRecursive(const vector<PResource> &pResources = {}) const;

  int num() const { return num_; }

 private:
  int num_;  // Id of the label
  RCNode *pNode_;  // current residing node of the label
  shared_ptr<RCArc> pInArc_;   // pointer to the arc trough which this label
  // was forward-expanded
  shared_ptr<RCArc> pOutArc_;   // pointer to the arc trough which this label
  // was backward-expanded
  shared_ptr<RCLabel> pPreviousLabel_;  // Label from which the current one
  // was expanded in forward propagation
  shared_ptr<RCLabel> pNextLabel_;  // Label from which the current one
  // was backward expanded
  double cost_;   // current cummulated cost of the label
  double baseCost_;  // cost without the dual values
  vector<ResourceValues> resourceValues_;   // vector containing all the
  // resource values

#ifdef NS_DEBUG
  double dualCost_;  // current part of the cost due to dual cost
  double consShiftCost_;  // current cumulated cost due to consecutive shifts
  double consWeekendShiftCost_;
  double totalShiftCost_;  // current cumulated cost due to total shifts
  double totalWeekendCost_;  // current cumulated cost due to worked weekends
  double preferencesCost_;
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

  void updateNumlabel(const PRCLabel& pL) { pL->num(numLabels_++); }

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
  bool operator()(const PRCLabel &pr1, const PRCLabel &pr2) {
    return pr1->cost() < pr2->cost();
  }
};

/**
* Class used to compare two labels costs in order to sort by increasing cost
*/
class LabelCostIncreasing{
 public:
  bool operator()(const PRCLabel &pr1, const PRCLabel &pr2) {
    return pr1->cost() < pr2->cost();
  }
};

/**
*  Class RCSolution
*  This class encodes a solution of the resource-constrained problem
*/

struct RCSolution : public Stretch {
  explicit RCSolution(int firstDay = -1,
                      std::vector<PShift> pShifts = {},
                      double cost = DBL_MAX,
                      double reducedCost = DBL_MAX) :
      Stretch(firstDay, std::move(pShifts)),
      cost_(cost),
      reducedCost_(reducedCost) {}

  explicit RCSolution(const Stretch& stretch,
                      double cost = DBL_MAX,
                      double reducedCost = DBL_MAX) :
      Stretch(stretch),
      cost_(cost),
      reducedCost_(reducedCost) {}

  RCSolution(const Stretch& stretch, const PRCLabel& pL) :
      Stretch(stretch),
      cost_(pL->baseCost()),
      reducedCost_(pL->cost())
#ifdef NS_DEBUG
      , pLabel_(pL)  // NOLINT
#endif
      {}

  std::string toString() const override;

  double cost() const { return cost_; }

  double reducedCost() const { return reducedCost_; }

  void addCost(double c, CostType t)  {
    cost_ += c;
    costs_[t] += c;
  }

  double costByType(CostType t) const {
    return costs_.at(t);
  }

  const std::map<CostType, double>& costs() const {
    return costs_;
  }

  std::string costsToString() const;

  void resetCosts() {
    cost_ = 0;
    costs_ = initCostPerType();
  }

  // Compare rotations on cost
  static bool compareCost(const RCSolution &sol1, const RCSolution &sol2);

  // Compare rotations on dual cost
  static bool compareReducedCost(
      const RCSolution &sol1, const RCSolution &sol2);

  static void sort(std::vector<RCSolution> *solutions) {
    std::stable_sort(solutions->begin(), solutions->end(),
                     [](const RCSolution &sol1, const RCSolution &sol2) {
                       return sol1.reducedCost() < sol2.reducedCost();
                     });
  }

#ifdef NS_DEBUG
  PRCLabel pLabel_ = nullptr;
#endif

 protected:
  double cost_, reducedCost_;
  std::map<CostType, double> costs_;
};

/**
 * Structure storing the information that is necessary to expand a label
 * through an arc
 */
struct Expander {
  explicit Expander(int rId, CostType type):
    indResource(rId), costType(type) {}
  virtual ~Expander() = default;

  // TODO(AL) : can we remove vChild as accessible from pLChild ?
  virtual bool expand(const PRCLabel &pLChild, ResourceValues *vChild) = 0;
  virtual bool expandBack(const PRCLabel &pLChild, ResourceValues *vChild) = 0;

  // !! BEWARE THAT THE RESOURCE ID IS FOR THE SPECIFIC SUBPROBLEM WHERE THE
  // EXPANDER BELONGS, IT CORRESPONDS TO THE POSITIONS OF THE RESOURCE IN THE
  // VECTOR OF RESOURCES OF THE RCSPP GRAPH
  const int indResource;
  const CostType costType;
};
typedef shared_ptr<Expander> PExpander;

/**  Base class that describes a resource of the RCGraph
*
*/
class Resource : public BaseResource {
 public:
  explicit Resource(std::string _name = "", int id = -1) :
      name(std::move(_name)), id_(id) {
    pFirstDay_ = std::make_shared<Day>(MONDAY, 0);
  }

  virtual ~Resource() = default;

  void setId(int id) { id_ = id; }

  int id() const { return id_; }

  PExpander initialize(const AbstractShift &prevAShift,
                       const Stretch &stretch,
                       const shared_ptr<RCArc> &pArc,
                       int indResource);

  // preprocess the input RCGraph to take the resource into consideration
  virtual void preprocess(const PRCGraph &pRCGraph) {}
  // preprocess the arc to take the resource into consideration
  virtual bool preprocess(const PRCArc &pA, double *cost) {
    cost = 0;
    return false;
  }

  //------------------------------------------------
  // Enumeration of sub paths
  //------------------------------------------------
  // This technique consists of getting rid of some resources by adding
  // supplement arcs in the rcGraph. Typically, we can delete the resources
  // corresponding to consecutive shifts types. Arcs with stretches
  // containing several shifts are added to represent a sequence of
  // successive shifts of the same type. To set up this technique, we have to
  // add arcs representing the different possibilities of sequences of shifts
  // types and we need to set the corresponding cost on these new arcs.  To
  // set up this technique, we have to  add arcs representing the different
  // possibilities of sequences of shifts types and we need to set the
  // corresponding costs of these new arcs (base cost, dual cost and cost due
  // to penalties on soft bounds for the number consecutive shifts types). Then,
  // we need to update the cost of the existing arcs (those which existed
  // before the enumeration process). We mainly add cost due to penalties on
  // soft bounds for the number of consecutive shifts types (See Antoine’s
  // internship report for more details).
  virtual void enumerate(const PRCGraph &pRCGraph, bool forceEnum) {}

  virtual DominationStatus dominates(
      RCLabel *pL1, RCLabel *pL2, double *cost) const;

  // return domination status of pL2
  DominationStatus dominates(RCLabel *pL1, RCLabel *pL2) const {
    return dominates(pL1, pL2, nullptr);
  }

  virtual bool merge(const ResourceValues &vForward,
                     const ResourceValues &vBack,
                     ResourceValues *vMerged,
                     const PRCLabel &pLMerged);

  virtual bool isDefaultDomination() const { return useDefaultDomination_; }

  virtual void useDefaultDomination() { useDefaultDomination_ = true; }

  // do nothing if not a SoftBoundedResource
  virtual void useAltenativeDomination() {}

  virtual bool isActive(int dayId, const PAbstractShift &pAShift) const {
    return true;
  }

  virtual bool isActive(int dssrLvl) const {
    return true;
  }

  // return true, if to be considered as a hard constraint
  virtual bool isHard() const { return false; }

  virtual bool isAnyWorkShiftResource() const { return false; }

  virtual int getConsumption(const State &initialState) const { return 0; }

  int totalNbDays() const { return totalNbDays_; }

  void totalNbDays(int totalNbDays) { totalNbDays_ = totalNbDays; }

  int firstDayId() const { return pFirstDay_->id; }

  void firstDayId(int id) {
    pFirstDay_ = pFirstDay_->addAndGet(id);
  }

  std::pair<int, int> getFirstLastDays(const Stretch &stretch) const {
    int firstDay = stretch.firstDayId(), lastDay = stretch.lastDayId();
    // if solving a cyclic problem, the first/last could be before the first day
    // of the horizon
    if (firstDay < pFirstDay_->id)
      firstDay += totalNbDays();
    if (lastDay < pFirstDay_->id)
      lastDay += totalNbDays();
    return {firstDay, lastDay};
  }

  bool isPreprocessed() const { return isPreprocessed_; }

  const std::string name;   // name of the resource

  CostType costType() const { return costType_; }

  // these two methods are used to split the resources between the master and
  // subproblem of each decomposition
  virtual bool isInRosterMaster() const = 0;
  virtual bool isInRotationMaster() const = 0;

 protected:
  int id_;  // id of the resource
  // allow to switch between two different domination functions is necessary
  bool useDefaultDomination_ = true;
  PDay pFirstDay_;  // first day of the horizon
  int totalNbDays_ = 0;  // Total number of days in the horizon
  bool isPreprocessed_ = false;  // true if the resource is completely dealt
  CostType costType_ = ANY_COST;

  // with by modifying the RCGraph during the preprocessing step
  virtual PExpander init(const AbstractShift &prevAShift,
                         const Stretch &stretch,
                         const shared_ptr<RCArc> &pArc,
                         int indResource);
};

/**
 * Class describing a resource which has both an upper and a lower bound
 */
class BoundedResource : public Resource {
 public:
  // Constructor
  BoundedResource(std::string _name, double lb, double ub) :
      Resource(std::move(_name)) {
    if (lb > ub)
      Tools::throwError("LB cannot be greater than the UB for "
                        "BoundedResource: %d > %d", lb, ub);
    setLb(lb);
    setUb(ub);
  }

  // true if rl1 dominates rl2, false otherwise
  DominationStatus dominates(
      RCLabel *pL1, RCLabel *pL2, double *cost) const override;

  int getLb() const {
    return lb_;
  }

  int getUb() const {
    return ub_;
  }

  virtual void setLb(int lb) {
    lb_ = lb;
  }

  virtual void setUb(int ub) {
    ub_ = ub;
  }

  virtual int maxConsumptionPerDay() const {
    return 1;
  }

 protected:
  int lb_;
  int ub_;
};
typedef shared_ptr<BoundedResource> PBoundedResource;

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
  DominationStatus dominates(
      RCLabel *pL1, RCLabel *pL2, double *cost) const override;

  void useAltenativeDomination() override {
    useDefaultDomination_ = false;
  }

  bool isHard() const override {
    return false;
  }

  double getLbCost() const {
    return lbCost_;
  }

  double getUbCost() const {
    return ubCost_;
  }

  void setLbCost(double cost) {
    lbCost_ = cost;
  }

  void setUbCost(double cost) {
    ubCost_ = cost;
  }

  double getLbCost(int consumption) const {
    if (consumption >= lb_) return 0;
    return lbCost_ * (lb_ - consumption);
  }

  double getUbCost(int consumption) const {
    if (consumption <= ub_) return 0;
    return ubCost_ * (consumption - ub_);
  }

  double getCost(int consumption) const {
    return getLbCost(consumption) + getUbCost(consumption);
  }

  virtual double getWorstLbCost(int consumption) const {
    return getLbCost(consumption);
  }

  virtual double getWorstUbCost(int consumption) const {
    return ubCost_ * std::min(consumption, ub_);
  }

  virtual double getWorstUbCost(int consumption, int nLeft) const {
    return getUbCost(consumption + nLeft);
  }

 protected:
  // costs
  double lbCost_ = 0.0;
  double ubCost_ = 0.0;

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;
};

/**
 * Bounded resource where both lower and upper bounds are hard
 */
class HardBoundedResource : public BoundedResource {
 public:
  // Constructor
  HardBoundedResource(std::string _name, int lb, int ub) :
      BoundedResource(std::move(_name), lb, ub) {
    costType_ = NO_COST;
  }

  // return true. If soft, should be overridden
  bool isHard() const override {
    return true;
  }

  bool merge(const ResourceValues &vForward,
             const ResourceValues &vBack,
             ResourceValues *vMerged,
             const PRCLabel &pLMerged) override;
};

#endif  // SRC_SOLVERS_MP_SP_RCSPP_RCLABEL_H_
