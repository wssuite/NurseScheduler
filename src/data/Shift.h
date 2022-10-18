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

#ifndef SRC_DATA_SHIFT_H_
#define SRC_DATA_SHIFT_H_

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "data/Day.h"
#include "tools/Tools.h"

using std::string;
using std::vector;
using std::shared_ptr;


struct AbstractShift {
  explicit AbstractShift(std::string  _name = "None") :
      name(std::move(_name)) {}
  virtual ~AbstractShift() = default;

  virtual bool isWork() const { return false; }
  virtual bool isRest() const { return false; }
  virtual bool isAnyWork() const {return false;}
  virtual bool isAnyRest() const {return false;}
  virtual bool isAnyOfType(int t) const {return false;}
  virtual bool isType(int t) const { return false; }
  virtual bool isSameType(const AbstractShift &s) const { return false; }
  virtual bool isShift(int i) const { return false; }
  virtual bool isSameShift(const AbstractShift &s) const { return false; }
  virtual int workTime() const { return 0; }

  virtual void print() const {
    std::cout << name << std::endl;
  }

  virtual bool includes(const AbstractShift &) const { return false; }

  bool equals(const AbstractShift & s) const {
    return includes(s) && s.includes(*this);
  }

  const string name;
};

typedef shared_ptr<AbstractShift> PAbstractShift;

struct ShiftType : public AbstractShift {
  // Default constructor builds a resting shift
  ShiftType() : AbstractShift(),
                type(0), work(false), rest(true) {}

  ShiftType(string str, int t) : AbstractShift(std::move(str)), type(t),
                                 work(type >= 1), rest(type == 0) {}

  ShiftType(const ShiftType &shift):
      AbstractShift(shift.name), type(shift.type),
      work(shift.work), rest(shift.rest) {}

  ShiftType(string str, int t, bool w, bool r) :
      AbstractShift(std::move(str)), type(t), work(w), rest(r) {}

  ~ShiftType() override = default;

  bool isRest() const override { return rest; }
  bool isWork() const override { return work; }
  bool isType(int t) const override { return t == type; }
  bool includes(const AbstractShift &s) const override {
    return s.isType(type);
  }
  bool isSameType(const AbstractShift &s) const override {
    return includes(s);
  }

  const int type;
  const bool work;
  const bool rest;
};

struct Shift : public ShiftType {
  // Default constructor builds a resting shift
  Shift() : ShiftType(),
            id(0),
            duration(0),
            timeStart(),
            timeEnd(),
            successors(vector<int>()),
            skills(vector<int>()) {}

  Shift(string str, int i, int t, int time, vector<int> list,
        vector<int> skList = vector<int>(),
        Tools::Time tStart = Tools::Time(),
        Tools::Time tEnd = Tools::Time()) :
      ShiftType(std::move(str), t),
      id(i),
      duration(time),
      timeStart(tStart),
      timeEnd(tEnd),
      successors(std::move(list)),
      skills(std::move(skList)) {}

  Shift(string str, int i, int t, int time, vector<int> list, bool w, bool r,
        vector<int> skList = vector<int>(),
        Tools::Time tStart = Tools::Time(),
        Tools::Time tEnd = Tools::Time()) :
      ShiftType(std::move(str), t, w, r),
      id(i),
      duration(time),
      timeStart(tStart),
      timeEnd(tEnd),
      successors(std::move(list)),
      skills(std::move(skList)) {}


  Shift(const Shift &shift):
      ShiftType(shift.name, shift.type),
      id(shift.id),
      duration(shift.duration),
      timeStart(shift.timeStart),
      timeEnd(shift.timeEnd),
      successors(vector<int>()),
      skills(vector<int>()) {
    for (auto successorId : shift.successors)
      successors.push_back(successorId);
    for (auto skillId : shift.skills)
      skills.push_back(skillId);
  }

  ~Shift() override = default;

  bool isShift(int i) const override { return i == id;}
  int workTime() const override { return duration; }
  bool includes(const AbstractShift &s) const override {
    return s.isShift(id);
  }
  bool canPrecede(const Shift &succShift) const {
    return std::any_of(successors.begin(), successors.end(),
                       [succShift](int s) { return s == succShift.id; });
  }
  bool canSucceed(const Shift &prevShift) const {
    return prevShift.canPrecede(*this);
  }
  // if the lists of skills is empty, the shift is for any skill
  bool hasSkills() const {return !skills.empty();}
  bool hasSkill(int skillId) const {
    if (!hasSkills()) return true;
    return std::any_of(skills.begin(), skills.end(),
                       [skillId](int sk) { return sk == skillId; });
  }

  bool isSameShift(const AbstractShift &s) const override {
    return includes(s);
  }

  // member variables
  //
  const int id;
  const int duration;
  Tools::Time timeStart;
  Tools::Time timeEnd;
  // ids of all the shift that can come after this one in a schedule
  vector<int> successors;
  // ids of the skills that are needed for this shift
  vector<int> skills;
};

typedef shared_ptr<Shift> PShift;


struct AnyRestShift : public AbstractShift {
  AnyRestShift(): AbstractShift("rest") {}

  bool isRest() const override { return true; }
  bool isWork() const override { return false; }
  bool isAnyRest() const override {return true;}
  bool includes(const AbstractShift &s) const override { return s.isRest(); }
};

struct AnyWorkShift : public AbstractShift {
  AnyWorkShift(): AbstractShift("work") {}

  bool isWork() const override { return true; }
  bool isType(int t) const override { return true; }
  bool isAnyWork() const override {return true;}
  bool includes(const AbstractShift &s) const override { return s.isWork(); }
};

struct AnyShift : public AbstractShift {
  AnyShift(): AbstractShift("any") {}

  bool isRest() const override { return true; }
  bool isWork() const override { return true; }
  bool isAnyRest() const override {return true;}
  bool isAnyWork() const override {return true;}
  bool includes(const AbstractShift &s) const override { return true; }
  bool isAnyOfType(int t) const override {return true;}
  bool isType(int t) const override { return true; }
  bool isSameType(const AbstractShift &s) const override { return true; }
  bool isShift(int i) const override { return true; }
  bool isSameShift(const AbstractShift &s) const override { return true; }
};

struct AnyOfTypeShift : public AbstractShift {
  explicit AnyOfTypeShift(int t, std::string _name = "") :
      AbstractShift(_name.empty() ?
                    "type_"+std::to_string(t) : std::move(_name)),
      type(t) {}

  bool isWork() const override { return type >= 1; }
  bool isRest() const override { return type == 0; }
  bool isAnyOfType(int t) const override { return type == t; }
  bool isType(int t) const override { return t == type; }
  bool includes(const AbstractShift &s) const override {
    return s.isType(type);
  }
  bool isSameType(const AbstractShift &s) const override {
    return includes(s);
  }

  const int type;
};

/** comparator for shifts based on id, type, rest, work **/
struct BaseShiftComparator {
  virtual bool equals(
      const AbstractShift &s1, const AbstractShift &s2) const = 0;
};
typedef shared_ptr<BaseShiftComparator> PBaseShiftComparator;

struct ShiftComparator : BaseShiftComparator {
  bool equals(
      const AbstractShift &s1, const AbstractShift &s2) const override {
    return s1.isSameShift(s2);
  }
};

struct ShiftTypeComparator : BaseShiftComparator {
  bool equals(
      const AbstractShift &s1, const AbstractShift &s2) const override {
    return s1.isSameType(s2);
  }
};

struct ShiftWorkComparator : BaseShiftComparator {
  // either both work or both rest
  bool equals(
      const AbstractShift &s1, const AbstractShift &s2) const override {
    return s1.isWork() ^ s2.isRest();
  }
};

/**
 * A pattern is a sequence of abstract shifts on abstract days. It aims at
 * characterizing sequences of assignments as generally as possible. For
 * instance one shift can be AnyWorkShift and the corresponding day can be
 * any saturday (which is an WeekDay with type 5).
 */
struct Pattern {
 public:
  Pattern(vector<PAbstractShift> pAShifts, vector<PAbstractDay> pADays) :
      size_(pAShifts.size()),
      pAShifts_(std::move(pAShifts)),
      pADays_(std::move(pADays)) {
    // there must be the same number of shifts and days in the pattern
    if (pAShifts_.size() != pADays_.size())
      Tools::throwError("Bad initialization of a pattern, there must be the "
                        "same number of shifts and days");
  }

  const int size_;
  const vector<PAbstractShift> pAShifts_;
  const vector<PAbstractDay> pADays_;
};

/**
 * A stretch is a sequence of shifts starting on a given day: the shifts may
 * be worked or rested
 */
class Stretch {
 public:
  explicit Stretch(int firstDayId = -1) :
      firstDayId_(firstDayId), duration_(0) {}

  Stretch(int firstDayId, const PShift& pShift) :
      firstDayId_(firstDayId),
      pShifts_({pShift}),
      duration_(pShift->duration) {
    pDays_.push_back(std::make_shared<Day>(firstDayId));
  }

  Stretch(const PDay& firstDay, const PShift& pShift) :
      pDays_({firstDay}),
      pShifts_({pShift}),
      firstDayId_(firstDay->id),
      duration_(pShift->duration) {}

  Stretch(int firstDay, vector<PShift> pShifts) :
      firstDayId_(firstDay),
      pShifts_(std::move(pShifts)),
      duration_(computeDuration()) {
    int dayId = firstDay;
    for (const auto &pS : pShifts_)
      pDays_.push_back(std::make_shared<Day>(dayId++));
  }

  Stretch(vector<PDay> pDays, vector<PShift> pShifts) :
      firstDayId_(pDays.front()->id),
      pDays_(std::move(pDays)),
      pShifts_(std::move(pShifts)),
      duration_(computeDuration()) {
#ifdef DBG
    if (pDays_.size() != pShifts_.size()) {
      Tools::throwError("stretches must have as many days as shifts");
    }
#endif
  }

  virtual ~Stretch() = default;

  void copy(const Stretch &st) {
    firstDayId_ = st.firstDayId_;
    pShifts_ = st.pShifts_;
    pDays_ = st.pDays_;
    duration_ = st.duration_;
  }

  void init(int firstDay, const vector<PShift>& pShifts) {
    pShifts_.clear();
    pDays_.clear();
    firstDayId_ = firstDay;
    int dayId = firstDayId_;
    for (const auto &pS : pShifts) {
      pShifts_.push_back(pS);
      pDays_.push_back(std::make_shared<Day>(dayId++));
    }
    duration_ = computeDuration();
  }

  void init(int firstDay, int nDays, const PShift &pSDefault) {
    vector<PShift> pShifts(nDays, pSDefault);
    init(firstDay, pShifts);
  }

  virtual int firstDayId() const { return firstDayId_; }
  virtual int nDays() const { return pShifts_.size(); }
  virtual int lastDayId() const { return firstDayId_+ nDays() - 1; }
  virtual const PShift& pShift(int dayId) const {
    return pShifts_.at(dayId- firstDayId());
  }
  virtual const vector<PDay> &pDays() const { return pDays_; }
  virtual const vector<PShift> &pShifts() const { return pShifts_; }
  virtual int shift(int dayId) const { return pShift(dayId)->id; }
  virtual PDay pDay(int ind) const { return pDays_[ind];}
  virtual int duration() const { return duration_; }
  int nWeekends(const Weekend& weekend) const {
    if (pDays_.empty()) return 0;
    return weekend.nWeekendsInInterval(*pDays_.front(), *pDays_.back());
  }

  virtual void pushFront(const Stretch &stretch) {
    firstDayId_ = stretch.firstDayId_;
    duration_ += stretch.duration_;
    pShifts_.insert(
        pShifts_.begin(), stretch.pShifts_.begin(), stretch.pShifts_.end());
    pDays_.insert(
        pDays_.begin(), stretch.pDays_.begin(), stretch.pDays_.end());
  }

  virtual void pushBack(const Stretch &stretch) {
    duration_ += stretch.duration_;
    pShifts_.insert(
        pShifts_.end(), stretch.pShifts_.begin(), stretch.pShifts_.end());
    pDays_.insert(
        pDays_.end(), stretch.pDays_.begin(), stretch.pDays_.end());
  }

  virtual void pushBack(const PShift &pS) {
    duration_ += pS->duration;
    pShifts_.push_back(pS);
    if (pDays_.empty())
      pDays_.push_back(std::make_shared<Day>(firstDayId_));
    else
      pDays_.push_back(pDays_.back()->next());
  }

  virtual void popBack() {
    assert(!pShifts_.empty());
    duration_ -= pShifts_.back()->duration;
    pShifts_.erase(--pShifts_.end());
    pDays_.erase(--pDays_.end());
  }

  void assignShift(int day, const PShift &pS) {
    duration_ += pS->duration - pShift(day)->duration;
    pShifts_[day] = pS;
  }

  // rotate of n days the stretch: put the n last shifts in front
  virtual void rotate(int n) {
    int length = pShifts_.size();
    auto start = pShifts_.end() - n;
    if (n < 0)  start -= length;
    vector<PShift> toInsert(start, pShifts_.end());
    pShifts_.insert(pShifts_.begin(), toInsert.begin(), toInsert.end());
    firstDayId_ -= n;
#ifdef DBG
    if (firstDayId_ < 0)
      Tools::throwError("stretch has a negative first day. "
                        "Too many shifts have been rotated.");
#endif
    pShifts_.resize(length);
  }

  bool operator!=(const Stretch &stretch) const {
    if (firstDayId() != stretch.firstDayId()) return true;
    if (nDays() != stretch.nDays()) return true;
    if (duration() != stretch.duration()) return true;
    for (auto it1=pShifts_.begin(), it2=stretch.pShifts_.begin();
         it1 != pShifts_.end(); it1++, it2++)
      if ((*it1)->id != (*it2)->id) return true;
    return false;
  }

  virtual std::string toString() const {
    std::stringstream buff;
    buff << "Stretch starting on day " << firstDayId_
         << " (length=" << nDays() << ", duration=" << duration() << "):"
         << std::endl;
    // display days
    for (int k = 0; k < firstDayId_; k++)
      buff << "|" << std::setw(SHIFT_PAD) << "";
    for (const PDay &pD : pDays_)
      buff << "|" << std::setw(SHIFT_PAD)
           << pD->toString().substr(0, SHIFT_PAD);
    buff << "|" << std::endl;
    // display shifts
    for (int k = 0; k < firstDayId_; k++)
      buff << "|" << std::setw(SHIFT_PAD) << "";
    for (const PShift &pS : pShifts_)
      if (pS->isRest())
        buff << "|" << std::setw(SHIFT_PAD) << REST_DISPLAY;
      else
        buff << "|" << std::setw(SHIFT_PAD) << pS->name.substr(0, SHIFT_PAD);
    buff << "|" << std::endl;
    return buff.str();
  }

 protected:
  int firstDayId_;
  vector<PDay> pDays_;
  vector<PShift> pShifts_;
  // Time duration (in a certain unit: day, hours, half-hours, ...)
  int duration_;

  int computeDuration() {
    int d = 0;
    for (const PShift &pS : pShifts_) d += pS->duration;
    return d;
  }
};

#endif  // SRC_DATA_SHIFT_H_
