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

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "tools/Tools.h"

using std::string;
using std::vector;
using std::shared_ptr;

struct AbstractShift {
  explicit AbstractShift(const std::string _name = "") : name(_name) {}
  virtual ~AbstractShift() = default;

  virtual bool isWork() const { return false; }
  virtual bool isRest() const { return false; }
  virtual bool isAnyWork() const {return false;}
  virtual bool isAnyOfType(int t) const {return false;}
  virtual bool isType(int t) const { return false; }
  virtual bool isShift(int i) const { return false; }
  virtual int workTime() const { return 0; }

  virtual void print() const {
    std::cout << name << std::endl;
  }

  virtual bool includes(const AbstractShift &) = 0;

  const string name;
};

typedef shared_ptr<AbstractShift> PAbstractShift;

struct Shift : public AbstractShift {
  // Default constructor builds a resting shift
  Shift() : AbstractShift(),
            id(0),
            type(0),
            duration(0),
            successors(vector<int>()),
            minCons(0),
            maxCons(99) {}
  Shift(int id, int type) : AbstractShift(),
                            id(id),
                            type(type),
                            duration(0),
                            successors(vector<int>()),
                            minCons(0),
                            maxCons(99) {}
  Shift(string str,
        int i,
        int t,
        int time = 0,
        vector<int> list = vector<int>(),
        int m = 0,
        int M = 99) : AbstractShift(std::move(str)),
                      id(i),
                      type(t),
                      duration(time),
                      successors(std::move(list)),
                      minCons(m),
                      maxCons(M) {}


  Shift(const Shift &shift):
      AbstractShift(shift.name),
      id(shift.id),
      type(shift.type),
      duration(shift.duration),
      successors(vector<int>()),
      minCons(shift.minCons),
      maxCons(shift.maxCons) {
    for (auto successorId : shift.successors)
      successors.push_back(successorId);
  }

  virtual ~Shift() = default;

  bool isRest() const override { return type == 0; }
  bool isWork() const override { return type >= 1; }
  bool isShift(int i) const override { return i == id; }
  bool isType(int t) const override { return t == type; }
  int workTime() const override { return duration; }
  bool includes(const AbstractShift &s) override { return s.isShift(id); }

  const int id;
  const int type;
  const int duration;
  vector<int> successors;
  const int minCons;
  const int maxCons;
};

typedef shared_ptr<Shift> PShift;

struct AnyWorkShift : public AbstractShift {
  AnyWorkShift(): AbstractShift("work") {}

  bool isWork() const override { return true; }
  bool isType(int t) const override { return true; }
  bool isAnyWork() const override {return true;}
  bool includes(const AbstractShift &s) override { return s.isWork(); }
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
  bool includes(const AbstractShift &s) override { return s.isType(type); }

  const int type;
};

/*struct Stint {
  Stint(int d, shared_ptr<AbstractShift> s) : day(d), shift(std::move(s)) {}

  int workTime() const { return shift->workTime(); }

  const int day;
  const shared_ptr<AbstractShift> shift;
};*/

/**
 * A stretch is a sequence of shifts starting on a given day: the shifts may
 * be worked or rested
 */
class Stretch {
 public:
  Stretch(const PShift& pShift, int firstDay) :
      firstDay_(firstDay),
      pShifts_({pShift}),
      duration_(pShift->duration) {}

  Stretch(vector<PShift> pShifts, int firstDay) :
      firstDay_(firstDay),
      pShifts_(std::move(pShifts)),
      duration_(computeDuration()) {}

  int firstDay() const { return firstDay_; }
  int nDays() const { return pShifts_.size(); }
  int lastDay() const { return firstDay_+ pShifts_.size() - 1; }
  const PShift& pShift(int ind) const { return pShifts_.at(ind); }
  const vector<PShift> &pShifts() const { return pShifts_; }
  int duration() const { return duration_; }

  bool operator!=(const Stretch &stretch) const {
    if (firstDay() != stretch.firstDay()) return true;
    if (nDays() != stretch.nDays()) return true;
    if (duration() != stretch.duration()) return true;
    for (auto it1=pShifts_.begin(), it2=stretch.pShifts_.begin();
         it1 != pShifts_.end(); it1++, it2++)
      if ((*it1)->id != (*it2)->id) return true;
    return false;
  }

 protected:
  const int firstDay_;
  const vector<PShift> pShifts_;
  // Time duration (in a certain unit: day, hours, half-hours, ...)
  const int duration_;

  int computeDuration() {
    int d = 0;
    for (const PShift &pS : pShifts_) d += pS->duration;
    return d;
  }
};

#endif  // SRC_DATA_SHIFT_H_
