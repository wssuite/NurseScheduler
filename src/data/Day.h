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

#ifndef SRC_DATA_DAY_H_
#define SRC_DATA_DAY_H_

#include <map>
#include <memory>
#include <string>
#include <utility>

#include "tools/Tools.h"

using std::string;
using std::shared_ptr;

/*
 * Below are all the methods related to days of the week
 * Beware that two concepts are involved:
 * - the day type is an integer index from 0 to 6 that characterizes  each
 * day of the week; Monday is 0, Tuesday is 1, etc.
 * - the day id characterizes a day of the scheduling horizon; the id of
 * the first day of the horizon 0, but it may not be a Monday
 */

enum DayOfWeek {
  MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY, NO_DAY_OF_WEEK
};

static const std::map<string, DayOfWeek> daysOfWeekByName = {
    {"Monday", MONDAY},
    {"Tuesday", TUESDAY},
    {"Wednesday", WEDNESDAY},
    {"Thursday", THURSDAY},
    {"Friday", FRIDAY},
    {"Saturday", SATURDAY},
    {"Sunday", SUNDAY},
    {"None", NO_DAY_OF_WEEK}
};
static std::map<DayOfWeek, string> namesByDayOfWeek =
    Tools::buildNamesByType(daysOfWeekByName);

static const std::map<string, DayOfWeek> daysOfWeekByShortName  = {
    {"Mon", MONDAY},
    {"Tue", TUESDAY},
    {"Wed", WEDNESDAY},
    {"Thu", THURSDAY},
    {"Fri", FRIDAY},
    {"Sat", SATURDAY},
    {"Sun", SUNDAY},
    {"None", NO_DAY_OF_WEEK}
};
static std::map<DayOfWeek, string> shortNamesByDayOfWeek =
    Tools::buildNamesByType(daysOfWeekByShortName);

// To get the day name from its type and vice-versa: day type 0 is for Mondays
const string& dayOfWeekToName(DayOfWeek dayOfWeek);
DayOfWeek nameToDayOfWeek(const string &dayName);
const string& dayOfWeekToShortName(DayOfWeek dayOfWeek);
DayOfWeek shortNameToDayOfWeek(const string &dayName);


struct AbstractDay {
  AbstractDay() = default;
  virtual ~AbstractDay() = default;

  virtual DayOfWeek getDayOfWeek() const { return NO_DAY_OF_WEEK; }
  virtual int getId() const {return -1;}
  virtual bool includes(const AbstractDay &) const = 0;
  virtual bool isWeekend(const AbstractDay& firstWeekendDay,
                         const AbstractDay& lastWeekendDay) const = 0;
  virtual string toString() const = 0;
};

typedef shared_ptr<AbstractDay> PAbstractDay;


struct AnyDay : public AbstractDay {
  AnyDay(): AbstractDay() {}

  bool includes(const AbstractDay &) const override { return true;}
  string toString() const override { return "Any day";}
  bool isWeekend(const AbstractDay& firstWeekendDay,
                 const AbstractDay& lastWeekendDay) const override {
    return true;
  }
};

struct WeekDay : public AbstractDay {
  explicit WeekDay(DayOfWeek t) :
      AbstractDay(), dayOfWeek(t) {
    if (t == NO_DAY_OF_WEEK)
      Tools::throwError("Cannot use NONE for a WeekDay.");
  }

  DayOfWeek getDayOfWeek() const override {return dayOfWeek;}

  bool includes(const AbstractDay &d) const override {
    return d.getDayOfWeek() == dayOfWeek;
  }

  bool isWeekend(const AbstractDay& firstWeekendDay,
                 const AbstractDay& lastWeekendDay) const override;

  string toString() const override {
    return string("Any" + dayOfWeekToShortName(dayOfWeek));
  }

  int daysOfWeekDifference(const WeekDay& day) const {
    return day.dayOfWeek >= dayOfWeek ?
           day.dayOfWeek - dayOfWeek + 1 : day.dayOfWeek + 7 - dayOfWeek + 1;
  }

  const DayOfWeek dayOfWeek;
};

struct Day : public WeekDay {
  static DayOfWeek firstDayOfWeek();
  static void setFirstDayOfWeek(DayOfWeek val);

  static bool isFirstDayOfWeek(int dayId, DayOfWeek first = MONDAY);
  static bool isLastDayOfWeek(int dayId, DayOfWeek first = SUNDAY);

  static DayOfWeek getDayOfWeek(int dayId);
  static int getDayId(DayOfWeek dayType, int weekId);

  static const string& toDayOfWeekName(int dayId);

  static const string& toDayOfWeekShortName(int dayId);

 private:
  static DayOfWeek firstDayOfWeek_;

 public:
  Day() : WeekDay(NO_DAY_OF_WEEK), id(-1) {}
  explicit Day(int dayId) :
      WeekDay(getDayOfWeek(dayId)),
      id(dayId) {
  }
  explicit Day(DayOfWeek type, int dayId) : WeekDay(type), id(dayId)  {}

  DayOfWeek getDayOfWeek() const override {return dayOfWeek;}
  int getId() const override {return id;}
  bool includes(const AbstractDay &d) const override {
    return d.getId() == id;
  }

  const string& toDayOfWeekName() const {
    return dayOfWeekToName(dayOfWeek);
  }

  const string& toDayOfWeekShortName() const {
    return dayOfWeekToShortName(dayOfWeek);
  }

  shared_ptr<Day> addAndGet(int i) const {
    auto newDayOfWeek = (DayOfWeek)((dayOfWeek + i)%7);
    return std::make_shared<Day>(newDayOfWeek >= 0 ? newDayOfWeek :
                                 (DayOfWeek)(newDayOfWeek+7), id+i);
  }
  shared_ptr<Day> previous() const {
    DayOfWeek newDayOfWeek =
        (dayOfWeek == MONDAY) ? SUNDAY : (DayOfWeek)(dayOfWeek-1);
    return std::make_shared<Day>(newDayOfWeek, id-1);
  }
  shared_ptr<Day> next() const {
    DayOfWeek newDayOfWeek =
        (dayOfWeek == SUNDAY) ? MONDAY : (DayOfWeek)(dayOfWeek+1);
    return std::make_shared<Day>(newDayOfWeek, id+1);
  }

  string toString() const override {
    return string(dayOfWeekToShortName(dayOfWeek) + " " + std::to_string(id));
  }

  const int id;
};

typedef shared_ptr<Day> PDay;

struct Weekend {
  static bool isDayOfWeekend(DayOfWeek dayOfWeek,
                             const AbstractDay& firstWeekendDay,
                             const AbstractDay& lastWeekendDay);

  static bool isWeekend(int dayId,
                        const AbstractDay& firstWeekendDay,
                        const AbstractDay& lastWeekendDay);

  explicit Weekend(DayOfWeek firstWeekendDayType = SATURDAY,
                   DayOfWeek lastWeekendDayType = SUNDAY):
      firstWeekendDay_(firstWeekendDayType),
      lastWeekendDay_(lastWeekendDayType),
      nDays_(firstWeekendDay_.daysOfWeekDifference(lastWeekendDay_)) {}

  Weekend(WeekDay firstWeekendDay, WeekDay lastWeekendDay):
      firstWeekendDay_(std::move(firstWeekendDay)),
      lastWeekendDay_(std::move(lastWeekendDay)),
      nDays_(firstWeekendDay_.daysOfWeekDifference(lastWeekendDay_)) {}

  const WeekDay& firstWeekendDay() const { return firstWeekendDay_; }
  const WeekDay& lastWeekendDay() const {return lastWeekendDay_; }

  bool isWeekend(int day) const {
    return isWeekend(day, firstWeekendDay_, lastWeekendDay_);
  }

  bool isWeekend(const AbstractDay &day) const {
    return day.isWeekend(firstWeekendDay_, lastWeekendDay_);
  }

  bool isWeekend(const PDay &pD) const {
    return isWeekend(*pD);
  }

  bool isFirstWeekendDay(int dayId) const {
    return Day::getDayOfWeek(dayId) == firstWeekendDay_.dayOfWeek;
  }

  bool isFirstWeekendDay(const PDay &pD) const {
    return pD->dayOfWeek == firstWeekendDay_.dayOfWeek;
  }

  bool isLastWeekendDay(int dayId) const {
    return Day::getDayOfWeek(dayId) == lastWeekendDay_.dayOfWeek;
  }

  bool isLastWeekendDay(const PDay &pD) const {
    return pD->dayOfWeek == lastWeekendDay_.dayOfWeek;
  }

  bool isWeekendDayButNotLastOne(int dayId) const {
    DayOfWeek dayOfWeek = Day::getDayOfWeek(dayId);
    if (firstWeekendDay_.dayOfWeek <= lastWeekendDay_.dayOfWeek)
      return (dayOfWeek >= firstWeekendDay_.dayOfWeek) &&
          (dayOfWeek < lastWeekendDay_.dayOfWeek);
    else
      return (dayOfWeek >= firstWeekendDay_.dayOfWeek) ||
          (dayOfWeek < lastWeekendDay_.dayOfWeek);
  }

  bool isWeekendDayButNotFirstOne(int dayId) const {
    DayOfWeek dayOfWeek = Day::getDayOfWeek(dayId);
    if (firstWeekendDay_.dayOfWeek <= lastWeekendDay_.dayOfWeek)
      return (dayOfWeek > firstWeekendDay_.dayOfWeek) &&
          (dayOfWeek <= lastWeekendDay_.dayOfWeek);
    else
      return (dayOfWeek > firstWeekendDay_.dayOfWeek) ||
          (dayOfWeek <= lastWeekendDay_.dayOfWeek);
  }

  int nWeekendsInInterval(const Day &firstDay, const Day &lastDay) const {
    return nWeekendsInInterval(firstDay.id, lastDay.id);
  }

  int nWeekendsInInterval(int firstDayId, int lastDayId) const;

  // if not in weekend, throw an error
  // otherwise return position in weekend starting at 0
  int positionInWeekend(const Day& day) const;

  int nDays() const { return nDays_; }

 protected:
  const WeekDay firstWeekendDay_;
  const WeekDay lastWeekendDay_;
  int nDays_;
};

#endif  // SRC_DATA_DAY_H_
