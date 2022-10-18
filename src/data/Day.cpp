/*
 * Copyright (C) 2022 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#include "Day.h"

// To get the day name from its type and vice-versa: day type 0 is for Mondays
const string& dayOfWeekToName(DayOfWeek dayOfWeek) {
  return namesByDayOfWeek.at(dayOfWeek);
}

DayOfWeek nameToDayOfWeek(const string &dayName) {
  return daysOfWeekByName.at(dayName);
}
const string& dayOfWeekToShortName(DayOfWeek dayOfWeek) {
  return shortNamesByDayOfWeek.at(dayOfWeek);
}

DayOfWeek shortNameToDayOfWeek(const string &dayName) {
  return daysOfWeekByShortName.at(dayName);
}

bool WeekDay::isWeekend(const AbstractDay& firstWeekendDay,
                        const AbstractDay& lastWeekendDay) const {
  return Weekend::isDayOfWeekend(dayOfWeek, firstWeekendDay, lastWeekendDay);
}

DayOfWeek Day::firstDayOfWeek_ = MONDAY;

DayOfWeek Day::firstDayOfWeek() { return firstDayOfWeek_; }
void Day::setFirstDayOfWeek(DayOfWeek val) { firstDayOfWeek_ = val; }

bool Day::isFirstDayOfWeek(int dayId, DayOfWeek first) {
  return getDayOfWeek(dayId) == first;
}
bool Day::isLastDayOfWeek(int dayId, DayOfWeek last) {
  return getDayOfWeek(dayId) == last;
}

int Day::getDayId(DayOfWeek dayType, int weekId) {
  return dayType - firstDayOfWeek() + 7 * weekId +
      (dayType >= firstDayOfWeek() ? 0 : 7);
}

DayOfWeek Day::getDayOfWeek(int dayId) {
  return (DayOfWeek)((dayId + firstDayOfWeek() + 7) % 7);
}

const string& Day::toDayOfWeekName(int dayId) {
  return dayOfWeekToName(getDayOfWeek(dayId));
}

const string& Day::toDayOfWeekShortName(int dayId) {
  return dayOfWeekToShortName(getDayOfWeek(dayId));
}

bool Weekend::isWeekend(int dayId,
                    const AbstractDay& firstWeekendDay,
                    const AbstractDay& lastWeekendDay) {
  return Weekend::isDayOfWeekend(
      Day::getDayOfWeek(dayId), firstWeekendDay, lastWeekendDay);
}

bool Weekend::isDayOfWeekend(DayOfWeek dayOfWeek,
                             const AbstractDay& firstWeekendDay,
                             const AbstractDay& lastWeekendDay) {
  if (firstWeekendDay.getDayOfWeek() <= lastWeekendDay.getDayOfWeek())
    return (dayOfWeek >= firstWeekendDay.getDayOfWeek()) &&
        (dayOfWeek <= lastWeekendDay.getDayOfWeek());
  else
    return (dayOfWeek >= firstWeekendDay.getDayOfWeek()) ||
        (dayOfWeek <= lastWeekendDay.getDayOfWeek());
}

int Weekend::nWeekendsInInterval(int firstDayId, int lastDayId) const {
  int nWeekends = 0;
  bool isCounted = false;
  for (int i = firstDayId; i <= lastDayId; i++) {
    if (isWeekend(i) && !isCounted) {
      ++nWeekends;
      isCounted = true;
    } else {
      isCounted = false;
    }
  }
  return nWeekends;
}

// if not in weekend, throw an error
// otherwise return position in weekend starting at 0
int Weekend::positionInWeekend(const Day& day) const {
  if (lastWeekendDay_.dayOfWeek >= firstWeekendDay_.dayOfWeek) {
    if ((day.dayOfWeek >= firstWeekendDay_.dayOfWeek) &&
        (day.dayOfWeek <= lastWeekendDay_.dayOfWeek))
      return day.dayOfWeek - firstWeekendDay_.dayOfWeek;
  } else {
    if ((day.dayOfWeek >= firstWeekendDay_.dayOfWeek) ||
        (day.dayOfWeek <= lastWeekendDay_.dayOfWeek))
      return day.dayOfWeek - firstWeekendDay_.dayOfWeek;
    else if (day.dayOfWeek <= lastWeekendDay_.dayOfWeek)
      return day.dayOfWeek + 7 - firstWeekendDay_.dayOfWeek;
  }
  Tools::throwError("Day (%s) is not in weekend (%s, %s)",
                    dayOfWeekToShortName(day.dayOfWeek).c_str(),
                    dayOfWeekToShortName(firstWeekendDay_.dayOfWeek).c_str(),
                    dayOfWeekToShortName(lastWeekendDay_.dayOfWeek).c_str());
  return -1;
}
