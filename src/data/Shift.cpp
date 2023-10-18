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

#include <set>

#include "data/Shift.h"

const PShift & AbstractShift::findIncludedShift(
    const std::vector<int> &shifts) const {
  for (const auto & pS : pShifts_) {
    for (int s : shifts)
      if (pS->isShift(s))
        return pS;
  }
  Tools::throwError(
      "There is no work shift defined within the given vector.");
  return pShifts_.front();
}

ShiftsFactory::ShiftsFactory(
    const vector<PShift> &pShifts,
    const std::vector<string> shiftTypeNames):
    pNoneShift_(std::shared_ptr<NoneShift>(new NoneShift())),
    pAnyShift_(std::shared_ptr<AnyShift>(new AnyShift())),
    pShifts_(pShifts),
    pEndShift_(std::make_shared<Shift>(-1, "End")) {
  // create None shift
  std::vector<int> successors;
  for (const PShift& pS : pShifts) successors.push_back(pS->id);
  pNoneShift_->addShift(std::make_shared<Shift>("None", -1, -1, 0, successors));

  // create work, rest and shift types abstract shifts
  vector2D<PShift> pShiftTypes;
  for (const PShift& pS : pShifts) {
    // add pShift as included for itself
    if (pS->pIncludedShifts().empty()) pS->addShift(pS);

    // check whether rest or work shift
    if (pS->isRest()) {
      // create rest if necessary
      if (pAnyRestShift_ == nullptr)
        pAnyRestShift_ = std::shared_ptr<AnyRestShift>(new AnyRestShift());
      pAnyRestShift_->addShift(pS);
    } else {
      // create work if necessary
      if (pAnyWorkShift_ == nullptr)
        pAnyWorkShift_ = std::shared_ptr<AnyWorkShift>(new AnyWorkShift());
      pAnyWorkShift_->addShift(pS);
    }

    // sort shift type
    if (pS->type >= pShiftTypes.size())
      pShiftTypes.resize(pS->type + 1);
    pShiftTypes[pS->type].push_back(pS);
  }

  // create shift type
  pAnyTypeShifts_.resize(pShiftTypes.size());
  for (int t=0; t < pShiftTypes.size(); t++) {
    string name;
    if (shiftTypeNames.empty()) {
      for (const PShift &pS : pShiftTypes[t]) {
        if (!name.empty()) name += ",";
        name += pS->name;
      }
    } else {
      name = shiftTypeNames.at(t);
    }
    pAnyTypeShifts_[t] =
        std::shared_ptr<AnyOfTypeShift>(new AnyOfTypeShift(t, name));
    for (const PShift& pS : pShiftTypes[t]) pAnyTypeShifts_[t]->addShift(pS);
  }
}

Shifts::Shifts(vector<PAbstractShift> pAShifts):
        // concatenate all names
        AbstractShift(std::accumulate(
                pAShifts.begin(), pAShifts.end(), std::string(),
                [](const std::string &a, const PAbstractShift &pS) {
                    return a + (a.length() > 0 ? "_" : "") + pS->name;
                })),
        pAShifts_(std::move(pAShifts)) {
  std::set<PShift> allShifts;
  for (const auto &pAS : pAShifts_)
    for (const auto &pS : pAS->pIncludedShifts())
      allShifts.insert(pS);
  // add all the included shifts
  pShifts_ = vector<PShift>(allShifts.begin(), allShifts.end());
}

bool Shifts::isWork() const {
  for (const auto &pS : pAShifts_)
    if (!pS->isWork()) return false;
  return true;
}

bool Shifts::isRest() const {
  for (const auto &pS : pAShifts_)
    if (!pS->isRest()) return false;
  return true;
}

bool Shifts::isType(int t) const {
  for (const auto &pS : pAShifts_)
    if (!pS->isType(t)) return false;
  return true;
}

bool Shifts::isSameType(const AbstractShift &s) const {
  for (const auto &pS : pAShifts_)
    if (!pS->isSameType(s)) return false;
  return true;
}

bool Shifts::includes(const AbstractShift &s) const {
  for (const auto &pS : pAShifts_)
    if (pS->includes(s)) return true;
  return false;
}
