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

