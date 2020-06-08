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

#ifndef SRC_SOLVERS_MP_SP_LONGROTATIONSP_H_
#define SRC_SOLVERS_MP_SP_LONGROTATIONSP_H_

#include <vector>

#include "solvers/mp/sp/RotationSP.h"

class LongRotationSP : public RotationSP {
 public:
  LongRotationSP() = default;
  virtual ~LongRotationSP();

  // Constructor that correctly sets the resource (time + bounds),
  // but NOT THE COST
  LongRotationSP(PScenario scenario,
                 int nbDays,
                 PConstContract contract,
                 std::vector<State> *pInitState);

  double startWorkCost(int a) const override;
  double historicalCost(int a) const override;

 protected:
  // Number of rotations found (that match the bound condition) at that
  // iteration
  int nVeryShortFound_ = 0;

  //-----------------------
  // THE BASE COSTS
  //-----------------------

  // SHORT SUCCESSIONS (computed when creating them)
  // For each size c \in [0,CDMin], for each short rotation of size c,
  // contains its base cost (independent from the date)
  vector2D<double>  baseArcCostOfShortSucc_;

  //-----------------------
  // THE SHORT SUCCESSIONS
  //-----------------------

  // SHORTSUCC -> OBJECTS
  //
  // Short successions (no starting date) -> those of all length
  // For each size c \in [0,CDMin], contains all allowed short successions of
  // that size (satisfies succession constraints)
  vector3D<int> allowedShortSuccBySize_;
  // contains the corresponding last shift performed
  vector2D<int> lastShiftOfShortSucc_;
  // For each size c \in [0,CDMin], for each short rotation of size c,
  // contains the number of consecutive days the last shift has been performed
  vector2D<int>  nLastShiftOfShortSucc_;
  // Minimum number of consecutive days worked for free
  // int CDMin_;
  // For each shift s, for each number of days n, contains the list of short
  // successions of size CDMin ending with n consecutive days of shift s
  vector3D<int> allShortSuccCDMinByLastShiftCons_;

  // SHORTSUCC -> FUNCTIONS
  //
  // Initializes all short successions, base costs, and corresponding vectors.
  // Should only be called ONCE.
  void initShortSuccessions();

  //----------------------------------------------------------------
  //
  // Cost computation of the "very" short rotations (< CD_min)
  //
  //----------------------------------------------------------------
  bool priceVeryShortRotationsFirstDay();
  bool priceVeryShortRotationsLastDay();
  bool priceVeryShortRotations();
  int priceVeryShortSameSizeRotations(int k, const vector2D<int> &succs);
  double costOfVeryShortRotation(int firstDay, const std::vector<int> &succ);


  //----------------------------------------------------------------
  //
  // Update of the costs / network for solve function
  //
  //----------------------------------------------------------------

  // FUNCTIONS -- SOLVE
  bool preprocess() override;
  bool solveShortRotations();

  // override creation of arcs source -> principal
  void createArcsSourceToPrincipal() override;

  // FUNCTIONS -- COSTS
  // Pricing of the short successions : only keep one of them, and the cost
  // of the corresponding arc
  void priceShortSucc();
  // Given a short succession and a start date, returns the cost of the
  // corresponding arc
  double costArcShortSucc(int size, int id, int startDate);

 public:
  // Some getters
  int nVeryShortFound() const { return nVeryShortFound_; }

  // Print functions.
  void printShortSucc() const;
  void printShortArcs() const;
};

#endif  // SRC_SOLVERS_MP_SP_LONGROTATIONSP_H_
