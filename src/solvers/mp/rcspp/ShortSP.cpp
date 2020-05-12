//
// Created by antoine legrain on 2020-04-12.
//

#include "ShortSP.h"

using std::string;
using std::vector;
using std::map;

//---------------------------------------------------------------------------
//
// C l a s s   S u b P r o b l e m
//
// Contains the shortest paths with resource constraints
//
//---------------------------------------------------------------------------

// Constructors and destructor
ShortSP::ShortSP() {}

ShortSP::ShortSP(PScenario scenario, int nbDays, PConstContract contract, vector<State>* pInitState):
    SubProblem(scenario, nbDays, contract,  pInitState) {

  minConsDays_ = contract->minConsDaysWork_;
  nLabels_ = 1;

  initShortSuccessions();
}

ShortSP::~ShortSP(){}

// Initializes the short successions. Should only be used ONCE (when creating the SubProblem).
void ShortSP::initShortSuccessions() {

  // Primary information needed
  //
  int nShiftsType = pScenario_->nbShiftsType_;
  int nShifts = pScenario_->nbShifts_;

  // Put an empty list of size 0 for all data because there exists no succession of length 0/
  //
  allowedShortSuccBySize_.emplace_back(vector2D<int>());
  lastShiftOfShortSucc_.emplace_back(vector<int>());
  nLastShiftOfShortSucc_.emplace_back(vector<int>());
  baseArcCostOfShortSucc_.emplace_back(vector<double>());

  // Initialize the other way round
  for(int s=0; s<nShiftsType; s++)
    allShortSuccCDMinByLastShiftCons_.emplace_back(vector2D<int>(maxCons(s)+1));

  // For all size from 1 to CD_min, compute all allowed shift successions.
  //
  for (int c = 1; c <= CDMin_; c++) {

    // Blank data
    vector2D<int> allSuccSizeC;
    vector<int> lastShiftSucc;
    vector<int> nLastShiftSucc;
    vector<double> arcCostSucc;

    // Compute new succession
    //

    // Size 1 -> special case, initialization -> add all single shift rotations
    //
    if (c == 1) {
      for (int s = 1; s < nShifts; s++) {
        allSuccSizeC.push_back({s});  // Add it to the possibilities
        lastShiftSucc.push_back(s);      // Record its last shift
        nLastShiftSucc.push_back(1);      // Only 1 successive performed so far
        arcCostSucc.push_back(0);            // No succession ended yet
      }
    }
    // Larger but not last -> Extend the previous one by extending each of size c-1 in all possible ways
    //
    else if (c < CDMin_) {
      // For each short rotation of size c-1
      for (unsigned int i = 0; i < allowedShortSuccBySize_[c - 1].size(); i++) {
        vector<int> succ = allowedShortSuccBySize_[c - 1][i];
        int lastSh = succ.back();
        int nLast = nLastShiftOfShortSucc_[c - 1][i];
        double cost = baseArcCostOfShortSucc_[c - 1][i];
        // For each possible new shift s.t. the succession is allowed
        for (int newSh = 1; newSh < nShifts; newSh++) {
          if (!pScenario_->isForbiddenSuccessorShift_Shift(newSh, lastSh)) {

            vector<int> newSucc = succ;
            newSucc.push_back(newSh);    // Create Succession
            allSuccSizeC.push_back(newSucc);              // Add it to the possibilities
            lastShiftSucc.push_back(newSh);                // Record its last shift

            int newShTypeID = pScenario_->shiftIDToShiftTypeID_[newSh];
            int lastShTypeID = pScenario_->shiftIDToShiftTypeID_[lastSh];

            // Depending on the previous one, update number of consecutive and cost
            if (newShTypeID == lastShTypeID) {
              nLastShiftSucc.push_back(nLast + 1);
              arcCostSucc.push_back(cost);
            } else {
              nLastShiftSucc.push_back(1);
              arcCostSucc.push_back(cost + pScenario_->consShiftCost(lastSh, nLast));
            }

          }
        }
      }
    }
    // Maximum allowed size -> more things to consider
    //
    else {
      // For each short rotation of size c-1
      for (unsigned int i = 0; i < allowedShortSuccBySize_[c - 1].size(); i++) {
        vector<int> succ = allowedShortSuccBySize_[c - 1][i];
        int lastSh = succ.back();
        int nLast = nLastShiftOfShortSucc_[c - 1][i];
        double cost = baseArcCostOfShortSucc_[c - 1][i];
        // For each possible new shift s.t. the succession is allowed
        for (int newSh = 1; newSh < nShifts; newSh++) {
          if (!pScenario_->isForbiddenSuccessorShift_Shift(newSh, lastSh)) {

            vector<int> newSucc(succ);
            newSucc.push_back(newSh);    // Create Succession
            allSuccSizeC.push_back(newSucc);              // Add it to the possibilities
            lastShiftSucc.push_back(newSh);                // Record its last shift
            int newNLast = 1;
            double newCost = cost;


            int newShTypeID = pScenario_->shiftIDToShiftTypeID_[newSh];
            int lastShTypeID = pScenario_->shiftIDToShiftTypeID_[lastSh];


            if (newShTypeID == lastShTypeID) {  // BUT : add the cost if longer than the maximum allowed
              newNLast += nLast;
              if (newNLast >= pScenario_->maxConsShiftsOfTypeOf(newSh)) {
                newCost += pScenario_->consShiftCost(lastSh, nLast);
              }
            } else {
              newCost += pScenario_->consShiftCost(lastSh, nLast);
            }
            nLastShiftSucc.push_back(newNLast);
            arcCostSucc.push_back(newCost);

            // Since it is the longest one, record the tables the other way round
            //
            int n = std::min(newNLast, maxCons(newShTypeID));
            allShortSuccCDMinByLastShiftCons_[newShTypeID][n].push_back(allSuccSizeC.size() - 1);

          }
        }
      }
    }

    // Store all vectors
    //
    allowedShortSuccBySize_.push_back(allSuccSizeC);
    lastShiftOfShortSucc_.push_back(lastShiftSucc);
    nLastShiftOfShortSucc_.push_back(nLastShiftSucc);
    baseArcCostOfShortSucc_.push_back(arcCostSucc);

  }

}


//--------------------------------------------
//
// Solve function
//
//--------------------------------------------

// Solve : Returns TRUE if negative reduced costs path were found; FALSE otherwise.
bool ShortSP::preprocess() {
  // find best start of rotation
  priceShortSucc();
  // find short rotations
  return solveShortRotations();
}


// For the short rotations, depends on the chosen option + on whether we want optimality (more important)
bool ShortSP::solveShortRotations() {
  nVeryShortFound_ = 0;
  bool ANS = false;
  // 0 -> no short rotations
  // 1 -> day-0 short rotations only
  // 2 -> day-0 and last-day short rotations only
  // 3 -> all short rotations
  // SOLVE_SHORT_NONE
  if (param_.shortRotationsStrategy_ == 0) {}
  // SOLVE_SHORT_DAY_0_ONLY
  else if (param_.shortRotationsStrategy_ == 1) {
    bool tmp_ANS = priceVeryShortRotationsFirstDay();
    ANS = ANS or tmp_ANS;
  }
  // SOLVE_SHORT_DAY_0_AND_LAST_ONLY
  else if (param_.shortRotationsStrategy_ == 2) {
    bool tmp_ANS = priceVeryShortRotationsFirstDay();
    bool tmp_ANS2 = priceVeryShortRotationsLastDay();
    ANS = tmp_ANS2 or ANS or tmp_ANS;
  }
  // SOLVE_SHORT_ALL
  else if (param_.shortRotationsStrategy_ == 3) {
    bool tmp_ANS = priceVeryShortRotations();
    ANS = ANS or tmp_ANS;
  } else {
    std::cout << "# INVALID / OBSOLETE OPTION FOR SHORT ROTATIONS" << std::endl;
    getchar();
    return false;
  }
  return ANS;
}

// Create all arcs whose origin is the source nodes (all go to short rotations nodes)
void ShortSP::createArcsSourceToPrincipal() {
  int origin = g_.source();
  for (PrincipalGraph &pg: principalGraphs_) {
    int s = pScenario_->shiftTypeIDToShiftID_[pg.shiftType()].front(); // any shift with the right shit type
    std::vector<int> shifts;
    Tools::initVector(shifts, minConsDays_, s);
    for (int k = CDMin_ - 1; k < nDays_; k++)
      for(int dest : pg.getDayNodes(k))
        // add one arc for each cons days available
        // the shifts will be updated based on their dual costs
        arcsFromSource_[pg.shiftType()][k].push_back(
            {addSingleArc(origin, dest, 0, startConsumption(k, shifts),
                          SOURCE_TO_PRINCIPAL, k - CDMin_ + 1)});
      }
}

double ShortSP::startWorkCost(int a) const {
  return g_.arcCost(a); // cost already updated;
}


// Pricing of the short successions : only keep one of them, and the cost of the corresponding arc
//
void ShortSP::priceShortSucc() {

  map<int, int> specialArcsSuccId;
  map<int, double> specialArcsCost;

  for (int s = 1; s < pScenario_->nbShiftsType_; s++) {
    for (int k = minConsDays_ - 1; k < nDays_; k++) {
      int max_cons = principalGraphs_[s].maxCons();
      for (int n = 0; n <= max_cons; n++) {

        int best_id = -1;
        double best_cost = MAX_COST;

        // CHECK THE ROTATIONS ONLY IF THE FIRST DAY IS ALLOWED
        if (startingDayStatus_[k - CDMin_ + 1]) {

          for (unsigned int i = 0; i < (allShortSuccCDMinByLastShiftCons_[s][n]).size(); i++) {
            int curSuccId = allShortSuccCDMinByLastShiftCons_[s][n][i];
            const vector<int>& succ = allowedShortSuccBySize_[CDMin_][curSuccId];

            // SUCCESSION IS TAKEN INTO ACCOUNT ONLY IF IT DOES NOT VIOLATE ANY FORBIDDEN DAY-SHIFT COUPLE
            if (canSuccStartHere(k - CDMin_ + 1, succ)) {
              double curCost = costArcShortSucc(CDMin_, curSuccId, k - CDMin_ + 1);

              // ONLY CASE WHEN THE DESTINATION NODE MAY HAVE TO CHANGE:
              // 1. Start date is 0
              // 2. Size of short succession is < than the number of levels maxValByShift[s]
              // 3. Number of last shifts cons in succession is CDMin_
              // 4. The shift is the same as the last one worked by the nurse at initial state
              if (k == CDMin_ - 1 and CDMin_ < max_cons and n == CDMin_ and s == pLiveNurse_->pStateIni_->shiftType_) {
                // a. Determine the destination of that arc
                int nConsWithPrev = CDMin_ + pLiveNurse_->pStateIni_->consShifts_;
                int nDestination = std::min(nConsWithPrev, max_cons);
                int a = arcsFromSource_[s][k][nDestination].front();
                // b. Store the succession ID + the special cost for that arc
                specialArcsSuccId[a] =curSuccId;
                specialArcsCost[a] = curCost;
              }

                // OTHER CASES ("REGULAR ONES")
              else if (curCost < best_cost) {
                best_id = curSuccId;
                best_cost = curCost;
              }
            }
          }
        }

        // if a best short succession has been found, update arc
        int a = arcsFromSource_[s][k][n].front();
        if (best_id >= 0) {
          g_.updateCost(a, best_cost);
          g_.updateShifts(a, allowedShortSuccBySize_[CDMin_][best_id]);
        }
          // otherwise, forbid it
        else g_.forbidArc(a);

//        if(!startingDayStatus_[k-CDMin_+1] || best_cost >= MAX_COST-1
      }
    }
  }

  // FOR THE SHIFTS ON THE FIRST DAY THAT EXTEND THE ONGOING WORK AT INITIAL STATE
  // (It's the only existing short succession that could be used by this arc)
  for (auto p: specialArcsSuccId) {
    int a = p.first;
    g_.authorizeArc(a);
    g_.updateCost(a, specialArcsCost.at(a));
    g_.updateShifts(a, allowedShortSuccBySize_[CDMin_][p.second]);
  }
}

// Given a short succession and a start date, returns the cost of the corresponding arc
//
double ShortSP::costArcShortSucc(int size, int succId, int startDate) {
  double ANS = 0;
  const vector<int>& succ = allowedShortSuccBySize_[size][succId];

  vector<int> vvv = {4,2,2};
  if(startDate == 0 && size == 3) {
    int i=0;
    while(i<size && vvv[i] == succ[i]) { i ++; }
    if(i==size)
      i = -1;
  }

  // A. COST: BASE COST
  //
  ANS += baseArcCostOfShortSucc_[size][succId];


  // B. COST: SPECIAL CASE FOR THE FIRST DAY
  //
  if (startDate == 0) {

    int shiftTypeIni = pLiveNurse_->pStateIni_->shiftType_;
    int nConsWorkIni = pLiveNurse_->pStateIni_->consDaysWorked_;
    int nConsShiftIni = pLiveNurse_->pStateIni_->consShifts_;

    int firstShift = succ.front();
    int firstShiftType = pScenario_->shiftIDToShiftTypeID_[firstShift];
    int nConsFirstShift = 0;
    for (int i = 0; i < size && pScenario_->shiftIDToShiftTypeID_[succ[i]] == firstShiftType; i++)
      nConsFirstShift++;

    // 1. The nurse was resting: pay more only if the rest is too short
    if (shiftTypeIni == 0) {
      int diffRest = pLiveNurse_->minConsDaysOff() - pLiveNurse_->pStateIni_->consDaysOff_;
      if(diffRest>0) ANS += diffRest * WEIGHT_CONS_DAYS_OFF;
    }

    // 2. The nurse was working
    else {
      // a. (i)   The nurse was working on a different shift: if too short, add the corresponding cost
      if (shiftTypeIni != firstShiftType) {
        int diff = pScenario_->minConsShiftsOf(shiftTypeIni) - nConsShiftIni;
        if(diff>0) ANS += diff * WEIGHT_CONS_SHIFTS;
      }
      // a. (ii)  The nurse was working on the same shift AND the short rotation contains other shifts (easy case for add/subtract)
      //            - Subtract the cost due to the consecutive beginning if > max
      //            - Subtract the cost due to the consecutive end of the initial state if > max
      //            - Add the consecutive cost for all shifts
      else if (nConsFirstShift < size) {
        int diffShift = nConsShiftIni - pScenario_->maxConsShiftsOf(shiftTypeIni);
        if(diffShift>0) ANS -= diffShift * WEIGHT_CONS_SHIFTS; // max penalty paid in the previous week that will be repaid
        ANS -= pScenario_->consShiftTypeCost(firstShiftType, nConsFirstShift); // penalty contained in the base cost
        ANS += pScenario_->consShiftTypeCost(firstShiftType, (nConsFirstShift + nConsShiftIni));
      }
      // a. (iii) The nurse was working on the same shift AND the short rotation only contains that shift:
      // recompute the cost -> easier, just the max cons shift
      else {
        ANS = 0;
        if (size + nConsShiftIni >= pScenario_->maxConsShiftsOf(shiftTypeIni))
          ANS += pScenario_->consShiftTypeCost(shiftTypeIni, nConsFirstShift + nConsShiftIni);
      }

      // b. If the number of consecutive days worked has already exceeded the max, subtract now the cost that will be readd later
      int diffWork = nConsWorkIni - pContract_->maxConsDaysWork_;
      if(diffWork>0) ANS -= diffWork * WEIGHT_CONS_DAYS_WORK;
    }
  }

  // C. COST: COMPLETE WEEKEND
  //
  ANS += startWeekendCosts_[startDate];

  // D. COST: PREFERENCES
  //
  for (int i = 0; i < size; i++) ANS += preferencesCosts_[startDate + i][succ[i]];

  // E. REDCOST: WEEKENDS
  //
  int nbWeekends = Tools::containsWeekend(startDate, startDate + size - 1);
  ANS -= nbWeekends * pCosts_->workedWeekendCost();

  // F. REDCOST: FIRST DAY (BACK TO WORK)
  //
  ANS -= pCosts_->startWorkCost(startDate);

  // G. REDCOST: EACH DAY/SHIFT REDUCED COST
  //
  for (int i = 0; i < size; i++)
    ANS -= pCosts_->workedDayShiftCost(startDate + i, succ[i]);

  return ANS;
}

//----------------------------------------------------------------
//
// Cost computation of the "very" short rotations (< CD_min)
//
//----------------------------------------------------------------

// Brutally try all possible short rotations from the very first day
bool ShortSP::priceVeryShortRotationsFirstDay(){
  int nFound = 0;
  for(int c=1; c<CDMin_; c++)
    nFound += priceVeryShortSameSizeRotations(0, allowedShortSuccBySize_[c]);
  return nFound > 0;
}

// Brutally try all possible short rotations that end on the last day
bool ShortSP::priceVeryShortRotationsLastDay(){
  int nFound = 0;
  for(int c=1; c<CDMin_; c++)
      nFound += priceVeryShortSameSizeRotations(nDays_-c, allowedShortSuccBySize_[c]);
  return nFound > 0;
}

// Brutally try all possible short rotations from every first day
bool ShortSP::priceVeryShortRotations() {
  int nFound = 0;
  for (int c = 1; c < CDMin_; c++) {
    const vector2D<int> &succs = allowedShortSuccBySize_[c];
    for (int k = 0; k <= nDays_ - c; k++)
      nFound += priceVeryShortSameSizeRotations(k, succs);
  }
  return nFound > 0;
}

int ShortSP::priceVeryShortSameSizeRotations(int k,  const vector2D<int>& succs) {
  int nFound = 0;
  if (startingDayStatus_[k])
    for (const vector<int> &succ: succs) {
      double redCost = costOfVeryShortRotation(k, succ);
      if (redCost < maxReducedCostBound_) {
        theSolutions_.emplace_back(k, succ, redCost);
        nPaths_++;
        nVeryShortFound_++;
        nFound++;
        bestReducedCost_ = std::min(bestReducedCost_, redCost);
      }
    }
  return nFound;
}


// Compute the cost of a single short rotation
double ShortSP::costOfVeryShortRotation(int startDate, const vector<int>& succ) {

  int endDate = startDate + succ.size() - 1;

  //check if any shift is forbidden
  for(int k=startDate; k<=endDate; k++)
    if(isDayShiftForbidden(k, succ[k-startDate]))
      return MAX_COST;

  // Regular costs
  double consDaysRegCost=0, consShiftsRegCost=0, completeWeekendRegCost=0, preferencesRegCost=0, shortRestRegCost=0;
  // Reduced costs
  double dayShiftsRedCost=0, startRedCost=0, endRedCost=0, weekendRedCost=0;

  // Initialize values
  int shift=0, consShifts=0, consDays=succ.size();

  // A. SPECIAL CASE OF THE FIRST DAY
  //
  if(startDate==0){
    // The nurse was working
    if(pLiveNurse_->pStateIni_->shift_ > 0){

      if(pScenario_->isForbiddenSuccessorShift_Shift( succ[0], pLiveNurse_->pStateIni_->shift_)){
        return MAX_COST;
      }

      // Change initial values
      shift = pLiveNurse_->pStateIni_->shift_;
      consShifts = pLiveNurse_->pStateIni_->consShifts_;
      consDays += pLiveNurse_->pStateIni_->consDaysWorked_;
      // If worked too much, subtract the already counted surplus
      consDaysRegCost -= std::max(0, pLiveNurse_->pStateIni_->consDaysWorked_ - pContract_->maxConsDaysWork_) * WEIGHT_CONS_DAYS_WORK;
      consShiftsRegCost -= std::max(0, pLiveNurse_->pStateIni_->consShifts_ - pScenario_->maxConsShiftsOfTypeOf(shift)) * WEIGHT_CONS_SHIFTS;
    }
      // The nurse was resting
    else {
      // Cost of a too short rest
      shortRestRegCost += std::max(0, pContract_->minConsDaysOff_ - pLiveNurse_->pStateIni_->consDaysOff_) * WEIGHT_CONS_DAYS_OFF;
    }
  }

  // B. REGULAR COST: CONSECUTIVE NUMBER OF DAYS (only if it does not end on last day)
  //
  if((int) (startDate+succ.size()) < nDays_ and consDays){
    consDaysRegCost += pContract_->consDaysCost(consDays);
  }

  // C. REGULAR COST: CONSECUTIVE SHIFTS
  //
  for(int k=startDate; k<=endDate; k++){
    int newShift = succ[k-startDate];
    int shiftType = pScenario_->shiftIDToShiftTypeID_[shift];
    int newShiftType = pScenario_->shiftIDToShiftTypeID_[newShift];

    if(newShiftType == shiftType){
      consShifts ++;
    } else {
      if(shiftType>0)
        consShiftsRegCost += pScenario_->consShiftCost(shift,consShifts);
      consShifts = 1;
      shift = newShift;
    }
    if(k==endDate and (k<nDays_-1 or consShifts > pScenario_->maxConsShiftsOfTypeOf(shift)))
      consShiftsRegCost +=pScenario_-> consShiftCost(shift, consShifts);
  }

  // D. REGULAR COST: COMPLETE WEEKENDS
  //
  completeWeekendRegCost = startWeekendCosts_[startDate] + endWeekendCosts_[endDate];

  // E. REGULAR COST: PREFERENCES
  //
  for(int k=startDate; k<=endDate; k++) preferencesRegCost += preferencesCosts_[k][ succ[k-startDate] ];


  // F. REDUCED COST: WEEKENDS
  //
  weekendRedCost -= Tools::containsWeekend(startDate, endDate) * pCosts_->workedWeekendCost();

  // F. REDUCED COST: FIRST DAY (BACK TO WORK)
  //
  startRedCost -= pCosts_->startWorkCost(startDate);

  // G. REDUCED COST: LAST DAY (BACK TO WORK)
  //
  endRedCost -= pCosts_->endWorkCost(endDate);

  // H. REDUCED COST: EACH DAY/SHIFT REDUCED COST
  //
  for(int k=startDate; k<=endDate; k++) dayShiftsRedCost -= pCosts_->workedDayShiftCost(k, succ[k - startDate]);


  // I. RETURN THE TOTAL COST
  //
  double regCost = consDaysRegCost + consShiftsRegCost + completeWeekendRegCost + preferencesRegCost + shortRestRegCost;
  double redCost = dayShiftsRedCost + startRedCost + endRedCost + weekendRedCost;
  double ANS = regCost + redCost;

  if(false)
  {
    std::cout << "# " << std::endl;
    std::cout << "#+---------------------------------------------+ (" << pLiveNurse_->name_ << ")" << std::endl;
    std::cout << "# " << startDate << "-";
    for(unsigned int i=0; i<succ.size(); i++) std::cout << pScenario_->intToShift_[succ[i]].at(0);
    std::cout << std::endl;
    std::cout << "# length " << consDays;
    if(startDate == 0 and pLiveNurse_->pStateIni_->shift_ > 0){
      std::cout << " (incl. " << pLiveNurse_->pStateIni_->consDaysWorked_ << " before planing horizon)";
    }
    std::cout << std::endl;

    std::cout << "# REG- Consecutive days cost   : " << consDaysRegCost << std::endl;
    std::cout << "# REG- Consecutive shifts cost : " << consShiftsRegCost << std::endl;
    std::cout << "# REG- Complete weekends cost  : " << completeWeekendRegCost << std::endl;
    std::cout << "# REG- Preferences cost        : " << preferencesRegCost << std::endl;
    std::cout << "# REG- Short rest before cost  : " << shortRestRegCost << std::endl;
    std::cout << "# REG-                 ~TOTAL~ : " << regCost << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# RED- Day-shifts cost         : " << dayShiftsRedCost << std::endl;
    std::cout << "# RED- Start work cost         : " << startRedCost << std::endl;
    std::cout << "# RED- End work cost           : " << endRedCost << std::endl;
    std::cout << "# RED- Weekend dual cost       : " << weekendRedCost << std::endl;
    std::cout << "# RED-                 ~TOTAL~ : " << redCost << std::endl;
    std::cout << "#+---------------------------------------------+" << std::endl;
    std::cout << "#                      ~TOTAL~ : " << ANS << std::endl;
    std::cout << "#+---------------------------------------------+" << std::endl;
    std::cout << "# " << std::endl;
    //getchar();
  }

  return ANS;
}

// Summary of the short successions generated
void ShortSP::printShortSucc() const {
  std::cout << "#   +------------+" << std::endl;
  std::cout << "#   | CD_min = " << CDMin_ << std::endl;
  std::cout << "#   +------------+" << std::endl;
  int nShortSucc = 0;
  int nShortRot = 0;
  for(unsigned int i=0; i<allowedShortSuccBySize_.size(); i++){
    vector2D<int> v2 = allowedShortSuccBySize_[i];
    std::cout << "#   | " << v2.size()  << " short successions of size " << i << std::endl;
    nShortSucc += v2.size();
    nShortRot += v2.size() * (nDays_-i+1);

  }
  std::cout << "#   +------------+" << std::endl;
  std::cout << "#   | " << nShortSucc << std::endl;
  std::cout << "#   +------------+" << std::endl;
}

// Prints all active pairs ( arcFromSource - corresponding short successions )
void ShortSP::printShortArcs() const {
  for(int s=1; s<pScenario_->nbShiftsType_; s++){
    for(int k=minConsDays_-1; k<nDays_; k++){
      for(int n=1; n<=principalGraphs_[s].maxCons(); n++){
        const Arc_Properties& arc_prop = arcsFromSource_[s][k][n].front();
        if(!arc_prop.forbidden){
          int v = principalGraphs_[s].getNode(k,n);
          std::cout << "# " << g_.shortNameNode(v) << " <- ";
          for(int s: arc_prop.shifts)
            std::cout << pScenario_->intToShift_[s].at(0);
          std::cout << std::endl;
        }
      }
    }
  }
}


