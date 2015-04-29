/*
* Greedy.cpp
*
*  Created on: 21 jan. 2015
*      Author: jeremy
*/

#include <algorithm>

#include "Greedy.h"

// method that insert a value and the associated index in vectros of values and
// indices so that the value vector is sorted in ascending order
void insertInVectors(double value, int index,
  vector<double>& valVect, vector<int>& indVect, int maxLength) {
  for (int i = 0; i < maxLength; i++) {
    if (value < valVect[i]) {
      valVect.insert(valVect.begin()+i, value);
      valVect.pop_back();
      indVect.insert(indVect.begin()+i, index);
      indVect.pop_back();
      break;
    }
  }
}

struct indexcost{
  int index;
  double cost;
};

bool compareCosts(indexcost ic1, indexcost ic2) {
  return ic1.cost < ic2.cost;
}


//-----------------------------------------------------------------------------
//
//  C l a s s   G r e e d y
//
//  Quick solution of the problem with a greedy
//
//-----------------------------------------------------------------------------

// Specific constructor
Greedy::Greedy(Scenario* pScenario, Demand* pDemand,
  Preferences* pPreferences, vector<State>* pInitState):
  Solver(pScenario, pDemand, pPreferences, pInitState) {

  // get the maximum rank
  //
  rankMax_ = 0;
  for (int i= 0; i < pScenario_->nbPositions(); i++)  {
    rankMax_ = std::max(rankMax_, pScenario_->pPositions()[i]->rank());
  }

  // initialize the sorted vectors with the input values
  //
  for (int i=0; i < pScenario_->nbNurses_; i++) {
    theNursesSorted_.push_back(theLiveNurses_[i]);
  }

  for (int sk=0; sk < pScenario_->nbSkills_; sk++) {
    skillsSorted_.push_back(sk);
  }

  for (int sh=1; sh < pScenario_->nbShifts_; sh++) {
    shiftsSorted_.push_back(sh);
  }

  // initialize weights that are used to penalize future violations
  //
  weightNbForbidden_ = (double) 1* WEIGHT_CONS_SHIFTS;
  weightRank_ =  (double) 1*WEIGHT_CONS_SHIFTS;
  weightCoverMin_ = (double)0*WEIGHT_OPTIMAL_DEMAND;

  // initialize the vector of excess in available nurses for a task with respect
  // to minimum demand
  //
  Tools::initVector3D(&shiftDemand_, pDemand_->nbDays_, pScenario_->nbShifts_, pScenario_->nbSkills_);
  for (int day = 0; day < pDemand_->nbDays_; day++){
    for (int sh = 1; sh < pScenario_->nbShifts_; sh++) {
      for (int sk = 0; sk < pScenario_->nbSkills_; sk++) {
        shiftDemand_[day][sh][sk] = -pDemand_->minDemand_[day][sh][sk];
        for (int n = 0; n < pScenario_->nbNurses_; n++) {
          LiveNurse* pNurse = theLiveNurses_[n];
          if (isFeasibleTask(*pNurse, day, sh, sk)) shiftDemand_[day][sh][sk]++;
        }
      }
    }
  }

}

//----------------------------------------------------------------------------
// Protected intermediate methods for the constructive greedy
//----------------------------------------------------------------------------

// Returns true if the input nurse will respect the hard constraints if she is
// assigned the input task
//
bool Greedy::isFeasibleTask(const LiveNurse &nurse, int day, int shift, int skill)  {
  // Check that the nurse has the assigned skill
  //
  if (!nurse.hasSkill(skill)) return false;

  // Check the forbidden successor constraint
  //
  int lastShift = nurse.states_[day].shift_;   // last shift assigned to the nurse
  if (pScenario_->isForbiddenSuccessor(shift, lastShift)) return false;

  return true;
}

// Method that return the cost for completing the input task with the input
// nurse
// The cost depends on the state of the nurse, but the method will not check
// the feasibility of the task
//
double Greedy::costTask(const LiveNurse &nurse, int day, int shift, int skill,
  vector<State>* states) {

  double cost = 0.0;

  // Get the current state
  //
  State state;
  if (states == NULL) {
    state = nurse.states_[day];
  }
  else {
    int iTmp = -1;
    for (int i = 0; i < (*states).size(); i++) {
      if ( (*states)[i].dayId_ == day) {
        state = (*states)[i];
        iTmp = i;
        break;
      }
    }
    if (iTmp < 0) {
      Tools::throwError("costTask: no state is chosen!");
    }

  }

  // Consecutive shifts
  //
  // penalize violation of maximum and minimum numbers of shifts
  int lastShift = state.shift_;
  int ecartShift, ecartOff, ecartDay;
  if (lastShift > 0 && shift == lastShift  && state.consShifts_ >= pScenario_->maxConsShifts_[lastShift]) {
    cost += WEIGHT_CONS_SHIFTS;
  }
  else if (lastShift > 0 && shift != lastShift && state.consShifts_ < pScenario_->minConsShifts_[lastShift]) {
    ecartShift = pScenario_->minConsShifts_[lastShift]-state.consShifts_;
    cost += WEIGHT_CONS_SHIFTS*ecartShift;
  }
  // penalize the action of taking a new shift whose number of minimum succesive
  // shifts will not be reached before the maximum number of worked days
  if (lastShift > 0 && shift != lastShift) {
    int maxDays = nurse.maxConsDaysWork()-state.consDaysWorked_;
    if (maxDays < pScenario_->minConsShifts_[shift]) {
      cost += WEIGHT_CONS_SHIFTS*(pScenario_->minConsShifts_[shift] - maxDays);
    }
  }
  // penalize the action to go from a shift with a lot of possible successors to
  // a shift with more forbidden successors
  if (lastShift > 0 && shift > 0) {
    cost += weightNbForbidden_*
      (pScenario_->nbForbiddenSuccessors_[shift] - pScenario_->nbForbiddenSuccessors_[lastShift]);
  }
  else if (lastShift == 0 && shift > 0) {
    cost += weightNbForbidden_*(pScenario_->nbForbiddenSuccessors_[shift]);

  }

  // Maximum and minimum working and resting days resting days
  //
  if (shift > 0) {
    if ( lastShift == 0 && state.consDaysOff_ < nurse.minConsDaysOff()) {
      ecartOff = nurse.minConsDaysOff()-state.consDaysOff_;
      cost += ecartOff * WEIGHT_CONS_DAYS_OFF;
    }
    else if (lastShift > 0 &&  state.consDaysWorked_ >= nurse.maxConsDaysWork()) {
      cost += WEIGHT_CONS_DAYS_WORK;
    }
  }
  else if (shift == 0) {
    if (lastShift == 0 && state.consDaysOff_ >= nurse.maxConsDaysOff() ) {
      cost += WEIGHT_CONS_DAYS_OFF;
    }
    else if (lastShift > 0 && state.consDaysWorked_ < nurse.minConsDaysWork()) {
      ecartShift = nurse.minConsDaysWork()-state.consDaysWorked_;
      cost += ecartShift * WEIGHT_CONS_DAYS_WORK;
    }
  }

  // treat the cases in which the nurse had no assigned task on the last day
  // this mentions that a new assignment before the minimum duration of a rest
  // and after the maximum duration of work will generate costs
  //
  if (lastShift < 0) {
    ecartOff = nurse.minConsDaysOff()-state.consDaysOff_;
    ecartDay = nurse.maxConsDaysWork()-state.consDaysWorked_;
    if (shift > 0) {
      if (-lastShift >= ecartDay && -lastShift < ecartOff) {
        cost += std::min((-lastShift-ecartDay+1)*WEIGHT_CONS_DAYS_WORK,
          (ecartOff+lastShift) * WEIGHT_CONS_DAYS_OFF);
      }
    }
  }


  // Preferences
  //
  if (nurse.wishesOff(day,shift)) {
    cost += WEIGHT_PREFERENCES;
  }

  // Complete week-ends
  //
  // penalize uncomplete week-ends when needed
  if (day%7 == 6 && nurse.needCompleteWeekends() ) {
    if ( (lastShift == 0 && shift > 0) || (lastShift > 0 && shift == 0) ) {
      cost += WEIGHT_COMPLETE_WEEKEND;
    }
    else if (lastShift < 0 && shift > 0) {
      int missingShifts = std::max(0, nurse.minConsDaysOff() - (-lastShift-1));
      int extraShifts = std::max(0, -lastShift+state.consDaysWorked_+1-nurse.maxConsDaysWork());
      if (missingShifts > 0 && extraShifts > 0) {
        cost += WEIGHT_COMPLETE_WEEKEND;
      }
    }
  }
  // when complete week-ends are needed and it is saturday, add a penalty if
  // working on sunday will create a violation in the maximum number of working
  // or resting days
  if ( day%7 == 5 && nurse.needCompleteWeekends() ) {
    if ( shift > 0 && state.consDaysWorked_+1 >= nurse.maxConsDaysWork() ) {
      cost += std::min(WEIGHT_COMPLETE_WEEKEND,WEIGHT_CONS_DAYS_WORK);
    }
    else if (shift == 0 && state.consDaysOff_+1 >= nurse.maxConsDaysOff()) {
      cost += std::min(WEIGHT_COMPLETE_WEEKEND,WEIGHT_CONS_DAYS_OFF);
    }
  }

  // Penalize the rank of the nurse : the higher the rank the less we want it
  // to work on tasks other nurses can do
  //
  Position position = *(nurse.pPosition_);
  cost += weightRank_ * (rankMax_-nurse.pPosition_->rank());

  // Penalize the violations of the limits in the total assignment and total
  // working week-ends
  // This must take into account the position of the day in the complete
  // scheduling process with an extra stochastic penalty if the nurse is
  // clearly below or above the average number of working days and working
  // week-ends (the penalty is not counted twice for the week-ends if the day
  // is sunday and saturday has been worked)
  //

  // get the average target number of working days for this week
  double avgMinDays, avgMaxDays, avgMaxWeekends, minDays, maxDays, maxWeekends;
  double factorWeek = (double)pScenario_->thisWeek()/(double) pScenario_->nbWeeks_;
  avgMinDays = factorWeek * (double) nurse.minTotalShifts();
  avgMaxDays = factorWeek * (double) nurse.maxTotalShifts();
  avgMaxWeekends = factorWeek*(double) nurse.maxTotalWeekends();
  minDays = (int) (avgMinDays-(1.0-factorWeek)*((double)nurse.minTotalShifts()/(double) pScenario_->nbWeeks_-(double)(day%7)));
  maxDays = (int) (avgMaxDays+(1.0-factorWeek)*((double)nurse.maxTotalShifts()/(double) pScenario_->nbWeeks_-(double)(day%7))) +1;
  maxWeekends = (int) avgMaxWeekends;

  if (state.totalDaysWorked_ >= maxDays) {
    cost += factorWeek*WEIGHT_TOTAL_SHIFTS;
  }
  if ( (day%7==5 || (day%7==6 && lastShift==0)) && state.totalWeekendsWorked_ >= maxWeekends) {
    cost += factorWeek*WEIGHT_TOTAL_WEEKENDS;
  }
  // if the nurse is below the number of total assignments, getting a shift is
  // rewarded with a negative cost
  // RqJO: not sure yet

  // Add a small penalty when the considered skill is not the rarest in the
  // skill list of the nurse
  // double maxRarity = skillRarity_[skill];
  // double arbitraryWeight = 5;
  // for (int sk = 0; sk < nurse.nbSkills_; sk++) {
  //   maxRarity = std::max(skillRarity_[sk],maxRarity);
  // }
  // cost += (maxRarity/skillRarity_[skill]-1)*arbitraryWeight;


  return cost;
}

// Necessary actions when assigning a task to a nurse
//
void Greedy::assignTaskToNurse(LiveNurse &nurse, int day, int shift, int skill) {

    nurse.roster_.assignTask(day, shift, skill);
    State nextState;
    nextState.addDayToState(nurse.states_[day], shift);
    nurse.states_[day+1] = nextState;

    // increment the vector of satisfied demand
    if (shift > 0) {
      satisfiedDemand_[day][shift][skill]++;
    }

    // fill the gaps in the roster if the last day has no assignment
    if (nurse.states_[day].shift_ < 0) {
      fillTheGaps(nurse, day);
      nextState.addDayToState(nurse.states_[day], shift);
      nurse.states_[day+1] = nextState;
    }
}

// Recursive function that tries to add rest or work to the input statesBlock
// and returns the best result if nbUnassigned days are treated
//
double Greedy::bestStatesBlock_rec(LiveNurse &nurse, vector<State> &states,
  vector<int> &shifts, vector<int> &skills, int day,int nbUnassigned, double costIni) {

    if (nbUnassigned <= 0) {
      Tools::throwError("bestStatesBlock_rec: The block is already complete!");
    }

    // copy the state, shift and skill vectors to get the two possible assignments
    // states1 will contain a rest on next day and states2 will contain a worked
    // day on next day
    //
    vector<State> states1 = states;
    vector<State> states2 = states;
    vector<int> shifts1 = shifts;
    vector<int> shifts2 = shifts;
    vector<int> skills1 = skills;
    vector<int> skills2 = skills;

    // A rest state is going to be appended to states1 and the best work state
    // will be appended to states2
    //
    State stateWork, stateRest;


    // Compute the cost of the rest, create and append the associated state
    //
    double costRest = costTask(nurse, day, 0, 0, &states);
    stateRest.addDayToState(states1.back(), 0);
    states1.push_back(stateRest);
    shifts1.push_back(0);
    skills1.push_back(-1);

    // Get the cost of the cheaper task that could be assigned to the nurse
    // Create the associated state and append it to states2
    //
    double costWork = 1.0e6;
    double costWorkTmp = 1.0e6;
    int shMin = 0, skMin = 0;
    for (int sh = 1; sh < pScenario_->nbShifts_; sh++) {
      // do not consider the shift if it creates a forbidden sequence with the
      // shift
      if (pScenario_->isForbiddenSuccessor(sh,states.back().shift_)) continue;

      // also reject shifts that will create forbidden sequence with the shift at
      // the end of the block of unassigned day if the end of the block is not
      // too far way
      if (std::max(0,nurse.minConsDaysWork()-states.back().consDaysWorked_)
          + nurse.minConsDaysOff() > nbUnassigned) {
        if (day+nbUnassigned < nurse.firstDay_+nurse.nbDays_) {
          if (pScenario_->isForbiddenSuccessor(nurse.roster_.shift(day+nbUnassigned),sh)) continue;
        }
      }
      for (int i = 0; i < nurse.nbSkills_; i++) {
        int sk = nurse.skills_[i];
        if( !isFeasibleTask(nurse, day, sh, sk) ) continue;
        double cost = costTask(nurse, day, sh, sk, &states);

        // only consider a low weight on demand to avoid taking a bad shift just
        // to satisfy a demand that may be covered by another nurse anyway
        // RqJO: this is a part of the code that would benefitiate from an
        // improvement
        double weightDemand = 10;//WEIGHT_OPTIMAL_DEMAND;
        cost -=
          (pDemand_->optDemand_[day][sh][sk]>satisfiedDemand_[day][sh][sk]) ?
          weightDemand: 0;

        // need to take into account right away the transition with the shift that
        // has been assigned at the end of the unassigned block
        State stateTmp;
        double costTmp= cost;
        if ( (day+nbUnassigned < nurse.firstDay_+nurse.nbDays_)
          && (sh != nurse.roster_.shift(day+nbUnassigned)) ) {
          stateTmp.addDayToState(states.back(), sh);
          int missingShifts = pScenario_->minConsShifts_[sh]-(nbUnassigned-1 + stateTmp.consShifts_);
          if (missingShifts > 0) {
            costTmp += missingShifts * WEIGHT_CONS_SHIFTS;
          }
        }
        if (costTmp < costWorkTmp) {
          shMin = sh;
          skMin = sk;
          costWork = cost;
          costWorkTmp = costTmp;
        }
      }
    }
    if (!shMin) {
      std::cout << "there is no possible task for the nurse!" << std::endl;
    }
    else {
      State stateWork;
      stateWork.addDayToState(states.back(), shMin);
      states2.push_back(stateWork);
      shifts2.push_back(shMin);
      skills2.push_back(skMin);
    }

    // To get the best total cost that can ensue from states1 and states2, call
    // recursively the function, unless all days have been assigned
    //
    double cost1, cost2;
    if (nbUnassigned-1 >= 1) {
      cost1 = bestStatesBlock_rec(nurse, states1, shifts1, skills1, day+1, nbUnassigned-1, costIni+costRest);
      cost2 = costWork == 1.0e6 ? 1.0e6 :
        bestStatesBlock_rec(nurse, states2, shifts2, skills2, day+1, nbUnassigned-1, costIni+costWork);
    }
    // If every day of the block has been assigned, compute the cost of the
    // shift just after the block and add it to the the cost
    else {
      cost1 = costIni+costRest;
      cost2 = costIni+costWork;
      // only add the cost of the next day if we did not reach the last day of
      // the considered period
      if (day+nbUnassigned < nurse.firstDay_+nurse.nbDays_) {
        int sh = nurse.roster_.skill(day+1);
        int sk = nurse.roster_.shift(day+1);
        cost1 += costTask(nurse, day+1, sh, sk, &states1);
        cost2 += costWork == 1.0e6 ? 0 : costTask(nurse, day+1, sh, sk, &states2);
      }
    }

    // Finally return update the states with best option Rest or Work
    // and return the best cost
    // Note: when the function reaches this part of the code, the vectors
    // states1 and states2 cover the whole block and cost1 and cost2 are the
    // costs for these complete roster block
    //
    if (cost1 <= cost2) {
      states = states1;
      shifts = shifts1;
      skills = skills1;
      return cost1;
    }
    else {
      states = states2;
      shifts = shifts2;
      skills = skills2;
      return cost2;
    }
  }

// When assigning a new task to a nurse that has an unassigned day just
// just before, find tasks/rest periods to assign in the preceeding block of
// unassinged days
// The input day shift are the shift about to be assigned to the nurse for the
// day
//
void Greedy::fillTheGaps(LiveNurse &nurse, int day) {

  int nbUnassigned = -nurse.states_[day].shift_;
  int dayFirst = day -nbUnassigned;// fist day of the block of unassigned

  // every possible sequence of assignments is tested and the best is kept
  // call a recursive function to perform the tests
  // the block of states is initialized with the first unassigned day
  vector<State> statesBlock;
  vector<int> shifts, skills;
  statesBlock.push_back(nurse.states_[dayFirst]);
  double costMin = bestStatesBlock_rec(nurse,statesBlock,shifts,skills,dayFirst,nbUnassigned, 0);

  // assign the best sequence of shifts and update the satisfied demand
  for (int i = 0; i < nbUnassigned; i++) {
    nurse.roster_.assignTask(dayFirst+i, shifts[i], skills[i]);
    nurse.states_[dayFirst+i+1] = statesBlock[i+1];
    if (shifts[i] > 0) {
      satisfiedDemand_[dayFirst+i][shifts[i]][skills[i]]++;
    }
  }
}

// Assign the unassigned nurses with best costs to the demand input tasks
// nbAssigned is the number of nurses that have actually obtained a new task
void Greedy::assignBestNursesToTask(int day, int sh, int sk, int demand,
  vector<LiveNurse*>& pNursesUnassigned, int &nbAssigned, bool isMinDemand) {

  // indices of the nurses with the minimum cost
  vector<int> nMin;

  int nbUnassigned = pNursesUnassigned.size();
  vector<indexcost> indexcostvect;

  for (int n=0; n < std::max(demand,nbUnassigned); n++) {
    indexcostvect.push_back({n,1.0e6});
  }

  for (int n = 0; n < nbUnassigned; n++)  {
    LiveNurse* pNurse = pNursesUnassigned[n];

    // consider the nurse only if the task respects the hard constraints
    if (isFeasibleTask(*pNurse, day, sh, sk)) {
      double cost = costTask(*pNurse, day, sh, sk);

      indexcostvect[n].cost = cost;
    }
  }

  // add a cost that takes into account the impossibility to cover the minimum
  // demand for the next day if we choose the wrong nurses to cover the demand
  if (isMinDemand && day+1 < pDemand_->nbDays_) {
    std::sort (indexcostvect.begin(), indexcostvect.end(), compareCosts);

    vector3D shiftDemandTmp = shiftDemand_;
    for (int n = 0; n < nbUnassigned; n++) {

      if (indexcostvect[n].cost >= 1.0e6) break;

      LiveNurse* pNurse = pNursesUnassigned[indexcostvect[n].index];
      bool createsInfeasibility = false;

      for (int shift = 1; shift < sh; shift++) {
        for (int i = 0; i < pNurse->nbSkills_; i++) {
          int skill = pNurse->skills_[i];
          if (shiftDemandTmp[day+1][shift][skill] <= 0) {
            createsInfeasibility = true;
            indexcostvect[n].cost += weightCoverMin_;
          }
        }
      }

      if (!createsInfeasibility) {
        for (int shift = 1; shift < sh; shift++) {
          for (int i = 0; i < pNurse->nbSkills_; i++) {
            int skill = pNurse->skills_[i];
            shiftDemandTmp[day+1][shift][skill]--;
          }
        }
      }
    }
  }


  // sort the cost vector to choose the least expensive nurses to cover the
  // demand
  std::sort (indexcostvect.begin(), indexcostvect.end(), compareCosts);

  // assign the task to the least expensive nurses
  for (int n=0; n < demand; n++) {
    if (indexcostvect[n].cost >= 1.0e6) break;

    LiveNurse* pNurse = pNursesUnassigned[indexcostvect[n].index];
    assignTaskToNurse(*pNurse, day, sh, sk);
    nbAssigned++;
    nMin.push_back(indexcostvect[n].index);

    // update the shift between the minimum demand and the available nurses
    if (isMinDemand && day+1 < pDemand_->nbDays_) {
      for (int shift = 1; shift < sh; shift++) {
        for (int i = 0; i < pNurse->nbSkills_; i++) {
          int skill = pNurse->skills_[i];
          shiftDemand_[day+1][shift][skill]--;
        }
      }
    }
  }


  // delete the nurses with a task
  // first sort nMin by ascending order to delete the higher indices
  // first (necessary to maintain the indices of the other nurses to
  // erase)
  std::sort(nMin.begin(), nMin.end());
  while (!nMin.empty()) {
    pNursesUnassigned.erase(pNursesUnassigned.begin()+nMin.back());
    nMin.pop_back();
  }
}


//----------------------------------------------------------------------------
// Constructive greedy algorithm
//
// Goes through the the demands in a chronoligcal order and assign the nurse
// that seems most appropriate to each task (shift/skill)
// Returns true if the minimum demand could be covered, and false otherwise
//
//----------------------------------------------------------------------------

bool Greedy::constructiveGreedy() {

  int nbDays = pDemand_->nbDays_, firstDay = pDemand_->firstDay_;
  int nbShifts = pScenario_->nbShifts_;
  int nbSkills = pScenario_->nbSkills_, nbNurses = pScenario_->nbNurses_;

  // initialize the algorithm by sorting the input data
  // the order of construction of the solution is essential in a constructive
  // greedy
  this->initializeConstructive();

  // First satisfy the minimum demand
  //
  for (int day = firstDay; day < firstDay+nbDays; day++) {
    // Initialize the set of nurses that are not assigned
    // we use this vector to avoid going through the nurses with a task
    vector<LiveNurse*> pNursesUnassigned;
    for (LiveNurse* pNurse:theNursesSorted_) {
      pNursesUnassigned.push_back(pNurse);
    }
    int nbUnassigned = nbNurses;

    // RqJO : l'ordre des skills/shifts pourrait etre change pour commencer par
    // ceux qui sont les plus critiques
    for (int sh:shiftsSorted_) { // recall that shift 0 is rest
      for (int sk:skillsSorted_) {

        // demand for the task
        int demand = pDemand_->minDemand_[day][sh][sk];
        if (demand <= 0) continue;

        int nbAssigned = 0;
        assignBestNursesToTask(day, sh, sk, demand, pNursesUnassigned, nbAssigned, true);
        if (nbAssigned <= 0) {
          std::cerr << "Day " << day << " ; shift " << pScenario_->intToShift_[sh];
          std::cerr <<" ; skill " << pScenario_->intToSkill_[sk];
          std::cerr << ": the demand cannot be covered.\n";
          return false;
        }
        nbUnassigned -= nbAssigned;
      }
    }

    // Once all the tasks of the day are assigned, treat the nurses with no task
    // After this step, the nurses with no assigned task are those that are
    // free to either take another working shift or a rest without penalty
    //
    for (LiveNurse* pNurse:pNursesUnassigned)  {
      State* pState = &pNurse->states_[day];

      // give a rest the nurses that need to rest to avoid penalties
      if ( pNurse->needRest(day) )  {
        assignTaskToNurse(*pNurse, day, 0, -1);
      }
      // assign the task that minimizes the penalties to the nurses that need
      // to work to avoid penalties
      else if ( false ) {//pNurse->needWork(day) ) {
        double costMin = 1.0e6;
        int shMin = 0, skMin = 0;

        for (int sh:shiftsSorted_) {
          // do not consider the shift if it creates a forbidden sequence
          if (pScenario_->isForbiddenSuccessor(sh,pState->shift_)) continue;
          for (int sk:skillsSorted_) {

            if( !isFeasibleTask(*pNurse, day, sh, sk) ) continue;
            double cost = costTask(*pNurse, day, sh, sk);
            cost -= WEIGHT_OPTIMAL_DEMAND
              *(pDemand_->minDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk]);
            if (cost < costMin) {
              shMin = sh;
              skMin = sk;
              costMin = cost;
            }
          }
        }
        if (!shMin) Tools::throwError("there is no possible task for the nurse!");
        else {
          assignTaskToNurse(*pNurse, day, shMin, skMin);
          pNurse->roster_.assignTask(day, shMin, skMin);
        }
      }
      // the other nurses simply get an unassigned state
      else {
        pNurse->roster_.assignTask(day, -1);
        State nextState;
        nextState.addDayToState(*pState, -1);
        pNurse->states_[day+1] = nextState;
      }
    }
  }

  // ------------------------------------------------------------------------
  // Try and satisfy the optimum demand with the unassigned nurses
  // At this stage the order of nurses/shifts/skills may not be relevant any
  // more
  // ------------------------------------------------------------------------

  // initialize weights that are used to penalize future violations
  weightNbForbidden_ = 0;
  weightRank_ = 0;
  weightCoverMin_ = 0;

  // go through the days
  for (int day = firstDay; day < firstDay+nbDays; day++) {
    // initialize the set of nurses that are not assigned on this day yet
    vector<LiveNurse*> pNursesUnassigned;
    int nbUnassigned = 0;
    for (LiveNurse* pNurse:theNursesSorted_) {
      State state = pNurse->states_[day+1];
      if (state.shift_ < 0) {
        pNursesUnassigned.push_back(pNurse);
        nbUnassigned++;
      }
    }

    // l'ordre des shifts/skills n'a peut-
    for (int sh : shiftsSorted_) { // recall that shift 0 is rest
      for (int sk: skillsSorted_) {
        // demand for the task
        int demand = pDemand_->optDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk];
        if (demand <= 0) continue;

        // initialize the indices and the costs of the nurses most interesting
        // for the task
        int nbAssigned = 0;
        assignBestNursesToTask(day, sh, sk, demand, pNursesUnassigned, nbAssigned, false);
        nbUnassigned -= nbAssigned;
      }
    }

    // Once all the tasks of the day are assigned, treat the nurses with no task
    // After this step, the nurses with no assigned task are those that are
    // free to either take another working shift or a rest without penalty
    //
    for (LiveNurse* pNurse: pNursesUnassigned)  {
      State* pState = &pNurse->states_[day];

      // give a rest the nurses that need to rest to avoid penalties
      if ( pNurse->needRest(day) )  {
        assignTaskToNurse(*pNurse, day, 0, -1);
      }
      // assign the task that minimizes the penalties to the nurses that need
      // to work to avoid penalties
      else if ( pNurse->needWork(day) ) {
        double costMin = 1.0e6;
        int shMin = 0, skMin = 0;
        for (int sh: shiftsSorted_) {
          // do not consider the shift if it creates a forbidden sequence
          if (pScenario_->isForbiddenSuccessor(sh,pState->shift_)) continue;
          for (int sk:skillsSorted_) {
            if( !isFeasibleTask(*pNurse, day, sh, sk) ) continue;
            double cost = costTask(*pNurse, day, sh, sk);
            cost -= WEIGHT_OPTIMAL_DEMAND
              *(pDemand_->minDemand_[day][sh][sk]-satisfiedDemand_[day][sh][sk]);
            if (cost < costMin) {
              shMin = sh;
              skMin = sk;
              costMin = cost;
            }
          }
        }
        if (!shMin) Tools::throwError("there is no possible task for the nurse!");
        else {
          assignTaskToNurse(*pNurse, day, shMin, skMin);
        }
      }

      // the other nurses simply get an unassigned state
      else {
        pNurse->roster_.assignTask(day, -1);
        State nextState;
        nextState.addDayToState(*pState, -1);
        pNurse->states_[day+1] = nextState;
      }
    }

    // If the last day of the demand is reached, fill the gaps in the end of the
    // rosters of the nurses that do not have an assignment every day in the
    // considered period
    //
    if (day == firstDay+nbDays-1) {
      for (LiveNurse* pNurse: theNursesSorted_) {
        if (pNurse->states_[day+1].shift_ < 0) {
          fillTheGaps(*pNurse, day+1);
        }
      }
    } // end fill the gaps

  }

  for(LiveNurse* pNurse: theLiveNurses_) {
     solution_.push_back(pNurse->roster_);
   }
   return true;
}



// Main method to solve the rostering problem for a given input
//
void Greedy::solve() {}

//------------------------------------------------------------------------------
// Create the vector of sorted nurses
// The nurses are ordered according to their position and the nurses that have
// the same position are shuffled
//------------------------------------------------------------------------------

void Greedy::sortShuffleTheNurses() {

  vector<Position*> positionsSorted;
  vector<vector<LiveNurse*>> nursePerPosition;

  // first, sort the position in the order in which the nurses should be treated
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    positionsSorted.push_back(pScenario_->pPosition(p));
  }
  std::sort(positionsSorted.begin(),positionsSorted.end(),comparePositions);

  // then organize the nurses depending on their position
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    vector<LiveNurse*> emptyVector;
    nursePerPosition.push_back(emptyVector);
  }
  for (int n=0; n < pScenario_->nbNurses_; n++) {
    LiveNurse* pNurse = theLiveNurses_[n];
    nursePerPosition[pNurse->pPosition()->id()].push_back(pNurse);
  }

  // shuffle the nurses that have the same position
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    std::random_shuffle(nursePerPosition[p].begin(),nursePerPosition[p].end());
  }

  // fill the sorted vector of live nurses
  theNursesSorted_.clear();
  for (int p=0; p < pScenario_->nbPositions(); p++) {
    int id = positionsSorted[p]->id();
    for (int n=0; n < nursePerPosition[id].size(); n++) {
      theNursesSorted_.push_back(nursePerPosition[id][n]);
    }
  }

  // todo: think about sorting the nurses according to their contract, we might
  // want to use the full-time nurses first
}

//------------------------------------------------------------------------------
// Initialize the greedy by preprocessing all the input attributes and sorting
// the shifts, skills, nurses
//------------------------------------------------------------------------------

void Greedy::initializeConstructive() {

  // Preprocess the attributes of the greedy solver
  // the result of the preprocessing will be very useful to sort the attributes
  // before greedily covering the demand
  //
  if (!pDemand_->isPreprocessed_) {
    pDemand_->preprocessDemand();
  }
  this->preprocessTheNurses();
  this->preprocessTheSkills();

  // sort the skills
  SkillSorter compareSkills(skillRarity_);
  std::stable_sort(skillsSorted_.begin(),skillsSorted_.end(),compareSkills);

  // sort the shifts (except the shift 0 which must always be rest)
  ShiftSorter compareShifts(pScenario_->nbForbiddenSuccessors_);
  std::stable_sort(shiftsSorted_.begin(),shiftsSorted_.end(),compareShifts);

  // sort the nurses
  sortShuffleTheNurses();
}
