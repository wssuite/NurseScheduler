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

  this->preprocessTheNurses();

  // get the maximum rank
  //
  rankMax_ = 0;
  for (int i= 0; i < pScenario_->nbPositions(); i++)  {
    std::cout << "Rank = " << pScenario_->pPositions()[i]->rank() << std::endl;
    rankMax_ = std::max(rankMax_, pScenario_->pPositions()[i]->rank());
  }

  // copy the live nurses in the sorted nurses vector
  //
  for (int i=0; i < pScenario_->nbNurses_; i++) {
    theNursesSorted_.push_back(theLiveNurses_[i]);
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
double Greedy::costTask(const LiveNurse &nurse, int day, int shift, int skill) {

  double cost = 0.0;

  // Count with a positive value the penalties associated with the assignement
  // of this task
  //
  State state = nurse.states_[day];

  // Consecutive shifts
  //
  // penalize violation of maximum and minimum numbers of shifts
  int lastShift = state.shift_;
  int ecartShift, ecartOff, ecartDay;
  if (lastShift > 0 && shift == lastShift  && state.consShifts_ >= pScenario_->maxConsShifts_[shift]) {
    cost += WEIGHT_CONS_SHIFTS;
  }
  else if (lastShift > 0 && shift != lastShift && state.consShifts_ < pScenario_->minConsShifts_[shift]) {
    ecartShift = pScenario_->minConsShifts_[shift]-state.consShifts_;
    cost += WEIGHT_CONS_SHIFTS*ecartShift;
  }
  // penalize the action of taking a new shift whose number of minimum succesive
  // shifts will nont be reached before the maximum number of worked days
  if (lastShift > 0 && shift != lastShift) {
    int maxDays = nurse.maxConsDaysWork()-state.consDaysWorked_;
    if (maxDays < pScenario_->minConsShifts_[shift]) {
      cost += WEIGHT_CONS_SHIFTS*(pScenario_->minConsShifts_[shift] - maxDays);
    }
  }
  // penalize the action to go from a shift with a lot of possible successors to
  // a shift with more forbidden successors
  double weightNbForbidden = 5;
  cost += weightNbForbidden*
    (pScenario_->nbForbiddenSuccessors_[shift] - pScenario_->nbForbiddenSuccessors_[lastShift]);

  // Maximum working days and minimum resting days
  //
  if (lastShift == 0 && state.consDaysOff_ < nurse.minConsDaysOff()) {
    ecartOff = nurse.minConsDaysOff()-state.consDaysOff_;
    cost += ecartOff * WEIGHT_CONS_DAYS_OFF;
  }
  else if ( state.consDaysWorked_ >= nurse.maxConsDaysWork()) {
    cost += WEIGHT_CONS_DAYS_WORK;
  }

  // treat the cases in which the nurse had no assigned task on the last day
  if (lastShift < 0) {
    ecartOff = nurse.minConsDaysOff()-state.consDaysOff_;
    ecartDay = nurse.maxConsDaysWork()-state.consDaysWorked_;

    if (-lastShift >= ecartDay && -lastShift < ecartOff) {
      cost += std::min((-lastShift-ecartDay+1)*WEIGHT_CONS_DAYS_WORK,
        (ecartOff+lastShift) * WEIGHT_CONS_DAYS_OFF);
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
  if (nurse.needCompleteWeekends() && day%7 == 6) {
    if (lastShift == 0) {
      cost += WEIGHT_COMPLETE_WEEKEND;
    }
  }
  // when complete week-ends are needed and it is saturday, add a penalty if
  // working on sunday will create a violation in the maximum number of worked
  // days
  if (nurse.needCompleteWeekends() && day%7 == 5) {
    if (state.consDaysWorked_+1 == nurse.maxConsDaysWork() ) {
      cost += std::min(WEIGHT_COMPLETE_WEEKEND,WEIGHT_CONS_DAYS_WORK);
    }
  }

  // Penalize the rank of the nurse : the higher the rank the less we want it
  // to work on tasks other nurses can do
  //
  double weightRank = 5; // valeur arbitraire à régler en fonction des autres pénalités
  Position position = *(nurse.pPosition_);
  cost += (rankMax_-nurse.pPosition_->rank())*weightRank;

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
  maxDays = (int) (avgMaxDays+(1.0-factorWeek)*((double)nurse.maxTotalShifts()/(double) pScenario_->nbWeeks_-(double)(day%7)));
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

// Assign the unassigned nurses with best costs to the demand input tasks
// nbAssigned is the number of nurses that have actually obtained a new task
void Greedy::assignBestNursesToTask(int day, int sh, int sk, int demand,
  vector<LiveNurse*>& pNursesUnassigned, int &nbAssigned) {

  // initialize the indices and the costs of the nurses most interesting
  // for the task
  vector<double> costMin;
  vector<int> nMin;
  for (int i = 0; i < demand; i++) {
    costMin.push_back(1.0e06);
    nMin.push_back(-1);
  }
  int nbUnassigned = pNursesUnassigned.size();
  for (int n = 0; n < nbUnassigned; n++)  {
    LiveNurse* pNurse = pNursesUnassigned[n];
    // consider the nurse only if the task respects the hard constraints
    if (isFeasibleTask(*pNurse, day, sh, sk)) {
      double cost = costTask(*pNurse, day, sh, sk);

      // insert the nurse in the vector of interesting choices if its cost
      // deserves it
      if (cost < costMin.back()) {
        insertInVectors(cost, n, costMin, nMin, demand);

        // the minimum achievable cost is 0, so stop the loop here if the
        // the cost is zero
        if (costMin.back() == 0) break;
      }
    }
  }

  // assign the task to each selected nurse
  nbAssigned = demand;
  while (!nMin.empty()) {
    if (nMin.back() < 0) {
      nMin.pop_back();
      nbAssigned--;
    }
    else break;
  }
  for (int i=0; i < nbAssigned; i++) {
    LiveNurse* pNurse = pNursesUnassigned[nMin[i]];
    pNurse->roster_.assignTask(day, sh, sk);
    State nextState;
    nextState.addDayToState(pNurse->states_[day], sh);
    pNurse->states_[day+1] = nextState;
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
//----------------------------------------------------------------------------

void Greedy::constructiveGreedy() {

  int nbDays = pDemand_->nbDays_, nbShifts = pScenario_->nbShifts_;
  int nbSkills = pScenario_->nbSkills_, nbNurses = pScenario_->nbNurses_;
  vector3D satisfiedDemand;
  Tools::initVector3D(&satisfiedDemand, nbDays, nbShifts,nbSkills);

  // First satisfy the minimum demand
  //
  for (int day = 0; day < nbDays; day++) {
    // Initialize the set of nurses that are not assigned
    // we use this vector to avoid going through the nurses with a task
    vector<LiveNurse*> pNursesUnassigned;
    for (int n = 0; n < nbNurses; n++) {
      pNursesUnassigned.push_back(theLiveNurses_[n]);
    }
    int nbUnassigned = nbNurses;

    // RqJO : l'ordre des skills/shifts pourrait être changé pour commencer par
    // ceux qui sont les plus critiques
    for (int sh = 1; sh < nbShifts; sh++) { // recall that shift 0 is rest
      for (int sk = 0; sk < nbSkills; sk++) {
        // demand for the task
        int demand = pDemand_->minDemand_[day][sh][sk];
        if (demand <= 0) continue;

        int nbAssigned = 0;
        assignBestNursesToTask(day, sh, sk, demand, pNursesUnassigned, nbAssigned);
        if (nbAssigned <= 0) Tools::throwError("there is not enough nurse for the task!");
        satisfiedDemand[day][sh][sk] = nbAssigned;
        nbUnassigned -= nbAssigned;
      }
    }

    // Once all the tasks of the day are assigned, treat the nurses with no task
    // After this step, the nurses with no assigned task are those that are
    // free to either take another working shift or a rest without penalty
    //
    for (int n = 0; n < nbUnassigned; n++)  {
      LiveNurse* pNurse = pNursesUnassigned[n];
      // No task has been assigned to the nurse if its most recent state is that
      // of the current day
      State* pState = &pNurse->states_[day];
      if (pState->dayId_ == day) {
        // give a rest the nurses that need to rest to avoid penalties
        if ( pNurse->needRest(day) )  {
          pNurse->roster_.assignTask(day, 0);
          State nextState;
          nextState.addDayToState(*pState, 0);
          pNurse->states_[day+1] = nextState;
        }
        // assign the task that minimizes the penalties to the nurses that need
        // to work to avoid penalties
        else if ( pNurse->needWork(day) ) {
          double costMin = 1.0e6;
          int shMin = 0, skMin = 0;
          for (int sh = 1; sh < nbShifts; sh++) {
            // do not consider the shift if it creates a forbidden sequence
            if (pScenario_->isForbiddenSuccessor(sh,pState->shift_)) continue;
            for (int sk = 0; sk < nbSkills; sk++) {
              if( !isFeasibleTask(*pNurse, day, sh, sk) ) continue;
              double cost = costTask(*pNurse, day, sh, sk);
              cost -= WEIGHT_OPTIMAL_DEMAND
                *(pDemand_->minDemand_[day][sh][sk]-satisfiedDemand[day][sh][sk]);
              if (cost < costMin) {
                shMin = sh;
                skMin = sk;
                costMin = cost;
              }
            }
          }
          if (!shMin) Tools::throwError("there is no possible task for the nurse!");
          else {
            pNurse->roster_.assignTask(day, shMin, skMin);
            State nextState;
            nextState.addDayToState(*pState, shMin);
            pNurse->states_[day+1] = nextState;
            satisfiedDemand[day][shMin][skMin] ++;
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
  }

  // Try and satisfy the optimum demand with the unassigned nurses
  //
  for (int day = 0; day < nbDays; day++) {
    // initialize the set of nurses that are not assigned on this day yet
    vector<LiveNurse*> pNursesUnassigned;
    int nbUnassigned = 0;
    for (int n = 0; n < nbNurses; n++) {
      State state = theLiveNurses_[n]->states_[day+1];
      if (state.shift_ < 0) {
        pNursesUnassigned.push_back(theLiveNurses_[n]);
        nbUnassigned++;
      }
    }

    // RqJO : l'ordre des skills/shifts pourrait être changé pour commencer par
    // ceux qui sont les plus critiques
    for (int sh = 1; sh < nbShifts; sh++) { // recall that shift 0 is rest
      for (int sk = 0; sk < nbSkills; sk++) {
        // demand for the task
        int demand = pDemand_->optDemand_[day][sh][sk]-satisfiedDemand[day][sh][sk];
        if (demand <= 0) continue;

        // initialize the indices and the costs of the nurses most interesting
        // for the task
        int nbAssigned = 0;
        assignBestNursesToTask(day, sh, sk, demand, pNursesUnassigned, nbAssigned);
        std::cout << "Shift " << sh << ", skill " << sk << " : ";
        std::cout << nbAssigned << " nurses are assigned for an extra demand of " << demand << std::endl;
        satisfiedDemand[day][sh][sk] += nbAssigned;
        nbUnassigned -= nbAssigned;
      }
    }

    // Once all the tasks of the day are assigned, treat the nurses with no task
    // After this step, the nurses with no assigned task are those that are
    // free to either take another working shift or a rest without penalty
    //
    for (int n = 0; n < nbUnassigned; n++)  {
      LiveNurse* pNurse = pNursesUnassigned[n];
      // No task has been assigned to the nurse if its most recent state is that
      // of the current day
      State* pState = &pNurse->states_[day];
      if (pState->dayId_ == day) {
        // give a rest the nurses that need to rest to avoid penalties
        if ( pNurse->needRest(day) )  {
          pNurse->roster_.assignTask(day, 0);
          State nextState;
          nextState.addDayToState(*pState, 0);
          pNurse->states_[day+1] = nextState;
        }
        // assign the task that minimizes the penalties to the nurses that need
        // to work to avoid penalties
        else if ( pNurse->needWork(day) ) {
          double costMin = 1.0e6;
          int shMin = 0, skMin = 0;
          for (int sh = 1; sh < nbShifts; sh++) {
            // do not consider the shift if it creates a forbidden sequence
            if (pScenario_->isForbiddenSuccessor(sh,pState->shift_)) continue;
            for (int sk = 0; sk < nbSkills; sk++) {
              if( !isFeasibleTask(*pNurse, day, sh, sk) ) continue;
              double cost = costTask(*pNurse, day, sh, sk);
              cost -= WEIGHT_OPTIMAL_DEMAND
                *(pDemand_->minDemand_[day][sh][sk]-satisfiedDemand[day][sh][sk]);
              if (cost < costMin) {
                shMin = sh;
                skMin = sk;
                costMin = cost;
              }
            }
          }
          if (!shMin) Tools::throwError("there is no possible task for the nurse!");
          else {
            pNurse->roster_.assignTask(day, shMin, skMin);
            State nextState;
            nextState.addDayToState(*pState, shMin);
            pNurse->states_[day+1] = nextState;
            satisfiedDemand[day][shMin][skMin] ++;
          }
        }
      }
    }
  }

  for(LiveNurse* pNurse: theLiveNurses_)
     solution_.push_back(pNurse->roster_);
}



// Main method to solve the rostering problem for a given input
//
void Greedy::solve() {}


// compare functions that can be used to sort the nurse before assigning them
// schedules in the greedy algorithm
//
bool compareNurses(const LiveNurse  &n1, const LiveNurse &n2) {

  // the first parameter for ordering the nurses is their positions
  // if they have different position, the position priority of the position
  // is sufficient to compare the nurses
  //

  return true;
}

// Build the sequence of positions reflecting the order in which the positions
// will be treated in the greedy
//
void Greedy::sortPositions() {

  vector<Position*> pPositions = pScenario_->pPositions();

  //---------------------------------------------------------------------------
  // A rank has been computed for each position
  // The greedy is going to treat the positions in ascending rank value
  // For a given rank, start with the position that contains the rarest skills
  // Rarity is defined by the number of nurse shifts available for this skill
  //---------------------------------------------------------------------------

  // initialize the vector of position order
  //
  vector<int> order;
  for (int i= 0; i < pScenario_->nbPositions(); i++)  order.push_back(0);

  // compute the maximum value of rarity among the skills of each position
  //
  for (int i= 0; i < pScenario_->nbPositions(); i++)  {
    int rarity = 0;
    // compute the highest rarity for the skills of the position
    // it has to be done in the preprocessing function of Solver, because it
    // depends on the history
    //

    // if rarity is higher
  }

  // get the maximum rank
  //
  int rankMax;
  for (int i= 0; i < pScenario_->nbPositions(); i++)  {
    rankMax = std::max(rankMax,pPositions[i]->rank());
  }

  // go through the positions for each rank value and set their order of
  // treatment.
  // the positions are already sorted in descending rarity, so they order is
  // set directly
  //
  int nextTreated = 0;
  for (int rank = 0; rank <= rankMax; rank++) {
    for (int i= 0; i < pScenario_->nbPositions(); i++)  {
      if (pPositions[i]->rank() == rank)  {
        order[pPositions[i]->id_] = nextTreated++;
      }
    }
  }
}
