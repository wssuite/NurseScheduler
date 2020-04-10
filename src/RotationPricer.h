/*
 * RotationPricer.h
 *
 * Allow to link the sub problems and the master problem trough scip
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#ifndef ROTATIONPRICER_H_
#define ROTATIONPRICER_H_

#include "MyTools.h"
#include "MasterProblem.h"
#include "SubProblem.h"
#include "Modeler.h"

/* namespace usage */
using namespace std;



//---------------------------------------------------------------------------
//
// C l a s s   R o t a t i o n P r i c e r
//
// Contains the pricer that generates new columns.
//
//---------------------------------------------------------------------------
class RotationPricer: public MyPricer
{
public:
   RotationPricer(MasterProblem* master, const char* name, SolverParam param);
   virtual ~RotationPricer();

   /* perform pricing */
   vector<MyVar*> pricing(double bound=0, bool before_fathom = true);

   // Initialize parameters
   void initPricerParameters(SolverParam param);

protected:

   // DATA - instance-related data
   //
   MasterProblem* pMaster_;
   Scenario* pScenario_;
   int nbDays_;
   Modeler* pModel_;
   vector<LiveNurse*> nursesToSolve_;
   // One subproblem per contract because the consecutive same shift constraints vary by contract.
   map<const Contract*, SubProblem*> subProblems_;

   // DATA - Solutions, rotations, etc.
   //
   vector<MyVar*> allNewColumns_;
   vector<Rotation> newRotationsForNurse_;
   int nbRotationsAddedToMasterSoFar_;

   // Stats on the number of subproblems solved and successfully solved
   int nbSPTried_;
   int nbSPSolvedWithSuccess_ ;

   // SETTINGS - Options for forbidden shifts, nurses, starting days, etc.
   //
   set<pair<int,int> > forbiddenShifts_;
   set<int> forbiddenNursesIds_;
   set<int> forbiddenStartingDays_;
   set<int> forbiddenEndingDays_;

   // SETTINGS - Options for the neighborhood. need of an original to reset at the end of each node when optimality has
   //            been reached. Here, we could have all the parameters as fields but they would be too many.
   //
   bool withSecondchance_ = false;
   bool shortSubproblem_ = true;
   int defaultSubprobemStrategy_ = 0;
   int secondchanceSubproblemStrategy_ = 0;
   int currentSubproblemStrategy_ = 0;

   // SETTINGS - Settings for the maximum number of problems to solve and of rotations to add to the master problem
   //
   int nbMaxRotationsToAdd_ = 0;
   int nbSubProblemsToSolve_ = 0;

public:

   // METHODS - Solutions, rotations, etc.
   //
   inline void resetSolutions(){
	   allNewColumns_.clear();
	   newRotationsForNurse_.clear();
	   forbiddenShifts_.clear();
	   nbRotationsAddedToMasterSoFar_ = 0;
	   nbSPSolvedWithSuccess_ = 0;
	   nbSPTried_ = 0;
   }

   // METHODS - Forbidden shifts, nurses, starting days, etc.
   //
   // !!! WARNING !!! : SOME METHODS ARE NOT YET IMPLEMENTED IN THE SUBPROBLEM (ALTHOUGH THE NECESSARY STRUCTURES MAY
   //                   ALREADY BE THERE !!!
   //
   // Shifts
   inline void forbidShift(int k, int s){forbiddenShifts_.insert(pair<int,int>(k,s));}
   inline void forbidShifts(set<pair<int,int> > shifts){ for(auto s : shifts) forbidShift(s.first, s.second);}
   inline void authorizeShift(int k, int s){forbiddenShifts_.erase(pair<int,int>(k,s));}
   inline void clearForbiddenShifts(){forbiddenShifts_.clear();}
   // Nurses
   inline void forbidNurse(int nurseId){forbiddenNursesIds_.insert(nurseId);}
   inline void forbidNurses(set<int,int> nurses){ for(auto n : nurses) forbidNurse(n);}
   inline void authorizeNurse(int nurseId){forbiddenNursesIds_.erase(nurseId);}
   inline void clearForbiddenNurses(){forbiddenNursesIds_.clear();}
   // Starting days
   inline void forbidStartingDay(int k){forbiddenStartingDays_.insert(k);}
   inline void forbidStartingDays(set<int> days){ for(auto d : days) forbidStartingDay(d);}
   inline void authorizeStartingDay(int k){forbiddenStartingDays_.erase(k);}
   inline void clearForbiddenStartingDays(){forbiddenStartingDays_.clear();}
   // Ending days
   inline void forbidEndingDay(int k){forbiddenEndingDays_.insert(k);}
   inline void forbidEndingDays(set<int> days){ for(auto d : days) forbidEndingDay(d);}
   inline void authorizeEndingDay(int k){forbiddenEndingDays_.erase(k);}
   inline void clearForbiddenEndingDays(){forbiddenEndingDays_.clear();}

   // Test functions
   inline bool isShiftForbidden(int k, int n){ return (forbiddenShifts_.find(pair<int,int>(k,n)) != forbiddenShifts_.end()); }
   inline bool isNurseForbidden(int n){ return (forbiddenNursesIds_.find(n) != forbiddenNursesIds_.end()); }
   inline bool isStartingDayForbidden(int k){ return (forbiddenStartingDays_.find(k) != forbiddenStartingDays_.end()); }
   inline bool isEndingDayForbidden(int k){ return (forbiddenEndingDays_.find(k) != forbiddenEndingDays_.end()); }



protected:
   /*
    * Methods
    */

   // Methods for the exhaustive / nonexhaustive search strategies
   //
//   inline void resetSearchParamToOriginal(){ currentPricerParam_ = originalPricerParam_; }
//   inline void authorizeExhaustiveSearch(){ currentPricerParam_.setToExhaustiveSearch();}



   //get the duals values per day and per shift for a nurse
   //
   vector< vector<double> > getWorkDualValues(LiveNurse* pNurse);
   vector<double> getStartWorkDualValues(LiveNurse* pNurse);
   vector<double> getEndWorkDualValues(LiveNurse* pNurse);
   double getWorkedWeekendDualValue(LiveNurse* pNurse);

   //compute some forbidden shifts from the lasts rotations and forbidden shifts
   void addForbiddenShifts();

   // Retrieve the right subproblem
   SubProblem* retriveSubproblem(LiveNurse*);

   // Add the rotations to the master problem
   void addRotationsToMaster();

   // Sort the rotations that just were generated for a nurse. Default option is sort by increasing reduced cost but we
   // could try something else (involving disjoint columns for ex.)
   void sortNewlyGeneratedRotations();

   // Set the subproblem options depending on the parameters
//   void setSubproblemOptions(vector<SolveOption>& options, int& maxRotationLengthForSubproblem, LiveNurse* pNurse);

   int nb_int_solutions_;

   // Print functions
   //
   void print_current_solution_();
   void printStatSPSolutions();




   // ----------------------------------------
   //
   // DBG - DEBUG & STATS DATA AND FUNCTIONS
   //
   // ----------------------------------------
   // DBG / Stats
   double timeInExSubproblems_ = 0;
   double timeForS_ = 0;
   double timeForN_ = 0;
   double timeForNL_ = 0;
   int nbSubproblems_ = 0;
   int nbExSubproblems_ = 0;
   int nbS_ = 0;
   int nbN_ = 0;
   int nbNL_ = 0;
   // DBG functions
   void recordSPStats(SubProblem* sp);
   void generateRandomForbiddenStartingDays();
   void checkForbiddenStartingDays();


};

#endif /* ROTATIONPRICER_H_ */
