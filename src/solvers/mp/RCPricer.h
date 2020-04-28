/*
 * RCPricer.h
 *
 * Allow to link the sub problems and the master problem trough scip
 *
 *  Created on: 2015-03-02
 *      Author: legraina
 */

#ifndef ROTATIONPRICER_H_
#define ROTATIONPRICER_H_

#include "tools/MyTools.h"
#include "MasterProblem.h"
#include "solvers/mp/rcspp/SubProblem.h"
#include "solvers/mp/modeler/Modeler.h"


//---------------------------------------------------------------------------
//
// C l a s s   R o t a t i o n P r i c e r
//
// Contains the pricer that generates new columns.
//
//---------------------------------------------------------------------------
class RCPricer: public MyPricer
{
public:
   RCPricer(MasterProblem* master, const char* name, const SolverParam& param);

   virtual ~RCPricer();

   /* perform pricing */
   std::vector<MyVar*> pricing(double bound=0, bool before_fathom = true);

   // Initialize parameters
   void initPricerParameters(const SolverParam& param);

protected:

   // DATA - instance-related data
   //
   MasterProblem* pMaster_;
   Scenario* pScenario_;
   int nbDays_;
   Modeler* pModel_;
    std::vector<LiveNurse*> nursesToSolve_;
   // One subproblem per contract because the consecutive same shift constraints vary by contract.
   std::map<const Contract*, SubProblem*> subProblems_;

   // DATA - Solutions, rotations, etc.
   //
   std::vector<MyVar*> allNewColumns_;
   std::vector<RCSolution> newSolutionsForNurse_;

   // Stats on the number of subproblems solved and successfully solved
   int nbSPTried_;
   int nbSPSolvedWithSuccess_ ;

   // SETTINGS - Options for forbidden shifts, nurses, starting days, etc.
   //
   std::set<std::pair<int,int> > forbiddenShifts_;
    std::set<int> forbiddenNursesIds_;
    std::set<int> forbiddenStartingDays_;
    std::set<int> forbiddenEndingDays_;

   // SETTINGS - Options for the neighborhood. need of an original to reset at the end of each node when optimality has
   //            been reached. Here, we could have all the parameters as fields but they would be too many.
   //
   bool withSecondchance_ = false;
   bool shortSubproblem_ = true;
   bool rosterSubproblem_ = false;
   int defaultSubprobemStrategy_ = 0;
   int secondchanceSubproblemStrategy_ = 0;
   int currentSubproblemStrategy_ = 0;

   // SETTINGS - Settings for the maximum number of problems to solve and of rotations to add to the master problem
   //
   int nbMaxColumnsToAdd_ = 0;
   int nbSubProblemsToSolve_ = 0;

public:

   // METHODS - Solutions, rotations, etc.
   //
   inline void resetSolutions(){
	   allNewColumns_.clear();
     newSolutionsForNurse_.clear();
	   forbiddenShifts_.clear();
	   nbSPSolvedWithSuccess_ = 0;
	   nbSPTried_ = 0;
   }

   // METHODS - Forbidden shifts, nurses, starting days, etc.
   //
   // !!! WARNING !!! : SOME METHODS ARE NOT YET IMPLEMENTED IN THE SUBPROBLEM (ALTHOUGH THE NECESSARY STRUCTURES MAY
   //                   ALREADY BE THERE !!!
   //
   // Shifts
   inline void forbidShift(int k, int s){forbiddenShifts_.insert(std::pair<int,int>(k,s));}
   inline void forbidShifts(const std::set<std::pair<int,int> > &shifts){ for(auto s : shifts) forbidShift(s.first, s.second);}
   inline void authorizeShift(int k, int s){forbiddenShifts_.erase(std::pair<int,int>(k,s));}
   inline void clearForbiddenShifts(){forbiddenShifts_.clear();}
   // Nurses
   inline void forbidNurse(int nurseId){forbiddenNursesIds_.insert(nurseId);}
   inline void forbidNurses(const std::set<int>& nurses){ for(auto n : nurses) forbidNurse(n);}
   inline void authorizeNurse(int nurseId){forbiddenNursesIds_.erase(nurseId);}
   inline void clearForbiddenNurses(){forbiddenNursesIds_.clear();}
   // Starting days
   inline void forbidStartingDay(int k){forbiddenStartingDays_.insert(k);}
   inline void forbidStartingDays(const std::set<int>& days){ for(int d : days) forbidStartingDay(d);}
   inline void authorizeStartingDay(int k){forbiddenStartingDays_.erase(k);}
   inline void clearForbiddenStartingDays(){forbiddenStartingDays_.clear();}
   // Ending days
   inline void forbidEndingDay(int k){forbiddenEndingDays_.insert(k);}
   inline void forbidEndingDays(const std::set<int> &days){ for(int d : days) forbidEndingDay(d);}
   inline void authorizeEndingDay(int k){forbiddenEndingDays_.erase(k);}
   inline void clearForbiddenEndingDays(){forbiddenEndingDays_.clear();}

   // Test functions
   inline bool isShiftForbidden(int k, int n){ return (forbiddenShifts_.find(std::pair<int,int>(k,n)) != forbiddenShifts_.end()); }
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
   vector2D<double> getShiftsDualValues(LiveNurse*  pNurse);
    std::vector<double> getStartWorkDualValues(LiveNurse* pNurse);
    std::vector<double> getEndWorkDualValues(LiveNurse* pNurse);
   double getWorkedWeekendDualValue(LiveNurse* pNurse);

   //compute some forbidden shifts from the lasts rotations and forbidden shifts
   void addForbiddenShifts();

   // Retrieve the right subproblem
   SubProblem* retriveSubproblem(LiveNurse* pNurse);

   // Add the rotations to the master problem
   int addColumnsToMaster(int nurseId);

   // Sort the rotations that just were generated for a nurse. Default option is sort by increasing reduced cost but we
   // could try something else (involving disjoint columns for ex.)
   void sortNewlyGeneratedSolutions();

   // Set the subproblem options depending on the parameters
//   void setSubproblemOptions(vector<SolveOption>& options, int& maxRotationLengthForSubproblem, pLiveNurse pNurse);

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
   std::minstd_rand rand_;
   // DBG functions
   void recordSPStats(SubProblem* sp);
   void generateRandomForbiddenStartingDays();
};

#endif /* ROTATIONPRICER_H_ */
