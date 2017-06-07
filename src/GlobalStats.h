//
// GlobalStats.h
//  RosterDesNurses
//
//  Created by Jeremy Omer on 16/02/2016.
//  Copyright (c) 2016 Jeremy Omer. All rights reserved.
//


#ifndef __GlobalStats__
#define __GlobalStats__

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "Solver.h"

//-----------------------------------------------------------------------------
//
//  S t r u c t   G l o b a l S t a t s
//
//  Set of statistics that are worth displaying in an article
//
//-----------------------------------------------------------------------------

struct GlobalStats{
public:
	// constructor and destructor
	//
	GlobalStats() {}
	~GlobalStats() {}

public:
	// status of the final solution: optimality can be proved only if the
	// optimality level is set to 3 in a solution phase that considers the
	// complete problem
	//
	Status status_=INFEASIBLE;

	// objective values for an algorithm in which one first phase is run to find
	// an initial solution and a second phase tries to improve it heuristically
	//
	double bestUB_=0.0;
	double bestUBInitial_=0.0;
	double rootLB_ = 0.0;
	double bestLB_ = 0.0;

	// global runtimes
	//
	double timeTotal_=0.0;
	double timeInitialSol_=0.0;
	double timeImproveSol_=0.0;

	// Details on Branch and price related runtimes
	//
	double timeGenColRoot_=0.0;
	double timeGenColMaster_=0.0;
	double timeGenSubProblems_=0.0;

	// details on Branch and price iterations
	//
	int itGenColInitial_=0;
	int itGenColImprove_=0;
	int nodesBBInitial_=0;
	int nodesBBImprove_=0;

	// number of problems solved in each phase
	//
	int itInitialSol_=0;
	int itImproveSol_=0;

	// STATISTICS ON THE LNS
	// Global stats
	double lnsImprovementValueTotal_=0.0;
	int lnsNbIterations_=0;
	int lnsNbIterationsWithImprovement_=0;

	// Counters of iterations with improvements depending on repair/destroy
	std::vector<int> nbImprovementsWithNursesSelection_;
	std::vector<int> nbImprovementsWithDaysSelection_;
	std::vector<int> nbImprovementsWithRepair_;

public:
	// write all the stats
	std::string toString();

	// Print to a string the statistics of the lns
	std::string lnsStatsToString();

	// add the information of a stat object to this one
	void add(const GlobalStats& stats);
};

#endif
