//
// GlobalStats.cpp
//  RosterDesNurses
//
//  Created by Jeremy Omer on 16/02/2016.
//  Copyright (c) 2016 Jeremy Omer. All rights reserved.
//


#include "GlobalStats.h"

//-----------------------------------------------------------------------------
//
//  S t r u c t   G l o b a l S t a t s
//
//  Set of statistics that are worth displaying in an article
//
//-----------------------------------------------------------------------------

// write all the stats
std::string GlobalStats::toString() {
	std::stringstream statStream;
	statStream.precision(1);
	statStream.setf( std::ios::fixed, std:: ios::floatfield );

	statStream << "final statistics= " << bestUBInitial_ << "\t" << bestUB_ << "\t";
	statStream << rootLB_ << "\t" << bestLB_ << "\t";
	statStream << (bestUBInitial_-bestLB_)/bestUBInitial_ << "\t" << (bestUB_-bestLB_)/bestUB_ << "\t";
	statStream << timeInitialSol_ << "\t" << timeImproveSol_ << "\t" << timeInitialSol_+timeImproveSol_ << "\t";
	statStream << timeGenColRoot_ << "\t" ;
	statStream << timeGenColMaster_ << "\t" ;
	statStream << timeGenSubProblems_ << "\t";
	statStream << itGenColInitial_ << "\t" << itGenColImprove_ << "\t";
	statStream << nodesBBInitial_ << "\t" << nodesBBImprove_ << "\t";

	return statStream.str();
}

// add the information of a stat object to this one
void GlobalStats::add(const GlobalStats& stats) {

	bestUBInitial_ += stats.bestUBInitial_;
	bestUB_ += stats.bestUB_;
	rootLB_ += stats.rootLB_;
	bestLB_ += stats.bestLB_;
	timeInitialSol_ += stats.timeInitialSol_;
	timeImproveSol_ += stats.timeImproveSol_;
	timeGenColRoot_ += stats.timeGenColRoot_;
	timeGenColMaster_ += stats.timeGenColMaster_;
	timeGenSubProblems_ += stats.timeGenSubProblems_;
	itGenColInitial_ += stats.itGenColInitial_;
	itGenColImprove_ += stats.itGenColImprove_;
	nodesBBInitial_ += stats.nodesBBInitial_;
	nodesBBImprove_ += stats.nodesBBImprove_;

	lnsImprovementValueTotal_+= stats.lnsImprovementValueTotal_;
	lnsNbIterations_+= stats.lnsNbIterations_;
	lnsNbIterationsWithImprovement_ += stats.lnsNbIterationsWithImprovement_;

	// Counters of iterations with improvements depending on repair/destroy
	if (!stats.nbImprovementsWithNursesSelection_.empty()) {
		if (nbImprovementsWithNursesSelection_.empty()) {
			nbImprovementsWithNursesSelection_.insert(nbImprovementsWithNursesSelection_.begin(),stats.nbImprovementsWithNursesSelection_.size(),0);
		}
		if (nbImprovementsWithDaysSelection_.empty()) {
			nbImprovementsWithDaysSelection_.insert(nbImprovementsWithDaysSelection_.begin(),stats.nbImprovementsWithDaysSelection_.size(),0);
		}
		if (nbImprovementsWithRepair_.empty()) {
			nbImprovementsWithRepair_.insert(nbImprovementsWithRepair_.begin(),stats.nbImprovementsWithRepair_.size(),0);
		}
		for (unsigned int i=0; i < stats.nbImprovementsWithNursesSelection_.size(); i++) {
			nbImprovementsWithNursesSelection_[i] += stats.nbImprovementsWithNursesSelection_[i];
		}
		for (unsigned int i=0; i < stats.nbImprovementsWithDaysSelection_.size(); i++) {
			nbImprovementsWithDaysSelection_[i] += stats.nbImprovementsWithDaysSelection_[i];
		}
		for (unsigned int i=0; i < stats.nbImprovementsWithRepair_.size(); i++) {
			nbImprovementsWithRepair_[i] += stats.nbImprovementsWithRepair_[i];
		}
	}
}

// Print to a string the statistics of the lns
//
std::string GlobalStats::lnsStatsToString() {
	std::stringstream statStream;
	statStream.precision(1);
	statStream.setf( std::ios::fixed, std:: ios::floatfield );

	statStream << "lns statistics: " << lnsImprovementValueTotal_ << " & ";
	statStream << lnsNbIterations_ << " & " << lnsNbIterationsWithImprovement_ << " &&";
	for (int it:nbImprovementsWithNursesSelection_) {
		statStream << " " << it << " &";
	}
	statStream << "&";
	for (int it:nbImprovementsWithDaysSelection_) {
		statStream << " " << it << " &";
	}
	statStream << "&";
	for (int it:nbImprovementsWithRepair_) {
		statStream << " " << it << " &";
	}
	statStream << "& ";

	return statStream.str();
}
