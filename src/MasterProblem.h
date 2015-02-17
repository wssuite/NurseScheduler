/*
 * MasterProblem.h
 *
 *  Created on: 3 févr. 2015
 *      Author: samuel
 */

#ifndef MASTERPROBLEM_H_
#define MASTERPROBLEM_H_

#include <cfloat>
#include "MyTools.h"


//-----------------------------------------------------------------------------
//
//  S t r u c t   R o t a t i o n
//
//  A rotation is a set of shifts for a set of consecutive days.
//  It has a cost and a dual cost (tbd).
//
//-----------------------------------------------------------------------------
struct Rotation {

	// Specific constructors and destructors
	//

	Rotation(map<int,int> shift, double cost = 999999, double dualCost = 999999) :
		cost_(cost), shifts_(shift), dualCost_(dualCost), length_(shift.size())
	{
		firstDay_ = 999;
		for(map<int,int>::iterator itS = shift.begin(); itS != shift.end(); ++itS)
			if(itS->first < firstDay_) firstDay_ = itS->first;
	};

	Rotation(int firstDay, vector<int> shiftSuccession, double cost = 999999, double dualCost = 999999) :
		cost_(cost), dualCost_(dualCost), firstDay_(firstDay), length_(shiftSuccession.size())
	{
		map<int,int> m;
		for(int k=0; k<shiftSuccession.size(); k++) m.insert(pair<int,int>( (firstDay+k) , shiftSuccession[k] ));
		shifts_ = m;
	}

	~Rotation(){};

	// Cost
	//
	double cost_;

	// Dual cost as found in the subproblem
	//
	double dualCost_;

	// Shifts to be performed
	//
	map<int,int> shifts_;

	// First worked day
	//
	int firstDay_;

	// Duration
	//
	int length_;

};

class MasterProblem : Solver{


public:


protected:


};



#endif /* MASTERPROBLEM_H_ */
