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

	Rotation(map<int,int> shift, double cost = 999999, double dualCost = 999999) : cost_(cost), shifts_(shift), dualCost_(dualCost) {};

	Rotation(int firstDay, vector<int> shiftSuccession){
		map<int,int> m;
		for(int k=0; k<shiftSuccession.size(); k++) m.insert(pair<int,int>( (firstDay+k) , shiftSuccession[k] ));
		shifts_ = m;
		cost_ = 999999;
		dualCost_ = 999999;
	}

	~Rotation(){};

	// Cost
	//
	double cost_;

	// Dual cost as found in the subproblem
	//
	const double dualCost_;

	// Shifts to be performed
	//
	const map<int,int> shifts_;

};

class MasterProblem : Solver{


public:


protected:


};



#endif /* MASTERPROBLEM_H_ */
