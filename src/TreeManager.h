/*
 * TreeManager.h
 *
 *  Created on: Dec 30, 2015
 *      Author: tonio
 */

#ifndef SRC_TREEMANAGER_H_
#define SRC_TREEMANAGER_H_

#include "MyTools.h"
#include "Modeler.h"
#include "MasterProblem.h"

// Node that correspond to branching on a set of columns (diving)
//
struct ColumnsNode: public MyNode{
	ColumnsNode(int index, MyNode* pParent, std::vector<Rotation>& rotations):
		MyNode(index, pParent), rotations_(rotations) {}

	string write() {
		stringstream out;
		out << "ColumnsNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",#Cols=" << rotations_.size() << ")";
		return out.str();
	}

	//vector of the columns on which we have branched.
	const vector<Rotation> rotations_;
};

// Node that correspond to branching on a resting day
//
struct RestNode: public MyNode{
	RestNode(int index, MyNode* pParent, LiveNurse* pNurse, int day, bool rest, vector<MyVar*>& restArcs):
		MyNode(index, pParent), pNurse_(pNurse), day_(day), rest_(rest), restArcs_(restArcs) {}

	string write() {
		stringstream out;
		out << "RestNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse=" << pNurse_->id_ << ",Day=" << day_ << ",";
		if(rest_) cout << "rest";
		else cout << "work";
		out << ")";
		return out.str();
	}

	//nurse and day on which we have branched for rest or work. pNurse_ can be 0
	const LiveNurse* pNurse_;
	const int day_;
	const bool rest_;
	vector<MyVar*> restArcs_;
};


// Node that correspond to branching on working shifts
//
struct ShiftNode: public MyNode{
	ShiftNode(int index, MyNode* pParent, LiveNurse* pNurse, int day, bool work, vector<int>& forbiddenShifts):
		MyNode(index, pParent), pNurse_(pNurse), day_(day), work_(work), forbiddenShifts_(forbiddenShifts) {}

	string write() {
		stringstream out;
		out << "ShiftNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse=" << pNurse_->id_ << ",Day=" << day_ << ",";
		if(work_) out << "work";
		else out << "N.A.";
		out << "). Forbidden shifts:";
		for(int s: forbiddenShifts_)
			out << " " << s;
		return out.str();
	}

	//nurse and day on which we have branched for rest or work. pNurse_ can be 0
	const LiveNurse* pNurse_;
	const int day_;
	const bool work_;
	vector<int> forbiddenShifts_;
};

// Node that corresponds to branching on a penalized original variable
// The variable can relate to cover constraints, total shifts or total week-ends
//
struct PenaltyNode: public MyNode{
	PenaltyNode(int index, MyNode* pParent, LiveNurse* pNurse, int day, bool work, vector<int>& forbiddenShifts):
		MyNode(index, pParent), pNurse_(pNurse), day_(day), work_(work), forbiddenShifts_(forbiddenShifts) {}

	string write() {
		stringstream out;
		out << "PenaltyNide: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse=" << pNurse_->id_ << ",Day=" << day_ << ",";
		if(work_) out << "work";
		else out << "N.A.";
		out << "). Forbidden shifts:";
		for(int s: forbiddenShifts_)
			out << " " << s;
		return out.str();
	}

	//nurse and day on which we have branched for rest or work. pNurse_ can be 0
	const LiveNurse* pNurse_;
	const int day_;
	const bool work_;
	vector<int> forbiddenShifts_;
};


struct NursesNumberNode: public MyNode{
	NursesNumberNode(int index, MyNode* pParent, MyVar* var, double lb, double ub):
		MyNode(index, pParent), pNumberOfNurses_(var), nursesLhs_(lb), nursesRhs_(ub) {}

	string write() {
		stringstream out;
		out << "NursesNumberNode: (depth=" << depth_ << ",LB=" << bestLB_;
		out << ",Var=" << pNumberOfNurses_->name_ << ",LB=" << nursesLhs_ << ",UB=" << nursesRhs_ << ")";
		return out.str();
	}

	//number of nurse on which we have branched. pNumberOfNurses_ can be 0
	const MyVar* pNumberOfNurses_;
	const double nursesLhs_, nursesRhs_;
};

struct RestTree: public MyTree{
	RestTree(Scenario* pScenario, Demand* pDemand);

	void logical_fixing();

	void addForbiddenShifts(LiveNurse* pNurse, set<pair<int,int> >& forbidenShifts);

	inline void pushBackNewNursesNumberNode(MyVar* var, double lb, double ub){
		NursesNumberNode* node = new NursesNumberNode(tree_.size(), currentNode_, var, lb, ub);
		pushBackNode(node);
	}

	inline void pushBackNewColumnsNode(std::vector<Rotation>& rotations){
		ColumnsNode* node = new ColumnsNode(tree_.size(), currentNode_, rotations);
		pushBackNode(node);
	}

	inline void pushBackNewRestNode(LiveNurse* pNurse, int day, bool rest, vector<MyVar*>& restArcs){
		RestNode* node = new RestNode(tree_.size(), currentNode_, pNurse, day, rest, restArcs);
		pushBackNode(node);
	}

	inline void pushBackNewShiftNode(LiveNurse* pNurse, int day, bool work, vector<int>& shifts){
		ShiftNode* node = new ShiftNode(tree_.size(), currentNode_, pNurse, day, work, shifts);
		pushBackNode(node);
	}

	inline bool is_columns_node() { return dynamic_cast<ColumnsNode*>(currentNode_)!=0; }

	inline bool continueDiving(){
		if(is_columns_node()) return true;
		return (diveDepth_ <= max(min_depth_+10, diveLength_));
	}

	void reset();

	void updateStats(MyNode* node);

	string writeBranchStats();

	string writeOneStat(string name, vector<pair<int,double>>& stats);

protected:
	Scenario* pScenario_;
	Demand* pDemand_;
	//Tree statistics
	//each pair represents first the number of times of the branching and second the increase of the LB
	vector<pair<int,double>> statsRestByDay_, statsWorkByDay_, statsRestByNurse_, statsWorkByNurse_;
	pair<int,double> statsCols_;
};

struct ColumnsComparator {
	virtual bool is_disjoint(Rotation& rot1, Rotation& rot2)=0;
};

struct DayDisjointComparator: public ColumnsComparator{
	bool is_disjoint(Rotation& rot1, Rotation& rot2) { return rot1.isDisjointWith(rot2); }
};

struct ShiftDisjointComparator: public ColumnsComparator{
	bool is_disjoint(Rotation& rot1, Rotation& rot2) { return rot1.nurseId_!=rot2.nurseId_ && rot1.isShiftDisjointWith(rot2); }
};

class DiveBranchingRule: public MyBranchingRule
{
public:
   DiveBranchingRule(MasterProblem* master, RestTree* tree, const char* name);
   virtual ~DiveBranchingRule() { }

   /* compute branching decisions */
   bool branching_candidates(MyBranchingCandidate& candidate);

   /* branch on a set of resting arcs */
   bool branchOnRestingArcs(MyBranchingCandidate& candidate);

   /* branch on a set of shifts */
   bool branchOnShifts(MyBranchingCandidate& candidate);

	/* branch on a set of original penalty variables */
   bool branchOnPenalty(MyBranchingCandidate& candidate);

   /* compute fixing decisions */
   bool column_candidates(MyBranchingCandidate& candidate);

   /* Choose columns */
   vector<MyVar*> chooseColumns(vector<pair<MyVar*,double>>& candidates, vector<Rotation>& rotations, double& maxValue, ColumnsComparator& comparator);

   /* compare columns */
   static bool compareColumnCloseToInt(pair<MyVar*, double> obj1, pair<MyVar*, double> obj2);

   static bool compareColumnCloseTo5(pair<MyVar*, double> obj1, pair<MyVar*, double> obj2);

protected:
   //Pointer to the master problem to link the master and the sub problems
   //
   MasterProblem* pMaster_;

   //Pointer to the rest tree
   //
   RestTree* tree_;

   //pointers to the data
   //
   Modeler* pModel_;
};

#endif /* SRC_TREEMANAGER_H_ */
