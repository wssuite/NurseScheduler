/*
 * TreeManager.h
 *
 *  Created on: Dec 30, 2015
 *      Author: tonio
 */

#ifndef SRC_TREEMANAGER_H_
#define SRC_TREEMANAGER_H_

#include "tools/MyTools.h"
#include "solvers/mp/modeler/Modeler.h"
#include "solvers/mp/MasterProblem.h"

// Node that correspond to branching on a set of columns (diving)
//
struct ColumnsNode: public MyNode{
	ColumnsNode(int index, MyNode* pParent, std::vector<PPattern>& patterns):
		MyNode(index, pParent), patterns_(patterns) {}

    std::string write() const {
      std::stringstream out;
		out << "ColumnsNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",#Cols=" << patterns_.size() << ")";
		return out.str();
	}

	//vector of the columns on which we have branched.
	const std::vector<PPattern> patterns_;
};

// Node that correspond to branching on a resting day
//
struct RestNode: public MyNode{
	RestNode(int index, MyNode* pParent, PLiveNurse pNurse, int day, bool rest, std::vector<MyVar*>& restArcs):
		MyNode(index, pParent), pNurse_(pNurse), day_(day), rest_(rest), restArcs_(restArcs) {}

    std::string write() const {
      std::stringstream out;
		out << "RestNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse=" << pNurse_->id_ << ",Day=" << day_ << ",";
		if(rest_) std::cout << "rest";
		else std::cout << "work";
		out << ")";
		return out.str();
	}

	//nurse and day on which we have branched for rest or work. pNurse_ can be 0
	const PLiveNurse pNurse_;
	const int day_;
	const bool rest_;
    std::vector<MyVar*> restArcs_;
};


// Node that correspond to branching on working shifts
//
struct ShiftNode: public MyNode{
	ShiftNode(int index, MyNode* pParent, PLiveNurse pNurse, int day, bool work, std::vector<int>& forbiddenShifts):
		MyNode(index, pParent), pNurse_(pNurse), day_(day), work_(work), forbiddenShifts_(forbiddenShifts) {}

    std::string write() const {
      std::stringstream out;
		out << "ShiftNode: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse=" << pNurse_->id_ << ",Day=" << day_ << ",";
		if(work_) out << "work";
		else out << "N.A.";
		out << "). Forbidden shifts:";
		for(int s: forbiddenShifts_)
			out << " " << s;
		return out.str();
	}

	//nurse and day on which we have branched for rest or work. pNurse_ can be 0
	const PLiveNurse pNurse_;
	const int day_;
	const bool work_;
    std::vector<int> forbiddenShifts_;
};

// Node that corresponds to branching on a penalized original variable
// The variable can relate to cover constraints, total shifts or total week-ends
//
struct PenaltyNode: public MyNode{
	PenaltyNode(int index, MyNode* pParent, PLiveNurse pNurse, int day, bool work, std::vector<int>& forbiddenShifts):
		MyNode(index, pParent), pNurse_(pNurse), day_(day), work_(work), forbiddenShifts_(forbiddenShifts) {}

    std::string write() const {
      std::stringstream out;
		out << "PenaltyNide: (depth=" << depth_ << ",LB=" << bestLB_ << ",Nurse=" << pNurse_->id_ << ",Day=" << day_ << ",";
		if(work_) out << "work";
		else out << "N.A.";
		out << "). Forbidden shifts:";
		for(int s: forbiddenShifts_)
			out << " " << s;
		return out.str();
	}

	//nurse and day on which we have branched for rest or work. pNurse_ can be 0
	const PLiveNurse pNurse_;
	const int day_;
	const bool work_;
    std::vector<int> forbiddenShifts_;
};


struct NursesNumberNode: public MyNode{
	NursesNumberNode(int index, MyNode* pParent, MyVar* var, double lb, double ub):
		MyNode(index, pParent), pNumberOfNurses_(var), nursesLhs_(lb), nursesRhs_(ub) {}

    std::string write() const {
      std::stringstream out;
		out << "NursesNumberNode: (depth=" << depth_ << ",LB=" << bestLB_;
		out << ",Var=" << pNumberOfNurses_->name_ << ",LB=" << nursesLhs_ << ",UB=" << nursesRhs_ << ")";
		return out.str();
	}

	//number of nurse on which we have branched. pNumberOfNurses_ can be 0
	const MyVar* pNumberOfNurses_;
	const double nursesLhs_, nursesRhs_;
};

struct RestTree: public MyTree{
	RestTree(PScenario pScenario, PDemand pDemand);

	void logical_fixing();

	void addForbiddenShifts(PLiveNurse pNurse, std::set<std::pair<int,int> >& forbidenShifts);

	inline void pushBackNewNursesNumberNode(MyVar* var, double lb, double ub){
		NursesNumberNode* node = new NursesNumberNode(tree_.size(), currentNode_, var, lb, ub);
		pushBackNode(node);
	}

	inline void pushBackNewColumnsNode(std::vector<PPattern>& patterns){
		ColumnsNode* node = new ColumnsNode(tree_.size(), currentNode_, patterns);
		pushBackNode(node);
	}

	inline void pushBackNewRestNode(PLiveNurse pNurse, int day, bool rest, std::vector<MyVar*>& restArcs){
		RestNode* node = new RestNode(tree_.size(), currentNode_, pNurse, day, rest, restArcs);
		pushBackNode(node);
	}

	inline void pushBackNewShiftNode(PLiveNurse pNurse, int day, bool work, std::vector<int>& shifts){
		ShiftNode* node = new ShiftNode(tree_.size(), currentNode_, pNurse, day, work, shifts);
		pushBackNode(node);
	}

	inline bool is_columns_node() const { return dynamic_cast<ColumnsNode*>(currentNode_)!=0; }

	inline bool continueDiving() const {
		if(is_columns_node()) return true;
		return (diveDepth_ <= std::max(min_depth_+10, diveLength_));
	}

	void reset();

	void updateStats(MyNode* node);

    std::string writeBranchStats() const;

    std::string writeOneStat(std::string name, const std::vector<std::pair<int,double>>& stats) const;

protected:
	PScenario pScenario_;
	PDemand pDemand_;
	//Tree statistics
	//each pair represents first the number of times of the branching and second the increase of the LB
  std::vector<std::pair<int,double>> statsRestByDay_, statsWorkByDay_, statsRestByNurse_, statsWorkByNurse_;
    std::pair<int,double> statsCols_;
};

struct ColumnsComparator {
	virtual bool is_disjoint(PPattern col1, PPattern col2)=0;
};

struct DayDisjointComparator: public ColumnsComparator{
	bool is_disjoint(PPattern col1, PPattern col2) { return col1->isDisjointWith(col2); }
};

struct ShiftDisjointComparator: public ColumnsComparator{
	bool is_disjoint(PPattern col1, PPattern col2) { return col1->nurseId_!=col2->nurseId_ && col1->isShiftDisjointWith(col2); }
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

   /* compute fixing decisions */
   bool column_candidates(MyBranchingCandidate& candidate);

   /* Choose columns */
   std::vector<MyVar*> chooseColumns(std::vector<std::pair<MyVar*,double>>& candidates,
       std::vector<PPattern>& patterns, double& maxValue, ColumnsComparator& comparator);

   /* compare columns */
   static bool compareColumnCloseToInt(std::pair<MyVar*, double> obj1, std::pair<MyVar*, double> obj2);

   static bool compareColumnCloseTo5(std::pair<MyVar*, double> obj1, std::pair<MyVar*, double> obj2);

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
