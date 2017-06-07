/*
 * TreeManager.cpp
 *
 *  Created on: Dec 30, 2015
 *      Author: tonio
 */

#include "TreeManager.h"

//////////////////////////////////////////////////////////////
//
// B C P T R E E S T A T S   M E T H O D S
//
//////////////////////////////////////////////////////////////

RestTree::RestTree(Scenario* pScenario, Demand* pDemand):
MyTree(), pScenario_(pScenario), pDemand_(pDemand),
statsRestByDay_(pDemand_->nbDays_), statsWorkByDay_(pDemand_->nbDays_),
statsRestByNurse_(pScenario->nbNurses()), statsWorkByNurse_(pScenario->nbNurses()) { }

// Forbid the shifts that are already fixed by the branching in the generation
// of new rotations
// This method is useful only for the rotation pricer
//
void RestTree::addForbiddenShifts(LiveNurse* pNurse, set<pair<int,int> >& forbidenShifts) {
	MyNode* node = currentNode_;
	vector<MyVar*> arcs;
	while(node->pParent_){
		// If on a rest node, we forbid shifts only if rest is forced
		// In that case, no rotation can inlude a shift on the fixed resting day
		//
		RestNode* restNode = dynamic_cast<RestNode*>(node);
		if(restNode != 0 && restNode->pNurse_ == pNurse)
		{
			if(restNode->rest_) {
				for(int i=1; i<pNurse->pScenario_->nbShifts_; ++i)	{
					forbidenShifts.insert(pair<int,int>(restNode->day_, i));
				}
			}
		}

		// If on a shift node, it is natural to forbid all the working forbidden
		// shifts in rotation pricing
		//
		ShiftNode* shiftNode = dynamic_cast<ShiftNode*>(node);
		if(shiftNode != 0 && shiftNode->pNurse_ == pNurse)
		{
			for(int s: shiftNode->forbiddenShifts_)	{
				if (s!=0) forbidenShifts.insert(pair<int,int>(shiftNode->day_, s));
			}
		}

		// If in a column node , forbid all the shifts that would be worked on a
		// day that is already covered by a fixed rotations
		// Moreover, there needs to be a resting day before and after each
		// rotation, so the shifts can also be forbidden on these two days (if
		// the rotation is not at an extremity of the horizon)
		ColumnsNode* columnsNode = dynamic_cast<ColumnsNode*>(node);
		if(columnsNode != 0 )
		{
			for(Rotation rot: columnsNode->rotations_) {
				if (rot.nurseId_ == pNurse->id_) {
					for (int day = rot.firstDay_-1; day <= rot.firstDay_+rot.length_; day++) {
						if (day<pScenario_->firstDay()) continue;
						if (day>=pScenario_->firstDay()+pScenario_->nbDays()) continue;
						for(int i=1; i<pNurse->pScenario_->nbShifts_; ++i)	{
							forbidenShifts.insert(pair<int,int>(day, i));
						}
					}
				}
			}
		}
		node = node->pParent_;

	}
}

void RestTree::updateStats(MyNode* node) {
	RestNode* restNode = dynamic_cast<RestNode*>(node);
	if(restNode != 0){
		double increaseLB = node->getBestLB() - node->pParent_->getBestLB();
		if(restNode->rest_){
			statsRestByDay_[restNode->day_].first ++;
			statsRestByDay_[restNode->day_].second += increaseLB;
			statsRestByNurse_[restNode->pNurse_->id_].first ++;
			statsRestByNurse_[restNode->pNurse_->id_].second += increaseLB;
		} else{
			statsWorkByDay_[restNode->day_].first ++;
			statsWorkByDay_[restNode->day_].second += increaseLB;
			statsWorkByNurse_[restNode->pNurse_->id_].first ++;
			statsWorkByNurse_[restNode->pNurse_->id_].second += increaseLB;
		}
	}

	ShiftNode* shiftNode = dynamic_cast<ShiftNode*>(node);
	if(shiftNode != 0){
		double increaseLB = node->getBestLB() - node->pParent_->getBestLB();
		if(!shiftNode->work_){
			statsRestByDay_[shiftNode->day_].first ++;
			statsRestByDay_[shiftNode->day_].second += increaseLB;
			statsRestByNurse_[shiftNode->pNurse_->id_].first ++;
			statsRestByNurse_[shiftNode->pNurse_->id_].second += increaseLB;
		} else{
			statsWorkByDay_[shiftNode->day_].first ++;
			statsWorkByDay_[shiftNode->day_].second += increaseLB;
			statsWorkByNurse_[shiftNode->pNurse_->id_].first ++;
			statsWorkByNurse_[shiftNode->pNurse_->id_].second += increaseLB;
		}
	}

	ColumnsNode* colsNode = dynamic_cast<ColumnsNode*>(node);
	if(colsNode!=0){
		statsCols_.first ++;
		statsCols_.second += node->getBestLB() - node->pParent_->getBestLB();
	}
}

string RestTree::writeBranchStats(){
	stringstream rep;
	rep << "";

	int nbNurses = statsRestByNurse_.size();
	int firstDay = pDemand_->firstDay_, nbDays = pDemand_->nbDays_;

	rep << "Stats on Columns" << std::endl;
	char buffer0[100];
	sprintf(buffer0, "Has branched %3d times with an average increased of the lower bound of %.2f", statsCols_.first, statsCols_.second / statsCols_.first);
	rep << buffer0 << endl;
	rep << "-------------------------------------"<< std::endl;


	rep << std::endl;
	rep << "Stats by Days" << std::endl;
	rep << "\t\t  ";
	for (int day = 0; day < nbDays; day++) rep << "|  " << Tools::intToDay(firstDay+day).at(0) << "  ";
	rep << "|" << std::endl;
	rep << writeOneStat("Rest", statsRestByDay_);
	rep << writeOneStat("Work", statsWorkByDay_);
	rep << "-------------------------------------"<< std::endl;


	rep << std::endl;
	rep << "Stats by Nurse" << std::endl;
	rep << "\t\t  ";
	for (int n = 0; n < nbNurses; n++){
		if(n<10) rep << "|" << pScenario_->theNurses_[n].name_ << " ";
		else rep << "|" << pScenario_->theNurses_[n].name_;
	}
	rep << "|" << std::endl;
	rep << writeOneStat("Rest", statsRestByNurse_);
	rep << writeOneStat("Work", statsWorkByNurse_);
	rep << "-------------------------------------"<< std::endl;

	return rep.str();
}

string RestTree::writeOneStat(string name, vector<pair<int,double>>& stats){
	stringstream rep;

	rep << name << "\t\t  ";
	for (int n = 0; n < stats.size(); n++) {
		pair<int,double>& p = stats[n];
		if(p.first){
			char buffer[100];
			sprintf(buffer, "|%3.2f", p.second / p.first);
			rep << buffer;
		}
		else rep << "| --- ";
	}
	rep << "|" << std::endl;

	return rep.str();
}

void RestTree::reset() {
	MyTree::reset();
	statsRestByDay_.clear();
	statsWorkByDay_.clear();
	statsRestByNurse_.clear();
	statsWorkByNurse_.clear();
	statsCols_.first = 0;
	statsCols_.second = 0;
}

//////////////////////////////////////////////////////////////
//
// B R A N C H I N G  R U L E
//
//////////////////////////////////////////////////////////////

template<typename T> bool compareObject(const pair<T,double>& p1, const pair<T,double>& p2){
	return (p1.second < p2.second);
}

/*************************************************************
 * Diving branching rule: dive then close to .5
 *************************************************************/

/* Constructs the branching rule object. */
DiveBranchingRule::DiveBranchingRule(MasterProblem* master, RestTree* tree, const char* name):
										MyBranchingRule(name), pMaster_(master), tree_(tree), pModel_(master->getModel())
{ }

//add all good candidates
bool DiveBranchingRule::column_candidates(MyBranchingCandidate& candidate){
	//look for fractional columns
	//Fix all column above BRANCH_LB
	//search the good candidates
	vector<pair<MyVar*,double> > candidates;
	vector<MyVar*> integerFixingCandidates;
	std::vector<MyVar*> fixingCandidates;
	std::vector<MyVar*> otherFixingCandidates;
	vector<Rotation> rotations;

	if ( (pModel_->getParameters().branchColumnUntilValue_&&pModel_->getParameters().branchColumnDisjoint_)
		|| (!pModel_->getParameters().branchColumnUntilValue_&&!pModel_->getParameters().branchColumnDisjoint_) )
	{
		Tools::throwError("DiveBranchingRule::column_candidates: the branching rule is not chosen properly!");
	}

	for(MyVar* var: pModel_->getActiveColumns()){

		// ROLLING: do not branch on relaxed columns
		//
		if ( pMaster_->isRelaxDay(var->getFirstDay()) || pMaster_->isFixDay(var->getFirstDay()) ) {
			continue;
		}

		// if we have already branched on this variable, continue
		if(var->getLB() > EPSILON) continue;

		double value = pModel_->getVarValue(var);
//		//if var fractional
//		if(value > EPSILON && value < 1 - EPSILON)
		//if var is non null
		if (value > 1-EPSILON) {
			integerFixingCandidates.push_back(var);
			Rotation rot(var->getPattern());
			rotations.push_back(rot);
			continue;
		}
		if(value > EPSILON)
			candidates.push_back(pair<MyVar*, double>(var, 1 - value));
	}

	stable_sort(candidates.begin(), candidates.end(), compareObject<MyVar*>);

	// Try the branching that mimics the heuristic
	//
	if (pModel_->getParameters().branchColumnUntilValue_) {
		if(candidates.size() > 0){
			double valueLeft = 1-EPSILON;
			for(pair<MyVar*,double>& p: candidates){
				if(p.second > valueLeft)
					break;
				else
					valueLeft -= p.second;
				fixingCandidates.push_back(p.first);
			}
		}
	}

	if (pModel_->getParameters().branchColumnDisjoint_) {
		double valueMax = pMaster_->pScenario_->nbWeeks_ / 2.0;
		DayDisjointComparator comp1 = DayDisjointComparator();
		fixingCandidates = chooseColumns(candidates, rotations, valueMax, comp1); // (pMaster_->pScenario_->nbWeeks_ + 1) / 2.0

		ShiftDisjointComparator comp2 = ShiftDisjointComparator();
		otherFixingCandidates = chooseColumns(candidates, rotations, valueMax, comp2); // (pMaster_->pScenario_->nbWeeks_ - 1)  / 2.0
	}


	if(fixingCandidates.empty()&&otherFixingCandidates.empty()&&integerFixingCandidates.empty()) return false;

	/* update candidate */
	MyBranchingNode& node = candidate.getChild(candidate.createNewChild());
	for(MyVar* var: fixingCandidates){
		int index = candidate.addBranchingVar(var);
		node.setLb(index, 1);
	}
	for(MyVar* var: otherFixingCandidates){
		int index = candidate.addBranchingVar(var);
		node.setLb(index, 1);
	}
	for(MyVar* var: integerFixingCandidates){
		int index = candidate.addBranchingVar(var);
		node.setLb(index, 1);
	}

	// Find the rotation to desactivate
	for(MyVar* var: pModel_->getActiveColumns()){
		if(var->getUB() == 0)
			continue;
		Rotation activeRot(var->getPattern());
		for(Rotation nodeRot: rotations) {

			if(activeRot==nodeRot || activeRot.nurseId_ != nodeRot.nurseId_ ||
				activeRot.firstDay_ > nodeRot.firstDay_+nodeRot.length_ ||
				activeRot.firstDay_+activeRot.length_ < nodeRot.firstDay_ ) continue;

			//add the variable to the candidate
			int index = candidate.addBranchingVar(var);

			//check if the shift is present in shifts
			//set the UB to 0 for the non-possible rotations
			node.setUb(index, 0);

			break;
		}
	}


	/* update tree */
	// append fixing candidates
	tree_->pushBackNewColumnsNode(rotations);

	return true;
}

bool DiveBranchingRule::branching_candidates(MyBranchingCandidate& candidate){
		return branchOnRestingArcs(candidate);
		//return branchOnShifts(candidate);
}


//-----------------------------------------------------------------------------
// Branch on a set of resting arcs
//-----------------------------------------------------------------------------

bool DiveBranchingRule::branchOnRestingArcs(MyBranchingCandidate& candidate){
	//set of rest variable closest to .5
	int bestDay = -1;
	double advantage = .2;
	LiveNurse* pBestNurse(0);
	double lowestScore = DBL_MAX;

	for(LiveNurse* pNurse: pMaster_->theLiveNurses_) {

		// ROLLING/LNS: do not branch on a nurse whose rotations are fixed
		//
		if (pMaster_->isFixNurse(pNurse->id_)) continue;

		for(int k=0; k<pMaster_->pDemand_->nbDays_; ++k){

			// ROLLING/LNS: do not branch on a day whose rotations are relaxed or fixed
			//
			if ( pMaster_->isRelaxDay(k) || pMaster_->isFixDay(k) ) {
				continue;
			}

			//choose the set of arcs the closest to .5
			double restValue = pModel_->getVarValue(pMaster_->getRestsPerDay(pNurse)[k]);

			//The value has to be not integer
			if(restValue < EPSILON || restValue > 1 - EPSILON)
				continue;

			double currentScore = 0.0;
			//look for the value the closest to .5
			double frac = restValue - floor(restValue);
			currentScore = abs(0.5-frac);
			if(Tools::isWeekend(k)) currentScore -= advantage;

			if(currentScore < lowestScore){
				bestDay = k;
				pBestNurse = pNurse;
				lowestScore = currentScore;
			}
		}
	}

	if(pBestNurse>0){
		//creating the branching cut
		char name[50];
		sprintf(name, "RestBranchingCons_N%d_%d", pBestNurse->id_, bestDay);
		vector<double> coeffs;
		vector<MyVar*> restingArcs;
		for(MyVar* var: pMaster_->getRestsPerDay(pBestNurse)[bestDay]){
			restingArcs.push_back(var);
			coeffs.push_back(1);
		}
		//create a new cons
		MyCons* cons;
		pModel_->createCutLinear(&cons, name, 0, 0, restingArcs, coeffs);

		/* update candidate */
		int index = candidate.addNewBranchingCons(cons);
		int index1 = candidate.createNewChild(), index2 = candidate.createNewChild();
		MyBranchingNode &restNode = candidate.getChild(index1), &workNode = candidate.getChild(index2);
		restNode.setLhs(index, 1);
		restNode.setRhs(index, 1);
		workNode.setLhs(index, 0);
		workNode.setRhs(index, 0);

		// Find the rotation to desactivate
		for(MyVar* var: pModel_->getActiveColumns()){
			if(var->getUB() == 0)
				continue;
			Rotation rot(var->getPattern());
			if(rot.nurseId_ != pBestNurse->id_ || rot.firstDay_ > bestDay || rot.firstDay_ + rot.length_ <= bestDay) continue;

			//add the variable to the candidate
			index = candidate.addBranchingVar(var);

			//check if the shift is present in shifts
			//set the UB to 0 for the non-possible rotations
			restNode.setUb(index, 0);
		}

		// Here : random choice to decide the order of the siblings
		if(Tools::randomInt(0, 1) == 0){
			tree_->pushBackNewRestNode(pBestNurse, bestDay, true, restingArcs);
			tree_->pushBackNewRestNode(pBestNurse, bestDay , false, restingArcs);
		} else {
			tree_->pushBackNewRestNode(pBestNurse, bestDay , false, restingArcs);
			tree_->pushBackNewRestNode(pBestNurse, bestDay , true, restingArcs);
			//put the workNode before the restNode
			candidate.swapChildren(index1, index2);
		}

		return true;
	}

	return branchOnShifts(candidate);
}


//-----------------------------------------------------------------------------
// Branch on a set of original penalty variables
//-----------------------------------------------------------------------------

bool DiveBranchingRule::branchOnPenalty(MyBranchingCandidate& candidate) {

	std::vector<std::pair<MyVar*,double> > candidates;
	std::vector<MyVar*> integerCandidates;

	// Check all the non-zero penalty variables and store them in a list of
	// candidates
	// The already integer ones are directly stored in a list of integer
	// candidates
	//
	for (MyVar* pVar: pMaster_->getMinWorkedDaysVars()) {
		double varValue = pModel_->getVarValue(pVar);
		if (varValue > EPSILON) {
			// if the variable is already integer and has not been branched on yet
			// record it in a specific list
			if (std::ceil(varValue)-varValue < EPSILON ) {
				if (varValue > pVar->getLB()+EPSILON) {
					integerCandidates.push_back(pVar);
				}
				continue;
			}
			candidates.push_back(std::pair<MyVar*,double>(pVar,varValue));
		}
	}
	for (MyVar* pVar: pMaster_->getMaxWorkedDaysVars()) {
		double varValue = pModel_->getVarValue(pVar);
		if (varValue > EPSILON) {
			// if the variable is already integer and has not been branched on yet
			// record it in a specific list
			if (std::ceil(varValue)-varValue < EPSILON ) {
				if (varValue > pVar->getLB()+EPSILON) {
					integerCandidates.push_back(pVar);
				}
				continue;
			}
			candidates.push_back(std::pair<MyVar*,double>(pVar,varValue));
		}
	}
	for (MyVar* pVar: pMaster_->getMaxWorkedWeekendVars()) {
		double varValue = pModel_->getVarValue(pVar);
		if (varValue > EPSILON) {
			candidates.push_back(std::pair<MyVar*,double>(pVar,varValue));
		}
	}
	for (int day= 0; day < pMaster_->getNbDays(); day++) {
		for (int shift = 0; shift < pMaster_->getNbShifts()-1; shift++) {
			for (MyVar* pVar: pMaster_->getOptDemandVars()[day][shift]) {
				double varValue = pModel_->getVarValue(pVar);
				if (varValue > EPSILON) {
					// if the variable is already integer and has not been branched on yet
					// record it in a specific list
					if (std::ceil(varValue)-varValue < EPSILON ) {
						if (varValue > pVar->getLB()+EPSILON) {
							integerCandidates.push_back(pVar);
						}
						continue;
					}
					candidates.push_back(std::pair<MyVar*,double>(pVar,varValue));
				}
			}
		}
	}

	// choose the candidate with value closest to integer and branch on the
	// farthest integer first to force a large difference in branch and bound
	if (!candidates.empty()) {

	}
	return true;
}


//-----------------------------------------------------------------------------
// Branch on a set of resting arcs
//-----------------------------------------------------------------------------

bool DiveBranchingRule::branchOnShifts(MyBranchingCandidate& candidate){
	//set of rest variable closest to .5
	int bestDay = -1;
	double advantage = .1;
	LiveNurse* pBestNurse(0);
	double lowestScore = DBL_MAX, value = 0;
	vector<int> forbiddenShifts;

	//compute the solution for each nurse, day, shift
	vector<MyVar*>& rotations = pModel_->getActiveColumns();
	vector<vector<vector<double>>> fractionalRoster = pMaster_->getFractionalRoster();

	//search for the best branching decision (set of shifts the closest to .5)
	for(LiveNurse* pNurse: pMaster_->theLiveNurses_) {

		// LNS: do not branch on a nurse whose rotations are fixed
		//
		if (pMaster_->isFixNurse(pNurse->id_)) continue;

		for(int k=0; k<pMaster_->pDemand_->nbDays_; ++k){

			// ROLLING/LNS:
			// do not branch on a day whose rotations are relaxed or fixed
			//
			if ( pMaster_->isRelaxDay(k) || pMaster_->isFixDay(k) ) {
				continue;
			}

			// Get the fraction that is spent on the resting shift and recover an
			// indexation of the working shifts from 1 to NbShifts_
			//
			vector<pair<int,double>> fractionalNurseDay;
			int s = 0;
			double valueLeft = 1;
			for(double val: fractionalRoster[pNurse->id_][k]){
				valueLeft -= val;
				fractionalNurseDay.push_back(pair<int,double>(++s, -val));
			}
			fractionalNurseDay.push_back(pair<int,double>(0, -valueLeft));

			sort(fractionalNurseDay.begin(), fractionalNurseDay.end(), compareObject<double>);

			//look for the balance between the shifts.
			double scoreNode1 = 0, scoreNode2 = 0;
			int nbShitfs1 = 0, nbShitfs2 = 0;
			vector<int> currentShifts;
			for(pair<int,double>& p: fractionalNurseDay){
				if(p.second == 0)
					if(nbShitfs1 >= nbShitfs2) ++ nbShitfs2;
					else{
						++ nbShitfs1;
						currentShifts.push_back(p.first);
					}
				else
					if(scoreNode1 > scoreNode2){
						++ nbShitfs2;
						scoreNode2 += -p.second;
					}
					else{
						++ nbShitfs1;
						scoreNode1 += -p.second;
						currentShifts.push_back(p.first);
					}
			}

			//look for the value the closest to .5
			double currentScore = abs(0.5-scoreNode1);
			if(Tools::isWeekend(k)) currentScore -= advantage;
			if(currentScore < lowestScore){
				bestDay = k;
				pBestNurse = pNurse;
				lowestScore = currentScore;
				forbiddenShifts = currentShifts;
			}
		}
	}

	if(pBestNurse>0){
		//creating the branching cut
		char name[50];
		sprintf(name, "ShiftBranchingCons_N%d_%d", pBestNurse->id_, bestDay);
		vector<double> coeffs;
		vector<MyVar*> restingArcs;
		for(MyVar* var: pMaster_->getRestsPerDay(pBestNurse)[bestDay]){
			restingArcs.push_back(var);
			coeffs.push_back(1);
		}
		//create a new cons
		MyCons* cons;
		pModel_->createCutLinear(&cons, name, 0, 1, restingArcs, coeffs);

		//compute the forbidden shifts
		vector<int> complementaryShifts;
		bool work = true;
		for(int i=0; i<pMaster_->pScenario_->nbShifts_; ++i){
			if(find(forbiddenShifts.begin(), forbiddenShifts.end(), i) == forbiddenShifts.end()){
				complementaryShifts.push_back(i);
				//if the forbidden shift 0 (rest) is in the complementary, it is possible to rest for the first node
				if(i==0) work = false;
			}
		}

		//add the node to the tree
		tree_->pushBackNewShiftNode(pBestNurse, bestDay, work, forbiddenShifts);
		tree_->pushBackNewShiftNode(pBestNurse, bestDay , !work, complementaryShifts);

		/* update candidate */
		int index = candidate.addNewBranchingCons(cons);
		int index1 = candidate.createNewChild(), index2 = candidate.createNewChild();
		MyBranchingNode &node1 = candidate.getChild(index1), &node2 = candidate.getChild(index2);
		// forbid to rest for the node where there are only working shifts
		if(work) node1.setRhs(index, 0);
		else node2.setRhs(index, 0);

		// Find the rotation to desactivate
		for(MyVar* var: pModel_->getActiveColumns()){
			if(var->getUB() == 0)
				continue;
			Rotation rot(var->getPattern());
			if(rot.nurseId_ != pBestNurse->id_ || rot.firstDay_ > bestDay || rot.firstDay_ + rot.length_ <= bestDay) continue;

			//add the variable to the candidate
			index = candidate.addBranchingVar(var);

			//check if the shift is present in shifts
			//set the UB to 0 for the non-possible rotations
			if( find(forbiddenShifts.begin(),forbiddenShifts.end(), rot.shifts_[bestDay]) != forbiddenShifts.end())
				node1.setUb(index, 0);
			else node2.setUb(index, 0);
		}
		return true;
	}

	return false;
}


//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

vector<MyVar*> DiveBranchingRule::chooseColumns(vector<pair<MyVar*,double>>& candidates, vector<Rotation>& rotations, double& maxValue, ColumnsComparator& comparator){
	vector<MyVar*> fixingCandidates;
	for(pair<MyVar*,double>& p: candidates){
		if(maxValue < p.second) continue;

		Rotation rot1(p.first->getPattern());
		//check if this rotation is totally disjoint with all the others
		//if not should be disjoint for the shift and the nurse
		bool isDisjoint = true;
		for(Rotation& rot2 : rotations)
			if( !comparator.is_disjoint(rot1, rot2) ){
				isDisjoint = false;
				break;
			}
		if(isDisjoint){
			fixingCandidates.push_back(p.first);
			maxValue -= p.second;
			rotations.push_back(rot1);
			//			cout << rot1.toString(7*pMaster_->pScenario_->nbWeeks_) << endl;
		}
	}

	return fixingCandidates;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool DiveBranchingRule::compareColumnCloseToInt(pair<MyVar*, double> obj1, pair<MyVar*, double> obj2){
	double frac1 = obj1.second - floor(obj1.second), frac2 = obj2.second - floor(obj2.second);
	double closeToInt1 = 1-frac1, closeToInt2 = 1-frac2;
	return (closeToInt1 < closeToInt2);
}


//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool DiveBranchingRule::compareColumnCloseTo5(pair<MyVar*, double> obj1, pair<MyVar*, double> obj2){
	double frac1 = obj1.second - floor(obj1.second), frac2 = obj2.second - floor(obj2.second);
	double closeTo5_1 = abs(0.5-frac1), closeTo5_2 = abs(0.5-frac2);
	return (closeTo5_1 < closeTo5_2);
}
