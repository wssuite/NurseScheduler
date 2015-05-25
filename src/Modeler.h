/*
 * ScipModeler.h:
 *    Tools to create a SCIP problem.
 *    In general, when you want to create a new scip object,
 *    you have to give a pointer to the pointer of the object.
 *
 *  Created on: 2015-02-23
 *      Author: legraina
 */

#ifndef SRC_MODELER_H_
#define SRC_MODELER_H_

/* standard library includes */
#include <cfloat>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>
#include "Solver.h"

#include "MyTools.h"

/* namespace usage */
using namespace std;

/*
 * Var types
 */
enum VarType {VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY};

/*
 * Rule Search Strategy
 */

enum SearchStrategy { BestFirstSearch, BreadthFirstSearch, DepthFirstSearch, HighestGapFirst };

/*
 * My Modeling objects
 * If the object is added to the vector objects_ of the Modeler, the modeler will also delete it at the end.
 */
struct MyObject {
   MyObject(const char* name):id_(s_count) {
      ++s_count;
      char* name2 = new char[255];
      strncpy(name2, name, 255);
      name_ = name2;
   }
   MyObject(const MyObject& myObject):id_(myObject.id_) {
      char* name2 = new char[255];
      strncpy(name2, myObject.name_, 255);
      name_ = name2;
   }
   virtual ~MyObject(){ delete[] name_; }
   //count object
   static unsigned int s_count;
   //for the map rotations_
   int operator < (const MyObject& m) const { return this->id_ < m.id_; }

   const char* name_;
private:
   const unsigned int id_;
};

struct MyVar: public MyObject{
   MyVar(const char* name, double cost, VarType type, double lb, double ub):
      MyObject(name), type_(type), cost_(cost), lb_(lb), ub_(ub)
   { }

   MyVar(const MyVar& var) :
      MyObject(var), type_(var.type_), cost_(var.cost_), lb_(var.lb_), ub_(var.ub_)
   { }

   virtual ~MyVar(){ }


   VarType getVarType() {return type_;}

   double getCost() { return cost_; }

   double getLB() { return lb_; }

   double getUB() { return ub_; }

   void setCost(double cost) { cost_=cost; }

   void setLB(double lb) { lb_=lb; }

   void setUB(double ub) { ub_=ub; }

   bool is_integer() { return type_ != VARTYPE_CONTINUOUS; }

protected:
   VarType type_; //type of the variable
   double cost_; //cost of the variable
   double lb_; //lower bound
   double ub_; //upper bound
};

struct MyCons: public MyObject{
   MyCons(const char* name, double lhs, double rhs):
      MyObject(name), lhs_(lhs), rhs_(rhs)
{ }

   MyCons(const MyCons& cons) :
      MyObject(cons.name_), lhs_(cons.lhs_), rhs_(cons.rhs_)
   { }

   virtual ~MyCons(){ }

   /*
    * Getters
    */

   double getLhs() { return lhs_; }

   double getRhs() { return rhs_; }

   void setLhs(double lhs) { lhs_=lhs; }

   void setRhs(double rhs) { rhs_=rhs; }

protected:
   double lhs_; //left hand side == lower bound
   double rhs_; //rihgt hand side == upper bound
};

/*
 * My pricer
 */
struct MyPricer{
   MyPricer(const char* name): name_(name){ }
   virtual ~MyPricer() { }

   //name of the pricer handler
   //
   const char* name_;

   /* perform pricing */
   //return true if optimal
   virtual bool pricing(double bound=0, bool before_fathom = true)=0;
};
/*
 * My branching rule
 */
struct MyBranchingRule{
   MyBranchingRule(const char* name): name_(name), searchStrategy_(BestFirstSearch) { }
   virtual ~MyBranchingRule() { }

   //name of the branching rule handler
   //
   const char* name_;

   /* compute logical fixing decisions */
   virtual void logical_fixing(vector<MyVar*>& fixingCandidates)=0; // stores the next candidates for logical fixing

   /* compute branching decisions */
   virtual void branching_candidates(vector<MyVar*>& branchingCandidates)=0; // stores the next candidates for branching

   void set_search_strategy(SearchStrategy searchStrategy){ searchStrategy_ = searchStrategy; }

protected:
   SearchStrategy searchStrategy_;
};

/* Exception to stop the solver */
struct FeasibleStop: public exception{
   FeasibleStop(string str){ cout << str << endl; }
};

struct InfeasibleStop: public exception{
   InfeasibleStop(string str){ cout << str << endl; }
};

struct OptimalStop: public exception{
   OptimalStop(string str){ cout << str << endl; }
};

class Modeler {
public:

   Modeler(): pPricer_(0), pBranchingRule_(0), best_ub(LARGE_SCORE) { }

   virtual ~Modeler(){
      for(MyObject* object: objects_)
         delete object;
   }

   //solve the model
   virtual int solve(bool relaxation = false)=0;

   //Reset and clear solving parameters
   virtual void reset() { best_ub = LARGE_SCORE; }

   //Add a pricer
   virtual int addObjPricer(MyPricer* pPricer){
      pPricer_ = pPricer;
      return 1;
   }

   //Add a branching rule
   virtual int addBranchingRule(MyBranchingRule* pBranchingRule){
      pBranchingRule_ = pBranchingRule;
      pBranchingRule_->set_search_strategy(searchStrategy_);
      return 1;
   }

   virtual void addForbidenShifts(LiveNurse* pNurse, set<pair<int,int> >& forbidenShifts) { }


   /*
    * Class methods for pricer and branching rule
    */

   //return true if optimal
   inline bool pricing(double bound=0, bool before_fathom = true){
      if(pPricer_)
         return pPricer_->pricing(bound, before_fathom);
      return true;
   }

   inline void branching_candidates(vector<MyVar*>& branchingCandidates){
      if(pBranchingRule_)
         pBranchingRule_->branching_candidates(branchingCandidates);
   }

   //remove all bad candidates from fixingCandidates
   inline void logical_fixing(vector<MyVar*>& fixingCandidates){
      if(pBranchingRule_)
         pBranchingRule_->logical_fixing(fixingCandidates);
   }
   //Set search strategy
   inline void set_search_strategy(SearchStrategy searchStrategy){
      if(pBranchingRule_)
         pBranchingRule_->set_search_strategy(searchStrategy);
   }

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
    */

   virtual int createVar(MyVar** var, const char* var_name, double objCoeff,
      double lb, double ub, VarType vartype, double score)=0;

   inline void createPositiveVar(MyVar** var, const char* var_name, double objCoeff, double score = 0, double ub = DBL_MAX){
      createVar(var, var_name, objCoeff, 0.0, ub, VARTYPE_CONTINUOUS, score);
   }

   inline void createIntVar(MyVar** var, const char* var_name, double objCoeff, double score = 0, double ub = DBL_MAX){
      createVar(var, var_name, objCoeff, 0, ub, VARTYPE_INTEGER, score);
      integerCoreVars_.push_back(*var);
   }

   inline void createBinaryVar(MyVar** var, const char* var_name, double objCoeff, double score = 0){
      createVar(var, var_name, objCoeff, 0.0, 1.0, VARTYPE_BINARY, score);
      binaryCoreVars_.push_back(*var);
   }

   virtual int createColumnVar(MyVar** var, const char* var_name, double objCoeff, double dualObj,
      double lb, double ub, VarType vartype, double score)=0;

   inline void createPositiveColumnVar(MyVar** var, const char* var_name, double objCoeff, double dualObj = 99999, double score = 0, double ub = DBL_MAX){
      createColumnVar(var, var_name, objCoeff, dualObj, 0.0, ub, VARTYPE_CONTINUOUS, score);
   }

   inline void createIntColumnVar(MyVar** var, const char* var_name, double objCoeff, double dualObj = 99999, double score = 0, double ub = DBL_MAX){
      createColumnVar(var, var_name, objCoeff, dualObj, 0, ub, VARTYPE_INTEGER, score);
   }

   inline void createBinaryColumnVar(MyVar** var, const char* var_name, double objCoeff, double dualObj = 99999, double score = 0){
      createColumnVar(var, var_name, objCoeff, dualObj, 0.0, 1.0, VARTYPE_BINARY, score);
   }

   /*
    * Create linear constraint:
    *    con is a pointer to the pointer of the constraint
    *    con_name is the name of the constraint
    *    lhs, rhs are the lower and upper bound of the constraint
    *    nonZeroVars is the number of non-zero coefficients to add to the constraint
    *    vars is an array of pointers to the variables to add to the constraints (with non-zero coefficient)
    *    coeffs is the array of coefficient to add to the constraints
    */

   virtual int createConsLinear(MyCons** cons, const char* con_name, double lhs, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {})=0;

   //Add a lower or equal constraint
   inline void createLEConsLinear(MyCons** cons, const char* con_name, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(cons, con_name, DBL_MIN, rhs, vars, coeffs);
   }

   //Add a greater or equal constraint
   inline void createGEConsLinear(MyCons** cons, const char* con_name, double lhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(cons, con_name, lhs, DBL_MAX, vars, coeffs);
   }

   //Add an equality constraint
   inline void createEQConsLinear(MyCons** cons, const char* con_name, double eq,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createConsLinear(cons, con_name, eq, eq, vars, coeffs);
   }

   //Add final linear constraints
   virtual int createFinalConsLinear(MyCons** cons, const char* con_name, double lhs, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {})=0;

   inline void createFinalLEConsLinear(MyCons** cons, const char* con_name, double rhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(cons, con_name, DBL_MIN, rhs, vars, coeffs);
   }

   inline void createFinalGEConsLinear(MyCons** cons, const char* con_name, double lhs,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(cons, con_name, lhs, DBL_MAX, vars, coeffs);
   }

   inline void createFinalEQConsLinear(MyCons** cons, const char* con_name, double eq,
      vector<MyVar*> vars = {}, vector<double> coeffs = {}){
      createFinalConsLinear(cons, con_name, eq, eq, vars, coeffs);
   }

   /*
    * Add variables to constraints
    */

   virtual int addCoefLinear(MyCons* cons, MyVar* var, double coeff, bool transformed=false)=0;

   /*
    * Add new Column to the problem
    */

   inline void createColumn(MyVar** var, const char* var_name, double objCoeff, double dualObj,  VarType vartype,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      switch(vartype){
      case VARTYPE_BINARY:
         createBinaryColumnVar(var, var_name, objCoeff, dualObj, score);
         break;
      case VARTYPE_INTEGER:
         createIntColumnVar(var, var_name, objCoeff, dualObj, score);
         break;
      default:
         createPositiveColumnVar(var, var_name, objCoeff, dualObj, score);
         break;
      }

      for(int i=0; i<cons.size(); i++)
         addCoefLinear(cons[i], *var, coeffs[i], transformed);
   }

   inline void createPositiveColumn(MyVar** var, const char* var_name, double objCoeff, double dualObj,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, dualObj, VARTYPE_CONTINUOUS, cons, coeffs, transformed, score);
   }

   inline void createBinaryColumn(MyVar** var, const char* var_name, double objCoeff, double dualObj,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, dualObj, VARTYPE_BINARY, cons, coeffs, transformed, score);
   }

   inline void createIntColumn(MyVar** var, const char* var_name, double objCoeff, double dualObj,
      vector<MyCons*> cons = {}, vector<double> coeffs = {}, bool transformed = false, double score = 0){
      createColumn(var, var_name, objCoeff, dualObj, VARTYPE_INTEGER, cons, coeffs, transformed, score);
   }

   /*
    * get the primal values
    */

   virtual bool isInteger(MyVar* var){
      double value = getVarValue(var);
      double fractionalPart = round(value) - value;
      if( abs(fractionalPart) < EPSILON )
         return true;
      return false;
   }

   virtual double getVarValue(MyVar* var)=0;

   inline vector<double> getVarValues(vector<MyVar*> vars){
      vector<double> values(vars.size());
      for(int i=0; i<vars.size(); ++i)
         values[i] = getVarValue(vars[i]);
      return values;
   }

   /*
    * Get the dual variables
    */

   virtual double getDual(MyCons* cons, bool transformed = false)=0;

   inline vector<double> getDuals(vector<MyCons*> cons, bool transformed = false){
      vector<double> dualValues(cons.size());
      for(int i=0; i<cons.size(); ++i)
         dualValues[i] = getDual(cons[i], transformed);
      return dualValues;
   }

   /*
    * Getters and setters
    */

   //compute the total cost of MyObject* in the solution
   virtual double getTotalCost(MyVar* var, bool print = false)=0;

   //compute the total cost of a vector of MyObject* in the solution
   template<typename T>  inline double getTotalCost(map<MyVar*, T> map0, bool print = false){
      double value = 0 ;
      for(pair<MyObject*, T> var: map0)
         value += getTotalCost(var.first, print);
      return value;
   }

   //compute the total cost of a multiple vectors of MyObject* in the solution
   template<typename V> inline double getTotalCost(vector<V> vector, bool print = false){
      double value = 0 ;
      for(V vect: vector)
         value += getTotalCost(vect, print);
      return value;
   }

   virtual double getObjective()=0;

   virtual double getRelaxedObjective() { return DBL_MAX; }

   /**************
    * Parameters *
    *************/
   virtual int setVerbosity(int v)=0;

   /**************
    * Outputs *
    *************/

   virtual int printStats()=0;

   virtual int printBestSol()=0;

   virtual int writeProblem(string fileName)=0;

   virtual int writeLP(string fileName)=0;

   virtual void toString(MyObject* obj){ cout << obj->name_ << endl; }

   /**************
    * Getters *
    *************/

   template<typename M> M getModel(){
      string error = "This template has not been implemented.";
      Tools::throwError(error.c_str());
   }

   inline vector<MyVar*>& getBinaryCoreVars(){ return binaryCoreVars_; }

   inline vector<MyVar*>& getIntegerCoreVars(){ return integerCoreVars_; }

   inline int getVerbosity() { return verbosity_; }

   inline double getBestUB() { return best_ub; }

   inline virtual void setBestUB(double ub) { if(ub < best_ub) best_ub = ub; }

   inline void setSearchStrategy(SearchStrategy searchStrategy){
      searchStrategy_ = searchStrategy;
      set_search_strategy(searchStrategy);
   }

   inline SearchStrategy getSearchStrategy(){ return searchStrategy_; }

   inline void setLastBranchingRest(pair<LiveNurse*, int> lastBranchingRest){ lastBranchingRest_ = lastBranchingRest; }

   inline pair<LiveNurse*, int> getLastBranchingRest() { return lastBranchingRest_; }

   inline void setParameters(SolverParam parameters){ 
    parameters_ = parameters;
    logfile_ = parameters.logfile_;
   }
   inline string logfile() {return logfile_;}

   inline SolverParam& getParameters() { return parameters_; }

   inline void setLogFile(string fileName) {logfile_ = fileName;}

protected:
   //store all MyObject*
   vector<MyObject*> objects_;
   vector<MyVar*> binaryCoreVars_;
   vector<MyVar*> integerCoreVars_;

   MyPricer* pPricer_;
   MyBranchingRule* pBranchingRule_;

   int verbosity_ = 0;
   SearchStrategy searchStrategy_ = BestFirstSearch;
   //best current upper bound found
   double best_ub;

   SolverParam parameters_;

   //strore the last branching decisions
   pair<LiveNurse*, int> lastBranchingRest_;

   // log file where outputs must be written
   string logfile_="";
};


#endif /* SRC_MODELER_H_ */
