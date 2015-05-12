/*
 * CoinModeler.h
 *
 *  Created on: 2015-04-01
 *      Author: legraina
 */

#ifndef SRC_COINMODELER_H_
#define SRC_COINMODELER_H_

#include "Modeler.h"

/* Coin includes */
#include <CoinPackedMatrix.hpp>

/*
 * My Variables
 */
//Coin var, just a virtual class
struct CoinVar: public MyObject {
   CoinVar(const char* name, int index, double cost, VarType type, double lb, double ub, double dualCost = 99999):
      MyObject(name), index_(index), type_(type), cost_(cost), lb_(lb), ub_(ub), dualCost_(dualCost)
   { }

   CoinVar(const CoinVar& var) :
      MyObject(var), index_(var.index_), type_(var.type_), cost_(var.cost_), dualCost_(var.dualCost_),
      lb_(var.lb_), ub_(var.ub_), indexRows_(var.indexRows_), coeffs_(var.coeffs_)
   { }

   virtual ~CoinVar(){ }

   /*
    * Setters and Getters
    */

   void addRow(int index, double coeff){
      indexRows_.push_back(index);
      coeffs_.push_back(coeff);
   }

   int getIndex() { return index_; }

   VarType getVarType() {return type_;}

   int getNbRows() { return indexRows_.size(); }

   vector<int>& getIndexRows() { return indexRows_; }

   int getIndexRow(int i) { return indexRows_[i]; }

   vector<double>& getCoeffRows() { return coeffs_; }

   double getCoeffRow(int i) { return coeffs_[i]; }

   double getCost() { return cost_; }

   double getLB() { return lb_; }

   double getUB() { return ub_; }

   bool is_integer() { return type_ != VARTYPE_CONTINUOUS; }

protected:
   int index_; //index of the column of the matrix here
   VarType type_; //type of the variable
   double cost_; //cost of the variable
   double dualCost_; //dualCost of the variable
   double lb_; //lower bound
   double ub_; //upper bound
   vector<int> indexRows_; //index of the rows of the matrix where the variable has non-zero coefficient
   vector<double> coeffs_; //value of these coefficients
};

/*
 * My Constraints
 */
//Coin cons, just a virtual class
struct CoinCons: public MyObject{
public:
   CoinCons(const char* name, int index, double lhs, double rhs):
      MyObject(name), index_(index), lhs_(lhs), rhs_(rhs)
{ }

   CoinCons(const CoinCons& cons) :
      MyObject(cons.name_), index_(cons.index_), lhs_(cons.lhs_), rhs_(cons.rhs_)
   { }

   virtual ~CoinCons(){ }

   /*
    * Getters
    */

   int getIndex() { return index_; }

   double getLhs() { return lhs_; }

   double getRhs() { return rhs_; }

protected:
   int index_; //index of the row of the matrix here
   double lhs_; //left hand side == lower bound
   double rhs_; //rihgt hand side == upper bound
};

class CoinModeler: public Modeler {
public:
   CoinModeler():
      Modeler()
{ }
   virtual ~CoinModeler() {}

   //solve the model
   virtual int solve(bool relaxation = false)=0;

   /*
    * Create variable:
    *    var is a pointer to the pointer of the variable
    *    var_name is the name of the variable
    *    lhs, rhs are the lower and upper bound of the variable
    *    vartype is the type of the variable: VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
    */

   //has to be implement to create the good var according to the coin modeler chosen (BCP, CBC ...)
   //WARNING: core vars have to all be created before creating column vars
   virtual int createCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub)=0;

   virtual int createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, double dualObj, VarType vartype, double lb, double ub)=0;

   int createVar(MyObject** var, const char* var_name, double objCoeff, double lb, double ub, VarType vartype, double score){
      if(lb==DBL_MIN)
         lb = -infinity;
      if(ub==DBL_MAX)
         ub = infinity;

      int index = coreVars_.size() + columnVars_.size();
      CoinVar* var2;
      createCoinVar(&var2, var_name, index, objCoeff, vartype, lb, ub);

      coreVars_.push_back(var2);
      *var = var2;

      return 1;
   }

   int createColumnVar(MyObject** var, const char* var_name, double objCoeff, double dualObj, double lb, double ub, VarType vartype, double score){
      if(lb==DBL_MIN)
         lb = -infinity;
      if(ub==DBL_MAX)
         ub = infinity;

      int index = coreVars_.size() + columnVars_.size();
      CoinVar* var2;
      createColumnCoinVar(&var2, var_name, index, objCoeff, dualObj, vartype, lb, ub);

      columnVars_.push_back(var2);
      *var = var2;

      return 1;
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
   //has to be implement to create the good cons according to the coin modeler chosen (BCP, CBC ...)
   virtual int createCoinConsLinear(CoinCons** con, const char* con_name, int index, double lhs, double rhs)=0;

   int createConsLinear(MyObject** con, const char* con_name, double lhs, double rhs,
      vector<MyObject*> vars, vector<double> coeffs){
      if(lhs==DBL_MIN)
         lhs = -infinity;
      if(rhs==DBL_MAX)
         rhs = infinity;

      int index = cons_.size();
      CoinCons* con2;
      createCoinConsLinear(&con2, con_name, index, lhs, rhs);

      for(int i=0; i<vars.size(); ++i)
         addCoefLinear(con2, vars[i], coeffs[i]);

      cons_.push_back(con2);
      *con = con2;

      return 1;
   }

   //Add there is no final linear constraints for BCP
   virtual int createFinalConsLinear(MyObject** con, const char* con_name, double lhs, double rhs,
      vector<MyObject*> vars = {}, vector<double> coeffs = {}){
      return createConsLinear(con, con_name, lhs, rhs, vars, coeffs);
   }

   /*
    * Add variables to constraints
    */

   int addCoefLinear(MyObject* cons, MyObject* var, double coeff, bool transformed=false){
      CoinVar* var2 = (CoinVar*) var;
      CoinCons* cons2 = (CoinCons*) cons;

      var2->addRow(cons2->getIndex(), coeff);

      return 1;
   }

   /*
    * Build the CoinedPackMatrix corresponding to the problem
    * If justCore = true, add just the core variables
    */

   /* build the problem */

   CoinPackedMatrix buildCoinMatrix(bool justCore = false){
      //define nb rows and col
      const int corenum = coreVars_.size();
      const int colnum = (justCore) ? corenum : corenum + columnVars_.size();

      //define a matrix as a vector of tuples (row_index, col_index, coeff)
      vector<int> row_indices, col_indices;
      vector<double> coeffs;

      for(int i=0; i<colnum; ++i){
         CoinVar* var(0);
         //copy of the core variables
         if(i<corenum)
            var = coreVars_[i];
         //copy of the column variables
         else
            var = columnVars_[i-corenum];

         //build the tuples of the matrix
         for(int j=0; j<var->getNbRows(); ++j){
            row_indices.push_back(var->getIndexRow(j));
            col_indices.push_back(var->getIndex());
            coeffs.push_back(var->getCoeffRow(j));
         }
      }

      //initialize the matrix
      CoinPackedMatrix matrix(false, &(row_indices[0]), &(col_indices[0]), &(coeffs[0]), row_indices.size());

      return matrix;
   }

   /*
    * Get the primal value
    */

   virtual double getVarValue(MyObject* var) { return 0; }

   /*
    * Get the dual variables
    */

   virtual double getDual(MyObject* cons, bool transformed = false) { return 0; }

   /**************
    * Parameters *
    *************/
   virtual int setVerbosity(int v) { return 0; }

   /**************
    * Outputs *
    *************/

   //compute the total cost of a var*
   double getTotalCost(MyObject* var){
      CoinVar* var2 = (CoinVar*) var;

      double value = getVarValue(var);
      return value *  var2->getCost();
   }

   virtual int printStats() { return 0; }

   virtual int printBestSol(){
      //print the value of the relaxation
      printf("%-30s %4.2f \n", "Relaxation:" , getRelaxedObjective());

      if(!loadBestSol())
         return 0;

      //print the objective value
      printf("%-30s %4.2f \n", "Objective:" , getObjective());

      if(verbosity_>=2){
         //print the value of the positive variables
         printf("%-30s \n", "Variables:");
         double tolerance = pow(.1, DECIMALS);
         //iterate on core variables
         for(CoinVar* var: coreVars_){
            double value = getVarValue(var);
            if( value > tolerance)
               printf("%-30s %4.2f (%6.0f) \n", var->name_ , value, var->getCost());
         }
         //iterate on column variables
         for(CoinVar* var: columnVars_){
            double value = getVarValue(var);
            if( value > tolerance)
               printf("%-30s %4.2f (%6.0f) \n", var->name_ , value, var->getCost());
         }

         printf("\n");
      }

      return 1;
   }

   virtual bool loadBestSol() { return false; }

   virtual double getObjective()=0;

   virtual double getRelaxedObjective() { return 0; }

   virtual int writeProblem(string fileName) { return 0; }

   virtual int writeLP(string fileName) { return 0; }


   //get the variables that are always present in the model
   vector<CoinVar*>& getCoreVars(){ return coreVars_; }

   //get the variables that are generating during the resolution (columns)
   vector<CoinVar*>& getColumns(){ return columnVars_; }

   int getNbColumns(){ return columnVars_.size(); }

   int getNbVars(){ return columnVars_.size()+coreVars_.size(); }

   vector<CoinCons*>& getCons(){ return cons_; }


protected:

   //Coin data
   double infinity = 1e40;

   vector<CoinVar*> coreVars_;
   vector<CoinVar*> columnVars_;
   vector<CoinCons*> cons_;
};

#endif /* SRC_COINMODELER_H_ */
