/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#include "solvers/mp/modeler/CbcModeler.h"


//-----------------------------------------------------------------------------
//
// C l a s s  C b c M o d e l e r
//
// All the information relative to a particular demand
//
//-----------------------------------------------------------------------------


// Constructor
//
CbcModeler::CbcModeler(const vector<CoinVar *> &coreVars,
                       const vector<CoinVar *> &columnVars,
                       const vector<CoinCons *> &cons) :
    CoinModeler(), primalValues_(0), objVal_(0), model_(NULL), pOsiSolver_(0) {
  initializeVectors(coreVars, columnVars, cons);
}

CbcModeler::CbcModeler(const vector<CoinVar *> &coreVars,
                       const vector<CoinVar *> &columnVars,
                       const vector<CoinCons *> &cons,
                       OsiSolverInterface *osiSolver_) :
    CoinModeler(),
    primalValues_(0),
    objVal_(0),
    model_(NULL),
    pOsiSolver_(osiSolver_) {
  initializeVectors(coreVars, columnVars, cons);
}

void CbcModeler::initializeVectors(const vector<CoinVar *> &coreVars,
                                   const vector<CoinVar *> &columnVars,
                                   const vector<CoinCons *> &cons) {
  int corenum = coreVars.size();
  int colnum = columnVars.size();
  int consnum = cons.size();

  for (int i = 0; i < corenum; i++) {
    CoinVar *coreVar = new CoinVar(*coreVars.at(i));
    coreVars_.push_back(coreVar);
    switch (coreVar->getVarType()) {
      case VARTYPE_BINARY:binaryCoreVars_.push_back(coreVar);
        break;
      case VARTYPE_INTEGER:integerCoreVars_.push_back(coreVar);
        break;
      default:
        break;
    }
    addObject(coreVar);
  }
  for (int i = 0; i < colnum; i++)
    addObject(new CoinVar(*columnVars.at(i)));

  for (int i = 0; i < consnum; i++) {
    coreCons_.push_back(new CoinCons(*cons.at(i)));
    addObject(coreCons_[i]);
  }
}

/*
 * Create variable:
 *    var is a pointer to the pointer of the variable
 *    var_name is the name of the variable
 *    lhs, rhs are the lower and upper bound of the variable
 *    vartype is the type of the variable: VARTYPE_CONTINUOUS, VARTYPE_INTEGER, VARTYPE_BINARY
 */
int CbcModeler::createVar(MyVar **var,
                          const char *var_name,
                          int index,
                          double objCoeff,
                          double lb,
                          double ub,
                          VarType vartype,
                          const std::vector<double> &pattern,
                          double score) {
  *var = new CoinVar(var_name, index, objCoeff, vartype, lb, ub, pattern);
  return 1;
}

int CbcModeler::createColumnVar(MyVar **var,
                                const char *var_name,
                                int index,
                                double objCoeff,
                                const std::vector<double> &pattern,
                                double dualObj,
                                double lb,
                                double ub,
                                VarType vartype,
                                double score) {
  *var = new CoinVar(var_name, index, objCoeff,
                     vartype, lb, ub, pattern, dualObj);
  return 1;
}

/*
* Create linear constraint:
*    con is a pointer to the pointer of the constraint
*    con_name is the name of the constraint
*    lhs, rhs are the lower and upper bound of the constraint
*/
int CbcModeler::createCoinConsLinear(CoinCons **con,
                                     const char *con_name,
                                     int index,
                                     double lhs,
                                     double rhs) {
  *con = new CoinCons(con_name, index, lhs, rhs);
  objects_.push_back(*con);
  return 1;
}

/*
 * Create the Clp solver and assign the result the Cbc model
 * The method can only be called after the creation of all the variables and
 * linear constraints
 *
*/
int CbcModeler::setModel() {
  if (model_ != NULL) delete model_;

  if (pOsiSolver_) {
    model_ = new CbcModel(*pOsiSolver_);
    return 1;
  }

  vector<double> collb, colub, obj, rowlb, rowub;

  CoinPackedMatrix *ctMatrix = buildCoinMatrix();

  const int corenum = coreVars_.size();
  const int colnum = coreVars_.size() + initialColumnVars_.size();
  const int rownum = coreCons_.size();

  // get the characteristics of the variables
  for (int i = 0; i < colnum; ++i) {
    MyVar *var = nullptr;
    // copy of the core variables
    if (i < corenum)
      var = coreVars_[i];
    else  // copy of the column variables
      var = initialColumnVars_[i - corenum];

    // build the vectors defining the LP
    collb.push_back(var->getLB());
    colub.push_back(var->getUB());
    obj.push_back(var->getCost());
  }

  // get the characteristics of the constraints
  for (int i = 0; i < rownum; i++) {
    rowlb.push_back(coreCons_[i]->getLhs());
    rowub.push_back(coreCons_[i]->getRhs());
  }

  OsiSolverInterface *solver = new OsiClpSolverInterface;
  solver->loadProblem(*ctMatrix,
                      &(collb[0]),
                      &(colub[0]),
                      &(obj[0]),
                      &(rowlb[0]),
                      &(rowub[0]));
  delete ctMatrix;

  // set the types of the variables
  for (int i = 0; i < colnum; i++) {
    MyVar *var(nullptr);
    // copy of the core variables
    if (i < corenum)
      var = coreVars_[i];
      // copy of the column variables
    else
      var = initialColumnVars_[i - corenum];

    switch (var->getVarType()) {
      case VARTYPE_BINARY:
      case VARTYPE_INTEGER:solver->setInteger(i);
        break;
      case VARTYPE_CONTINUOUS:
      default:solver->setContinuous(i);
        break;
    }
  }

  model_ = new CbcModel(*solver);

  delete solver;

  return 1;
}

///*
// * Set a high priority on all the variable returned by the branching rule
// */
// void CbcModeler::setBranchingRule() {
//  vector<MyObject*> integerVariables;
//  // put all the integer variables
//  for (MyVar *coreVar : coreVars_)
//    if (coreVar->getVarType() != VARTYPE_CONTINUOUS)
//      integerVariables.push_back(coreVar);
//  for (MyVar *col : initialColumnVars_)
//    if (col->getVarType() != VARTYPE_CONTINUOUS)
//      integerVariables.push_back(col);
//
//  // remove the worst/best candidates, keep the medium ones
//  vector<MyObject*> branchingCandidates = integerVariables;
//  branching_candidates(branchingCandidates);
//
//  // remove the bad candidates, keep the best
//  vector<MyObject*> fixingCandidates = integerVariables;
//  logical_fixing(fixingCandidates);
//
//  // set the priorities: 1 highest and 100 lowest
//  int priorities[integerVariables.size()];  // NOLINT
//  int index = 0;
//  auto it = integerVariables.begin(),
//      itMedium = branchingCandidates.begin(),
//      itBest = fixingCandidates.begin();
//  while (it != integerVariables.end()) {
//    if (*it == *itBest) {
//      priorities[index] = 1;
//      ++itBest;
//    } else if (*it == *itMedium) {
//      priorities[index] = 50;
//      ++itMedium;
//    } else {
//      priorities[index] = 100;
//    }
//
//    ++index;
//    ++it;
//  }
//
//  // set to the highest priority (1) the var in branchingCandidates
//  model_->passInPriorities(priorities, false);
//}

/*
* get the primal values
*/
double CbcModeler::getVarValue(MyVar *var) {
  if (primalValues_ == nullptr) {
    Tools::throwError("Primal solution has not been initialized.");
  }
  return primalValues_[var->getIndex()];
}

// solve the model
int CbcModeler::solve(bool relaxation) {
  this->setModel();
//  this->setBranchingRule();
//  // set an upper bound
//  if (best_ub < DBL_MAX)
//    model_->setCutoff(best_ub);
  // set verbosity
  model_->setLogLevel(verbosity_);
//  // set solving time
//  if (max_solving_time < DBL_MAX)
//    model_->setMaximumSeconds(max_solving_time);
  model_->branchAndBound();
  this->setSolution();

  return model_->status();
}

// Set the value of the solution obtained after solving the MILP
void CbcModeler::setSolution() {
  objVal_ = model_->getCurrentObjValue();
  primalValues_ = model_->bestSolution();
}

/**************
 * Outputs *
 *************/

int CbcModeler::printStats() {
  std::cout << "Status of the solution = " << model_->status();
  if (model_->isProvenOptimal()) {
    std::cout << "The current solution is optimal." << std::endl;
  } else if (model_->isProvenInfeasible()) {
    std::cout << "The problem is infeasible." << std::endl;
  }
//  else if (model_->isSolutionLimitReached()) {
//
//  } else if (model_->isNodeLimitReached()) {
//
//  } else if (model_->isAbandoned()) {
//
//  }

  return model_->status();
}

/* Print the solution.  CbcModel clones the solver so we
   need to get current copy from the CbcModel */
int CbcModeler::printBestSol() {
  FILE *pFile;
  pFile = logfile_.empty() ? stdout : fopen(logfile().c_str(), "a");
  if (primalValues_ == 0)
    Tools::throwError("Primal solution has not been initialized.");

  int numberColumns = model_->solver()->getNumCols();

  // print the objective value
  fprintf(pFile, "%-30s %4.2f \n", "Objective:", objVal_);

  // print the value of the positive variables
  fprintf(pFile, "%-30s \n", "Variables:");
  double tolerance = pow(.1, DECIMALS);
  // iterate on core variables
  for (MyVar *var : coreVars_) {
    double value = getVarValue(var);
    if (fabs(value) > tolerance)
      fprintf(pFile,
              "%-30s %4.2f (%6.0f) \n",
              var->name_,
              value,
              var->getCost());
  }

  // iterate on column variables
  for (MyVar *var : initialColumnVars_) {
    double value = getVarValue(var);
    if (fabs(value) > tolerance)
      fprintf(pFile,
              "%-30s %4.2f (%6.0f) \n",
              var->name_,
              value,
              var->getCost());
  }

  fprintf(pFile, "\n");
  if (!logfile_.empty()) fclose(pFile);

  return 1;
}

// Write the MILP in the format implicity required by the filename
int CbcModeler::writeProblem(std::string filename) {
  if (model_ == NULL) return 1;

  // get the extension of the file
  std::string extension = filename.substr(filename.find_last_of(".") + 1);

  // use the relevant method depending on the extension
  if (!extension.compare("mps")) {
    model_->solver()->writeMps(filename.c_str());
  } else if (!extension.compare("lp")) {
    model_->solver()->writeLp(filename.c_str());
  } else {
    Tools::throwError("CbcModeler::writeLP: the extension of the file does "
                      "not match any available method. Use.mps or .lp.");
  }

  return 0;
}

int CbcModeler::writeLP(std::string filename) {
  model_->solver()->writeLp(filename.c_str());
  return 1;
}
