/*
* CoinModeler.h
*
*  Created on: 2015-04-01
*      Author: legraina
*/

#ifndef SRC_COINMODELER_H_
#define SRC_COINMODELER_H_

#include "solvers/mp/modeler/Modeler.h"

/* Coin includes */
#include <CoinPackedMatrix.hpp>

/*
* My Constraints
*/
//Coin cons, just a virtual class
struct CoinCons: public MyCons{
public:
	CoinCons(const char* name, int index, double lhs, double rhs):
	MyCons(name, index, lhs, rhs)
	{ }

	CoinCons(const CoinCons& cons) :
	MyCons(cons)
	{ }

	virtual ~CoinCons(){ }
};

/*
* My Variables
*/
//Coin var, just a virtual class
struct CoinVar: public MyVar {
	CoinVar(const char* name, int index, double cost, VarType type, double lb, double ub,
	const std::vector<double>& pattern = DEFAULT_PATTERN, double dualCost = 99999,
	const std::vector<int>& indexRows = Tools::EMPTY_INT_VECTOR, const std::vector<double>& coeffs = Tools::EMPTY_DOUBLE_VECTOR):
	MyVar(name, index, cost, type, lb, ub, pattern), dualCost_(dualCost), indexRows_(indexRows), coeffs_(coeffs)
	{ }

	CoinVar(const CoinVar& var) :
	MyVar(var), dualCost_(var.dualCost_), indexRows_(var.indexRows_), coeffs_(var.coeffs_)
	{ }

	virtual ~CoinVar(){ }

	/*
	* Setters and Getters
	*/

	void addRow(int index, double coeff){
		indexRows_.push_back(index);
		coeffs_.push_back(coeff);
	}

	void toString(std::vector<CoinCons*>& cons) {
    std::cout << name_ << ":";
		for(unsigned int i=0; i<indexRows_.size(); ++i)
      std::cout << " " << cons[indexRows_[i]]->name_ << ":" << coeffs_[i];
    std::cout << std::endl;
	}

	int getNbRows() { return indexRows_.size(); }

    std::vector<int>& getIndexRows() { return indexRows_; }

	int getIndexRow(int i) { return indexRows_[i]; }

    std::vector<double>& getCoeffRows() { return coeffs_; }

	double getCoeffRow(int i) { return coeffs_[i]; }

protected:
	double dualCost_; //dualCost of the variable
    std::vector<int> indexRows_; //index of the rows of the matrix where the variable has non-zero coefficient
    std::vector<double> coeffs_; //value of these coefficients
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

protected:
	//has to be implement to create the good var according to the coin modeler chosen (BCP, CBC ...)
	//WARNING: core vars have to all be created before creating column vars
	virtual int createCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, VarType vartype, double lb, double ub, const std::vector<double>& pattern = DEFAULT_PATTERN)=0;

	virtual int createColumnCoinVar(CoinVar** var, const char* var_name, int index, double objCoeff, const std::vector<double>& pattern, double dualObj, VarType vartype, double lb, double ub)=0;

	int createVar(MyVar** var, const char* var_name, int index, double objCoeff, double lb, double ub, VarType vartype, const std::vector<double>& pattern, double score){

		CoinVar* var2;
		createCoinVar(&var2, var_name, index, objCoeff, vartype, lb, ub, pattern);
		coreVars_.push_back(var2);
		*var = var2;

		return 1;
	}

	int createColumnVar(MyVar** var, const char* var_name, int index, double objCoeff, const std::vector<double>& pattern, double dualObj, double lb, double ub, VarType vartype, double score){

		CoinVar* var2;
		createColumnCoinVar(&var2, var_name, index, objCoeff, pattern, dualObj, vartype, lb, ub);
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

	int createConsLinear(MyCons** con, const char* con_name, int index, double lhs, double rhs,
                       std::vector<MyVar*> vars, std::vector<double> coeffs){

		CoinCons* con2;
		createCoinConsLinear(&con2, con_name, index, lhs, rhs);

		for(unsigned int i=0; i<vars.size(); ++i)
		addCoefLinear(con2, vars[i], coeffs[i]);

		cons_.push_back(con2);
		*con = con2;

		return 1;
	}

	//Add there is no final linear constraints for BCP
	virtual int createFinalConsLinear(MyCons** con, const char* con_name, int index, double lhs, double rhs,
                                    std::vector<MyVar*> vars = {}, std::vector<double> coeffs = {}){
		return createConsLinear(con, con_name, index, lhs, rhs, vars, coeffs);
	}

public:
	/*
	* Add variables to constraints
	*/

	int addCoefLinear(MyCons* cons, MyVar* var, double coeff, bool transformed=false){
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
	CoinPackedMatrix buildCoinMatrix(const std::vector<MyVar*>& columns = EMPTY_VARS){
		//define nb rows and col
		const int corenum = coreVars_.size(), colnum = corenum+columns.size();
		//define a matrix as a vector of tuples (row_index, col_index, coeff)
    std::vector<int> row_indices, col_indices;
    std::vector<double> coeffs;

		for(int i=0; i<colnum; ++i){
			CoinVar* var(0);
			//copy of the core variables
			if(i<corenum)
			var = (CoinVar*) coreVars_[i];
			//copy of the column variables
			else
			var = (CoinVar*) columns[i-corenum];

			//build the tuples of the matrix
			int col_index = var->getIndex();
			for(int j=0; j<var->getNbRows(); ++j){
				row_indices.push_back(var->getIndexRow(j));
				col_indices.push_back(col_index);
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

	virtual double getVarValue(MyVar* var) { return 0; }

	/*
	* Get the dual variables
	*/

	virtual double getDual(MyCons* cons, bool transformed = false) { return 0; }

	/*
	* Get the reduced cost
	*/

	virtual double getReducedCost(MyVar* var) { return 0; }

	/**************
	* Parameters *
	*************/
	virtual int setVerbosity(int v) { return 0; }

	/**************
	* Outputs *
	*************/

	//compute the total cost of a var*
	double getTotalCost(MyVar* var, bool print = false){
		CoinVar* var2 = (CoinVar*) var;

		double value = getVarValue(var);
		if(print && value>EPSILON)
      std::cout << var->name_ << ": " << value << "*" << var2->getCost() << std::endl;
		return value *  var2->getCost();
	}

	virtual int nbSolutions() { return 0; }

	inline virtual double getObjective(){ return pTree_->getBestUB(); }
	virtual double getObjective(int index) { return LARGE_SCORE; }

	virtual int printBestSol(){
		FILE * pFile;
		pFile = logfile_.empty() ? stdout : fopen (logfile_.c_str(),"a");
		//print the value of the relaxation
		fprintf(pFile,"%-30s %4.2f \n", "Relaxation:" , getRelaxedObjective());

		if(!loadBestSol())
		return 0;

		//print the objective value
		fprintf(pFile,"%-30s %4.2f \n", "Objective:" , Modeler::getObjective());

		//display all objective solution found
		fprintf(pFile,"%-30s \n", "All Objective:");
		for(int n=0; n<nbSolutions(); ++n)
		fprintf(pFile,"Solution number: %-5d: %4.2f \n", n , getObjective(n));

		if(verbosity_>=2){
			//print the value of the positive variables
			fprintf(pFile,"%-30s \n", "Variables:");
			double tolerance = pow(.1, DECIMALS);
			//iterate on core variables
			for(CoinVar* var: coreVars_){
				double value = getVarValue(var);
				if( value > tolerance)
				fprintf(pFile,"%-30s %4.2f (%6.0f) \n", var->name_ , value, var->getCost());
			}
			//iterate on column variables
			for(MyVar* var: activeColumnVars_){
				double value = getVarValue(var);
				if( value > tolerance)
				fprintf(pFile,"%-30s %4.2f (%6.0f) \n", var->name_ , value, var->getCost());
			}

			fprintf(pFile,"\n");
		}
		if (!logfile_.empty()) fclose(pFile);

		return 1;
	}

	virtual bool loadBestSol() { return false; }

	virtual int writeProblem(std::string fileName) { return 0; }

	virtual int writeLP(std::string fileName) { return 0; }

	virtual void toString(MyObject* obj){
		CoinVar* var = dynamic_cast<CoinVar*>(obj);
		if(var) var->toString(cons_);
		else Modeler::toString(obj);
	}

	//get the variables that are always present in the model
  std::vector<CoinVar*>& getCoreVars(){ return coreVars_; }

    std::vector<CoinCons*>& getCons(){ return cons_; }

protected:

    std::vector<CoinVar*> coreVars_;
    std::vector<CoinCons*> cons_;
};

#endif /* SRC_COINMODELER_H_ */
