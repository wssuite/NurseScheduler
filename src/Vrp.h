/*
 * Vrp.h
 *
 *  Created on: 2015-02-19
 *      Author: legraina
 */

#ifndef SRC_VRP_H_
#define SRC_VRP_H_

/* standard library includes */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* user defined includes */
#include "Pricer_vrp.h"

class Vrp{
public:
	// Constructor and destructor
	//
	Vrp(string dataFile);
	~Vrp();

public:
	//solve the vrp with scip
	int solve(string dataFile);
	//read the data problem
	int read_problem(
	   const char*           filename,           /**< filename */
	   int&                  num_nodes,          /**< number of nodes in instance */
	   int&                  capacity,           /**< capacity in instance */
	   std::vector<int>&          demand,             /**< array of demands of instance */
	   std::vector<std::vector<int> >& distance            /**< distances between nodes */
	   );
};

#endif /* SRC_VRP_H_ */
