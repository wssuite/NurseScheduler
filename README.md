# Static nurse scheduler

[![Build Status](https://travis-ci.org/legraina/StaticNurseScheduler.svg?branch=master)](https://travis-ci.org/legraina/StaticNurseScheduler)

The C++ code and the instances shared on this repository allow to solve a static variant of the nurse scheduling  problem described in the second international nurse rostering competition (INRC-II). The difference with the original formulation is that the competition deals with a dynamic revealing of the data, where the demand and the nurses' preferences on a given week are only revealed once the schedules of the previous weeks are computed. Here, the demands and the preferences of the complete horizon are known beforehand.
Every information about INRC-II can be found on their website http://mobiz.vives.be/inrc2/, and the initial description of the problem is given in:
[1] S. Ceschia, N. Thi, T. Dang, and P. De Causmaecker, "Second International Nurse Rostering Competition (INRC-II): Problem Description and Rules." p. 1–18, 2015.
The methods implemented in this code are all described in the following manuscript (still under revision). Please cite this reference in any use of our code.
[2] A. Legrain, J. Omer and S. Rosat, "A rotation-based branch-and-price approach for the nurse scheduling problem", p. 1-28, 2017, submitted manuscript.
These two references can also be found in the directory ./references.

Guide
------------------

The following describes how to handle our code.

1) To install the required libraries and build the code, please follow the instructions detailed in the INSTALL.md file.

2) Content of the project.

	a. ./src : source code (executables are stored in ./bin after building and object files are stored in ./obj)

	b. ./datasets: We provide the benchmark used by the organizers of the INRC2 [1] in the ./datasets directory with the format nXXXwY, where XXX refers to the number of nurses and Y is the number of weeks in the planning horizon. For each number of nurses and planning horizon, several history and demand files are provided thus allowing to test a very large number of different instances.

	c. ./paramfiles : Directory where all the parameters of the solution methods are stored. This is where the particular method executed when running the executable is chosen. The parameters files initially present in the directory are those used for the tests in [2].

	d. validator.jar is the java executable provided by the organizers of INRC2 to check the validity of a solution and compute its cost independently.

	e. ./scripts : Directory including several scripts for the execution of multiple runs at once

	f. ./bashfiles : Directory containing bash files that execute the solver on an instance (the files are created by the scripts)

3) Global structure of the code:

	Every source file is in the ./src directory, where header files are used to declare the classes, and methods.

	a. The main is in "DeterministicMain.cpp" and "DeterministicMain_test.cpp" is for the definition of some unitary tests.

	b. Input data and basic preprocessing methods are in "Nurse.h/.cpp", "Roster.h", "Scenario.h/.cpp" and stored in an instance of the class defined in "SolverInput.h"

	c. "Solver.h" stores an abstract solver class and the definitions of several other classes of objects manipulated by the algorithm. In particular a LiveNurse has the constant attributes of a Nurse and other attributes that will be modified during the execution of the solution algorithm. The files "InitializeSolver.h/.cpp" runs some preprocessing actions before actually solving the problem.

	d. "DeterministicSolver.h/.cpp" contains the declaration and the structure of every algorithm described in [2] as well as all the methods involved in the large neighborhood search. Every method that BCP needs redefined for the branch-and-price algorithm are in "BcpModeler.h/.cpp", "CoinModeler.h" and in "TreeManager.h/.cpp". The global structure of the column generation subproblem, including the construction of the constrained shortest path network, is in "Subproblem.h/.cpp", and the dynamic programming algorithm that solves the subproblems is implemented in "RotationPricer.h/.cpp" (it is adapted from an algorithm found in the Boost library).

	e. The postprocessing/display/parsing methods are in "ReadWrite.h/.cpp", "GlobalStats.h/.cpp", and "InputPaths.h/.cpp".

	f. The files "MyTools.h/.cpp" contain intermediary methods frequently used in the code.

4) Execution of the code:

	a. A typical execution of our code is done from the root directory of the project with the following list of arguments:
	
	````bash
	./bin/staticscheduler --dir datasets/ --instance n030w4 --weeks 6-2-9-1 --his 1 --param paramfiles/default.txt --sol outfiles/default/n030w4_1_6-2-9-1 --timeout 780
	
	--dir is followed by the directory where the instance is stored
	--instance is the name of the subdirectory of where the specific instance is stored
	--weeks is the sequence of week files numbers that form the complete horizon
	--his is the number of the history file among those published
	--param is followed by the name of the parameter file used in this run
	--sol is the directory where the solution will be stored
	--timeout is the total execution time
	````
	
	The validator can then be run by:
	````bash
	java -jar validator.jar --sce datasets/n030w4/Sc-n030w4.txt --his datasets/n030w4/H0-n030w4-1.txt --weeks datasets/n030w4/WD-n030w4-6.txt datasets/n030w4/WD-n030w4-2.txt datasets/n030w4/WD-n030w4-9.txt datasets/n030w4/WD-n030w4-1.txt --sols outfiles/default/n030w4_1_6-2-9-1/sol-week0.txt outfiles/default/n030w4_1_6-2-9-1/sol-week1.txt outfiles/default/n030w4_1_6-2-9-1/sol-week2.txt outfiles/default/n030w4_1_6-2-9-1/sol-week3.txt > outfiles/default/n030w4_1_6-2-9-1/validator.txt
	````
	
	or:
	````bash
	./validator.sh n030w4 6-2-9-1 1 outfiles/default/n030w4_1_6-2-9-1
	````
	
	All the results can then be found in the "outfiles/default/n030w4_1_6-2-9-1" directory (replace default with the name of the parameter file you used)

	b. Other options for a quicker run of the code are:
	
		- run the solver with default options on the instance n005w4_1_1-6-2-9-1:
		````bash
		./bin/staticscheduler
		````
		
		- run the solver on the instance n005w4_0_2-0-2-1 with options defined in paramfiles/default.txt:
		````bash
		./bin/staticscheduler --dir datasets/ --instance n005w4 --his 0 --weeks 2-0-2-1 --param paramfiles/default.txt
		````
		
		- run the solver on the instance n005w4_0_2-0-2-1 with default options:
		````bash
		./bin/deterministicroster --his testdatasets/n005w4/H0-n005w4-0.txt --sce testdatasets/n005w4/Sc-n005w4.txt --week testdatasets/n005w4/WD-n005w4-2.txt  --week testdatasets/n005w4/WD-n005w4-0.txt --week testdatasets/n005w4/WD-n005w4-2.txt --week testdatasets/n005w4/WD-n005w4-1.txt
		````
		
		- run a test with name testname:
		````bash
		./bin/staticscheduler --test testname
		````

	c. Scritps for running several instances at once :
	
		- to write the bash files that run the solver on all the instances with a specific set of parameters defined in the file "paramfiles/param.txt" :
		````bash
		./scripts/writeAllRuns.sh param
		````
		
		- to run all these bashfiles one after the other :
		````bash
		./scripts/runDir.sh param
		````
		The outputs will then be written in "outfiles/param/"


5) There are some random aspects in our solver (in the large neighborhood search for instance) and in the third party libraries that are called by our solver. For instance, the perturbations added by CLP to avoid degeneracy will not impact the objective value, but they can impact the specific optimal solution, and hence the dual solution, which can lead to differences in the subproblem. As a consequence, the solution values can be slightly different from those reported in [2]. In our tests on several different machines, this has not impacted the interpretations and comparisons discussed in [2] though.

6) Description of some notations that appear in the code/comment:

	- rotation: sequence of working days for a nurse. The shifts that are covered and the skills that are used can be different on each day. A rotation starts at the beginning of a week or after a resting day, and it ends at the end of a week or before a resting day.
	
	- break/pause/holiday/rest period: sequence of resting days. This period starts at the beginning of a week or after a working day, and it ends at the end of a week or before a working day.
	
	- stint: one rotation followed by a break
