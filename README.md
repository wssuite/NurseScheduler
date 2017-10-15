# Dynamic Nurse Scheduler
The C++ code shared on this repository has been submitted by the Polytechnique team to the International Nurse Rostering Competition II. Every information about INRC-II can be found on their website http://mobiz.vives.be/inrc2/, and the initial description of the problem is given in:

[1] S. Ceschia, N. Thi, T. Dang, and P. De Causmaecker, "Second International Nurse Rostering Competition (INRC-II): Problem Description and Rules." p. 1–18, 2015.

[2] A. Legrain, J. Omer and S. Rosat, "An Online Stochastic Algorithm for a Dynamic Nurse Scheduling Problem", p. 1-26, 2017, submitted manuscript. The methods implemented in this code are all described in this manuscript (still under revision). Please cite it in any use of our code.

[3] A. Legrain, J. Omer and S. Rosat, "A rotation-based branch-and-price approach for the nurse scheduling problem", p. 1-28, 2017, submitted manuscript. This reference describes the offline solver used for solving the nurse scheduling problem.

These three references can also be found in the directory ./references.

Guide
------------------

The following describes how to handle our code.

1) To install the required libraries and build the code, please follow the instructions detailed in the INSTALL.md file. You can also read DOCKER.md if you want to run the code within a container.

2) Content of the project.

	a. ./src : source code (executables are stored in ./bin after building and object files are stored in ./obj)
	
	b. ./datasets: We provide the benchmark used by the organizers of the INRC2 [1] in the ./datasets directory with the format nXXXwY, where XXX refers to the number of nurses and Y is the number of weeks in the planning horizon. For each number of nurses and planning horizon, several history and demand files are provided thus allowing to test a very large number of different instances.
	
	c. ./Simulator.jar is the java executable provided by the organizers of INRC2 to run the solver week by week on a sequence of instances.
	
	d. ./validator.jar is the java executable provided by the organizers of INRC2 to check the validity of a solution and compute its cost independently.
	
	e. ./generateScript.sh is a bash script that can generate a specific bash script to run the Simulator and the validator on a given instance.

3) Global structure of the code:
	Every source file is in the ./src directory, where header files are used to declare the classes, and methods.
	
	a. The main is in "main.cpp".
	
	b. Input data and basic preprocessing methods are in "Nurse.h/.cpp", "Roster.h", "Scenario.h/.cpp" and stored in an instance of the class defined in "SolverInput.h"
	
	c. "Solver.h" stores an abstract solver class and the definitions of several other classes of objects manipulated by the algorithm. In particular a LiveNurse has the constant attributes of a Nurse and other attributes that will be modified during the execution of the solution algorithm. The files "InitializeSolver.h/.cpp" runs some preprocessing actions before actually solving the problem.
	
	d. "StochasticSolver.h/.cpp" contains the declaration and the structure of every algorithm described in [2]. Every method that BCP needs redefined for the branch-and-price algorithm are in "BcpModeler.h/.cpp", "CoinModeler.h". The global structure of the column generation subproblem, including the construction of the constrained shortest path network, is in "Subproblem.h/.cpp", and the dynamic programming algorithm that solves the subproblems is implemented in "RotationPricer.h/.cpp" (it is adapted from an algorithm found in the Boost library).
	
	e. The IO methods are in "ReadWrite.h/.cpp".
	
	f. The files "MyTools.h/.cpp" contain intermediary methods frequently used in the code.

4) Execution of the code:

	a. Generate a script for a given instance and seeds (not compulsory):
	````bash
	./generateScript.sh n005w4_1-2-3-3_0 22-36-96-5
	````
  
  
	b. Then, run it:
	````bash
	./n005w4_1-2-3-3_0_22-36-96-5.sh
	````


5) There are some random aspects in our solver (in the sampling of the demand) and in the third party libraries that are called by our solver. For instance, the perturbations added by CLP to avoid degeneracy will not impact the objective value, but they can impact the specific optimal solution, and hence the dual solution, which can lead to differences in the subproblem. As a consequence, the solution values can be slightly different from those reported in [2]. In our tests on several different machines, this has not impacted the interpretations and comparisons discussed in [2] though.

6) Description of some notations that appear in the code/comment

	- rotation: sequence of working days for a nurse. The shifts that are covered and the skills that are used can be different on each day. A rotation starts at the beginning of a week or after a resting day, and it ends at the end of a week or before a resting day.
	
	- break/pause/holiday/rest period: sequence of resting days. This period starts at the beginning of a week or after a working day, and it ends at the end of a week or before a working day.
	
	- stint: one rotation followed by a break
