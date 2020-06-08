# Nurse scheduler

[![Docker tests](https://github.com/wssuite/NurseScheduler/workflows/Docker%20tests/badge.svg)](https://github.com/wssuite/NurseScheduler/actions?query=workflow%3A%22Docker+tests%22+branch%3Amaster+event%3Apush)
[![Docker Automated build](https://img.shields.io/docker/automated/legraina/nurse-scheduler)](https://hub.docker.com/repository/docker/legraina/nurse-scheduler/)
[![DOI](https://zenodo.org/badge/150300357.svg)](https://zenodo.org/badge/latestdoi/150300357)


The C++ code and the instances shared on this repository allow to solve a static variant of the nurse scheduling problem described in the second international nurse rostering competition (INRC-II) as well as the original dynamic problem.
Every information about INRC-II can be found on their website http://mobiz.vives.be/inrc2/, and the initial description of the problem as well as the methodology corresponding to this code are given in:
[1] S. Ceschia, N. Thi, T. Dang, and P. De Causmaecker, "Second International Nurse Rostering Competition (INRC-II): Problem Description and Rules.", CoRR, p. 1–18, 2015.

[2] A. Legrain, J. Omer and S. Rosat, "An Online Stochastic Algorithm for a Dynamic Nurse Scheduling Problem", European Journal of Operational Research, 2018. The methods implemented in this code are all described in this manuscript. Please cite it in any use of our dynamic code. You can also use this doi: [10.5281/zenodo.3460633](https://doi.org/10.5281/zenodo.3460633) to make a direct reference to the dynamic code.

[3] A. Legrain, J. Omer and S. Rosat, "A rotation-based branch-and-price approach for the nurse scheduling problem", p. 1-28, 2017, submitted manuscript. This reference describes the offline solver used for solving the nurse scheduling problem. Please cite it in any use of our static code. You can also use this doi: [10.5281/zenodo.3460634](https://doi.org/10.5281/zenodo.3460634) to make a direct reference to the static code.

These two references can also be found in the directory ./references.

Guide
------------------

The following describes how to handle our code.

1. To install the required libraries and build the code, please follow the instructions detailed in the INSTALL.md file. You can also read DOCKER.md if you want to run the code within a container.

2. Content of the project.

	a. ./src : source code (executables are stored in ./bin after building and object files are stored in ./obj)

	b. ./datasets: We provide the benchmark used by the organizers of the INRC2 [1] in the ./datasets directory with the format nXXXwY, where XXX refers to the number of nurses and Y is the number of weeks in the planning horizon. For each number of nurses and planning horizon, several history and demand files are provided thus allowing to test a very large number of different instances.

	c. ./paramfiles : Directory where all the parameters of the solution methods are stored. This is where the particular method executed when running the executable is chosen. The parameters files initially present in the directory are those used for the tests in [2].

	d. ./Simulator.jar is the java executable provided by the organizers of INRC2 to run the solver week by week on a sequence of instances.

	e. validator.jar is the java executable provided by the organizers of INRC2 to check the validity of a solution and compute its cost independently.

	f. ./scripts : Directory including several scripts for the execution of multiple runs at once

	g. ./bashfiles : Directory containing bash files that execute the solver on an instance (the files are created by the scripts)
	
	h. ./benchmark : Directory containing the different benchmark on which the code can be run in a static and dynamic version. The best obtained results are also stored there. You can edit this file and run our code on a benchmark with the python script "run-benchmark.py".
	
3. Global structure of the code:

	Every source file is in the ./src directory, where header files are used to declare the classes, and methods.

	a. The main is in "DeterministicMain.cpp" and "DynamicMain.cpp". "DeterministicMain_test.cpp" is for the definition of some unitary tests.

	b. Input data and basic preprocessing methods are stored in directory "data/"
	and correspond to the files "Nurse.h/.cpp", "Roster.h", "Scenario.h/.cpp". 

	c. The directory "solvers/" contains all the algorithms. In particular:
	- "Solver.h" stores an abstract solver class and the definitions of several other classes of objects manipulated by the algorithm. In particular a LiveNurse has the constant attributes of a Nurse and other attributes that will be modified during the execution of the solution algorithm. The files "InitializeSolver.h/.cpp" runs some preprocessing actions before actually solving the problem.
    - "DeterministicSolver.h/.cpp" contains the declaration and the structure of every algorithm described in [2] as well as all the methods involved in the large neighborhood search. 
    - Every method that BCP needs are redefined for the branch-and-price algorithm are in "mp/modeler/BcpModeler.h/.cpp", "mp/modeler/CoinModeler.h" and in "mp/TreeManager.h/.cpp". The global structure of the column generation subproblem, including the construction of the constrained shortest path network, is in "mp/sp/Subproblem.h/.cpp", and the dynamic programming algorithm that solves the subproblems is implemented in "mp/sp/rcspp/BoostRCGraph.h/.cpp" (it is adapted from an algorithm found in the Boost library).
    - "StochasticSolver.h/.cpp" contains the declaration and the structure of the algorithm described in [3] and really close to the one submitted to INRCII.

	d. The postprocessing/display/parsing methods are stored in the directory "tools/".
	   Especially, the files "MyTools.h/.cpp" contain intermediary methods frequently used in the code.

4. Execution of the deterministic solver:

	a. A typical execution of our code is done from the root directory of the project with the following list of arguments:

	```bash
	./bin/staticscheduler --dir datasets/ --instance n030w4 --weeks 6-2-9-1 --his 1 --param paramfiles/default.txt --sol outfiles/default/n030w4_1_6-2-9-1 --timeout 780

	--dir is followed by the directory where the instance is stored
	--instance is the name of the subdirectory of where the specific instance is stored
	--weeks is the sequence of week files numbers that form the complete horizon
	--his is the number of the history file among those published
	--param is followed by the name of the parameter file used in this run
	--sol is the directory where the solution will be stored
	--timeout is the total execution time
	```

	The validator can then be run by:
	```bash
	java -jar validator.jar --sce datasets/n030w4/Sc-n030w4.txt --his datasets/n030w4/H0-n030w4-1.txt --weeks datasets/n030w4/WD-n030w4-6.txt datasets/n030w4/WD-n030w4-2.txt datasets/n030w4/WD-n030w4-9.txt datasets/n030w4/WD-n030w4-1.txt --sols outfiles/default/n030w4_1_6-2-9-1/sol-week0.txt outfiles/default/n030w4_1_6-2-9-1/sol-week1.txt outfiles/default/n030w4_1_6-2-9-1/sol-week2.txt outfiles/default/n030w4_1_6-2-9-1/sol-week3.txt > outfiles/default/n030w4_1_6-2-9-1/validator.txt
	```
	or:
	```bash
	./validator.sh n030w4 6-2-9-1 1 outfiles/default/n030w4_1_6-2-9-1
	```

	All the results can then be found in the "outfiles/default/n030w4_1_6-2-9-1" directory (replace default with the name of the parameter file you used)

  b. Other options for a quicker run of the code are:
  
   - run the solver with default options on the instance n005w4_1_1-6-2-9-1:
   ````bash
   ./bin/staticscheduler
   ````
	
   - run the solver on the instance n005w4_0_2-0-2-1 with options defined in paramfiles/default.txt:
   ```bash
   ./bin/staticscheduler --dir datasets/ --instance n005w4 --his 0 --weeks 2-0-2-1 --param paramfiles/default.txt
   ```
	
   - run the solver on the instance n005w4_0_2-0-2-1 with default options:
   ```bash
   ./bin/staticscheduler --his testdatasets/n005w4/H0-n005w4-0.txt --sce testdatasets/n005w4/Sc-n005w4.txt --week testdatasets/n005w4/WD-n005w4-2.txt  --week testdatasets/n005w4/WD-n005w4-0.txt --week testdatasets/n005w4/WD-n005w4-2.txt --week testdatasets/n005w4/WD-n005w4-1.txt
   ```
	
   - run a test with name testname:
   ```bash
   ./bin/staticscheduler --test testname
   ```

  c. Scripts located in folder "scripts/" to generate new scripts that run the determistic solver. Note that the outputs will then be written in "outfiles/param/".

   - writeRun.sh writes a bash file that runs the solver on a specific instance with a specific set of parameters defined in the folder "paramfiles/". For example:
   ````bash
   ./scripts/writeRun.sh --instance n005w4_0_2-0-2-1 --param default.txt
   ````
   You can use the option "-h" to see all the available flags.
   
   - writeAllRuns.sh writes bash files that run the solver on all the instances with a specific set of parameters defined in the folder "paramfiles/". The instances to run are hard-written in the script. Example of use:
   ````bash
   ./scripts/writeAllRuns.sh --param lns_feas.txt
   ````
   You can use the option "-h" to see all the available flags (the same than writeRun.sh except --instance which is useless).
   
   - runDir.sh runs all the bashfiles one after the other within the folder associated to a param folder name. Example for the bashfiles within "bashfiles/lns_repeat/":
   ````bash
   ./scripts/runDir.sh lns_repeat
   ````

5. Execution of the stochastic solver:

	a. Generate a script for a given instance and seeds for example:
	````bash
	./scripts/writeDynamicRun.sh -i n005w4_1-2-3-3_0 -s 22-36-96-5
	````
	You can use the option "-h" to see all the available flags.

	b. Then, run it:
	````bash
	./scripts/n005w4_1-2-3-3_0_22-36-96-5.sh
	````

6. There are some random aspects in our solver (in the large neighborhood search for instance) and in the third party libraries that are called by our solver. For instance, the perturbations added by CLP to avoid degeneracy will not impact the objective value, but they can impact the specific optimal solution, and hence the dual solution, which can lead to differences in the subproblem. As a consequence, the solution values can be slightly different from those reported in [2] and [3]. In our tests on several different machines, this has not impacted the interpretations and comparisons discussed in [2] and [3] though.

7. Description of some notations that appear in the code/comment:

	- rotation: sequence of working days for a nurse. The shifts that are covered and the skills that are used can be different on each day. A rotation starts at the beginning of a week or after a resting day, and it ends at the end of a week or before a resting day.

	- break/pause/holiday/rest period: sequence of resting days. This period starts at the beginning of a week or after a working day, and it ends at the end of a week or before a working day.

	- stint: one rotation followed by a break
