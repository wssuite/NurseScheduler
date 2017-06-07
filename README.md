# Static nurse scheduler
# ======================

The code and the instances shared on this repository allow to solve a static variant of the nurse scheduling  problem as described in the INRC2 competition. The difference with the original formulation is that the competition deals with a dynamic revealing of the data, where the demand and the nurses' preferences on a given week are only revealed once the schedules of the previous weeks are computed. Here, the demands and the preferences of the complete horizon are known beforehand. 
Every information about INRC2 can be found on their website http://mobiz.vives.be/inrc2/, and the initial description of the problem is given in:
[1] S. Ceschia, N. Thi, T. Dang, and P. De Causmaecker, "Second International Nurse Rostering Competition ( INRC-II ): Problem Description and Rules." p. 1–18, 2015.
The method implemented in this code are all described in the following manuscript (still under revision). Please cite this reference in any use of our code.
[2] A. Legrain, J. Omer and S. Rosat, "A rotation-based branch-and-price approach for the nurse scheduling problem", p. 1-28, 2017, submitted manuscript.
These two references can also be found in the directory ./references.

The following decribes how to handle our code.

1) To install the required libraries and build the code, please follow the instructions detailed in the INSTALL.md file.

2) Content of the projet.
	a) ./src : source code (executables are stored in ./bin after building and object files are stored in ./obj)
	b) ./datasets: We provide the benchmark used by the organizers of the INRC2 [1] in the ./datasets directory with the format nXXXwY, where XXX refers to the nnumber of nurses and Y is the number of weeks in the planning horizon. For each number of nurses and planning horizon, several history and demand files are provided thus allowing to test a very large number of different instances.
	c) ./paramfiles : Directory where all the parameters of the solution methods are stored. This is where the particular method executed when running the executable is chosen. The parameters files initially present in the directory are those used for the tests in [2].
	d) validator.jar is the java executable provided by the organizers of INRC2 to check the validity of a solution and compute its cost independently.


2) A typical execution of our code is done from the root directory of the project with the following list of arguments

	./bin/staticscheduler --dir datasets/ --instance n030w4 --weeks 6-2-9-1 --his 1 --param paramfiles/default.txt --sol outfiles/default/n030w4_1_6-2-9-1 --timeout 780

	--dir is followed by the directory where the instance is stored
	--instance is the name of the subdirectory of where the specific instance is stored
	--weeks is the sequence of week files numbers that form the complete horizon
	--his is the number of the history file among those published
	--param is followed by the name of the parameter file used in this run
	--sol is the directory where the solution will be stored
	--timeout is the total execution time

	The validator can then be run by:
	java -jar validator.jar --sce datasets/n030w4/Sc-n030w4.txt --his datasets/n030w4/H0-n030w4-1.txt --weeks datasets/n030w4/WD-n030w4-6.txt datasets/n030w4/WD-n030w4-2.txt datasets/n030w4/WD-n030w4-9.txt datasets/n030w4/WD-n030w4-1.txt --sols outfiles/default/n030w4_1_6-2-9-1/sol-week0.txt outfiles/default/n030w4_1_6-2-9-1/sol-week1.txt outfiles/default/n030w4_1_6-2-9-1/sol-week2.txt outfiles/default/n030w4_1_6-2-9-1/sol-week3.txt > outfiles/default/n030w4_1_6-2-9-1/validator.txt

3) Other options for a quicker run of the code are

	a) run the solver with default options on the instance n005w4_1_1-6-2-9-1
	./bin/staticscheduler
	b) run the solver on the instance n005w4_0_2-0-2-1 with options defined in paramfiles/default.txt
	./bin/staticscheduler --dir datasets/ --instance n005w4 --his 0 --weeks 2-0-2-1 --param paramfiles/parameters.txt
	c) run the solver on the instance n005w4_0_2-0-2-1 with default options
	./bin/deterministicroster --his testdatasets/n005w4/H0-n005w4-0.txt --sce testdatasets/n005w4/Sc-n005w4.txt --week testdatasets/n005w4/WD-n005w4-2.txt  --week testdatasets/n005w4/WD-n005w4-0.txt --week testdatasets/n005w4/WD-n005w4-2.txt --week testdatasets/n005w4/WD-n005w4-1.txt
	d) run a test with name testname
	./bin/staticscheduler --test testname

4) Description of some notations that appear in the code/comment
	- rotation: sequence of working days for a nurse. The shifts that are covered and the skills that are used can be different on each day. A rotation starts at the beginning of a week or after a resting day, and it ends at the end of a week or before a resting day.
	- break/pause/holiday/rest period: sequence of resting days. This period starts at the beginning of a week or after a working day, and it ends at the end of a week or before a working day.
	- stint: one rotation followed by a break
