/*
 * OptimalSolver.cpp
 *
 *  Created on: 2015-05-16
 *      Author: legraina
 */

#include "main_test.h"
#include "MyTools.h"


// Function for solving the optimal solution
int main(int argc, char** argv)
{
   int locINT;

   char* inst = argv[1];

   std::istringstream(argv[2]) >> locINT;
   int historyID = locINT;
   vector<int> numberWeek;
   for(int i=3; i<argc-1; ++i){
      std::istringstream(argv[i]) >> locINT;
      numberWeek.push_back(locINT);
   }
//   if( (numberWeek.size()%4) != 0 )
//      Tools::throwError("Bad instance.");
   string outdir = argv[argc-1];

   // Time the complete execution of the algorithm
   Tools::Timer* timertotal = new Tools::Timer();
   timertotal->init();
   timertotal->start();

   string data = "datasets/";
   string scenarPath = data + inst + "/Sc-" + inst + ".txt";

   SolverParam optParam;
   optParam.nbDiveIfMinGap_ = 2;
   optParam.nbDiveIfRelGap_ = 4;
//   testMultipleWeeksDeterministic(data, inst, historyID, numberWeek, GENCOL, "outfiles/Competition/"+outdir+"/Opt", optParam);
   testMultipleWeeksStochastic(data, inst, historyID, numberWeek, GENCOL, "outfiles/Competition/"+outdir+"/Opt");

   // Display the total time spent in the algorithm
   timertotal->stop();
   std::cout << "Total time spent in the algorithm : " << timertotal->dSinceInit() << std::endl;

   // free the allocated pointers
   //
   delete timertotal;

//vector<string> instances =
//{
//   "n030w4_1_6-2-9-1", "n030w4_1_6-7-5-3", "n030w8_1_2-7-0-9-3-6-0-6", "n030w8_1_6-7-5-3-5-6-2-9",
//   "n040w4_0_2-0-6-1", "n040w4_2_6-1-0-6", "n040w8_0_0-6-8-9-2-6-6-4", "n040w8_2_5-0-4-8-7-1-7-2",
//   "n050w4_0_0-4-8-7", "n050w4_0_7-2-7-2", "n050w8_1_1-7-8-5-7-4-1-8", "n050w8_1_9-7-5-3-8-8-3-1",
//   "n060w4_1_6-1-1-5", "n060w4_1_9-6-3-8", "n060w8_0_6-2-9-9-0-8-1-3", "n060w8_2_1-0-3-4-0-3-9-1",
//   "n080w4_2_4-3-3-3", "n080w4_2_6-0-4-8", "n080w8_1_4-4-9-9-3-6-0-5", "n080w8_2_0-4-0-9-1-9-6-2",
//   "n100w4_0_1-1-0-8", "n100w4_2_0-6-4-6", "n100w8_0_0-1-7-8-9-1-5-4", "n100w8_1_2-4-7-9-3-9-2-8",
//   "n120w4_1_4-6-2-6", "n120w4_1_5-6-9-8", "n120w8_0_0-9-9-4-5-1-0-3", "n120w8_1_7-2-6-4-5-2-0-2"
//};

}
