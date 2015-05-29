/*
 * OptimalSolver.cpp
 *
 *  Created on: 2015-05-16
 *      Author: legraina
 */

#include "main_test.h"
#include "MyTools.h"

// n030w4 1 4 6 2 9 1 n030w4_1_6-2-9-1
// Function for solving the optimal solution
int main(int argc, char** argv)
{
	std::cout << "# Launching Optimal Roster..." << std::endl;

   int locINT;

   string inst = argv[1];

   std::istringstream(argv[2]) >> locINT;
   int historyID = locINT;

   std::istringstream(argv[3]) >> locINT;
   int nbWeeks = locINT;
   if(5+nbWeeks > argc)
      Tools::throwError("Bad input format. Should be: instance_name historyID numberWeeks vector<int>weekIndices outdir prefix numberTest solver_options generation_options evaluations_options. After outdir, the arguments are optional.");

   vector<int> numberWeek;
   for(int i=4; i<4+nbWeeks; ++i){
      std::istringstream(argv[i]) >> locINT;
      numberWeek.push_back(locINT);
   }

   string outdir = argv[4+nbWeeks];

   string prefix = "";
   if(5+nbWeeks < argc)
        prefix = argv[5+nbWeeks];

   int givenSeed;
   if(6+nbWeeks < argc){
	      std::istringstream(argv[6+nbWeeks]) >> locINT;
	      givenSeed = locINT;
   }

   string data = "datasets/";
   string scenarPath = data + inst + "/Sc-" + inst + ".txt";
   string outpath = "outfiles/Competition/"+outdir+"/"+prefix;
   string outfile = outpath+"sol-week";
   string logfile = outpath+"Log.txt";
   string optionspath = data+"optionFiles/";

   string stoOptionsFile = optionspath+"stochastic_solver.txt";
   string geneOptionsFile = optionspath+"generation_solver.txt";
   string evaOptionsFile = optionspath+"evaluation_solver.txt";

   // Time the complete execution of the algorithm
   Tools::Timer* timertotal = new Tools::Timer();
   timertotal->init();
   timertotal->start();

//   SolverParam optParam;
//   optParam.printEverySolution_ = true;
//   optParam.weekIndices_ = numberWeek;
//   optParam.outfile_ = outfile;
//   optParam.logfile_ = logfile;
//   optParam.solveToOptimality_ = true;
//   optParam.nbDiveIfMinGap_ = 2;
//   optParam.nbDiveIfRelGap_ = 8;
//   testMultipleWeeksDeterministic(data, inst, historyID, numberWeek, GENCOL, "outfiles/Competition/"+outdir+"/"+prefix, optParam);


   StochasticSolverOptions stochasticSolverOptionsScore;
   setStochasticSolverOptions(stochasticSolverOptionsScore, SUNGRID, inst, outfile, outpath,
		   stoOptionsFile, geneOptionsFile, evaOptionsFile);

   srand(givenSeed);
   int seed = rand();

   pair<double, int> p = testMultipleWeeksStochastic(data, inst, historyID, numberWeek, stochasticSolverOptionsScore, outpath+"score_", seed);

   char results[250];
   sprintf(results, "Seed %d; Cost %.2f; NbGene %d; NbEval %d; WeightStrat %d; RankStrat %d;  nbDaysGeneration %d; nbDaysEvaluation %d;",
		   seed, p.first, p.second, stochasticSolverOptionsScore.nEvaluationDemands_,
		   stochasticSolverOptionsScore.generationParameters_.weightStrategy_, stochasticSolverOptionsScore.rankingStrategy_,
		   7+stochasticSolverOptionsScore.nExtraDaysGenerationDemands_, stochasticSolverOptionsScore.nDaysEvaluation_);
   string sensibilityOutfile = outpath+"score_sensibility.txt";
   Tools::LogOutput sensibilityStream(sensibilityOutfile, true);
   sensibilityStream << results << std::endl;

   StochasticSolverOptions stochasticSolverOptionsMean;
   setStochasticSolverOptions(stochasticSolverOptionsMean, SUNGRID, inst, outfile, outpath,
		   stoOptionsFile, geneOptionsFile, evaOptionsFile);

   stochasticSolverOptionsMean.rankingStrategy_ = RK_MEAN;

   p = testMultipleWeeksStochastic(data, inst, historyID, numberWeek, stochasticSolverOptionsMean, outpath+"mean_", seed);

   char results2[250];
   sprintf(results2, "Seed %d; Cost %.2f; NbGene %d; NbEval %d; WeightStrat %d; RankStrat %d;  nbDaysGeneration %d; nbDaysEvaluation %d;",
		   seed, p.first, p.second, stochasticSolverOptionsMean.nEvaluationDemands_,
		   stochasticSolverOptionsMean.generationParameters_.weightStrategy_, stochasticSolverOptionsMean.rankingStrategy_,
		   7+stochasticSolverOptionsMean.nExtraDaysGenerationDemands_, stochasticSolverOptionsMean.nDaysEvaluation_);
   string sensibilityOutfile2 = outpath+"mean_sensibility.txt";
   Tools::LogOutput sensibilityStream2(sensibilityOutfile2, true);
   sensibilityStream2 << results2 << std::endl;

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
