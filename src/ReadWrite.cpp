#include "ReadWrite.h"
#include "MyTools.h"
#include <iostream>
#include <streambuf>
#include <fstream>
#include <math.h>
#include <time.h>


//--------------------------------------------------------------------------
// Methods that read all the input files and store the content in the
// input scenario instance
//

// Read the scneario file and store the content in a Scenario instance
//
void readScenario(std::string strWeekFile, Scenario* pScenario) {


	// open the file
	std::fstream file;
	std::cout << "Reading " << strWeekFile << std::endl;
	file.open(strWeekFile.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "While trying to read " << strWeekFile << std::endl;
		ToolsThrow("The input file was not opened properly!");
	}

	std::string title;
	char   charTmp[256];
	int tmp;

	// fill the attributes of the scenario structure
	//
	readUntilChar(&file, '=', &title);
	while ( file.good() ){
		if (!strcmp(title.c_str(), "SCENARIO")) {
			file >>charTmp;
			pScenario->name(charTmp);
		}
		else if (!strcmp(title.c_str(), "WEEKS"))	{
			file >> intTmp;
			pScenario->nbWeeks(intTmp);
		}
		else if (!strcmp(title.c_str(), "SKILLS"))	{
			file >> intTmp;
			pScenario->nbSkills(intTmp);
			for (int i = 0; i < intTmp; i++) 
		}

		file >> pParam_->tIni;
		else if (!strcmp(title.c_str(), "tFin"))
		file >> pParam_->tFin;
		else if (!strcmp(title.c_str(), "coCordesV"))
		file >> pParam_->coCordesV;
		else if (!strcmp(title.c_str(), "dCV"))
		file >> pParam_->dCV;
		else if (!strcmp(title.c_str(), "coCordesA"))
		file >> pParam_->coCordesA;
		else if (!strcmp(title.c_str(), "coCostLines"))
		file >> pParam_->coCostLines;
		else if (!strcmp(title.c_str(), "fixedCostPercent"))
		file >> pParam_->fixedCostPercent;
		else
			ToolsThrow("The title of the line is unknown!");

				file.getline(charTmp, 256);
				if (!ToolsReadUntilChar(&file, '=', &title)) continue;
			}


}

// Read the Week file and store the content in a Scenario instance
//
void readWeek(std::string strWeekFile, Scenario* pScenario);


// Read the history file
//
void readHistory(std::string strHistoryFile, Scenario* pScenario);

// Read the input custom file and store the content in a Scenario instance
//
void readCustom(std::string strCustomInputFile, Scenario* pScenario);

//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// Methods that write the ouputs of the solver

// Write the solution file for the current week
//
void writeSolution(std::string strCustomOutputFile, Solution* pSolution);

// Write the output custom file from values in the scenario and the solution
// instances
//
void writeCustom(std::string strCustomOutputFile, Scenario* pScenario, Solution* pSolution);

//--------------------------------------------------------------------------

// Set the parameters values by reading them in a file
//
void CftSolver::readParam(std::string fileName)
{

	std::string title;
	char   charTmp[256];
	std::fstream file;

	std::cout << fileName << std::endl;
	file.open(fileName.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "In displaySummary, fileName = " << fileName << std::endl;
		ToolsThrow("The parameter file was not opened properly!");
	}

	// fill the attributes of the parameters structure
	//
	readUntilChar(&file, '=', &title);
	while ( file.good() ){
		if (!strcmp(title.c_str(), "tMan"))
			file >> pParam_->tMan;
		else if (!strcmp(title.c_str(), "tIni"))
			file >> pParam_->tIni;
		else if (!strcmp(title.c_str(), "tFin"))
			file >> pParam_->tFin;
		else if (!strcmp(title.c_str(), "coCordesV"))
			file >> pParam_->coCordesV;
		else if (!strcmp(title.c_str(), "dCV"))
			file >> pParam_->dCV;
		else if (!strcmp(title.c_str(), "coCordesA"))
			file >> pParam_->coCordesA;
		else if (!strcmp(title.c_str(), "coCostLines"))
			file >> pParam_->coCostLines;
		else if (!strcmp(title.c_str(), "fixedCostPercent"))
			file >> pParam_->fixedCostPercent;
		else
			ToolsThrow("The title of the line is unknown!");

		file.getline(charTmp, 256);
		if (!ToolsReadUntilChar(&file, '=', &title)) continue;
	}
	file.close();
}

// Set the flights by reading data from files
//
void CftSolver::readFlights(std::string fileName)
{

	std::string title;
	char   charTmp[256];
	std::fstream file;

	std::cout << fileName << std::endl;
	file.open(fileName.c_str(), std::fstream::in);
	if (!file.is_open()) {
		std::cout << "In displaySummary, fileName = " << fileName << std::endl;
		ToolsThrow("The flights file was not opened properly!");
	}

	// set the horizontal separation norm
	//
	ToolsReadUntilChar(&file, '=', &title);
	file >> pParam_->hNorm;

	// set the number of flights
	//
	file.getline(charTmp, 256);
	ToolsReadUntilChar(&file, '=', &title);
	file >> coFlights_;

	// initialize the flight structures
	//
	for (int i = 0; i < coFlights_ ; i++)	{
		//CftFlight* pFlight = new CftFlight();
		vFlights_.push_back(CftFlight());
	}

	// fill the attributes of the flights structures
	//
	file.getline(charTmp, 256);
	ToolsReadUntilChar(&file, '=', &title);
	while (strcmp(title.c_str(), "nConflits")){
		for (int i = 0; i < coFlights_ ; i++) {
			if (!strcmp(title.c_str(), "xIni"))
				file >> vFlights_[i].xIni;
			else if (!strcmp(title.c_str(), "yIni"))
				file >> vFlights_[i].yIni;
			else if (!strcmp(title.c_str(), "xFin"))
				file >> vFlights_[i].xFin;
			else if (!strcmp(title.c_str(), "yFin"))
				file >> vFlights_[i].yFin;
			else if (!strcmp(title.c_str(), "tIni"))
				file >> vFlights_[i].tIni;
			else if (!strcmp(title.c_str(), "tFin"))
				file >> vFlights_[i].tFin;
			else if (!strcmp(title.c_str(), "vxIni"))
				file >> vFlights_[i].vxIni;
			else if (!strcmp(title.c_str(), "vxFin"))
				file >> vFlights_[i].vxFin;
			else if (!strcmp(title.c_str(), "vyIni"))
				file >> vFlights_[i].vyIni;
			else if (!strcmp(title.c_str(), "vyFin"))
				file >> vFlights_[i].vyFin;
			else if (!strcmp(title.c_str(), "vyFin"))
				file >> vFlights_[i].vyFin;
			else if (!strcmp(title.c_str(), "vMin"))
				file >> vFlights_[i].vMin;
			else if (!strcmp(title.c_str(), "vMax"))
				file >> vFlights_[i].vMax;
			else if (!strcmp(title.c_str(), "aMax"))
				file >> vFlights_[i].aMax;
			else if (!strcmp(title.c_str(), "omegaMax"))
				file >> vFlights_[i].omegaMax;
			else if (!strcmp(title.c_str(), "costAlpha"))
				file >> vFlights_[i].alphaCost;
			else if (!strcmp(title.c_str(), "costBeta"))
				file >> vFlights_[i].betaCost;
			else if (!strcmp(title.c_str(), "costIndex"))
				file >> vFlights_[i].costIndex;
			else
				ToolsThrow("The title of the line is unknown!");
		}

		file.getline(charTmp, 256);
		ToolsReadUntilChar(&file, '=', &title);
	}

	// convert the speed values to get NM/min and compute speed norms
	for (int i = 0; i < coFlights_; i++)	{
		vFlights_[i].vxIni = vFlights_[i].vxIni/60.0;
		vFlights_[i].vyIni = vFlights_[i].vyIni/60.0;
		vFlights_[i].vxFin = vFlights_[i].vxFin/60.0;
		vFlights_[i].vyFin = vFlights_[i].vyFin/60.0;
		vFlights_[i].vMin = vFlights_[i].vMin/60.0;
		vFlights_[i].vMax = vFlights_[i].vMax/60.0;

		vFlights_[i].vIni = sqrt( pow(vFlights_[i].vxIni, 2.0)+ pow(vFlights_[i].vyIni, 2.0) );
		vFlights_[i].vFin = sqrt( pow(vFlights_[i].vxFin, 2.0)+ pow(vFlights_[i].vyFin, 2.0) );

		// test the validity of some parameters
		if (vFlights_[i].tFin < vFlights_[i].tIni)
			ToolsThrow("The final time is smaller than the initial one!");
	}

	// detect the flights with nonconstant altitudes and store them in a different
	// vector
	int coEvol = 0, coFl = 0;
	pFlightsInEvol_ = (int*) malloc(coFlights_*sizeof(int));
	pFlightsPhase_  = (int*) malloc((coFlights_*sizeof(int)));
	for (int i = 0; i < coFlights_; i++)  {
		pFlightsInEvol_[i] = 0;
		if (vFlights_[i].tIni > 0)  {
			pFlightsInEvol_[coEvol] = i;
			pFlightsPhase_[i] = coEvol;
			vEvol_.push_back(vFlights_[i]);
			coEvol_++; coEvol++;
		}
		else {
			pFlightsPhase_[i] = -coFl-1;
			coFl++;
		}
	}
	for (int i = 0; i < coEvol_; i++) {
		vFlights_.erase(vFlights_.begin()+pFlightsInEvol_[i]-i);
		coFlights_--;
	}

	// adjust the final and initial times of each flight in order to get multiples
	// of the maneuvers duration
	// the time interval is reduced by assuming the initial and final speed
	// vectors are respected at the beginning and the end
	for (int i = 0; i < coFlights_; i++)  {
		double tExcess = fmod(vFlights_[i].tFin, pParam_->tMan);
		vFlights_[i].tFin = vFlights_[i].tFin - tExcess;
		if (tExcess > 1.0e-6) {
			vFlights_[i].xFin = vFlights_[i].xFin - tExcess*vFlights_[i].vxFin;
			vFlights_[i].yFin = vFlights_[i].yFin - tExcess*vFlights_[i].vyFin;
		}
		pParam_->tFin = std::max(pParam_->tFin, vFlights_[i].tFin);

		tExcess = fmod(vFlights_[i].tIni, pParam_->tMan);
		double tSlack = (tExcess <= 1.0e-06) ? 0 : pParam_->tMan - tExcess;
		vFlights_[i].tIni = vFlights_[i].tIni + tSlack;
		if (tSlack > 0) {
			vFlights_[i].xIni = vFlights_[i].xIni + tSlack*vFlights_[i].vxIni;
			vFlights_[i].yIni = vFlights_[i].yIni + tSlack*vFlights_[i].vyIni;
		}
	}

	// in contrast, the time interval is enlarged for flights in evolution without
	// real problem since the speed vector must be constant
	for (int i = 0; i < coEvol_; i++)  {
		double tExcess, tSlack;

		// set the initial time to a multiple of tMan and ensure it is greater than
		// the global initial time of the instance
		if (vEvol_[i].tIni >= pParam_->tIni) {
			tExcess = fmod(vEvol_[i].tIni, pParam_->tMan);
			vEvol_[i].tIni = vEvol_[i].tIni - tExcess;
			if (tExcess > 1.0e-6) {
				vEvol_[i].xIni = vEvol_[i].xIni - tExcess*vEvol_[i].vxIni;
				vEvol_[i].yIni = vEvol_[i].yIni - tExcess*vEvol_[i].vyIni;
			}
		}
		else  {
			std::cout << "The flight in evolution number " << i << " has a too small initial time !" << std::endl;
			tSlack = pParam_->tIni-vEvol_[i].tIni;
			vEvol_[i].tIni = vEvol_[i].tIni + tSlack;
			if (tSlack > 0) {
				vEvol_[i].xIni = vEvol_[i].xIni + tSlack*vEvol_[i].vxIni;
				vEvol_[i].yIni = vEvol_[i].yIni + tSlack*vEvol_[i].vyIni;
			}
		}

		// set the final time to a multiple of tMan and ensure it is smaller than
		// the global initial time of the instance
		// no need to touch the final position it is not used in what follows
		// (because the speed of the flights in evolution is not modified)
		if (vEvol_[i].tFin <= pParam_->tFin) {
			tExcess = fmod(vEvol_[i].tFin, pParam_->tMan);
			tSlack = (tExcess <= 1.0e-06) ? 0 : pParam_->tMan - tExcess;

			vEvol_[i].tFin = vEvol_[i].tFin + tSlack;
		}
		else  {
			std::cout << "The flight in evolution number " << i << " has a too large final time !" << std::endl;
			tExcess = vEvol_[i].tFin-pParam_->tFin;
			vEvol_[i].tFin = vEvol_[i].tFin - tExcess;
		}
	}

	// fill the table of conflicts
	//
	int coCfts;
	file >> coCfts;
	std::cout << "Number of conflicts: " << coCfts << std::endl;
	if (coCfts) {
		pConflicts_ = (cft*) malloc(coCfts*sizeof(cft));
		pConflictsEvol_ = (cft*) malloc(coCfts*sizeof(cft));
	}
	else {
		file.getline(charTmp, 256);
		ToolsReadUntilChar(&file, '=', &title);
	}
	for (int i = 0; i < coCfts; i++)  {
		int fl1, fl2;
		file.getline(charTmp, 256);
		ToolsReadUntilChar(&file, '=', &title);
		file >> fl1 >> fl2;
		if (pFlightsPhase_[fl1] >= 0 && pFlightsPhase_[fl2] >= 0) {
			ToolsThrow("Two flights in conflict are in evolution!");
		}
		else if (pFlightsPhase_[fl1] >= 0)  {
			coConflictsEvol_++;
			pConflictsEvol_[i].fl1 = -pFlightsPhase_[fl2]-1;
			pConflictsEvol_[i].fl2 = pFlightsPhase_[fl1];
			std::cout << "Conflict in evolution" << coConflictsEvol_ << ": ";
			std::cout << pConflictsEvol_[i].fl1 << " vs " << pConflictsEvol_[i].fl2 << std::endl;
		}
		else if (pFlightsPhase_[fl2] >= 0)  {
			coConflictsEvol_++;
			pConflictsEvol_[i].fl1 = -pFlightsPhase_[fl1]-1;
			pConflictsEvol_[i].fl2 = pFlightsPhase_[fl2];
			std::cout << "Conflict in evolution" << coConflictsEvol_ << ": ";
			std::cout << pConflictsEvol_[i].fl1 << " vs " << pConflictsEvol_[i].fl2 << std::endl;
		}
		else {
			coConflicts_++;
			pConflicts_[i].fl1 = -pFlightsPhase_[fl1]-1;
			pConflicts_[i].fl2 = -pFlightsPhase_[fl2]-1;
			std::cout << "Conflict " << coConflictsEvol_ << ": ";
			std::cout << pConflicts_[i].fl1 << " vs " << pConflicts_[i].fl2 << std::endl;
		}
	}

	file.getline(charTmp, 256);
	ToolsReadUntilChar(&file, '=', &title);
	file >> coTrails_;
	std::cout << "Number of trails: " << coTrails_ << std::endl;
	if (coTrails_)
		pTrails_ = (cft*) malloc(coTrails_*sizeof(cft));
	for (int i = 0; i < coTrails_; i++)  {
		file.getline(charTmp, 256);
		ToolsReadUntilChar(&file, '=', &title);
		file >> pTrails_[i].fl1;
		file >> pTrails_[i].fl2;
		std::cout << "Trail " << i << ": ";
		std::cout << pTrails_[i].fl1 << " vs " << pTrails_[i].fl2 << std::endl;
	}

	file.getline(charTmp, 256);
	ToolsReadUntilChar(&file, '=', &title);
	file >> this->coParallel_;
	std::cout << "Number of parallel trajectories: " << coParallel_ << std::endl;
	if (coParallel_)
		pParallel_ = (cft*) malloc(coParallel_*sizeof(cft));
	for (int i = 0; i < coParallel_; i++)  {
		file.getline(charTmp, 256);
		ToolsReadUntilChar(&file, '=', &title);
		file >> pParallel_[i].fl1;
		file >> pParallel_[i].fl2;
		std::cout << "Parallel trajectories " << i << ": ";
		std::cout << pParallel_[i].fl1 << " vs " << pParallel_[i].fl2 << std::endl;
	}

	file.close();
}


//--------------------------------------------------------------------------
// Useful parsing functions
//

// Read a file stream until the separating character is met
// Store the characters read until the separating character in pStrRead
//
bool readUntilChar(std::fstream *pFile, char separater, std::string *pStrRead) {
	char cTmp = 'A';

	// empty the title string if it is not
	//
	if (!pStrRead->empty())
		pStrRead->erase();

	// go through the file until the delimiter is met
	//
	if (file->good()) {
		cTmp = file->get();
	}
	while (cTmp != separater && file->good() )  {
		pStrRead->push_back(cTmp);
		cTmp = file->get();
	}

	if (!file->good())
		return false;

	return true;
}

//--------------------------------------------------------------------------
