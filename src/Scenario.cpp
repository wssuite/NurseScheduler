//
//  Scenario.cpp
//  Project:RosterDesNurses
//
//  Created by Jérémy Omer on 18/12/2013.
//  Copyright (c) 2014 Jérémy Omer. All rights reserved.
//


#include "Scenario.h"
#include "MyTools.h"
#include "Nurse.h"
#include "Solution.h"

#include <iostream>
#include <streambuf>
#include <fstream>
#include <math.h>
#include <time.h>

// Read the scneario file and store the content in a Scenario instance
//
// void readScenario(std::string strWeekFile, Scenario* pScenario) {
//
//
//         // open the file
//         std::fstream file;
//         std::cout << "Reading " << strWeekFile << std::endl;
//         file.open(strWeekFile.c_str(), std::fstream::in);
//         if (!file.is_open()) {
//                 std::cout << "While trying to read " << strWeekFile << std::endl;
//                 ToolsThrow("The input file was not opened properly!");
//         }
//
//         std::string title;
//         char charTmp[256];
//         int tmp;
//
//         // fill the attributes of the scenario structure
//         //
//         readUntilChar(&file, '=', &title);
//         while ( file.good() ) {
//             if (!strcmp(title.c_str(), "SCENARIO")) {
//                         file >>charTmp;
//                         pScenario->name(charTmp);
//                 }
//                 else if (!strcmp(title.c_str(), "WEEKS")) {
//                         file >> intTmp;
//                         pScenario->nbWeeks(intTmp);
//                 }
//                 else if (!strcmp(title.c_str(), "SKILLS"))   {
//                         file >> intTmp;
//                         pScenario->nbSkills(intTmp);
//                         for (int i = 0; i < intTmp; i++) {
//
//                         }
//                 }
//
//                 file >> pParam_->tIni;
//                 else if (!strcmp(title.c_str(), "tFin"))
//                         file >> pParam_->tFin;
//                 else if (!strcmp(title.c_str(), "coCordesV"))
//                         file >> pParam_->coCordesV;
//                 else if (!strcmp(title.c_str(), "dCV"))
//                         file >> pParam_->dCV;
//                 else if (!strcmp(title.c_str(), "coCordesA"))
//                         file >> pParam_->coCordesA;
//                 else if (!strcmp(title.c_str(), "coCostLines"))
//                         file >> pParam_->coCostLines;
//                 else if (!strcmp(title.c_str(), "fixedCostPercent"))
//                         file >> pParam_->fixedCostPercent;
//                 else
//                         ToolsThrow("The title of the line is unknown!");
//
//                 file.getline(charTmp, 256);
//                 if (!ToolsReadUntilChar(&file, '=', &title)) continue;
//         }
// }
// }
