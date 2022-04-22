/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */


#include "ParseArguments.h"

#include <string>
#include <vector>

#include "tools/InputPaths.h"
#include "tools/Tools.h"

using std::string;

/******************************************************************************
* Read the arguments in non compact format
******************************************************************************/
InputPaths *readNonCompactArguments(int argc, char **argv) {
  InputPaths *pInputPaths = new InputPaths();

  // Read the arguments and store them in inputPaths
  //
  int narg = 1;
  while (narg < argc) {
    std::cout << "arg = " << argv[narg] << " " << argv[narg + 1] << std::endl;
    // problem with the java simulator that add quote
    // marks in the arguments, which messes with the open file methods
    // the code below is here to remove these quote marks
    //
    string str(argv[narg + 1]);
    std::size_t found = str.find("\"");
    while (found != std::string::npos) {
      str.erase(found, 1);
      found = str.find("\"");
    }

    if (!strcmp(argv[narg], "--sce")) {
      pInputPaths->scenario(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--his")) {
      pInputPaths->history(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--week")) {
      pInputPaths->addWeek(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--sol")) {
      pInputPaths->solutionPath(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--param")) {
      pInputPaths->paramFile(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--timeout")) {
      pInputPaths->timeOut(std::stoi(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--cusIn")) {
      pInputPaths->customInputFile(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--cusOut")) {
      pInputPaths->customOutputFile(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--verbose")) {
      pInputPaths->verbose(std::stoi(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--rand")) {
      pInputPaths->randSeed(std::stoi(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--sp-type")) {
      pInputPaths->SPType(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--n-threads")) {
      pInputPaths->nThreads(std::stoi(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--sp-strategy")) {
      pInputPaths->SPStrategy(std::stoi(str));
      narg += 2;
    } else if (!strcmp(argv[narg], "--rcspp-type")) {
      pInputPaths->RCSPPType(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--n-candidates")) {
      pInputPaths->nCandidates(std::stoi(str));
      narg += 2;
    } else {
      Tools::throwError(
          "main: the argument (%s) does not match the expected list!",
          argv[narg]);
    }
  }
  // Throw an error if a necessary input file is missing
  if (pInputPaths->scenario().empty() || pInputPaths->history().empty()
      || pInputPaths->weeks().empty()) {
    throw Tools::myException(
        "readNonCompactArguments: A necessary file name is missing!",
        __LINE__);
  }

  // Non compact format is only for release versions, so no log file is required
  pInputPaths->logPath("");

  return pInputPaths;
}

/******************************************************************************
* Read the arguments in compact format
******************************************************************************/

InputPaths *readCompactArguments(int argc, char **argv) {
  // Default arguments are set to enable simple call to the function without
  // argument.
  std::string dataDir, instanceName, solutionPath, logPath, paramFile,
      SPType, RCSPPType;
  int historyIndex = -1, randSeed = 0, nTreads = -1,
      SPStrategy = -1, verbose = -1, nCandidates = -1, timeOut = -1;
  std::vector<int> weekIndices;
  bool cyclic = false;

  // Read the arguments and store them in inputPaths
  //
  int narg = 1;
  while (narg < argc) {
    const char *arg = argv[narg];
    std::cout << "arg = " << arg << " " << argv[narg + 1] << std::endl;
    // WARNING: problem with the java simulator that add quote
    // marks in the arguments, which messes with the open file methods
    // the code below is here to remove these quote marks
    //
    string str(argv[narg + 1]);
    std::size_t found = str.find("\"");
    while (found != std::string::npos) {
      str.erase(found, 1);
      found = str.find("\"");
    }

    if (!strcmp(arg, "--dir")) {
      dataDir = str;
      narg += 2;
    } else if (!strcmp(arg, "--instance")) {
      instanceName = str;
      narg += 2;
    } else if (!strcmp(arg, "--his")) {
      historyIndex = std::stoi(str);
      narg += 2;
    } else if (!strcmp(arg, "--weeks")) {
      weekIndices = Tools::tokenize<int>(str, '-');
      narg += 2;
    } else if (!strcmp(arg, "--sol")) {
      solutionPath = str;
      narg += 2;
    } else if (!strcmp(arg, "--log")) {
      logPath = str;
      narg += 2;
    } else if (!strcmp(arg, "--param")) {
      paramFile = str;
      narg += 2;
    } else if (!strcmp(arg, "--verbose")) {
      verbose = std::stoi(str);
      narg += 2;
    } else if (!strcmp(arg, "--timeout")) {
      timeOut = std::stoi(str);
      narg += 2;
    } else if (!strcmp(arg, "--rand")) {
      randSeed = std::stoi(str);
      narg += 2;
    } else if (!strcmp(arg, "--sp-type")) {
      SPType = str;
      narg += 2;
    } else if (!strcmp(argv[narg], "--n-threads")) {
      nTreads = std::stoi(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--sp-strategy")) {
      SPStrategy = std::stoi(str);
      narg += 2;
    } else if (!strcmp(argv[narg], "--rcspp-type")) {
      RCSPPType = str;
      narg += 2;
    } else if (!strcmp(argv[narg], "--cyclic")) {
      cyclic = true;
      narg += 1;
    } else if (!strcmp(argv[narg], "--n-candidates")) {
      nCandidates = std::stoi(str);
      narg += 2;
    } else {
      Tools::throwError(
          "main: the argument (%s) does not match the expected list!",
          argv[narg]);
    }
  }

  if (cyclic && historyIndex >= 0)
    Tools::throwError("If cyclic option is enable, you can't use "
                      "an historical state for the nurses.");

  // Initialize the input paths
  //
  // if cyclic, do not enter any history
  if (cyclic)
    return new InputPaths(dataDir,
                          instanceName,
                          weekIndices,
                          solutionPath,
                          logPath,
                          paramFile,
                          timeOut,
                          verbose,
                          randSeed,
                          SPType,
                          SPStrategy,
                          RCSPPType,
                          nTreads,
                          nCandidates);

  return new InputPaths(dataDir,
                        instanceName,
                        historyIndex,
                        weekIndices,
                        solutionPath,
                        logPath,
                        paramFile,
                        timeOut,
                        verbose,
                        randSeed,
                        SPType,
                        SPStrategy,
                        RCSPPType,
                        nTreads,
                        nCandidates);
}

/******************************************************************************
* Read the arguments with the right method for the format
******************************************************************************/
InputPaths *readArguments(int argc, char **argv) {
  // check if the arg --dir is present
  int narg = 1;
  while (narg < argc && strcmp(argv[narg], "--dir"))
    narg++;
  if (narg == argc)
    return readNonCompactArguments(argc, argv);
  return readCompactArguments(argc, argv);
}
