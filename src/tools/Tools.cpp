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

#include "Tools.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <stdexcept>

// necessary because OS X does not have clock_gettime, using clock_get_time
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

using std::minstd_rand;
using std::mutex;
using std::thread;

namespace Tools {

// Compare functions to sort
bool compareDecreasing(int i, int j) { return (i > j); }

// Throw an exception with the input message
//
void throwException(const char *exceptionMsg) {
  fprintf(stderr, "Exception caught: %s\n", exceptionMsg);
  throw NSException(exceptionMsg);
}

void throwException(const std::string &exceptionMsg) {
  throwException(exceptionMsg.c_str());
}

void throwError(const char *exceptionMsg) {
  fprintf(stderr, "Error caught: %s\n", exceptionMsg);
  throw exceptionMsg;
}

void throwError(const std::string &exceptionMsg) {
  throwError(exceptionMsg.c_str());
}

// Display a debug message
//
void debugMsg(const char *debugMsg, int debugLevel) {
  if (debugLevel >= DEBUG)
    printf("%s\n", debugMsg);
}

// Store the characters read until the separating character in pStrRead
//
bool readUntilOneOfTwoChar(std::fstream *pFile,
                           char separator1,
                           char separator2,
                           std::string *pStrRead) {
  char cTmp = 'A';

  // empty the title string if it is not
  //
  if (!pStrRead->empty())
    pStrRead->erase();

  // go through the file until the delimiter is met
  //
  if (pFile->good()) {
    cTmp = pFile->get();
  }
  while (cTmp != separator1 && cTmp != separator2 && pFile->good()) {
    pStrRead->push_back(cTmp);
    cTmp = pFile->get();
  }

  return pFile->good();
}

// Read a file stream until the separating character is met
//
bool readUntilChar(std::fstream *file, char separator, std::string *pTitle) {
  char cTmp = 'A';

  // empty the title string if it is not
  //
  if (!pTitle->empty())
    pTitle->erase();

  // go through the file until the delimiter is met
  //
  if (file->good()) {
    cTmp = file->get();
  }
  while (cTmp != separator && file->good()) {
    pTitle->push_back(cTmp);
    cTmp = file->get();
  }

  if (!file->good())
    return false;

  return true;
}

// Checks if the string (sentence) ends with the given substring (word)
//
bool strEndsWith(std::string sentence, std::string word) {
  int lWord = word.length();
  int lSentence = sentence.length();
  if (lWord > lSentence)
    return false;

  std::string endOfSentence = sentence.substr(lSentence - lWord, lWord);
  return (!strcmp(word.c_str(), endOfSentence.c_str()));
}

// random generator of tools
std::minstd_rand rdm0(0);

// Initialize the random generator with a given seed
void initializeRandomGenerator() {
  rdm0 = getANewRandomGenerator();
}

void initializeRandomGenerator(int rdmSeed) {
  rdm0 = getANewRandomGenerator(rdmSeed);
}

// Create a random generator
// the objective is to be sure to have always the same sequence of number
//
minstd_rand getANewRandomGenerator() {
  return getANewRandomGenerator(rdm0());
}

minstd_rand getANewRandomGenerator(int rdmSeed) {
  std::cout << "The new random seed of random generator is " << rdmSeed
            << std::endl;
  minstd_rand rdm(rdmSeed);
  return rdm;
}

// round with probability
int roundWithProbability(double number) {
  // round with a certain probability to the floor or the ceil
  double probFactor = number - floor(number);
  if (randomDouble(0, 1) < probFactor)
    return static_cast<int>(floor(number));
  else
    return static_cast<int>(ceil(number));
}

// Returns an integer with random value (uniform) within [minVal, maxVal]
//
int randomInt(int minVal, int maxVal) {
  if (minVal > maxVal) {
    throwError("Tools::randomInt: minVal must be smaller than maxVal");
  }
  return minVal + rdm0() % (maxVal - minVal + 1);
}

// Returns a double with random value (uniform) within [minVal, maxVal]
//
double randomDouble(double minVal, double maxVal) {
  return ((maxVal - minVal) * (rdm0() * 1.0 / (rdm0.max()))  + minVal);
  // RAND_MAX) + minVal);
}

// Initializes a vector< double > of size m
// with random values (uniform) within [minVal, maxVal]
//
std::vector<double> randomDoubleVector(int m, double minVal, double maxVal) {
  std::vector<double> v1D;
  for (int i = 0; i < m; i++) {
    double a = randomDouble(minVal, maxVal);
    v1D.push_back(a);
  }
  return v1D;
}

// Returns a vector< vector< double > > of size m x n
// with random values (uniform) within [minVal, maxVal]
//
std::vector<std::vector<double> > randomDoubleVector2D(int m,
                                                       int n,
                                                       double minVal,
                                                       double maxVal) {
  std::vector<std::vector<double> > v2D;
  for (int i = 0; i < m; i++) {
    std::vector<double> v1D = randomDoubleVector(n, minVal, maxVal);
    v2D.push_back(v1D);
  }
  return v2D;
}

// For an input vector of weights w, draw the index of one element assuming
// that the probability of each index is w_i/sum{w_i}
//
int drawRandomWithWeights(std::vector<double> weights) {
  if (weights.empty()) {
    throwError(
        "Tools::drawRandomWithWeights: the weight vector cannot be empty!");
  }

  // compute the sum of the weights
  double sumWeights = 0.0;
  for (double w : weights) {
    sumWeights += w;
  }

  // draw the index
  double randNumber = randomDouble(0, sumWeights);
  int index = 0;
  double partialSum = weights[0];
  for (unsigned int i = 0; i < weights.size() - 1; i++) {
    if (randNumber <= partialSum) {
      return index;
    } else {
      partialSum += weights[i + 1];
      index++;
    }
  }
  return index;
}

// Draw randomly nbInd different indices between indMin and indMax
//
std::vector<int> drawRandomIndices(int nbInd, int indMin, int indMax) {
  // return value
  std::vector<int> randIndVector;

  // number and vector of indices that have not been drawn
  int indWidth = indMax - indMin + 1;
  std::vector<int> indVector;
  for (int i = indMin; i <= indMax; i++) {
    indVector.push_back(i);
  }

  // return all the indices if they are less numerours than the number of
  // required indices
  if (nbInd >= indWidth) {
    randIndVector = indVector;
  } else {
    // draw sequentially the indices
    for (int i = 0; i < nbInd; i++) {
      int randint = rdm0() % indWidth;
      randIndVector.push_back(indVector[randint]);
      indVector.erase(indVector.begin() + randint);
      indWidth--;
    }
  }
  return randIndVector;
}

// To get the day from its id and vice-versa
// First day is always supposed to be a Monday
//
std::string intToDay(int dayId) {
  if ((dayId % 7) == 0) return "Mon";
  else if ((dayId % 7) == 1) return "Tue";
  else if ((dayId % 7) == 2) return "Wed";
  else if ((dayId % 7) == 3) return "Thu";
  else if ((dayId % 7) == 4) return "Fri";
  else if ((dayId % 7) == 5) return "Sat";
  else
    return "Sun";
}

int dayToInt(std::string day) {
  if (day == "Mon") return 0;
  if (day == "Tue") return 1;
  if (day == "Wed") return 2;
  if (day == "Thu") return 3;
  if (day == "Fri") return 4;
  if (day == "Sat") return 5;
  if (day == "Sun") return 6;
  else
    return -1;
}

bool isSaturday(int dayId) {
  return (dayId % 7) == 5;
}

bool isSunday(int dayId) {
  return (dayId % 7) == 6;
}

bool isWeekend(int dayId) {
  return isSaturday(dayId) || isSunday(dayId);
}

bool isFirstWeekDay(int dayId) {
  return (dayId % 7 == 0);
}

bool isFirstWeekendDay(int dayId) {
  return isSaturday(dayId);
}

bool isLastWeekendDay(int dayId) {
  return isSunday(dayId);
}

bool isWeekendDayButNotLastOne(int dayId) {
  return isSaturday(dayId);
}

int nWeekendsInInterval(int start, int end) {
  int nWeekends = 0;
  bool isCounted = false;
  for (int i = start; i <= end; i++) {
    if (isWeekend(i) && !isCounted) {
      ++nWeekends;
      isCounted = true;
    } else {
      isCounted = false;
    }
  }
  return nWeekends;
}

// constructor of Timer
//
Timer::Timer(bool start) : cpuInit_(std::chrono::system_clock::now()),
                           cpuSinceStart_(0),
                           cpuSinceInit_(0),
                           coStop_(0),
                           isStarted_(false),
                           isStopped_(true) {
  if (start) this->start();
}

// start the timer
//
void Timer::start() {
  if (isStarted_)
    throwError("Trying to start an already started timer!");

  cpuInit_ = std::chrono::system_clock::now();
  isStarted_ = true;
  isStopped_ = false;
}

// Stop the time and update the times spent since the last start and since the
// initialization
//
void Timer::stop() {
  if (isStopped_)
    throwError("Trying to stop an already stopped timer!");

  cpuSinceStart_ = std::chrono::system_clock::now() - cpuInit_;
  cpuSinceInit_ += cpuSinceStart_;

  isStarted_ = false;
  isStopped_ = true;
  coStop_++;
}  // end stop

// Get the time spent since the initialization of the timer without stopping it
//
double Timer::dSinceInit() const {
  std::chrono::duration<double> cpuCurrent = cpuSinceInit_;
  if (isStarted_)
    cpuCurrent += std::chrono::system_clock::now() - cpuInit_;
  return getSeconds(cpuCurrent);
}  // end dSinceInit


// Get the time spent since the last start of the timer without stopping it
//
double Timer::dSinceStart() {
  if (isStarted_) {
    cpuSinceStart_ = std::chrono::system_clock::now() - cpuInit_;
  }
  return getSeconds(cpuSinceStart_);
}  // end dSinceStart

double Timer::getSeconds(const std::chrono::duration<double>& d) const {
  return d.count();
//  return std::chrono::duration_cast<std::chrono::seconds>(d).count();
}

//
// ThreadsPool implementation
//

int ThreadsPool::maxGlobalThreads_ = thread::hardware_concurrency();
int ThreadsPool::nGlobalThreadsAvailable_ = ThreadsPool::maxGlobalThreads_;
bool ThreadsPool::settingMaxGlobalThreads_ = false;
std::mutex ThreadsPool::mThreadMutex_;
std::condition_variable ThreadsPool::cThreadReleased_;

int ThreadsPool::getNGlobalThreadsAvailable() {
  return nGlobalThreadsAvailable_;
}

int ThreadsPool::getMaxGlobalThreads() {
  return maxGlobalThreads_;
}

int ThreadsPool::setMaxGlobalThreads(int maxThreads, bool wait) {
  std::unique_lock<mutex> lock(mThreadMutex_);

  // set to true the flag settingMaxGlobalThreads_
  // if already setting the max number of threads, print a warning and return
  if (settingMaxGlobalThreads_) {
    std::cerr
        << "WARNING: another thread is already trying to "
           "change the max number of available threads"
        << std::endl;
    return maxGlobalThreads_;
  }
  settingMaxGlobalThreads_ = true;

  // check max number of cores on the system
  int sysThreads = thread::hardware_concurrency();
  if (maxThreads > sysThreads) {
    std::cerr << "WARNING: the system has " << sysThreads
              << " cores, and we are trying to set the maximum threads to "
              << maxThreads
              << ": the maximum is automatically blocked at " << sysThreads
              << " threads." << std::endl;
    maxThreads = sysThreads;
  }

  // try to update the max
  int diff = maxThreads - maxGlobalThreads_;
  // ensure that the maximum can be decreased
  while (nGlobalThreadsAvailable_ + diff < 0) {
    // wait to be able to decrease totally the number of global threads
    if (wait) {
      cThreadReleased_.wait(lock);  // wait until notification
    } else {  // just decrease what is currently available
      diff = -nGlobalThreadsAvailable_;
      std::cerr
          << "WARNING: the maximum global number of threads "
             "cannot be decreased to "
          << maxThreads
          << " as " << (maxGlobalThreads_ - nGlobalThreadsAvailable_)
          << " of them are currently used."
          << "The max has been instead set to " << (maxGlobalThreads_ + diff)
          << "." << std::endl;
    }
  }

  // update the number of threads
  maxGlobalThreads_ += diff;
  nGlobalThreadsAvailable_ += diff;
  settingMaxGlobalThreads_ = false;

  return maxGlobalThreads_;
}

ThreadsPool::ThreadsPool()
    : maxThreads_(maxGlobalThreads_), nThreadsAvailable_(maxThreads_) {}

ThreadsPool::ThreadsPool(int nThreads) {
  if (nThreads <= maxGlobalThreads_) {
    maxThreads_ = nThreads;
  } else {
    maxThreads_ = maxGlobalThreads_;
    std::cerr << "WARNING: the local pool cannot use " << nThreads
              << " threads, instead it uses the maximum allowed: "
              << maxGlobalThreads_ << std::endl;
  }
  nThreadsAvailable_ = maxThreads_;
}

void ThreadsPool::run(Job job) {
  if (maxThreads_ <= 1) {  // no concurrency
    job();
  } else {
    reserve();  // reserve a thread
    thread t([this, job]() {
      job();  // run the function
      release();  // release the thread
    });
    t.detach();
  }
}

// wait until all threads of the pool are available.
// returns true if some threads were still running
bool ThreadsPool::wait() {
  return wait(maxThreads_);
}

bool ThreadsPool::wait(int nThreadsToWait) {
  std::unique_lock<mutex> lock(mThreadMutex_);  // create a lock on the mutex
  // all threads of the local pool have ended
  if (nThreadsAvailable_ == maxThreads_)
    return false;
  // wait for all threads to finish
  int threadsAvailableTarget =
      std::min(maxThreads_, nThreadsAvailable_ + nThreadsToWait);
  cThreadReleased_.wait(lock, [&]() {
    return !settingMaxGlobalThreads_
        && nThreadsAvailable_ >= threadsAvailableTarget;
  });
  return true;
}

void ThreadsPool::reserve() {
  std::unique_lock<mutex> lock(mThreadMutex_);  // create a lock on the mutex
  // reserve ont thread if available, otherwise wait for one
  bool reserved = false;
  while (!reserved) {
    // wait until one thread is available and settingMaxGlobalThreads_ is false
    cThreadReleased_.wait(lock, [this]() {
      return nGlobalThreadsAvailable_ > 0 && nThreadsAvailable_ > 0
          && !settingMaxGlobalThreads_;
    });
    // reserved one thread
    --nGlobalThreadsAvailable_;  // reserve one global thread
    --nThreadsAvailable_;  // reserve one local thread
    reserved = true;
  }
}

void ThreadsPool::release() {
  std::lock_guard<mutex> lock(mThreadMutex_);  // wait until can lock the mutex
  // release one thread
  ++nGlobalThreadsAvailable_;
  ++nThreadsAvailable_;
  cThreadReleased_.notify_all();  // notify any thread that were waiting
}

}  // namespace Tools
