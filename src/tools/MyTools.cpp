//
//  MyTools.cxx
//  IDSReseau
//
//  Created by J������r������my Omer on 19/11/2013.
//  Copyright (c) 2013 J������r������my Omer. All rights reserved.
//

#include "tools/MyTools.h"
#include <stdlib.h>     /* rand */
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <stdexcept>
#include <math.h>       /* pow */

// necessary because OS X does not have clock_gettime, using clock_get_time
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

using namespace std;

namespace Tools {

// Compare functions to sort
  bool compareDecreasing(int i, int j) { return (i > j); }


// Throw an exception with the input message
//
  void throwException(const char *exceptionMsg) {
    printf("Exception caught: %s\n", exceptionMsg);
    throw NSException(exceptionMsg);
  }

  void throwException(const std::string& exceptionMsg) {
    throwException(exceptionMsg.c_str());
  }

  void throwError(const char *exceptionMsg) {
    printf("Error caught: %s\n", exceptionMsg);
    throw exceptionMsg;
  }

  void throwError(const std::string& exceptionMsg) {
    throwError(exceptionMsg.c_str());
  }

// Display a debug message
//
  void debugMsg(const char *debugMsg, int debugLevel) {

    if (debugLevel >= DEBUG)
      printf("%s\n", debugMsg);
  }

// Read a file stream until the separating character is met
//
  bool readUntilChar(fstream &file, char separateur, string &pTitle) {
    char cTmp = 'A';

    // empty the title string if it is not
    //
    if (!pTitle.empty())
      pTitle.erase();

    // go through the file until the delimiter is met
    //
    if (file.good()) {
      cTmp = file.get();
    }
    while (cTmp != separateur && file.good()) {
      pTitle.push_back(cTmp);
      cTmp = file.get();
    }

    if (!file.good())
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
    else {
      string endOfSentence = sentence.substr(lSentence - lWord, lWord);
      return (!strcmp(word.c_str(), endOfSentence.c_str()));
    }
  }

// Parse an int list written as string with a char delimiter
//
  std::vector<int> parseList(std::string strList, char delimiter) {

    std::vector<int> intList;
    std::stringstream ss(strList);
    int i;

    while (ss >> i) {
      intList.push_back(i);

      if (ss.peek() == delimiter) {
        ss.ignore();
      }
    }
    return intList;
  }

//random generator of tools
  std::minstd_rand rdm0(0);

// Initialize the random generator with a given seed
  void initializeRandomGenerator() {
    rdm0 = getANewRandomGenerator();
  }

  void initializeRandomGenerator(int rdmSeed) {
    rdm0 = getANewRandomGenerator(rdmSeed);
  }

//Create a random generator
//the objective is to be sure to have always the same sequence of number
//
  minstd_rand getANewRandomGenerator() {
    return getANewRandomGenerator(rdm0());
  }

  minstd_rand getANewRandomGenerator(int rdmSeed) {
    std::cout << "The new random seed of random generator is " << rdmSeed << std::endl;
    minstd_rand rdm(rdmSeed);
    return rdm;
  }

//round with probability
  int roundWithProbability(double number) {
    //round with a certain probability to the floor or the ceil
    double probFactor = number - floor(number);
    if (randomDouble(0, 1) < probFactor) return (int) floor(number);
    else return (int) ceil(number);
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
    return ((maxVal - minVal) * ((double) rdm0() / (double) RAND_MAX) + minVal);
  }

// Initializes a vector< double > of size m with random values (uniform) within [minVal, maxVal]
//
  std::vector<double> randomDoubleVector(int m, double minVal, double maxVal) {
    std::vector<double> v1D;
    for (int i = 0; i < m; i++) {
      double a = randomDouble(minVal, maxVal);
      v1D.push_back(a);
    }
    return v1D;
  }

// Returns a vector< vector< double > > of size m x n with random values (uniform) within [minVal, maxVal]
//
  std::vector<std::vector<double> > randomDoubleVector2D(int m, int n, double minVal, double maxVal) {
    std::vector<std::vector<double> > v2D;
    for (int i = 0; i < m; i++) {
      std::vector<double> v1D = randomDoubleVector(n, minVal, maxVal);
      v2D.push_back(v1D);
    }
    return v2D;
  }

// For an input vector of weights w, draw the index of one element assuming that
// the probability of each index is w_i/sum{w_i}
//
  int drawRandomWithWeights(std::vector<double> weights) {
    if (weights.empty()) {
      throwError("Tools::drawRandomWithWeights: the weight vector cannot be empty!");
    }

    // compute the sum of the weights
    double sumWeights = 0.0;
    for (double w:weights) {
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
    //return value
    std::vector<int> randIndVector;

    // number and vector of indices that have not been drawn
    int indWidth = indMax - indMin + 1;
    std::vector<int> indVector;
    for (int i = indMin; i <= indMax; i++) {
      indVector.push_back(i);
    }

    //return all the indices if they are less numerours than the number of
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
    else return "Sun";
  }

  int dayToInt(std::string day) {
    if (day == "Mon") return 0;
    if (day == "Tue") return 1;
    if (day == "Wed") return 2;
    if (day == "Thu") return 3;
    if (day == "Fri") return 4;
    if (day == "Sat") return 5;
    if (day == "Sun") return 6;
    else return -1;
  }

  bool isSaturday(int dayId) {
    if ((dayId % 7) == 5) return true;
    return false;
  }

  bool isSunday(int dayId) {
    if ((dayId % 7) == 6) return true;
    return false;
  }

  bool isWeekend(int dayId) {
    if (isSaturday(dayId) || isSunday(dayId)) return true;
    return false;
  }

  int containsWeekend(int start, int end) {
    int nbWeekend = 0;
    for (int i = start; i <= end; i++)
      if (isWeekend(i)) {
        ++nbWeekend;
        ++i;
      }
    return nbWeekend;
  }

// constructor of Timer
//
  Timer::Timer() : isInit_(0), isStarted_(0), isStopped_(0) {
    this->init();
  }

// initialize the timer
//
  void Timer::init() {
    coStop_ = 0;
    cpuInit_.tv_sec = 0;
    cpuInit_.tv_nsec = 0;
    cpuSinceStart_.tv_sec = 0;
    cpuSinceStart_.tv_nsec = 0;
    cpuSinceInit_.tv_sec = 0;
    cpuSinceInit_.tv_nsec = 0;
    isInit_ = 1;
    isStarted_ = 0;
    isStopped_ = 1;
  }

// start the timer
//
  void Timer::start() {

    if (isStarted_)
      throwError("Trying to start an already started timer!");

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    cpuInit_.tv_sec = mts.tv_sec;
    cpuInit_.tv_nsec = mts.tv_nsec;

#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuInit_);
#endif

    cpuSinceStart_.tv_sec = 0;
    cpuSinceStart_.tv_nsec = 0;
    isStarted_ = 1;
    isStopped_ = 0;

  }

// Stop the time and update the times spent since the last start and since the
// initialization
//
  void Timer::stop() {

    if (isStopped_)
      throwError("Trying to stop an already stopped timer!");

    timespec cpuNow;

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    cpuNow.tv_sec = mts.tv_sec;
    cpuNow.tv_nsec = mts.tv_nsec;

#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuNow);
#endif

    if (cpuNow.tv_nsec - cpuInit_.tv_nsec < 0) {
      cpuSinceStart_.tv_sec = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
      cpuSinceStart_.tv_nsec = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
    } else {
      cpuSinceStart_.tv_sec = cpuNow.tv_sec - cpuInit_.tv_sec;
      cpuSinceStart_.tv_nsec = cpuNow.tv_nsec - cpuInit_.tv_nsec;
    }

    if (cpuSinceStart_.tv_nsec + cpuSinceInit_.tv_nsec >= 1e09) {
      cpuSinceInit_.tv_sec = cpuSinceInit_.tv_sec + cpuSinceStart_.tv_sec + 1;
      cpuSinceInit_.tv_nsec = cpuSinceInit_.tv_nsec + cpuSinceStart_.tv_nsec - 1e09;
    } else {
      cpuSinceInit_.tv_sec = cpuSinceInit_.tv_sec + cpuSinceStart_.tv_sec;
      cpuSinceInit_.tv_nsec = cpuSinceInit_.tv_nsec + cpuSinceStart_.tv_nsec;
    }

    isStarted_ = 0;
    isStopped_ = 1;
    coStop_++;

  } // end stop

// Get the time spent since the initialization of the timer without stopping it
//
  const double Timer::dSinceInit() {

    if (isStarted_) {
      timespec cpuNow;

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
      clock_serv_t cclock;
      mach_timespec_t mts;
      host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
      clock_get_time(cclock, &mts);
      mach_port_deallocate(mach_task_self(), cclock);
      cpuNow.tv_sec = mts.tv_sec;
      cpuNow.tv_nsec = mts.tv_nsec;

#else
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuNow);
#endif

      timespec cpuTmp;
      if (cpuNow.tv_nsec - cpuInit_.tv_nsec < 0) {
        cpuTmp.tv_sec = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
        cpuTmp.tv_nsec = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
      } else {
        cpuTmp.tv_sec = cpuNow.tv_sec - cpuInit_.tv_sec;
        cpuTmp.tv_nsec = cpuNow.tv_nsec - cpuInit_.tv_nsec;
      }

      timespec cpuCurrent;
      if (cpuTmp.tv_nsec + cpuSinceInit_.tv_nsec >= 1e09) {
        cpuCurrent.tv_sec = cpuSinceInit_.tv_sec + cpuTmp.tv_sec + 1;
        cpuCurrent.tv_nsec = cpuSinceInit_.tv_nsec + cpuTmp.tv_nsec - 1e09;
      } else {
        cpuCurrent.tv_sec = cpuSinceInit_.tv_sec + cpuTmp.tv_sec;
        cpuCurrent.tv_nsec = cpuSinceInit_.tv_nsec + cpuTmp.tv_nsec;
      }

      return (double) cpuCurrent.tv_sec + (double) cpuCurrent.tv_nsec / 1.0e09;
    } else if (isStopped_) {
      return (double) cpuSinceInit_.tv_sec + (double) cpuSinceInit_.tv_nsec / 1.0e09;
    } else
      throwError("Trying to get the value of an unitialized timer!");

    return 0.0;

  } //end dSinceInit


// Get the time spent since the last start of the timer without stopping it
//
  const double Timer::dSinceStart() {

    if (!isStarted_ && !isStopped_)
      throwError("Trying to get the value of an unitialized timer!");
    else if (isStarted_) {
      timespec cpuNow;

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
      clock_serv_t cclock;
      mach_timespec_t mts;
      host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
      clock_get_time(cclock, &mts);
      mach_port_deallocate(mach_task_self(), cclock);
      cpuNow.tv_sec = mts.tv_sec;
      cpuNow.tv_nsec = mts.tv_nsec;

#else
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuNow);
#endif

      if (cpuNow.tv_nsec - cpuInit_.tv_nsec < 0) {
        cpuSinceStart_.tv_sec = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
        cpuSinceStart_.tv_nsec = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
      } else {
        cpuSinceStart_.tv_sec = cpuNow.tv_sec - cpuInit_.tv_sec;
        cpuSinceStart_.tv_nsec = cpuNow.tv_nsec - cpuInit_.tv_nsec;
      }
    }

    return (double) cpuSinceStart_.tv_sec + (double) cpuSinceStart_.tv_nsec / 1.0e09;

  } //end dSinceStart


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
    unique_lock<mutex> lock(mThreadMutex_);

    // set to true the flag settingMaxGlobalThreads_
    if(settingMaxGlobalThreads_) {// if already setting the max number of threads, print a warning and return
      std::cerr << "WARNING: another thread is already trying to change the max number of available threads" << std::endl;
      return maxGlobalThreads_;
    }
    settingMaxGlobalThreads_ = true;

    // check max number of cores on the system
    int sysThreads = thread::hardware_concurrency();
    if(maxThreads > sysThreads) {
      std::cerr << "WARNING: the system has " << sysThreads
                << " cores, and we are trying to set the maximum threads to " << maxThreads
                << ": the maximum is automatically blocked at " << sysThreads << " threads." << std::endl;
      maxThreads = sysThreads;
    }

    // try to update the max
    int diff = maxThreads - maxGlobalThreads_;
    while(nGlobalThreadsAvailable_ + diff < 0) { // ensure that the maximum can be decreased
      if (wait) { // wait to be able to decrease totally the number of global threads
        cThreadReleased_.wait(lock); // wait until notification
      } else { // just decrease what is currently available
        diff = -nGlobalThreadsAvailable_;
        std::cerr << "WARNING: the maximum global number of threads cannot be decreased to " << maxThreads
                  << " as " << (maxGlobalThreads_ - nGlobalThreadsAvailable_) << " of them are currently used."
                  << "The max has been instead set to " << (maxGlobalThreads_ + diff) << "." << std::endl;
      }
    }

    // update the number of threads
    maxGlobalThreads_ += diff;
    nGlobalThreadsAvailable_ += diff;
    settingMaxGlobalThreads_ = false;

    return maxGlobalThreads_;
  }

  ThreadsPool::ThreadsPool(): maxThreads_(maxGlobalThreads_), nThreadsAvailable_(maxThreads_) {}

  ThreadsPool::ThreadsPool(int nThreads) {
    if(nThreads <= maxGlobalThreads_) {
      maxThreads_ = nThreads;
    } else {
      maxThreads_ = maxGlobalThreads_;
      std::cerr << "WARNING: the local pool cannot use " << nThreads
                << " threads, instead it uses the maximum allowed: " << maxGlobalThreads_ << std::endl;
    }
    nThreadsAvailable_= maxThreads_;
  }

  void ThreadsPool::run(Job job) {
    if(maxThreads_ <= 1) // no concurrency
      job();
    else {
      reserve(); // reserve a thread
      thread t([this, job]() {
          job(); // run the function
          release(); // release the thread
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
    unique_lock<mutex> lock(mThreadMutex_); // create a lock on the mutex
    // all threads of the local pool have ended
    if (nThreadsAvailable_ == maxThreads_)
      return false;
    // wait for all threads to finish
    int threadsAvailableTarget = std::min(maxThreads_, nThreadsAvailable_ + nThreadsToWait);
    cThreadReleased_.wait(lock, [&]() {
        return !settingMaxGlobalThreads_ && nThreadsAvailable_ >= threadsAvailableTarget;
    });
    return true;
  }

  void ThreadsPool::reserve() {
    unique_lock<mutex> lock(mThreadMutex_); // create a lock on the mutex
    // reserve ont thread if available, otherwise wait for one
    bool reserved = false;
    while (!reserved) {
      // wait until one thread is available and settingMaxGlobalThreads_ is false
      cThreadReleased_.wait(lock, [this]() {
        return nGlobalThreadsAvailable_ > 0 && nThreadsAvailable_ > 0 && !settingMaxGlobalThreads_;
      });
      // reserved one thread
      --nGlobalThreadsAvailable_; // reserve one global thread
      --nThreadsAvailable_; // reserve one local thread
      reserved = true;
    }
  }

  void ThreadsPool::release() {
    lock_guard<mutex> lock(mThreadMutex_); // wait until can lock the mutex
    ++nGlobalThreadsAvailable_; // release a thread for the global number of threads
    ++nThreadsAvailable_; // release a thread for the local number of threads
    cThreadReleased_.notify_all(); // notify any thread that were waiting
  }

} // end namespace
