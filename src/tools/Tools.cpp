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

#include <chrono>  // NOLINT
#include <ctime>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <streambuf>
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

const char COMMENT_KEY = '#';

// Compare functions to sort
bool compareDecreasing(int i, int j) { return (i > j); }

// Throw an exception with the input message
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

//  struct rusage {
//    struct timeval ru_utime; /* user time used */
//    struct timeval ru_stime; /* system time used */
//
//    long   ru_maxrss;        /* maximum resident set size */
//    The maximum resident set size used, in kilobytes.
//    That is, the maximum number of kilobytes of physical memory
//    that processes used simultaneously.
//
//    long   ru_ixrss;         /* integral shared memory size */
//    long   ru_idrss;         /* integral unshared data size */
//    long   ru_isrss;         /* integral unshared stack size */
//    long   ru_minflt;        /* page reclaims */
//    long   ru_majflt;        /* page faults */
//    long   ru_nswap;         /* swaps */
//    long   ru_inblock;       /* block input operations */
//    long   ru_oublock;       /* block output operations */
//    long   ru_msgsnd;        /* messages sent */
//    long   ru_msgrcv;        /* messages received */
//    long   ru_nsignals;      /* signals received */
//    long   ru_nvcsw;         /* voluntary context switches */
//    long   ru_nivcsw;        /* involuntary context switches */
//  };

rusage getRUsage() {
  struct rusage use;
  getrusage(RUSAGE_SELF, &use);
  return use;
}

float getResidentMemoryGB() {
  return getRUsage().ru_maxrss * 10e-10;
}

tm *dateForDay(const tm *const startDate, const int &dayId) {
  tm* newTime = new tm;
  memcpy(newTime, startDate, sizeof(tm));

  newTime->tm_mday += dayId;
  time_t nt_seconds = mktime(newTime) - timezone;
  delete newTime;

  return gmtime(&nt_seconds);
}

// Store the characters read until the separating character in pStrRead
//
char readUntilOneOfTwoChar(std::fstream *pFile,
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

  return cTmp;
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

bool readUntilAndWhileChar(std::fstream *file, char separator, std::string
*pTitle) {
  char cTmp = 'A';

  // empty the title string if it is not
  //
  if (!pTitle->empty())
    pTitle->erase();

  // go through the file until the delimiter is met
  //
  if (file->good()) {
    cTmp = file->get();
    pTitle->push_back(cTmp);
  }
  while (cTmp != separator && file->good()) {
    cTmp = file->get();
    pTitle->push_back(cTmp);
  }
  // keep going through the file while the delimiter goes on
  while (cTmp == separator && file->good()) {
    cTmp = file->get();
  }

  if (!file->good())
    return false;

  return true;
}

tm * readDateFromStr(const std::string& dateStr) {
  std::stringstream ss(dateStr);
  char charTmp;
  int yyyy, mm, dd;
  ss >> yyyy;
  ss >> charTmp;
  ss >> mm;
  ss >> charTmp;
  ss >> dd;
  std::tm time_in = { 0, 0, 0,  // second, minute, hour
                      dd, mm-1, yyyy - 1900 };  // 1-based day, 0-based month,
  // year since 1900
  std::time_t time_temp = std::mktime(&time_in);
  // Note: Return value of localtime is not threadsafe, because it might be
  // (and will be) reused in subsequent calls to std::localtime!
  return std::localtime(&time_temp);
}

Tools::Time readHourFromStr(const std::string& hourStr) {
  std::stringstream sstream(hourStr);
  char charTmp;
  int hh, mm, ss;
  sstream >> hh;
  sstream >> charTmp;
  sstream >> mm;
  sstream >> charTmp;
  sstream >> ss;
  return Time(hh, mm);
}

int readBoundedResource(std::fstream *file,
                        int* lbOn, int* lb, int* lbCost,
                        int* ubOn, int* ub, int* ubCost) {
  std::string strTmp;
  Tools::readUntilChar(file, '(', &strTmp);
  *file >> *ubOn;
  Tools::readUntilChar(file, '|', &strTmp);
  *file >> *ubCost;
  Tools::readUntilChar(file, '|', &strTmp);
  *file >> *ub;
  Tools::readUntilChar(file, '(', &strTmp);
  *file >> *lbOn;
  Tools::readUntilChar(file, '|', &strTmp);
  *file >> *lbCost;
  Tools::readUntilChar(file, '|', &strTmp);
  *file >> *lb;
  return 0;
}

int readUbResource(std::fstream *file,
                   int* ubOn, int* ub, int* ubCost) {
  std::string strTmp;
  Tools::readUntilChar(file, '(', &strTmp);
  *file >> *ubOn;
  Tools::readUntilChar(file, '|', &strTmp);
  *file >> *ubCost;
  Tools::readUntilChar(file, '|', &strTmp);
  *file >> *ub;
  return 0;
}

// Checks if the string (sentence) starts with the given substring (word)
//
bool strStartsWith(std::string sentence, std::string word) {
  int lWord = word.length();
  int lSentence = sentence.length();
  if (lWord > lSentence)
    return false;

  std::string startOfSentence = sentence.substr(0, lWord);
  return (!strcmp(word.c_str(), startOfSentence.c_str()));
}

// Checks if the string (sentence) starts with the comment character
bool strStartsWithComment(std::string sentence) {
  if (sentence.empty()) return false;
  return sentence[0] == COMMENT_KEY;
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

std::vector<std::string> split(std::string sentence,
                               std::string delimiter) {
  std::vector<std::string> words;
  size_t pos = 0;
  std::string token;
  while ((pos = sentence.find(delimiter)) != std::string::npos) {
    words.emplace_back(sentence.substr(0, pos));
    sentence.erase(0, pos + delimiter.length());
  }
  words.push_back(sentence);  // add last word
  return words;
}

// convert to UPPER CASES
std::string toUpperCase(std::string str) {
  for (int i = 0; i < str.length(); i++)
    str[i] = toupper(str[i]);
  return str;
}

// convert to lower case
std::string toLowerCase(std::string str) {
  for (int i = 0; i < str.length(); i++)
    str[i] = tolower(str[i]);
  return str;
}

/************************************************************************
* Read the options
*************************************************************************/
void loadOptions(
    const std::string &strOptionFile,
    std::function<bool(const std::string&, std::fstream *file)> load) {
  // open the file
  std::fstream file;
  openFile(strOptionFile, &file);

  // fill the attributes of the options structure
  std::string field;
  while (file.good()) {
    char sep = Tools::readUntilOneOfTwoChar(&file, '\n', '=', &field);
    // ignore line
    if (field.empty() || Tools::strEndsWith(field, "\n")) {
      continue;
    } else if (Tools::strStartsWithComment(field)) {
      // read line
      if (sep == '=')
        Tools::readUntilChar(&file, '\n', &field);
    } else if (load(field, &file)) {
      // read line
      Tools::readUntilChar(&file, '\n', &field);
    } else {
      Tools::throwError(
          "Field not recognized (%s) when reading the options from file %s",
          field.c_str(), strOptionFile.c_str());
    }
  }
}

void openFile(const std::string &fileName, std::fstream *file) {
  // open the file
#ifdef DBG
  std::cout << "Reading " << fileName << std::endl;
#endif
  file->open(fileName.c_str(), std::fstream::in);
  if (!file->is_open()) {
    std::cout << "While trying to read the file " << fileName << std::endl;
    std::cout << "The input file was not opened properly!" << std::endl;
    Tools::throwError("The input file (%s) was not opened properly!",
                      fileName.c_str());
  }
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

bool mkdirs(const std::string& dirPath) {
  // create directory if needed
#if __cplusplus >= 201703L
  return std::filesystem::create_directories(dirPath);
#else
  return false;
#endif
}

LogOutput::LogOutput(std::string logName, bool append) :
    LogOutput(std::move(logName), 0, append) {}

LogOutput::LogOutput(std::string logName, int width, bool append):
     width_(width), precision_(5), logName_(std::move(logName)) {
  if (logName_.empty()) {
    pLogStream_ = &(std::cout);
    isStdOut_ = true;
  } else {
    try {
      auto pos = logName_.rfind('/');
      if (pos != std::string::npos) {
        std::string outdir = logName_.substr(0, pos);
        mkdirs(outdir);
      }
      if (append) {
        pLogStream_ = new std::ofstream(logName_.c_str(), std::fstream::app);
      } else {
        pLogStream_ = new std::ofstream(logName_.c_str(), std::fstream::out);
      }
      isStdOut_ = false;
    } catch (const std::exception &e) {
      std::cout << "LogOutput::LogOutput() caught an exception=: "
                << e.what() << std::endl;
    }
  }
}

LogOutput::~LogOutput() {
  close();
  delete pLogStream_;
}

void LogOutput::addCurrentTime() {
  auto now = std::chrono::system_clock::now();
  std::time_t current_time = std::chrono::system_clock::to_time_t(now);
  std::tm* time_info = std::localtime(&current_time);
  char buffer[128];
  strftime(buffer, sizeof(buffer), "%F %T", time_info);
  (*pLogStream_) << "[" << buffer << "]    ";
}

void LogOutput::close() {
  try {
    auto *pStream = dynamic_cast<std::ofstream *>(pLogStream_);
    if (pStream) {
      if (pStream->is_open()) pStream->close();
    } else {
      pLogStream_ = nullptr;
    }
  } catch (const std::exception &e) {
    std::cout << "LogOutput::close() caught an exception=: "
              << e.what() << std::endl;
  }
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

void PThreadsPool::store() {
  if (get())
    get()->store(*this);
}

PThreadsPool::~PThreadsPool() {
  if (get())
    get()->removePtr();
}

int ThreadsPool::maxGlobalThreads_ = 1;
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

int ThreadsPool::setMaxGlobalThreadsToMax(bool wait) {
  return ThreadsPool::setMaxGlobalThreads(
      thread::hardware_concurrency(), wait);
}

int ThreadsPool::setMaxGlobalThreads(int maxThreads, bool wait) {
  std::unique_lock<mutex> l(mThreadMutex_);

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
      cThreadReleased_.wait(l);  // wait until notification
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

GlobalThreadsPool globalThreadsPool;

PThreadsPool ThreadsPool::newThreadsPool() {
  return newThreadsPool(maxGlobalThreads_);
}

PThreadsPool ThreadsPool::newThreadsPool(int nThreads) {
  return PThreadsPool(new ThreadsPool(nThreads));
}

void ThreadsPool::runOneJob(Job job, bool force_detach) {
  newThreadsPool(job.pTask_->nThreads())->run(job, force_detach);
}

ThreadsPool::ThreadsPool() : ThreadsPool(maxGlobalThreads_) {}

ThreadsPool::ThreadsPool(int nThreads):
    pThreadsPool_(nullptr), nActivePtr_(0) {
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

void ThreadsPool::run(Job job, bool force_detach) {
  if (pThreadsPool_ == nullptr)
    throwError("PThreadsPool has been deleted, "
               "cannot run a new job in this pool.");
  PTask pTask = job.pTask_;
  if (pTask->nThreads() > maxThreads_) {
    std::cerr << "Cannot run on thread pool of " << maxThreads_
              << " threads a job needing " << pTask->nThreads()
              << " threads." << std::endl;
    std:: cerr << "The job number of threads is thus decreased." << std::endl;
    pTask->nThreads(maxThreads_);
  }
  if (pTask->pThreadsPool_ != nullptr) {
    pTask->resume();  // if already attached, resume the job
  } else if (maxGlobalThreads_ <= 1 && !force_detach) {  // no concurrency
    pTask->run();
  } else {
    reserve(pTask->nThreads());  // reserve a thread
    addPtr(pTask);  // attach the task to the pool
    globalThreadsPool.run([this, pTask]() {
      pTask->run();  // run the function
      release(pTask->nThreads());  // release the thread
      pTask->detach();  // detach from the pool
    });
  }
}

void Job::run() {
  if (pTask_ != nullptr)
    pTask_->run();
}

void Job::pause() {
  if (pTask_ != nullptr)
    pTask_->pause();
}

void Job::resume() {
  if (pTask_ != nullptr)
    pTask_->resume();
}

void Job::wait() {
  if (pTask_ != nullptr)
    pTask_->wait();
}

void Task::run() {
  try {
    func_();
  } catch (const std::exception &e) {
    std::cout << "Task::run() caught an exception=: "
              << e.what() << std::endl;
  }
}

void Task::attach(const PThreadsPool& pThreadsPool) {
  std::unique_lock<mutex> l(mutex_);
  pThreadsPool_ = pThreadsPool;
}

void Task::detach() {
  std::unique_lock<mutex> l(mutex_);
  pThreadsPool_ = nullptr;
  cResume_.notify_all();
}

bool Task::active() {
  std::unique_lock<mutex> l(mutex_);
  return pThreadsPool_ != nullptr;
}

void Task::askStop() {
  std::unique_lock<mutex> l(mutex_);
  needStop_ = true;
  needPause_ = false;
  if (paused_) {
    // reserve again the threads
    pThreadsPool_->reserve(nThreads_);
    // and notify the paused thread
    paused_ = false;
    cResume_.notify_all();
  }
}

bool Task::shouldStop() {
  std::unique_lock<mutex> l(mutex_);
  return needStop_;
}

void Task::askPause() {
  std::unique_lock<mutex> l(mutex_);
  needPause_ = true;
}

bool Task::shouldPause() {
  std::unique_lock<mutex> l(mutex_);
  return needPause_;
}

void Task::pause() {
  std::unique_lock<mutex> l(mutex_);
  if (pThreadsPool_ == nullptr) return;
  // release the threads
  pThreadsPool_->release(nThreads_);
  // then wait
  needPause_ = false;
  paused_ = true;
  cResume_.wait(l);
}

void Task::resume() {
  std::unique_lock<mutex> l(mutex_);
  if (pThreadsPool_ == nullptr || !paused_) return;
  // reserve again the threads
  pThreadsPool_->reserve(nThreads_);
  // and notify the paused thread
  paused_ = false;
  cResume_.notify_all();
}

void Task::wait() {
  std::unique_lock<mutex> l(mutex_);
  while (pThreadsPool_ != nullptr)
    cResume_.wait(l);
}

Thread::Thread(): f_(nullptr), needStop_(false) {
  t_ = std::thread([this]() {
    // run an infinite loop while not finished
    std::unique_lock<std::mutex> l(mutex_);
    while (!needStop_) {
      // if nothing to do, wait
      if (f_ == nullptr)
        cWaiting_.wait(l);
      // if need to stop
      if (needStop_)
        break;
      // else run
      l.unlock();
      f_();
      // once finished, remove f_
      l.lock();
      f_ = nullptr;
      globalThreadsPool.makeAvailable(this);
    }
    delete this;
  });
  t_.detach();
}

void Thread::run(std::function<void(void)> f) {
  std::lock_guard<std::mutex> l(mutex_);
  if (f_ != nullptr)
    Tools::throwError("The thread is already associated to a job.");
  f_ = std::move(f);
  cWaiting_.notify_one();
}

void Thread::stop() {
  std::lock_guard<std::mutex> l(mutex_);
  needStop_ = true;
  cWaiting_.notify_one();
}

GlobalThreadsPool::~GlobalThreadsPool() {
  for (Thread *pThread : activeThreads)
    pThread->stop();
}

void GlobalThreadsPool::run(std::function<void(void)> f) {
  // retrieve an available thread
  std::unique_lock<std::mutex> l(mutex_);
  Thread *pThread;
  if (threadsAvailable_.empty()) {
    pThread = new Thread();
    activeThreads.push_back(pThread);
  } else {
    pThread = threadsAvailable_.front();
    threadsAvailable_.pop_front();
  }
  l.unlock();

  // attach the job
  pThread->run(std::move(f));
}

void GlobalThreadsPool::makeAvailable(Thread *pThread) {
  std::lock_guard<std::mutex> l(mutex_);
  threadsAvailable_.push_back(pThread);
}

// wait until all threads of the pool are available.
// returns true if some threads were still running
bool ThreadsPool::wait() {
  return wait(maxThreads_);
}

bool ThreadsPool::wait(int nThreadsToWait) {
  std::unique_lock<mutex> l(mThreadMutex_);  // create a lock on the mutex
  // all threads of the local pool have ended
  if (nThreadsAvailable_ == maxThreads_)
    return false;
  // wait for all threads to finish
  int threadsAvailableTarget =
      std::min(maxThreads_, nThreadsAvailable_ + nThreadsToWait);
  cThreadReleased_.wait(l, [&]() {
    return !settingMaxGlobalThreads_
        && nThreadsAvailable_ >= threadsAvailableTarget;
  });
  return true;
}

void ThreadsPool::store(const PThreadsPool& pT) {
  std::lock_guard<std::recursive_mutex> l(mutex_);
  nActivePtr_++;
  pThreadsPool_ = pT;
}

void ThreadsPool::removePtr() {
  std::lock_guard<std::recursive_mutex> l(mutex_);
  nActivePtr_--;
  // only delete the pointer if there is no active pointer anymore
  if (nActivePtr_ == 0)
    pThreadsPool_ = nullptr;
}

void ThreadsPool::addPtr(PTask pTask) {
  std::lock_guard<std::recursive_mutex> l(mutex_);
  nActivePtr_++;
  pTask->attach(pThreadsPool_);  // attach to pool
}

void ThreadsPool::reserve(int nThreads) {
  std::unique_lock<mutex> l(mThreadMutex_);
  // reserve ont thread if available, otherwise wait for one
  bool reserved = false;
  while (!reserved) {
    // wait until one thread is available and settingMaxGlobalThreads_ is false
    cThreadReleased_.wait(l, [this]() {
      return nGlobalThreadsAvailable_ > 0 && nThreadsAvailable_ > 0
          && !settingMaxGlobalThreads_;
    });
    // reserved one thread
    int minThreads = std::min(
        nThreads, std::min(nGlobalThreadsAvailable_, nThreadsAvailable_));
    nGlobalThreadsAvailable_ -= minThreads;  // reserve global threads
    nThreadsAvailable_ -= minThreads;  // reserve local threads
    nThreads -= minThreads;
    if (nThreads == 0)
      reserved = true;
  }
}

void ThreadsPool::release(int nThreads) {
  std::lock_guard<mutex> l(mThreadMutex_);
  // release threads
  nGlobalThreadsAvailable_ += nThreads;
  nThreadsAvailable_ += nThreads;
  cThreadReleased_.notify_all();  // notify any thread that were waiting
}

}  // namespace Tools
