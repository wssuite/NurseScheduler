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

#ifndef SRC_TOOLS_MYTOOLS_H_
#define SRC_TOOLS_MYTOOLS_H_

#define _USE_MATH_DEFINES  // needed for the constant M_PI

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include <chrono>  // NOLINT (suppress cpplint error)
#include <memory>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <streambuf>
#include <string>
#include <vector>
#include <exception>
#include <algorithm>
#include <cfloat>
#include <random>
#include <functional>
#include <utility>
#include <thread>  // NOLINT (suppress cpplint error)
#include <mutex>  // NOLINT (suppress cpplint error)
#include <condition_variable>  // NOLINT (suppress cpplint error)

static const int DEBUG = 1;
static const char REST_SHIFT[] = "None";
static const int DECIMALS = 3;  // precision when printing floats
static const int NB_SHIFT_UNLIMITED = 28;

static const int LARGE_SCORE = 9999999;
static const int LARGE_TIME = 9999999;


// definitions of multi-dimensional int vector types
//
template<class T> using vector2D = std::vector<std::vector<T>>;
template<class T> using vector3D = std::vector<vector2D<T>>;
template<class T> using vector4D = std::vector<vector3D<T>>;

namespace Tools {

// Compare functions to sort
bool compareDecreasing(int i, int j);

template<typename T>
bool compareObject(const std::pair<T, double> &p1,
                   const std::pair<T, double> &p2) {
  return (p1.second < p2.second);
}

// class defining my own type of exceptions
//
class myException : public std::exception {
 public:
  myException(const char *Msg, int Line) {
    std::ostringstream oss;
    oss << "Error line " << Line << " : " << Msg;
    this->msg = oss.str();
  }

  virtual ~myException() throw() {}

  virtual const char *what() const throw() { return this->msg.c_str(); }

 private:
  std::string msg;
};

// Throw an exception with the input message
//
struct NSException : std::exception {
  explicit NSException(const char *what) : std::exception(), what_(what) {}
  explicit NSException(std::string what) :
      std::exception(), what_(what) {}
  template<typename ...Args>
  NSException(const char *str, Args... args) {
    char buff[999];
    snprintf(buff, sizeof(buff), str, args...);
    what_ = buff;
  }

  const char *what() const throw() {
    return what_.c_str();
  }

 private:
  std::string what_;
};

void throwException(const char *exceptionMsg);
void throwException(const std::string &exceptionMsg);

template<typename ...Args>
void throwException(const char *str, Args... args) {
  char buff[999];
  snprintf(buff, sizeof(buff), str, args...);
  throwException(buff);
}

void throwError(const char *exceptionMsg);
void throwError(const std::string &exceptionMsg);

template<typename ...Args>
void throwError(const char *str, Args... args) {
  char buff[999];
  snprintf(buff, sizeof(buff), str, args...);
  throwError(buff);
}

// Display a debug message
//
void debugMsg(const char *debugMsg, int debugLevel);

// Read a file stream until the separating character (or one of them) is met
// Store the characters read until the separating character in pStrRead
//
bool readUntilOneOfTwoChar(std::fstream *pFile,
                           char separater1,
                           char separater2,
                           std::string *pStrRead);

// Read a file stream until the separating character is met
//
bool readUntilChar(std::fstream *file, char separator, std::string *pTitle);

// Checks if the string (sentence) ends with the given substring (word)
//
bool strEndsWith(std::string sentence, std::string word);

// Parse a T list written as string with a char delimiter
//
template<typename T>
std::vector<T> tokenize(std::string str, char delim) {
  std::vector<T> Tlist;
  size_t start;
  size_t end = 0;
  T i;
  while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
    end = str.find(delim, start);
    std::stringstream ss(str.substr(start, end - start));
    ss >> i;
    if (ss.rdbuf()->in_avail() >0)
      Tools::throwError("%s cannot be tokenized with the type %s",
          str.c_str(), typeid(T).name());
    Tlist.push_back(i);
  }
  return Tlist;
}

// Create a random generator
// the objective is to be sure to have always the same sequence of number
//
std::minstd_rand getANewRandomGenerator();

std::minstd_rand getANewRandomGenerator(int rdmSeed);

// Initialize the random generator with a given seed
void initializeRandomGenerator();

void initializeRandomGenerator(int rdmSeed);

// round with probability
int roundWithProbability(double number);

// convert a number to a string
//
template<typename T>
std::string itoa(T num) {
  std::stringstream stream;
  stream << std::setprecision(1) << std::fixed << num;
  return stream.str();
}

// initializes a 1D, 2D, 3D or 4D Vector of the given size
// (filled only with val)
//
template<class T>
void initVector(std::vector<T> *v, int m, T val) {
  v->clear();
  v->resize(m, val);
}

template<class T>
void initVector2D(vector2D<T> *v2D, int m, int n, T val) {
  v2D->clear();
  v2D->resize(m);
  for (std::vector<T> &v : *v2D) initVector(&v, n, val);
}

template<class T>
void initVector3D(vector3D<T> *v3D, int m, int n, int p, T val) {
  v3D->clear();
  v3D->resize(m);
  for (vector2D<T> &v2 : *v3D) initVector2D(&v2, n, p, val);
}

template<class T>
void initVector4D(vector4D<T> *v4D, int m, int n, int p, int q, T val) {
  v4D->clear();
  v4D->resize(m);
  for (vector3D<T> &v3 : *v4D) initVector3D(&v3, n, p, q, val);
}

// return the object at position i (support negative index)
template<class T>
const T &get(const std::vector<T> &v, int i) {
  if (i >= 0) return v[i];
  return v[v.size() + i];
}

// Returns an integer with random value (uniform) within [minVal, maxVal]
//
int randomInt(int minVal, int maxVal);

// Returns a double with random value (uniform) within [minVal, maxVal]
//
double randomDouble(double minVal, double maxVal);

// Creates 1D/2D vectors with random values (uniform) within a given range.
//
std::vector<double> randomDoubleVector(int m, double minVal, double maxVal);

std::vector<std::vector<double> > randomDoubleVector2D(int m,
                                                       int n,
                                                       double minVal,
                                                       double maxVal);

// For an input vector of weights w, draw the index of one element assuming that
// the probability of each index is w_i/sum{w_i}
//
int drawRandomWithWeights(std::vector<double> weights);

// Draw randomly nbInd different indices between indMin and indMax
//
std::vector<int> drawRandomIndices(int nbInd, int indMin, int indMax);

// Appends the values of v2 to at the end of v1
//
template<typename T>
std::vector<T> appendVectors(
    const std::vector<T> &v1, const std::vector<T> &v2) {
  std::vector<T> ANS;
  for (int i = 0; i < v1.size(); i++) ANS.push_back(v1[i]);
  for (int i = 0; i < v2.size(); i++) ANS.push_back(v2[i]);
  return ANS;
}

// To get the day from its id and vice-versa
// First day is always supposed to be a Monday
//
std::string intToDay(int dayId);

int dayToInt(std::string day);

bool isSaturday(int dayId);

bool isSunday(int dayId);

bool isWeekend(int dayId);
bool isFirstWeekDay(int dayId);
bool isLastWeekendDay(int dayId);
int containsWeekend(int startDate, int endDate);

// High resolution timer class to profile the performance of the algorithms
// Warning : the timer class of the stl to not seem to be portable I observed
// problems on windows for instance and it requires some precompiler
// instructions to work on mac
//

  typedef std::chrono::duration<int, std::nano> nanoseconds_type;


class Timer {
 public:
  // constructor and destructor
  //
  explicit Timer(bool start = false);

  ~Timer() {}

 private:
  std::chrono::time_point<std::chrono::system_clock> cpuInit_;
  std::chrono::duration<double> cpuSinceStart_, cpuSinceInit_;

  int coStop_;  // number of times the timer was stopped
  bool isStarted_;
  bool isStopped_;

  double getSeconds(const std::chrono::duration<double>& d) const;

 public:
  void start();

  void stop();

  void reset() {
    stop();
    start();
  }

  // get the time spent since the initialization of the timer and since the last
  // time it was started
  //
  double dSinceInit() const;

  double dSinceStart();
};

// Instantiate an obect of this class to write directly in the attribute log
// file.
// The class can be initialized with an arbitrary width if all the outputs must
// have the same minimum width. Precision can be set to force a maximum width
// as far as floating numbers are concerned.
//
class LogOutput {
 private:
  std::ostream *pLogStream_;
  int width_;
  int precision_;
  std::string logName_ = "";

 public:
  explicit LogOutput(std::string logName, bool append = false) :
      width_(0), precision_(5), logName_(logName) {
    if (logName.empty()) {
      pLogStream_ = &(std::cout);
    } else if (append) {
      pLogStream_ = new std::ofstream(logName.c_str(), std::fstream::app);
    } else {
      pLogStream_ = new std::ofstream(logName.c_str(), std::fstream::out);
    }
    // logStream_.open(logName.c_str(), std::fstream::out);
  }

  LogOutput(std::string logName, int width, bool append = false)
      : width_(width), precision_(5), logName_(logName) {
    if (append) {
      pLogStream_ = new std::ofstream(logName.c_str(), std::fstream::app);
    } else {
      pLogStream_ = new std::ofstream(logName.c_str(), std::fstream::out);
    }
  }

  ~LogOutput() {
    if (!logName_.empty() && pLogStream_)
      delete pLogStream_;
    pLogStream_ = nullptr;
  }

  void close() {
    std::ofstream *pStream = dynamic_cast<std::ofstream *>(pLogStream_);
    if (pStream) {
      if (pStream->is_open()) pStream->close();
    } else {
      pLogStream_ = nullptr;
    }
  }

  // switch from unformatted to formatted inputs and reversely
  //
  void switchToFormatted(int width) {
    width_ = width;
  }

  void switchToUnformatted() {
    width_ = 0;
  }

  // set flags in the log stream
  //
  void setFlags(std::ios_base::fmtflags flag) { pLogStream_->flags(flag); }

  // modify the precision used to write in the stream
  //
  void setPrecision(int precision) { precision_ = precision; }

  // modify the width of the fields
  //
  void setWidth(int width) { width_ = width; }

  void endl() { (*pLogStream_) << std::endl; }

  // redefine the output function
  //
  template<typename T>
  LogOutput &operator<<(const T &output) {
    pLogStream_->width(width_);
    pLogStream_->unsetf(std::ios::floatfield);
    pLogStream_->precision(precision_);
    (*pLogStream_) << std::left << std::setprecision(5) << output;

    return *this;
  }

  // write in the command window as well as in the log file
  //
  template<typename T>
  void print(const T &output) {
    (*pLogStream_) << output;
    std::cout << output;
  }

  LogOutput &operator<<(std::ostream &(*func)(std::ostream &)) {
    func(*pLogStream_);
    return *this;
  }
};

// Can be used to create an output stream that writes with a
// constant width in all future calls of << operator
//
class FormattedOutput {
 private:
  int width;
  std::ostream &stream_obj;

 public:
  FormattedOutput(std::ostream &obj, int w) : width(w), stream_obj(obj) {}

  template<typename T>
  FormattedOutput &operator<<(const T &output) {
    stream_obj.width(width);
    stream_obj << output;

    return *this;
  }

  FormattedOutput &operator<<(std::ostream &(*func)(std::ostream &)) {
    func(stream_obj);
    return *this;
  }
};

// Create a pool of threads that can be used to run functions in parallel.
// Each pool has a certain number of threads available,
// but there is also a global limit on the number of threads used.
typedef std::function<void(void)> Job;

class ThreadsPool {
 public:
  static int getNGlobalThreadsAvailable();

  static int getMaxGlobalThreads();

  // if maximum number of threads decreases,
  // wait for some to be released if wait is true.
  // otherwise decrease of the number of available threads and print a warning.
  // return the new true maxGlobalThreads_
  static int setMaxGlobalThreads(int maxThreads, bool wait = true);

 private:
  // global maximum number of threads
  static int maxGlobalThreads_;
  // global number of available threads  for the local pool
  static int nGlobalThreadsAvailable_;
  // true if currently trying to set the maxGlobalThreads_
  static bool settingMaxGlobalThreads_;
  static std::mutex mThreadMutex_;
  static std::condition_variable cThreadReleased_;

 public:
  ThreadsPool();

  explicit ThreadsPool(int nThreads);

  void run(Job job);

  // wait until all threads of the pool are available.
  // returns true if some threads were still running
  bool wait();

  bool wait(int nThreadsToWait);

 private:
  int maxThreads_;  // local maximum number of threads
  // local number of available threads for the local pool
  int nThreadsAvailable_;

  void reserve();

  void release();
};

}  // namespace Tools
#endif  // SRC_TOOLS_MYTOOLS_H_
