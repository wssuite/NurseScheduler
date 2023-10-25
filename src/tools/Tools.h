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

#ifndef SRC_TOOLS_TOOLS_H_
#define SRC_TOOLS_TOOLS_H_

#define _USE_MATH_DEFINES  // needed for the constant M_PI

#include <sys/time.h>
#include <sys/resource.h>

#include <cmath>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <chrono>  // NOLINT (suppress cpplint error)
#include <memory>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <list>
#include <map>
#include <set>
#include <streambuf>
#include <string>
#include <utility>
#include <algorithm>
#include <vector>
#include <exception>
#include <cfloat>
#include <random>
#include <functional>
#include <thread>  // NOLINT (suppress cpplint error)
#include <mutex>  // NOLINT (suppress cpplint error)
#include <condition_variable>  // NOLINT (suppress cpplint error)
#include <cassert>

static const char COMMENT_KEY = '#';
static const int SHIFT_PAD = 3;
static const char REST_SHIFT[] = "Rest";
static const char REST_DISPLAY[] = " - ";  // should be of the size of pad
static const int DECIMALS = 3;  // precision when printing floats
static const int NB_SHIFT_UNLIMITED = 28;

static const int LARGE_INT = 9999;
static const int INFEAS_COST = 9999;
static const int HARD_COST = 999999;
static const int LARGE_TIME = 999999;

bool isLargeNumber(double c);
bool isInfeasibleCost(double c);
bool isHardCost(double c);
bool isSoftCost(double c);


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

template<typename ...Args>
std::string string(const char *str, Args... args) {
  char buff[999];
  snprintf(buff, sizeof(buff), str, args...);
  return buff;
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
      std::exception(), what_(std::move(what)) {}
  template<typename ...Args>
  NSException(const char *str, Args... args) {
    char buff[999];
    snprintf(buff, sizeof(buff), str, args...);
    what_ = buff;
  }

  const char *what() const throw() override {
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

//  struct rusage {
//    struct timeval ru_utime; /* user time used */
//    struct timeval ru_stime; /* system time used */
//    long   ru_maxrss;        /* maximum resident set size */
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

rusage getRUsage();

float getResidentMemoryGB();

// Struct for very simple storage of an hour of the day
//
struct Time {
  Time() : hh(0), mm(0) {}
  Time(double hour, double min) : hh(hour), mm(min) {
    assert(hour >= 0 && hour <= 23);
    assert(min >= 0 && min <= 59);
  }

  // number of hours elapsed since time
  double diff(const Time &time) const {
    // if time is smaller, we assume it is on the previous day
    if (hh > time.hh || (hh == time.hh && mm >= time.mm))
      return (hh - time.hh) + (mm - time.mm) / 60.0;
    else
      return 24.0 + (hh - time.hh) + (mm - time.mm) / 60.0;
  }
  // number of hours elapsed since time
  bool isAfter(const Time &time) const {
    // if time is smaller, we assume it is on the previous day
    return (hh > time.hh || (hh == time.hh && mm >= time.mm));
  }

  bool equals(const Time &time) const {
    return (time.hh == hh) && (time.mm == mm);
  }

  const double hh;
  const double mm;
};

tm *dateForDay(const tm *startDate, const int &dayId);

// Read a file stream until the separating character (or one of them) is met
// Store the characters read until the separating character in pStrRead
char readUntilOneOfTwoChar(std::fstream *pFile,
                           char separater1,
                           char separater2,
                           std::string *pStrRead);

// Store the next line which is not a comment in pStrRead
bool readLinesUntilNotAComment(
    std::fstream *pFile, std::string *pStrRead, bool acceptEmptyLine = true);

// read until end of line
bool readLine(std::fstream *pFile, std::string *pStrRead);

// Read a file stream until the separating character is met
bool readUntilChar(std::fstream *file, char separator, std::string *pTitle);

// Read a file stream until the separating character is met and keep reading
// while the delimiter repeats
bool readUntilAndWhileChar(std::fstream *file, char separator, std::string
*pTitle);

// Checks if the string (sentence) starts with the given substring (word)
bool strStartsWith(std::string sentence, std::string word);

// Checks if the string (sentence) starts with the comment character
bool strStartsWithComment(std::string sentence);

// Checks if the string (sentence) ends with the given substring (word)
//
bool strEndsWith(std::string sentence, std::string word);

// split string
std::vector<std::string> split(std::string sentence, std::string delimiter);

// Read a dash separated date in format yyyy-mm-dd and return a tm* object
//
tm *readDateFromStr(const std::string &dateStr);
int readNbWeeksFromDate(const std::string &dateStrBeg,
                        const std::string &dateStrEnd);
// Read a colon separated hour in format hh:mm:ss and return a tm* object
//
Tools::Time readHourFromStr(const std::string &hourStr);

// Read the parameters of a bounded resource
//
int readBoundedResource(std::fstream *file,
                        int *lbOn, int *lb, int *lbCost,
                        int *ubOn, int *ub, int *ubCost);
int readUbResource(std::fstream *file,
                   int *ubOn, int *ub, int *ubCost);

// Parse a T list written as string with a char delimiter
//
template<typename T>
std::vector<T> tokenize(std::string str, char delim) {
  std::vector<T> Tlist;
  size_t start;
  size_t end = 0;
  while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
    T i;
    end = str.find(delim, start);
    std::stringstream ss(str.substr(start, end - start));
    ss >> i;
    if (ss.rdbuf()->in_avail() > 0)
      Tools::throwError("%s cannot be tokenized with the type %s",
                        str.c_str(), typeid(T).name());
    Tlist.push_back(i);
  }
  if (str.back() == delim) {
    T i;
    Tlist.push_back(i);
  }
  return Tlist;
}

// convert to UPPER CASES
std::string toUpperCase(std::string str);

// convert to lower case
std::string toLowerCase(std::string str);

// trim from start (in place)
void ltrim(std::string *s);

// trim from end (in place)
void rtrim(std::string *s);

// trim from both ends (in place)
void trim(std::string *s);

std::string loadOptions(
    const std::string &strOptionFile,
    std::function<bool(const std::string &, std::fstream *file)>);

void openFile(const std::string &fileName, std::fstream *file);

// Create a random generator
// the objective is to be sure to have always the same sequence of number
//
std::minstd_rand getANewRandomGenerator(bool printSeed = false);

std::minstd_rand getANewRandomGenerator(int rdmSeed, bool printSeed = false);

// Initialize the random generator with a given seed
void initializeRandomGenerator();

void initializeRandomGenerator(int rdmSeed);

// round with probability
int roundWithProbability(double number);

// great common divider
int gcd(int a, int b);
int gcd(const std::vector<int> &numbers);
int gcd(const std::set<int> &numbers);

template <class _InputIterator>
int gcd(_InputIterator it, _InputIterator end) {
  if (it == end) return 0;
  int d = *(it++);
  for (; it != end; it++)
    d = gcd(d, *it);
  return d;
}

// convert a number to a string
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

template<class T>
void initVector(std::vector<T> *v, int m) {
  v->clear();
  v->resize(m);
}

template<class T>
void initVector2D(vector2D<T> *v2D, int m, int n) {
  v2D->clear();
  v2D->resize(m);
  for (std::vector<T> &v : *v2D) initVector(&v, n);
}

template<class T>
void initVector3D(vector3D<T> *v3D, int m, int n, int p) {
  v3D->clear();
  v3D->resize(m);
  for (vector2D<T> &v2 : *v3D) initVector2D(&v2, n, p);
}

template<class T>
void initVector4D(vector4D<T> *v4D, int m, int n, int p, int q) {
  v4D->clear();
  v4D->resize(m);
  for (vector3D<T> &v3 : *v4D) initVector3D(&v3, n, p, q);
}

// return the object at position i (support negative index)
template<class T>
const T &get(const std::vector<T> &v, int i) {
  if (i >= 0) return v[i];
  return v[v.size() + i];
}

template<class T>
void insert_back(std::vector<T> *v1, const std::vector<T> &v2) {
  v1->insert(v1->end(), v2.begin(), v2.end());
}

// remove an element from a list
template<class T>
bool erase(std::vector<T> *vec, const T &el, bool throwErrorNotFound = false) {
  auto it = std::find(vec->begin(), vec->end(), el);
  if (it != vec->end()) {
    vec->erase(it);
    return true;
  }
#ifdef NS_DEBUG
  if (throwErrorNotFound)
    throwError("Element to delete not found in vector.");
#endif
  return false;
}

template<typename T>
static inline bool includes(const std::vector<T> &vec1,
                            const std::vector<T> &vec2) {
  auto it1 = vec1.begin();
  for (auto it2 = vec2.begin(); it2 != vec2.end(); it2++) {
    while (*it1 < *it2) {
      if (it1 == vec1.end() - 1) return false;
      else
        it1++;
    }
    if (*it1 != *it2) return false;
  }
  return true;
}

// return the name for the enum from a given map of name
template<typename T>
static const std::string &getNameForEnum(
    const std::map<std::string, T> &typesByName, T type) {
  for (const auto &p : typesByName)
    if (p.second == type) return p.first;
  Tools::throwError("No name found in the map for the given type");
  return typesByName.begin()->first;
}

template<typename T>
static std::map<T, std::string> buildNamesByType(
    const std::map<std::string, T> &typesByName) {
  std::map<T, std::string> namesByType;
  for (const auto &p : typesByName)
    namesByType[p.second] = p.first;
  return namesByType;
}

template<typename T>
static std::map<T, std::string> buildPrettyNamesByType(
    const std::map<std::string, T> &typesByName,
    const std::string &delimiter = "_") {
  std::map<T, std::string> prettyNamesByType;
  for (const auto &p : typesByName) {
    std::vector<std::string> words = split(p.first, delimiter);
    std::string sentence;
    for (const std::string &word : words) {
      std::string s = toLowerCase(word);
      if (sentence.empty())  // First case will be an upper case
        s[0] = toupper(s[0]);
      else  // add a space between words
        sentence += " ";
      sentence += s;
    }
    prettyNamesByType[p.second] = sentence;
  }
  return prettyNamesByType;
}

template<class T>
class FixedSizeList {
 public:
  explicit FixedSizeList(int size) : fixedSize_(size) {}

  typedef typename std::list<T>::iterator iterator;

  T *insert(iterator it, T v) {
    T *v2 = &(*list_.insert(it, std::move(v)));
    if (list_.size() > fixedSize_) list_.resize(fixedSize_);
    return v2;
  }

  T *push_back(T v) {
    if (list_.size() >= fixedSize_)
      throwException("FixedSizeList has already reached its size of %s",
                     fixedSize_);
    list_.push_back(v);
    return &list_.back();
  }

  const T &front() const { return list_.front(); }

  const std::list<T> &list() const { return list_; }

  iterator begin() { return list_.begin(); }

  iterator end() { return list_.end(); }

  const int fixedSize_;

 private:
  std::list<T> list_;
};

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
  std::vector<T> ANS = v1;
  ANS.insert(ANS.end(), v2.begin(), v2.end());
  return ANS;
}

template<typename T>
std::vector<T> appendVectors(const vector2D<T> &vectors) {
  std::vector<T> ANS;
  for (const std::vector<T>& v : vectors)
    ANS.insert(ANS.end(), v.begin(), v.end());
  return ANS;
}

template<typename T>
std::vector<T> reccursiveAppendVectors(const vector2D<T> &vectors) {
  return appendVectors(vectors);
}

template<typename T>
std::vector<T> reccursiveAppendVectors(const vector3D<T> &vectors) {
  return reccursiveAppendVectors(appendVectors(vectors));
}

template<typename T>
std::vector<T> reccursiveAppendVectors(const vector4D<T> &vectors) {
  return reccursiveAppendVectors(appendVectors(vectors));
}

// High resolution timer class to profile the performance of the algorithms
// Warning : the timer class of the stl to not seem to be portable I observed
// problems on windows for instance and it requires some precompiler
// instructions to work on mac
typedef std::chrono::duration<int, std::nano> nanoseconds_type;

class Timer {
 public:
  // constructor and destructor
  //
  explicit Timer(std::string name, bool start = false);

  ~Timer() {}

 private:
  std::string name_;

  std::chrono::time_point<std::chrono::system_clock> cpuInit_;
  std::chrono::duration<double> cpuSinceStart_, cpuSinceInit_;

  int nStop_;  // number of times the timer was stopped
  bool isStarted_;
  bool isStopped_;

  double getSeconds(const std::chrono::duration<double> &d) const;

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

bool mkdirs(const std::string &directoryPath);

class LogOutput {
 private:
  std::ostream *pLogStream_;
  int width_;
  int precision_;
  std::string logName_ = "";
  bool isStdOut_;
  LogOutput *pLogCout;

 public:
  explicit LogOutput(std::string logName = "",
                     bool append = true,
                     bool alwaysPrint = false);

  LogOutput(std::string logName, int width,
            bool append = false, bool alwaysPrint = false);

  ~LogOutput();

  void close();

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
    if (pLogCout) (*pLogCout) << output;
    return *this;
  }

  // write in the command window as well as in the log file
  //
  template<typename T>
  void print(const T &output) {
    (*pLogStream_) << output;
    if (!isStdOut_) std::cout << output;
  }

  template<typename T, typename ...Args>
  void print(const char *str, const T &arg0, Args... args) {
    char buff[999];
    snprintf(buff, sizeof(buff), str, arg0, args...);
    print(buff);
  }

  template<typename T>
  void printnl(const T &output) {
    (*pLogStream_) << output << std::endl;
    if (!isStdOut_) std::cout << output << std::endl;
  }

  template<typename T, typename ...Args>
  void printnl(const char *str, const T &arg0, Args... args) {
    char buff[999];
    snprintf(buff, sizeof(buff), str, arg0, args...);
    printnl(buff);
  }

  LogOutput &operator<<(std::ostream &(*func)(std::ostream &)) {
    func(*pLogStream_);
    if (pLogCout) (*pLogCout) << func;
    return *this;
  }

  LogOutput &addCurrentTime();
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
class ThreadsPool;

struct PThreadsPool : public std::shared_ptr<ThreadsPool> {
  template<typename ...Args>
  PThreadsPool(Args... args):  // NOLINT
      std::shared_ptr<ThreadsPool>(args...) {
    store();
  }

  void store();

  ~PThreadsPool();
};

class Job;

class Task {
  explicit Task(std::function<void(Job)> f, int _nThreads) :
      func_(std::move(f)), nThreads_(_nThreads) {}

  // run the task
  void run(const Job &job);

  // attach to the pool
  void attach(const PThreadsPool &pThreadsPool);

  // detach from the pool
  void detach();

  // check if the task is still active for the pool (i.e. attached)
  bool active();

  // ask the task to stop
  void askStop();

  // check if the task has been asked to stop
  bool shouldStop();

  // ask the task to pause
  void askPause();

  // check if the task has been asked to pause
  bool shouldPause();

  // pause the current thread and wait to be resumed by another thread
  void pause();

  // resume the paused thread
  void resume();

  // wait for the task to finish
  void wait();

  int nThreads() const { return nThreads_; }

  void nThreads(int n) { nThreads_ = n; }

  int nThreads_;
  std::function<void(Job)> func_;
  PThreadsPool pThreadsPool_ = nullptr;

  std::mutex mutex_;
  std::condition_variable cResume_;

  bool needStop_ = false;
  bool needPause_ = false;
  bool paused_ = false;

  friend class Job;
  friend class ThreadsPool;
};
typedef std::shared_ptr<Task> PTask;

/* A wrapper around PTask. Also check if a job is associated to a task, i.e.
 * the class Job has not been initialized with the default constructor */
class Job {
 public:
  Job() : pTask_(nullptr) {}
  explicit Job(std::function<void(Job)> f, int nThreads = 1) :
      pTask_(new Task(std::move(f), nThreads)) {}

  void copy(const Job &job) {
    pTask_ = job.pTask_;
  }

  // true if a task is attached to the job
  bool defined() const {
    return pTask_ != nullptr;
  }

  // check if any active tsk is attached
  bool active() {
    return pTask_ != nullptr && pTask_->active();
  }

  // true if not active
  bool finish() {
    return !active();
  }

  // ask to the attached task to stop if any
  void askStop() {
    if (pTask_) pTask_->askStop();
  }

  // check if the attached task has been asked to stop if any
  bool shouldStop() {
    if (pTask_) return pTask_->shouldStop();
    return false;
  }

  // ask to the attached task to pause if any
  void askPause() {
    if (pTask_) pTask_->askPause();
  }

  // check if the attached task has been asked to pause if any
  bool shouldPause() {
    if (pTask_) return pTask_->shouldPause();
    return false;
  }

  // run the task
  void run();

  // pause the current thread and wait to be resumed by another thread
  // if any attached task
  void pause();

  // resume the paused thread if any attached task
  void resume();

  // wait for the task to finish
  void wait();

 private:
  PTask pTask_ = nullptr;

  friend class ThreadsPool;
};

class Thread {
 public:
  Thread();
  ~Thread() = default;

  void stop();

  void run(std::function<void(void)> f);

 private:
  std::function<void(void)> f_;
  bool needStop_;
  std::mutex mutex_;
  std::condition_variable cWaiting_;
  std::thread t_;
};

class GlobalThreadsPool {
 public:
  GlobalThreadsPool() = default;
  ~GlobalThreadsPool();

  void run(std::function<void(void)> f);

  void makeAvailable(Thread *pThread);

 private:
  std::list<Thread *> threadsAvailable_, activeThreads;
  std::mutex mutex_;
};

class ThreadsPool {
 public:
  static int getNGlobalThreadsAvailable();

  static int getMaxGlobalThreads();

  // if maximum number of threads decreases,
  // wait for some to be released if wait is true.
  // otherwise decrease of the number of available threads and print a warning.
  // return the new true maxGlobalThreads_
  static int setMaxGlobalThreads(int maxThreads, bool wait = true);
  static int setMaxGlobalThreadsToMax(bool wait = true);

  // method to create a new ThreadsPool
  static PThreadsPool newThreadsPool();
  static PThreadsPool newThreadsPool(int nThreads);

  // create a pool to run only one job
  static void runOneJob(Job job, bool force_detach = false);

 private:
  // global maximum number of threads
  static int maxGlobalThreads_;
  // global number of available threads  for the local pool
  static int nGlobalThreadsAvailable_;
  // true if currently trying to set the maxGlobalThreads_
  static bool resizingMaxGlobalThreads_;
  static std::mutex mThreadMutex_;
  static std::condition_variable cThreadReleased_;

 public:
  void run(Job job, bool force_detach = false);

  // wait until all threads of the pool are available.
  // returns true if some threads were still running
  bool wait();

  bool wait(int nThreadsToWait);

  bool available() const { return nThreadsAvailable_ > 0; }

  bool busy() const { return nThreadsAvailable_ < maxThreads_; }

  bool idle() const { return nThreadsAvailable_ == maxThreads_; }

 private:
  ThreadsPool();

  explicit ThreadsPool(int nThreads);

  int maxThreads_;  // local maximum number of threads
  // local number of available threads for the local pool
  int nThreadsAvailable_;

  void reserve(int nThreads);

  void release(int nThreads);

  // store the shared_pointer to ensure that a copy is always stored and
  // thus we can control when to destroy each pool. It is necessary to keep
  // at least one copy of the shared pointer until each job running
  // in the pool finishes
  PThreadsPool pThreadsPool_;
  // count the number of jobs that has been run and are still not deleted
  int nActivePtr_;
  std::recursive_mutex mutex_;

  void store(const PThreadsPool &pT);

  void addPtr(PTask pTask);
  void removePtr();

  friend class PThreadsPool;
  friend class Task;
};

}  // namespace Tools

#endif  // SRC_TOOLS_TOOLS_H_
