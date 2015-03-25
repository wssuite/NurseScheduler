//
//  MyTools.h
//
//  Created by J��r��my Omer on 19/11/2013.
//  Copyright (c) 2013 J��r��my Omer. All rights reserved.
//

#ifndef __MyTools__
#define __MyTools__

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <streambuf>
#include <string>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <exception>
#include <algorithm>
#include <cfloat>

#define _USE_MATH_DEFINES // needed for the constant M_PI
#include <math.h>

using std::map;
using std::string;
using std::vector;
using std::set;

static const int DEBUG = 1;
static const string REST_SHIFT = "None";


// definitions of multi-dimensional int vector types
//
typedef vector<vector<int> > vector2D;
typedef std::vector<std::vector<std::vector<int> > > vector3D;

namespace Tools{


// class defining my own type of exceptions
//
class myException : public std::exception
{
public:
  myException( const char * Msg, int Line )
  {
      std::ostringstream oss;
      oss << "Erreur ligne " << Line << " : "
          << Msg;
      this->msg = oss.str();
  }

  virtual ~myException() throw() {}

  virtual const char * what() const throw() { return this->msg.c_str();}

private:
    std::string msg;
};

// Throw an exception with the input message
//
void throwError(const char* exceptionMsg);


// Display a debug message
//
void debugMsg(const char* debugMsg, int debugLevel);

// Read a file stream until the separating character is met
//
bool readUntilChar(std::fstream *file, char separateur, std::string *pTitle);

// initializes a 1D, 2D or 3D Vector of the given size (filled only with zeroes)
//
void initVector(vector<int>* v1D, int m);
void initVector2D(vector2D* v2D, int m, int n);
void initVector3D(vector3D* v3D, int m, int n, int p);
void initDoubleVector(vector<double>* v1D, int m);
void initDoubleVector2D(vector< vector< double > >* v2D, int m, int n);

// To get the day from its id and vice-versa
// First day is always supposed to be a Monday
//
string intToDay(int dayId);
int dayToInt(string day);
bool isSaturday(int dayId);
bool isSunday(int dayId);
bool isWeekend(int dayId);
bool containsWeekend(int startDate, int endDate);

// High resolution timer class to profile the performance of the algorithms
// Warning : the timer class of the stl to not seem to be portable I observed
// problems on windows for instance and it requires some precompiler
// instructions to work on mac
//
class Timer
{
public:

	//constructor and destructor
	//
	Timer();
	~Timer() {}

private:
	timespec cpuInit_;
	timespec cpuSinceStart_;
	timespec cpuSinceInit_;

	int coStop_;	//number of times the timer was stopped
	bool isInit_;
	bool isStarted_;
	bool isStopped_;

public:
	void init();
	bool isInit() {return isInit_;}
	void start();
	void stop();

	// get the time spent since the initialization of the timer and since the last
	// time it was started
	//
	const double dSinceInit();
	const double dSinceStart();

};

// Instantiate an obect of this class to write directly in the attribute log
// file.
// The class can be initialized with an arbitrary width if all the outputs must
// have the same minimum width. Precision can be set to force a maximum width
// as far as floating numbers are concerned.
//
class LogOutput
{
private:
	std::fstream logStream_;
  bool isFormatted_;
	int width_;

public:
	LogOutput(string logName):width_(0) {
		logStream_.open(logName.c_str(), std::fstream::out);
	}
	LogOutput(string logName, int width):width_(width){
		logStream_.open(logName.c_str(), std::fstream::out);
	}

	~LogOutput() {logStream_.close();}

	// switch from unformatted to formatted inputs and reversely
	//
	void switchToFormatted(int width) {
		width_ = width;
	}
	void switchToUnformatted() {
		width_=0;
	}

	// set flags in the log stream
	//
	void setFlags(std::ios_base::fmtflags flag) {logStream_.flags(flag);}

	// modify the precision used to write in the stream
	//
	void setPrecision(int precision) {logStream_.precision(precision);}

	// redefine the output function
	//
	template<typename T>
	LogOutput& operator<<(const T& output)
	{
		logStream_.width(width_);
		logStream_ << output;

		return *this;
	}

	// write in the command window as well as in the log file
	//
	template<typename T>
	void print(const T& output) {
		logStream_ << output;
		std::cout << output;
	}

	LogOutput& operator<<(std::ostream& (*func)(std::ostream&)) {

		func(logStream_);
		return *this;
	}
};

// Can be used to create an output stream that writes with a
// constant width in all future calls of << operator
//
class FormattedOutput
{
private:
	int width;
	std::ostream& stream_obj;

public:
	FormattedOutput(std::ostream& obj, int w): width(w), stream_obj(obj) {}

	template<typename T>
	FormattedOutput& operator<<(const T& output)
	{
		stream_obj.width(width);
		stream_obj << output;

		return *this;
	}

	FormattedOutput& operator<<(std::ostream& (*func)(std::ostream&))
	{
		func(stream_obj);
		return *this;
	}
};

}
#endif /* defined(__IDSReseau__MyTools__) */
