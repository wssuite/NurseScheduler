//
//  MyTools.h
//
//  Created by Jérémy Omer on 19/11/2013.
//  Copyright (c) 2013 Jérémy Omer. All rights reserved.
//

#ifndef __MyTools__
#define __MyTools__

#include <iostream>
#include <fstream>
#include <streambuf>
#include <time.h>
#include <string>
#include <vector>

#define _USE_MATH_DEFINES // needed for the constant M_PI
#include <math.h>

using std::vector;

static const int DEBUG = 1;
// definitions of multi-dimensional int vector types
//
typedef vector<vector<int> > vector2D;
typedef std::vector<std::vector<std::vector<int> > > vector3D;

namespace Tools{
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

	int coStop_;
	bool isInit_;
	bool isStarted_;
	bool isStopped_;

public:
	void init();
	bool isInit() {return isInit_;}
	void start();
	void stop();

	const double dSinceInit();
	const double dSinceStart();

};

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

// structure for a segment of R
//
struct Seg {
  double ini, end;

  // constructor
  Seg(): ini(1.0), end(-1.0) {};
  Seg (double d1, double d2): ini(d1), end(d2) {}
  Seg (Seg const &seg): ini(seg.ini), end(seg.end) {}

  // set the value of a segment
  void setSeg(Seg seg) {
    ini = seg.ini;
    end = seg.end;
  }
  void set(double dIni, double dEnd)  {
    ini = dIni;
    end = dEnd;
  }

  // union and intersection with another segment
  //
  void interSeg(Seg const& seg, Seg &segInter);
  void unionSeg(Seg const& seg, Seg &segUnion1, Seg &segUnion2);

  // test whether the segment in attribute is empty
  bool empty() const {return ini > end;}
};



// structure for an interval formed by a union of segments of R
//
struct Interval {
  int coSeg;
  std::vector<Seg> vSeg;

  // constructor/destructor
  Interval (): coSeg(0) {}
  ~Interval () {
    vSeg.clear();
  }

  // add a segment at the end of the interval
  //
  void addSeg(Seg const &seg) {
    vSeg.push_back(seg);
    coSeg++;
  }

  // remove the last segment of the interval
  //
  void removeLastSeg() {
    vSeg.pop_back();
    coSeg--;
  }

  // copy the content of an interval
  //
  void copy(Interval const &interval);

  // clear the content of the interval
  //
  void clear()  {
    coSeg = 0;
    vSeg.clear();
  }

  // union and intersection of a segment with the interval
  //
  void interSeg(Seg const &seg, Interval &intervalInter);
  void unionSeg(Seg const &seg);

  // intersection with another segment (modifies the input interval)
  //
  void interInterval(Interval const &interval, Interval &intervalInter);

  // Wrap the interval to [-pi;pi[
  //
  void wrapToPi();

};

// norm of a two dimensions vector
//
inline double norm(double x, double y)  {return sqrt(pow(x,2) + pow(y,2));}

// Wrap an angle value to the interval [-pi;pi[
//
double wrapAnglePi(const double angle);

// Throw an exception with the input message
//
void throwError(const char* exceptionMsg);

// Display a debug message
//
void debugMsg(const char* debugMsg, int debugLevel);

// Read a file stream until the separating character is met
//
bool readUntilChar(std::fstream *file, char separateur, std::string *pTitle);

// Solve the second degree equation ax^2+bx+c=0
// Returns true when there is a solution in R and false otherwise
//
bool secondDegree(const double &a, const double &b, const double &c,
                       double& x1, double& x2);

}
#endif /* defined(__IDSReseau__MyTools__) */
