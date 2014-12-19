//
//  MyTools.cxx
//  IDSReseau
//
//  Created by Jérémy Omer on 19/11/2013.
//  Copyright (c) 2013 Jérémy Omer. All rights reserved.
//

#include "MyTools.h"
#include <stdlib.h>     /* rand */
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <math.h>       /* pow */

using namespace std;

namespace Tools{
Timer::Timer():isInit_(0), isStarted_(0), isStopped_(0) {
	this->init();
}

void Timer::init()	{
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
void Timer::start()  {

	if (isStarted_)
		throwError("Trying to start an already started timer!");

	//clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &this->cpuInit_);
	cpuSinceStart_.tv_sec   = 0;
	cpuSinceStart_.tv_nsec  = 0;
	isStarted_ = 1;
	isStopped_ = 0;

}

void Timer::stop() {

	if ( isStopped_ )
		throwError("Trying to stop an already stopped timer!");

	/*timespec cpuNow;
	//clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuNow);

	if ( cpuNow.tv_nsec-cpuInit_.tv_nsec < 0 ) {
		cpuSinceStart_.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
		cpuSinceStart_.tv_nsec  = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
	} else {
		cpuSinceStart_.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec;
		cpuSinceStart_.tv_nsec  = cpuNow.tv_nsec - cpuInit_.tv_nsec;
	}

	if ( cpuSinceStart_.tv_nsec + cpuSinceInit_.tv_nsec >= 1e09 ) {
		cpuSinceInit_.tv_sec  = cpuSinceInit_.tv_sec + cpuSinceStart_.tv_sec + 1;
		cpuSinceInit_.tv_nsec = cpuSinceInit_.tv_nsec + cpuSinceStart_.tv_nsec - 1e09;
	} else  {
		cpuSinceInit_.tv_sec  = cpuSinceInit_.tv_sec + cpuSinceStart_.tv_sec;
		cpuSinceInit_.tv_nsec = cpuSinceInit_.tv_nsec + cpuSinceStart_.tv_nsec;
	}
*/
	isStarted_ = 0;
	isStopped_ = 1;
	coStop_ ++;
}

const double Timer::dSinceInit() {

	if (isStarted_) {
		//timespec cpuNow;
		//clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuNow);

		//timespec cpuTmp;
		/*if ( cpuNow.tv_nsec-cpuInit_.tv_nsec < 0 ) {
			cpuTmp.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
			cpuTmp.tv_nsec  = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
		} else {
			cpuTmp.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec;
			cpuTmp.tv_nsec  = cpuNow.tv_nsec - cpuInit_.tv_nsec;
		}*/

		//timespec cpuCurrent;
		/*if ( cpuTmp.tv_nsec + cpuSinceInit_.tv_nsec >= 1e09 ) {
			cpuCurrent.tv_sec  = cpuSinceInit_.tv_sec + cpuTmp.tv_sec + 1;
			cpuCurrent.tv_nsec = cpuSinceInit_.tv_nsec + cpuTmp.tv_nsec - 1e09;
		} else  {
			cpuCurrent.tv_sec  = cpuSinceInit_.tv_sec + cpuTmp.tv_sec;
			cpuCurrent.tv_nsec = cpuSinceInit_.tv_nsec + cpuTmp.tv_nsec;
		}*/

		return  0.0;
		//return (double) cpuCurrent.tv_sec + (double)cpuCurrent.tv_nsec/1.0e09;
	}
	else if (isStopped_)  {
		return (double) cpuSinceInit_.tv_sec + (double)cpuSinceInit_.tv_nsec/1.0e09;
	}
	else
		throwError("Trying to get the value of an unitialized timer!");

  return 0.0;

}

const double Timer::dSinceStart() {

	if (!isStarted_ && !isStopped_)
		throwError("Trying to get the value of an unitialized timer!");
	else if (isStarted_) {
		//timespec cpuNow;
		//clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpuNow);

		/*if ( cpuNow.tv_nsec-cpuInit_.tv_nsec < 0 ) {
			cpuSinceStart_.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
			cpuSinceStart_.tv_nsec  = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
		} else {
			cpuSinceStart_.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec;
			cpuSinceStart_.tv_nsec  = cpuNow.tv_nsec - cpuInit_.tv_nsec;
		}*/
	}

	return (double) cpuSinceStart_.tv_sec + (double)cpuSinceStart_.tv_nsec/1.0e09;

}

// intersection of two segments of R, sInter is set to NULL if the intersection
// is empty
//
void Seg::interSeg(Seg const& seg, Seg &segInter)  {

  segInter.ini = std::max(this->ini, seg.ini);
  segInter.end = std::min(this->end, seg.end);

}

void Seg::unionSeg(Seg const& seg, Seg &segUnion1,
                        Seg &segUnion2)  {
  segUnion1.set(1.0,-1.0);
  segUnion2.set(1.0,-1.0);


  if (this->empty())   {
    if (seg.empty()) segUnion1.setSeg(seg);
  }
  else if (seg.empty()) {
    segUnion1.setSeg(*this);
  }
  else if ((this->ini > seg.end+1.0e-09) || (seg.ini > this->end+1.0e-09)) {
    segUnion1.setSeg(*this);
    segUnion2.setSeg(seg);
  }
  else {
    segUnion1.set(std::min(this->ini, seg.ini), std::max(this->end, seg.end));
  }
}

// copy the content of an interval
//
void Interval::copy(Interval const &interval) {

  this->clear();

  for (int i = 0; i < interval.coSeg; i++) {
    this->addSeg(interval.vSeg[i]);
  }
}

// intersection of a segment with the interval
//
void Interval::interSeg(const Seg &seg, Interval &intervalInter) {

  for (int i = 0; i < this->coSeg; i++) {
    Seg segInter;

    vSeg[i].interSeg(seg, segInter);

    if (!segInter.empty()) intervalInter.addSeg(segInter);
  }
}

// union of a segment with the interval
//
void Interval::unionSeg(Seg const &seg){

  if (seg.empty()) return;

  int count = 0;
  Seg segUnion(seg);
  while (count < this->coSeg) {
    Seg segTemp(vSeg[count]);

    if ((segTemp.ini > segUnion.end+1.0e-09) || (segUnion.ini > segTemp.end+1.0e-09)) {
      //intervalUnion.addSeg(segTemp);
      count++;
    }
    else {
      segUnion.set(std::min(segUnion.ini, segTemp.ini), std::max(segUnion.end, segTemp.end));
      vSeg.erase(vSeg.begin()+count);
      coSeg--;
    }
  }

  this->addSeg(segUnion);
}

// intersection with another segment
//
void Interval::interInterval(Interval const &interval,
                                  Interval & intervalInter)  {

  for (int i = 0; i < interval.coSeg; i++)  {
    Seg segTemp(interval.vSeg[i]);

    // the intersection is added to the input interval, it is ok if the segments
    // included in this and in interval are disjoint
    this->interSeg(segTemp, intervalInter);
  }

}

// Wrap the interval to [-pi;pi[
//
void Interval::wrapToPi() {

  int count = 0;
  int coSegIni = this->coSeg;
  while (count < coSegIni) {

    if (vSeg[count].end-vSeg[count].ini >= 2*M_PI)
      vSeg[count].set(-M_PI,M_PI);
    else {
      double dIni = wrapAnglePi(vSeg[count].ini);
      double dEnd = dIni + vSeg[count].end-vSeg[count].ini;

      if (dEnd > M_PI)  {
        vSeg[count].set(dIni, M_PI-1.0e-09);

        dEnd = wrapAnglePi(dEnd);
        Seg newSeg(-M_PI, dEnd);
        this->addSeg(newSeg);
      }
      else
        vSeg[count].set(dIni, dEnd);
    }

    count++;
  }

  Interval intervalTemp;
  intervalTemp.copy(*this);

  while (intervalTemp.coSeg > 0) {
    this->unionSeg(intervalTemp.vSeg.back());
    intervalTemp.removeLastSeg();
  }

}

// Wrap an angle value to the interval [-pi;pi[
//
double wrapAnglePi(const double angle) {

  double remainder = (double) fmod(angle, 2*M_PI) + ((angle >= 0) ? 0:2*M_PI);
  return  (remainder < M_PI) ? remainder:remainder-2*M_PI;

}

// Throw an exception with the input message
//
void throwError(const char* exceptionMsg)  {
	try {
		throw std::string(exceptionMsg);
	} catch (const std::string str) {
		printf("Exception caught: %s\n", str.c_str());
		throw;

	}
}

// Display a debug message
//
void debugMsg(const char* debugMsg, int debugLevel)	{

	if (debugLevel >= DEBUG)
		printf("%s\n", debugMsg);
}

// Read a file stream until the separating character is met
//
bool readUntilChar(fstream *file, char separateur, string *pTitle) {
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
  while (cTmp != separateur && file->good() )  {
    pTitle->push_back(cTmp);
    cTmp = file->get();
  }

  if (!file->good())
		return false;

	return true;
}

// Solve the second degree equation ax^2+bx+c=0
// Returns true when there is a solution in R and false otherwise
//
bool secondDegree(const double &a, const double &b, const double &c,
                       double& x1, double& x2)  {
  double delta = pow(b, 2) - 4*a*c;

  if (delta < -1.0e-9)
    return false;
  else if (delta >= -1.0e-9 && delta <=0)
    delta = 0;

  double sqrtDelta = sqrt(delta);
  x1 = (-b-sqrtDelta)/(2*a);
  x2 = (-b+sqrtDelta)/(2*a);

  return true;
}
}
