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

// necessary because OS X does not have clock_gettime, using clock_get_time
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

using namespace std;

namespace Tools{

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


// Returns a vector2D of the given size (filled only with zeroes)
//
vector2D initVector2D(int m, int n){
	vector2D rep;
	for (int i=0; i<m; i++){
		vector<int> v;
		for(int j=0; j<n; j++){
			v.push_back(0);
		}
		rep.push_back(v);
	}
	return rep;
}

// Returns a vector3D of the given size (filled only with zeroes)
//
vector3D initVector3D(int m, int n, int p){
	vector3D rep;
	for (int i=0; i<m; i++){
		vector2D vv = initVector2D(n,p);
		rep.push_back(vv);
	}
	return rep;
}

// To get the day from its id and vice-versa
// First day is always supposed to be a Monday
//
string intToDay(int dayId){
	if((dayId%7)==0) return "Mon";
	else if((dayId%7)==1) return "Tue";
	else if((dayId%7)==2) return "Wed";
	else if((dayId%7)==3) return "Thu";
	else if((dayId%7)==4) return "Fri";
	else if((dayId%7)==5) return "Sat";
	else return "Sun";
}
int dayToInt(string day){
	if(day=="Mon") return 0;
	if(day=="Tue") return 1;
	if(day=="Wed") return 2;
	if(day=="Thu") return 3;
	if(day=="Fri") return 4;
	if(day=="Sat") return 5;
	if(day=="Sun") return 6;
	else return -1;
}


// constructor of Timer
//
Timer::Timer():isInit_(0), isStarted_(0), isStopped_(0) {
	this->init();
}

// initialize the timer
//
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

// start the timer
//
void Timer::start()  {

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

		cpuSinceStart_.tv_sec   = 0;
		cpuSinceStart_.tv_nsec  = 0;
		isStarted_ = 1;
		isStopped_ = 0;

	}

	// Stop the time and update the times spent since the last start and since the
	// initialization
	//
	void Timer::stop() {

if ( isStopped_ )
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

	isStarted_ = 0;
	isStopped_ = 1;
	coStop_ ++;

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
		if ( cpuNow.tv_nsec-cpuInit_.tv_nsec < 0 ) {
			cpuTmp.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
			cpuTmp.tv_nsec  = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
		} else {
			cpuTmp.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec;
			cpuTmp.tv_nsec  = cpuNow.tv_nsec - cpuInit_.tv_nsec;
		}

		timespec cpuCurrent;
		if ( cpuTmp.tv_nsec + cpuSinceInit_.tv_nsec >= 1e09 ) {
			cpuCurrent.tv_sec  = cpuSinceInit_.tv_sec + cpuTmp.tv_sec + 1;
			cpuCurrent.tv_nsec = cpuSinceInit_.tv_nsec + cpuTmp.tv_nsec - 1e09;
		} else  {
			cpuCurrent.tv_sec  = cpuSinceInit_.tv_sec + cpuTmp.tv_sec;
			cpuCurrent.tv_nsec = cpuSinceInit_.tv_nsec + cpuTmp.tv_nsec;
		}

		return (double) cpuCurrent.tv_sec + (double)cpuCurrent.tv_nsec/1.0e09;
	}
	else if (isStopped_)  {
		return (double) cpuSinceInit_.tv_sec + (double)cpuSinceInit_.tv_nsec/1.0e09;
	}
	else
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

		if ( cpuNow.tv_nsec-cpuInit_.tv_nsec < 0 ) {
			cpuSinceStart_.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec - 1;
			cpuSinceStart_.tv_nsec  = 1e09 + cpuNow.tv_nsec - cpuInit_.tv_nsec;
		} else {
			cpuSinceStart_.tv_sec   = cpuNow.tv_sec - cpuInit_.tv_sec;
			cpuSinceStart_.tv_nsec  = cpuNow.tv_nsec - cpuInit_.tv_nsec;
		}
	}

	return (double) cpuSinceStart_.tv_sec + (double)cpuSinceStart_.tv_nsec/1.0e09;

} //end dSinceStart


} // end namespace
