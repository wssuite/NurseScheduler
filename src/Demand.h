//
//  Demand.h
//  RosterDesNurses
//

#ifndef __Demand__
#define __Demand__

#include <map>
#include <string>
#include <vector>

#include "MyTools.h"

using std::map;
using std::pair;
using std::string;
using std::vector;


//-----------------------------------------------------------------------------
//
//	C l a s s  D e m a n d
//
// All the information relative to a particular demand
//
//-----------------------------------------------------------------------------

class Demand {

public:

  // generic constructor and destructor
  Demand(int nbDays, int firstDay, int nbShifts, int nbSkills, std::string name,
  vector3D minDemand, vector3D optDemand);
  ~Demand();

  // constant attributes of the demand
  //
public:

  // name of the demand
  //
  std::string name_;

  // number of days covered by the demand and index of the first day
  //
  int nbDays_;
  const int firstDay_;

  // number of shifts per day and number of skills to cover
  const int nbShifts_, nbSkills_;

  // minimum and optimal demand for each day, shift and skill
  //
  vector3D minDemand_;
  vector3D optDemand_;

  // preprocessed attributes aggregating the information of the demand
  //
public:
  // total demand in the minimal and optimal demands
  //
  int minTotal_, optTotal_;

  // total demand per skill in the minimal and optimal demands
  //
  vector<int> minPerSkill_, optPerSkill_;

  // total demand per shift in the minimal and optimal demands
  //
  vector<int> minPerShift_, optPerShift_;

  // total demand per day in the minimal and optimal demands
  //
  vector<int> minPerDay_, optPerDay_;

  // highest demands per skill over the considered period
  //
  vector<int> minHighestPerSkill_, optHighestPerSkill_;

private:
  bool isPreprocessed_;

public:

  // compute all the potentially helpful attributes of a demand
  // this includes the total demand per skill, per shift,
  void preprocessDemand();

  // add another week demand at the end of the current one and create a new one
  // update all the parameters
  // return the new demand
  void push_back(Demand* pDemand);

  // display the demand, and include the preprocessed information if the input
  // boolean is set to true
  //
  string toString(bool withPreprocessedInfo);

};

#endif
