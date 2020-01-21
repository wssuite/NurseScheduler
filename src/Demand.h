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

/* namespace usage */
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
	Demand(): name_(""), nbDays_(0), firstDay_(0), nbShifts_(0), nbSkills_(0){}
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

public:

  // total demand in the minimal and optimal demands
  //
  int minTotal_, optTotal_;

  // preprocessed attributes aggregating the information of the demand
  //
  bool isPreprocessed_;

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

protected:

  // modify the demand by randomly swapping the demand of nnSwaps days
  //
  void swapDays(int nbSwaps);

  // modify the demand by randomly swapping the demand of nbSwaps shifts
  // the swapped shifts necessarily correspond to the same skill
  //
  void swapShifts(int nbSwaps);

  // perturb the demand by adding demand in a number of shifts randomly chosen
  // in the interval [minPerturb,maxPerturb]
  // the perturbed shifts are also randomly chosen
  // for a given skill the demand on a shift cannot become greater than the
  // largest demand observed on the week
  //
  void perturbShifts(int minPerturb, int maxPerturb);

public:

  // Index of the last day covered by the demand
  int lastDay() {return firstDay_+nbDays_-1;}

  // compute all the potentially helpful attributes of a demand
  // this includes the total demand per skill, per shift,
  void preprocessDemand();

  // add another week demand at the end of the current one
  // update all the parameters
  void push_back(Demand* pDemand);

  // Returns a new demand that appends pDemand to the current one
  Demand* append(Demand* pDemand);

  // display the demand, and include the preprocessed information if the input
  // boolean is set to true
  //
  string toString(bool withPreprocessedInfo);

  // copy the input demand and apply a perturbation to generate random demand
  //
  Demand* randomPerturbation();

  // shorten the demand by keeping only the nbDays in [begin, end)
  // return a new demand
  Demand* keep(int begin, int end);

  // shorten the demand by keeping only the nbDays first days
  //
  void keepFirstNDays(int nbDays);

  // shorten the demand by removing the nbDays first days
  //
  void removeFirstNDays(int nbDays);

	// remove a list of skills from the demand
	//
	void removeSkills(vector<int> skills);

};

#endif
