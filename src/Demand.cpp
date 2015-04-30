#include "MyTools.h"
#include "Demand.h"
#include "Scenario.h"

#include <sstream>
#include <math.h>


//-----------------------------------------------------------------------------
//
//	C l a s s  D e m a n d
//
// All the information relative to a particular demand
//
//-----------------------------------------------------------------------------

// constructor and destructor
//
Demand::Demand(int nbDays, int firstDay, int nbShifts, int nbSkills, std::string name,
   vector3D minDemand, vector3D optDemand): name_(name),
   nbDays_(nbDays), firstDay_(firstDay), nbShifts_(nbShifts), nbSkills_(nbSkills),
   minDemand_(minDemand), optDemand_(optDemand),
   minTotal_(0), optTotal_(0), isPreprocessed_(false)
{
   // run the preprocessing
   this->preprocessDemand();
}

Demand::~Demand()
{}

// compute all the potentially helpful attributes of a demand
// this includes the total demand per skill, per shift,
//
void Demand::preprocessDemand() {
   // initialize the preprocessed vectors
   Tools::initVector(&minPerDay_, nbDays_);
   Tools::initVector(&optPerDay_, nbDays_);
   Tools::initVector(&minPerShift_, nbShifts_);
   Tools::initVector(&optPerShift_, nbShifts_);
   Tools::initVector(&minPerSkill_, nbSkills_);
   Tools::initVector(&optPerSkill_, nbSkills_);
   Tools::initVector(&minHighestPerSkill_, nbSkills_);
   Tools::initVector(&optHighestPerSkill_, nbSkills_);

   for (int day = 0; day < nbDays_; day++)	{
      for (int shift = 1; shift < nbShifts_; shift++) {
         for (int skill = 0; skill < nbSkills_; skill++)	{
            // update the total demand
            minTotal_ += minDemand_[day][shift][skill];
            optTotal_ += optDemand_[day][shift][skill];

            // update the demand per day
            minPerDay_[day] += minDemand_[day][shift][skill];
            optPerDay_[day] += optDemand_[day][shift][skill];

            // update the demand per shift
            minPerShift_[shift] += minDemand_[day][shift][skill];
            optPerShift_[shift] += optDemand_[day][shift][skill];

            // update the demand per skill
            minPerSkill_[skill] += minDemand_[day][shift][skill];
            optPerSkill_[skill] += optDemand_[day][shift][skill];

            // update the demand per day
            minHighestPerSkill_[skill] =
               std::max(minDemand_[day][shift][skill],minHighestPerSkill_[skill]);
            optHighestPerSkill_[skill] =
               std::max(optDemand_[day][shift][skill],optHighestPerSkill_[skill]);
         }
      }
   }
   isPreprocessed_ = true;
}

void Demand::push_back(Demand* pDemand){
   // check if same scenario
   if( (nbShifts_ != pDemand->nbShifts_) || (nbSkills_ != pDemand->nbSkills_) ){
      string error = "Demands are not compatible";
      Tools::throwError(error.c_str());
   }

   /*
    * Build new demand
    */

   // number of days covered by the demand and index of the first day
   //
   nbDays_ += pDemand->nbDays_;

   //pushes back the second demand on the first
   for(vector2D vector: pDemand->minDemand_)
      minDemand_.push_back(vector);
   for(vector2D vector: pDemand->optDemand_)
      optDemand_.push_back(vector);

   // run the preprocessing
   this->preprocessDemand();
}

// modify the demand by randomly swapping the demand of nnSwaps days
//
void Demand::swapDays(int nbSwaps) {

  for (int i=0; i < nbSwaps; i++) {
    int day1 = rand()%nbDays_;
    int day2 = rand()%nbDays_;

    // save the demand on day 1
    vector2D minDemandTmp, optDemandTmp;
    Tools::initVector2D(&minDemandTmp,nbShifts_,nbSkills_,0);
    Tools::initVector2D(&optDemandTmp,nbShifts_,nbSkills_,0);

    for (int sh = 1; sh < nbShifts_; sh++) {
      for (int sk = 0; sk < nbSkills_; sk++) {
        minDemandTmp[sh][sk] = minDemand_[day1][sh][sk];
        optDemandTmp[sh][sk] = optDemand_[day1][sh][sk];
      }
    }

    // make the modification in the demand
    for (int sh = 1; sh < nbShifts_; sh++) {
      for (int sk = 0; sk < nbSkills_; sk++) {
        minDemand_[day1][sh][sk] = minDemand_[day2][sh][sk];
        optDemand_[day1][sh][sk] = optDemand_[day2][sh][sk];

        minDemand_[day2][sh][sk] = minDemandTmp[sh][sk];
        optDemand_[day2][sh][sk] = optDemandTmp[sh][sk];
      }
    }
  }
}

// modify the demand by randomly swapping the demand of nbSwaps shifts
// the swapped shifts necessarily correspond to the same skill
//
void Demand::swapShifts(int nbSwaps) {
  for (int i=0; i < nbSwaps; i++) {
    int sk = rand()%nbSkills_;
    int day1 = rand()%nbDays_;
    int sh1 = 1+rand()%(nbShifts_-1); // make sure shift 0 is not taken
    int day2 = rand()%nbDays_;
    int sh2 = 1+rand()%(nbShifts_-1); // make sure shift 0 is not taken

    // save the demand on day1/shift1
    int minDemandTmp, optDemandTmp;
    minDemandTmp = minDemand_[day1][sh1][sk];
    optDemandTmp = optDemand_[day1][sh1][sk];

    // make the modification in the demand
    minDemand_[day1][sh1][sk] = minDemand_[day2][sh2][sk];
    optDemand_[day1][sh1][sk] = optDemand_[day2][sh2][sk];

    minDemand_[day2][sh2][sk] = minDemandTmp;
    optDemand_[day2][sh2

    ][sk] = optDemandTmp;
  }
}

// perturb the demand by adding demand in a number of shifts randomly chosen
// in the interval [minPerturb,maxPerturb]
// if the generate number is negative, then shifts are removed
// the perturbed shifts are also randomly chosen
// for a given skill the demand on a shift cannot become greater than the
// largest demand observed on the week
//
void Demand::perturbShifts(int minPerturb, int maxPerturb) {
  // generate the number of perturbations
  int nbPerturb = rand()%(maxPerturb-minPerturb) + minPerturb;
  int valPerturb = (nbPerturb>=0)? 1:-1;
  nbPerturb = fabs(nbPerturb);

  // preprocess the demand to find the highest demand per skill
  if (!isPreprocessed_) this->preprocessDemand();

  int coTrials=0; // used to avoid infinite loop
  for (int i=0; i < fabs(nbPerturb); i++) {
    // draw the particular shift whose demand should be perturbed
    // select only a demand that is below the highest demand of the week for
    // this shift if demand should be added (this constraint is added to avoid
    // non feasibility due to the perturbation)
    bool isAtUpperBound = true;
    int day, sh, sk;
    while (isAtUpperBound && coTrials < 10*nbPerturb) {
      day = rand()%nbDays_;
      sh = 1+rand()%(nbShifts_-1); // make sure shift 0 is not taken
      sk = rand()%nbShifts_;
      isAtUpperBound = (valPerturb >= 0)? (minDemand_[day][sh][sk] >= minHighestPerSkill_[sk]):false;
      coTrials++;
    }
    if (coTrials >= 10*nbPerturb) break;

    // perturb the demand
    minDemand_[day][sh][sk] += valPerturb;
    optDemand_[day][sh][sk] += valPerturb;
  }
}

// copy the input demand and apply a perturbation to generate random demand
//
Demand* Demand::randomPerturbation() {
  Demand* pDemand = new Demand(*this);

  // three different types of perturbations are made
  // the order does not seem to be important
  pDemand->swapDays(nbDays_/2);
  pDemand->swapShifts(nbDays_*nbSkills_);
  pDemand->perturbShifts(-nbDays_,nbDays_);

  // get the main characteristics of the new demand
  pDemand->preprocessDemand();

  return pDemand;
}


// display the demand, and include the preprocessed information if the input
// boolean is set to true
//
string Demand::toString(bool withPreprocessedInfo) {

   std::stringstream rep;

   rep << "# " << std::endl;
   rep << "# DEMAND" << std::endl;

   // describe the demand being written
   //
   rep << "# " << std::endl;
   rep << "# Name of the demand: " << name_ << std::endl;

   rep << "# The demand refers to " << nbSkills_ << " skills for " ;
   rep << nbShifts_-1 << " shifts per day on "<< nbDays_ << " days" << std::endl;
   rep << std::endl;

   // write the number of staff required per shift for each skill
   //
   rep << "#\t\t\t\t\t";
   for(int dayId=0; dayId<7; dayId++){
      rep << " " << Tools::intToDay(dayId) << "\t";
   }
   rep << "# " << std::endl;
   for(int sh = 0; sh < nbShifts_; sh ++){
      for (int sk = 0; sk < nbSkills_; sk++){
         // string str = "#   " + Tools::intToShift_[sh] + " " + Scenario::intToSkill_[sk] + " ";
         // rep << str;
         // if(str.length() < 16) rep << "\t";
         rep <<  "#\tShift " <<  sh  <<  " Skill " <<  sk << " " << "\t";
         for(int day = 0; day < 7; day ++){
            rep << "\t(" << minDemand_[day][sh][sk] << "," << optDemand_[day][sh][sk] << ")";
         }
         rep << std::endl;
      }
      rep << "# " << std::endl;
   }

   // write the preprocessed indicators if preprocess was run
   //
   if (withPreprocessedInfo)	{

      if (!isPreprocessed_)
         Tools::throwError("Trying to write the preprocessed information of a demand that was not preprocessed!");
      // enumerate the global indicators
      //
      rep << "# " << std::endl;
      rep << "# PREPROCESSED DATA ON THE DEMAND" << std::endl;

      rep << "# Total minimum demand = " << minTotal_ << std::endl;
      rep << "# Total optimal demand = " << optTotal_ << std::endl;

      rep << "# " << std::endl;
      rep  << "# Demand per day" << std::endl;
      for (int i = 0; i < nbDays_; i++)	{
         rep << "#\t\t" << Tools::intToDay(i+firstDay_) << " (" << i+firstDay_ << ")" << ": ";
         rep << "minimum = " << minPerDay_[i] << " ; optimal = " << optPerDay_[i];
         rep << std::endl;
      }

      rep << "# " << std::endl;
      rep  << "# Demand per shift" << std::endl;
      rep << "#\t\tShift 0 is rest;" << std::endl;
      for (int i = 1; i < nbShifts_; i++)	{
         rep << "#\t\tShift " << i << ": ";
         rep << "minimum = " << minPerShift_[i] << " ; optimal = " << optPerShift_[i];
         rep << std::endl;
      }

      rep << "# " << std::endl;
      rep  << "# Demand per skill" << std::endl;
      for (int i = 0; i < nbSkills_; i++)	{
         rep << "#\t\tSkill " << i << ": ";
         rep << "minimum = " << minPerSkill_[i] << " ; optimal = " << optPerSkill_[i];
         rep << std::endl;
      }

      rep << "# " << std::endl;
      rep <<  "# Highest demand per skill for one shift" << std::endl;
      for (int i = 0; i < nbSkills_; i++)	{
         rep << "#\t\tSkill " << i << ": ";
         rep << "minimum = " << minHighestPerSkill_[i] << " ; optimal = " << optHighestPerSkill_[i];
         rep << std::endl;
      }
   }

   return rep.str();
}
