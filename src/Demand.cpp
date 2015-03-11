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
      for (int shift = 0; shift < nbShifts_; shift++) {
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
