/*
* Greedy.cpp
*
*  Created on: 21 jan. 2015
*      Author: jeremy
*/


#include "Greedy.h"

//-----------------------------------------------------------------------------
//
//  C l a s s   G r e e d y
//
//  Quick solution of the problem with a greedy
//
//-----------------------------------------------------------------------------


// Main method to solve the rostering problem for a given input
//
void Greedy::solve() {}


// Build the sequence of positions reflecting the order in which the positions
// will be treated in the greedy
//
vector<Position*> Greedy::sortPositions() {

  vector<Position*> pPositions = pScenario_->pPositions();

  //---------------------------------------------------------------------------
  // A rank has been computed for each position
  // The greedy is going to treat the positions in ascending rank value
  // For a given rank, start with the position that contains the rarest skills
  // Rarity is defined by the number of nurse shifts available for this skill
  //---------------------------------------------------------------------------

  // initialize the vector of position order
  //
  vector<int> order;
  for (int i= 0; i < pScenario_->nbPositions(); i++)  order.push_back(0);

  // compute the maximum value of rarity among the skills of each position
  //
  for (int i= 0; i < pScenario_->nbPositions(); i++)  {
    int rarity = 0;
    // compute the highest rarity for the skills of the position
    // it has to be done in the preprocessing function of Solver, because it
    // depends on the history
    //

    // if rarity is higher
  }

  // get the maximum rank
  //
  int rankMax;
  for (int i= 0; i < pScenario_->nbPositions(); i++)  {
    rankMax = std::max(rankMax,pPositions[i]->rank());
  }

  // go through the positions for each rank value and set their order of
  // treatment.
  // the positions are already sorted in descending rarity, so they order is
  // set directly
  //
  int nextTreated = 0;
  for (int rank = 0; rank <= rankMax; rank++) {
    for (int i= 0; i < pScenario_->nbPositions(); i++)  {
      if (pPositions[i]->rank() == rank)  {
        order[pPositions[i]->id_] = nextTreated++;
      }
    }
  }



}
