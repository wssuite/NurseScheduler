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

  // First, organize the positions in ranks
  // rank i contains all the positions that are dominated only by positions
  // with a rank smaller than i
  //


}
