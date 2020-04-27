//
// Created by antoine legrain on 2020-04-19.
//

#ifndef NURSESCHEDULER_ROSTERSP_H
#define NURSESCHEDULER_ROSTERSP_H

#include "SubProblem.h"

class RosterSP : public SubProblem {
  public:

    RosterSP();
    virtual ~RosterSP();

    // Constructor that correctly sets the resource (time + bounds), but NOT THE COST
    //
    RosterSP(Scenario* scenario, int nbDays, const Contract* contract, std::vector<State>* pInitState);

//    double startWorkCost(int a) const override;


  protected:

    //----------------------------------------------------------------
    //
    // Update of the costs / network for solve function
    //
    //----------------------------------------------------------------

    // FUNCTIONS -- SOLVE
    //

    // return a function that will post process any path found by the RC graph
    std::function<void (spp_res_cont&)> postProcessResCont() const override ;

    // Creates all nodes of the rcspp (including resource window)
    void createNodes() override;

    // override creation of arcs source -> principal
    void createArcsSourceToPrincipal() override;
    void createArcsAllPriceLabels() override;
};


#endif //NURSESCHEDULER_ROSTERSP_H
