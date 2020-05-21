//
// Created by antoine legrain on 2020-05-21.
//

#ifndef NURSESCHEDULER_ROTATIONSP_H
#define NURSESCHEDULER_ROTATIONSP_H

#include "SubProblem.h"


class RotationSP: public SubProblem {
  public:
    RotationSP() = default;

    RotationSP(PScenario scenario, int nbDays, PConstContract contract, std::vector<State>* pInitState);

    virtual ~RotationSP();

  protected:

    // Creates all nodes of the rcspp (including resource window)
    virtual void createNodes() override;

    // override creation of arcs source -> principal
    virtual void createArcsSourceToPrincipal() override;
    virtual void createArcsAllPriceLabels() override;
};


#endif //NURSESCHEDULER_ROTATIONSP_H
