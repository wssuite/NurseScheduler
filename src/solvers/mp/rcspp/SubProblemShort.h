//
// Created by antoine legrain on 2020-04-12.
//

#ifndef NURSESCHEDULER_SUBPROBLEMSHORT_H
#define NURSESCHEDULER_SUBPROBLEMSHORT_H

#include "SubProblem.h"


class SubProblemShort : public SubProblem {

  public:

    SubProblemShort();
    virtual ~SubProblemShort();

    // Constructor that correctly sets the resource (time + bounds), but NOT THE COST
    //
    SubProblemShort(Scenario* scenario, int nbDays, const Contract* contract, std::vector<State>* pInitState);

    double startWorkCost(int a) const override;


  protected:

    // Number of rotations found (that match the bound condition) at that iteration
    //
    int nVeryShortFound_ = 0;

    //-----------------------
    // THE BASE COSTS
    //-----------------------

    // SHORT SUCCESSIONS (computed when creating them)
    vector2D<double> baseArcCostOfShortSucc_;										// For each size c \in [0,CDMin], for each short rotation of size c, contains its base cost (independent from the date)

    //-----------------------
    // THE SHORT SUCCESSIONS
    //-----------------------

    // SHORTSUCC -> OBJECTS
    //
    // Short successions (no starting date) -> those of all length
    vector3D<int> allowedShortSuccBySize_;														// For each size c \in [0,CDMin], contains all allowed short successions of that size (satisfies succession constraints)
    vector2D<int> lastShiftOfShortSucc_;															// For each size c \in [0,CDMin], for each short rotation of size c, contains the corresponding last shift performed
    vector2D<int> nLastShiftOfShortSucc_;														// For each size c \in [0,CDMin], for each short rotation of size c, contains the number of consecutive days the last shift has been performed
    // Objects for short successions of maximal size CDMin
    //	int CDMin_;																				// Minimum number of consecutive days worked for free

    vector3D<int> allShortSuccCDMinByLastShiftCons_;												// For each shift s, for each number of days n, contains the list of short successions of size CDMin ending with n consecutive days of shift s

    // SHORTSUCC -> FUNCTIONS
    //
    // Initializes all short successions, base costs, and corresponding vectors. Should only be called ONCE.
    void initShortSuccessions();

    //----------------------------------------------------------------
    //
    // Cost computation of the "very" short rotations (< CD_min)
    //
    //----------------------------------------------------------------
    bool priceVeryShortRotationsFirstDay();
    bool priceVeryShortRotationsLastDay();
    bool priceVeryShortRotations();
    double costOfVeryShortRotation(int firstDay, const std::vector<int>& succ);


    //----------------------------------------------------------------
    //
    // Update of the costs / network for solve function
    //
    //----------------------------------------------------------------

    // FUNCTIONS -- SOLVE
    //
    bool preprocess() override;
    bool solveShortRotations();

    // override creation of arcs source -> principal
    void createArcsSourceToPrincipal() override;

    // FUNCTIONS -- COSTS
    //
    // Pricing of the short successions : only keep one of them, and the cost of the corresponding arc
    void priceShortSucc();
    // Given a short succession and a start date, returns the cost of the corresponding arc
    double costArcShortSucc(int size, int id, int startDate);


  public:

    // // Some getters
    // //
    inline int nVeryShortFound() const {return nVeryShortFound_;}

    // Print functions.
    //
    void printShortSucc() const;
    void printShortArcs() const;

};


#endif //NURSESCHEDULER_SUBPROBLEMSHORT_H
