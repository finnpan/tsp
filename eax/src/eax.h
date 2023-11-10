/*
License: see tsp.h
 */

#pragma once

#include "tsp.h"

/* TODO:
 * 1. Diversity preservation: Entropy
 * 2. Eset Type: Block2
 */

class Cross {
public:
    Cross (const Evaluator* e, int nPop);
    ~Cross ();
    void SetParents (const Indi& tPa1, const Indi& tPa2, int numKids);
    void DoIt (Indi& tKid, Indi& tPa2, int numKids, int flagP);

    int fNumGeneratedCh;

private:
    /* Step 2 of EAX */
    void SetABcycle (const Indi& parent1, const Indi& parent2, int numKids);
    void FormABcycle (); /* Store an AB-cycle found */
    /* Apply an AB-cycle to an intermediate solution */
    void ChangeSol (Indi& tKid, int ABnum, int type);
    void MakeCompleteSol (Indi& tKid); /* Step 5 of EAX */
    void MakeUnit (); /* Step 5-1 of EAX */
    void BackToPa1 (Indi& tKid); /* Undo the parent p_A */
    void GoToBest (Indi& tKid); /* Modify tKid to the best offspring solution */

    void CheckValid (Indi& indi); /* For debug */

private:
    int **_nearData;
    int *_koritsu, *_bunki, *_koriInv, *_bunInv;
    int _koritsuMany, _bunkiMany;
    int _stAppear;
    int *_checkKoritsu;
    int *_route;
    int _flagSt, _flagCircle, _prType;
    int **_abCycle;
    int *_permu;
    int _numABcycle;
    int _posiCurr;

    int *_c;

    // Speed Up Start
    int *_order;
    int *_inv;
    int **_segment;
    int *_segUnit;

    int _numUnit;
    int _numSeg;
    int *_segPosiList;
    int _numSPL;
    int *_linkAPosi;
    int **_linkBPosi;
    int *_posiSeg;
    int *_numElementInUnit;
    int *_centerUnit;
    int _numElementInCU;
    int *_listCenterUnit;

    EvalType *_gainAB;
    EvalType _gainModi;
    int _numModiEdge;
    int _numBestModiEdge;
    int **_modiEdge;
    int **_bestModiEdge;
    int _numAppliedCycle;
    int _numBestAppliedCycle;
    int *_appliedCylce;
    int *_bestAppliedCylce;
    // Speed Up End

    int _numABcycleInEset;
    int *_abCycleInEset;

    const int _numCity;
    const int _maxNumABcycle;
    const Evaluator* const _eval;
    const int _numPop;
};

class EAX {
public:
    EAX ();
    ~EAX ();
    void Define ();

    void DoIt ();
    void Init ();
    bool TerminationCondition ();
    void SetAverageBest ();
    void InitPop ();
    void SelectForMating ();
    void GenerateKids (int s);

    Evaluator* const _eval;     /* Distance of the edges */
    int* _indexForMating;       /* Mating list (r[] in the paper) */
    Indi* _curPop;             /* Current population members */
    Kopt* _kopt;               /* Local search with the 2-opt neighborhood */
    Cross* _cross;             /* Eede assembly crossover */

    const char *_fileNameTSP;  /* File name of an TSP instance */
    bool _silent{false};

    int _numPop;       /* Number of population members (N_pop in the paper) */
    int _numKids;      /* Number of offspring solutions (N_ch in the paper) */

    Indi _best;                /* Best solution in the current population */
    int _curNumGen;          /* The current number of generations */
    long int _accumurateNumCh;
    int _stagBest;
    double _averageValue;
    EvalType _bestValue;
    int _bestIndex;
};
