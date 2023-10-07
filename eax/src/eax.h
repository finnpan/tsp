/*
Original License:
  MIT License

  Copyright (c) 2021 Yuichi Nagata

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

Repo:
  https://github.com/nagata-yuichi/GA-EAX

Commit:
  015dfbe9f267230f78787bd244af393ffc018900
 */

#pragma once

#include "indi.h"

class Cross {
public:
    Cross (int N);
    ~Cross ();
    void SetParents (const Indi& tPa1, const Indi& tPa2, int numOfKids);
    void DoIt (Indi& tKid, Indi& tPa2, int numOfKids, int flagP);

    int fNumOfGeneratedCh;
    const Evaluator* eval;
    int _numOfPop;

private:
    /* Step 2 of EAX */
    void SetABcycle (const Indi& parent1, const Indi& parent2, int numOfKids);
    void FormABcycle (); /* Store an AB-cycle found */
    /* Apply an AB-cycle to an intermediate solution */
    void ChangeSol (Indi& tKid, int ABnum, int type);
    void MakeCompleteSol (Indi& tKid); /* Step 5 of EAX */
    void MakeUnit (); /* Step 5-1 of EAX */
    void BackToPa1 (Indi& tKid); /* Undo the parent p_A */
    void GoToBest (Indi& tKid); /* Modify tKid to the best offspring solution */

    void CheckValid (Indi& indi); /* For debug */

private:
    int _n;
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
    int _maxNumABcycle;

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
    int *_listOfCenterUnit;

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
};

class Kopt;
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

    Evaluator* _eval;     /* Distance of the edges */
    Cross* _cross;             /* Eede assembly crossover */
    Kopt* _kopt;               /* Local search with the 2-opt neighborhood */
    const char *_fileNameTSP;  /* File name of an TSP instance */

    int _numOfPop;       /* Number of population members (N_pop in the paper) */
    int _numOfKids;      /* Number of offspring solutions (N_ch in the paper) */
    Indi* _curPop;             /* Current population members */
    Indi _best;                /* Best solution in the current population */
    int _curNumOfGen;          /* The current number of generations */
    long int _accumurateNumCh;

    int _stagBest;
    double _averageValue;
    EvalType _bestValue;
    int _bestIndex;
    int* _indexForMating;       /* Mating list (r[] in the paper) */
};
