/*
License: see util.h
 */

#pragma once

#include "util.h"

/* TODO:
 * 1. Diversity preservation: Entropy
 * 2. E-sets Type: Block2
 */

class Cross {
public:
    Cross(const Evaluator* e, int nPop);
    ~Cross();

    void DoIt(Indi& kid, Indi& pa2, int nKids);

private:
    void BuildABcycle(const Indi& pa1, const Indi& pa2, int nKids);
    void BuildABcycle_0(int stAppear, int& posiCurr);
    void ChangeSol(Indi& kid, int idx, bool reverse, bool updateSeg = true);
    void MakeUnit();
    EvalType MakeCompleteSol(Indi& kid);
    void BackToPa1(Indi& kid, int appliedCycle);
    void GoToBest(Indi& kid, int bestAppliedCycle);

private:
    const Evaluator* const _eval;
    const int _numCity;
    const int _numPop;
    const int _maxNumABcycle;

    int *_pa1Route, *_pa1RouteInv;

    int _numABcycle;
    int** _ABcycleList;
    int* _permuABCycle;
    EvalType* _gainABcycle;

    int** _overlapEdges;
    int *_cycBuf1, *_cycBuf1Inv;
    int *_cycBuf2, *_cycBuf2Inv;
    int _cycBuf1Num, _cycBuf2Num;
    int* _cycRoute;
    int* _ABCycle;

    int _numModiEdge;
    int** _modiEdge;
    int _numBestModiEdge;
    int** _bestModiEdge;

    int** _segment;
    int* _segUnit;
    int _numUnit;
    int _numSeg;
    int* _segPosiList;
    int _numSPL;
    int* _linkAPosi;
    int** _linkBPosi;
    int* _posiSeg;
    int* _numElementInUnit;
    int* _centerUnit;
    int* _listCenterUnit;
};

class EAXGA {
public:
    EAXGA(int nPop = 100, int nKid = 30);
    ~EAXGA();

    bool Define(const char* tspFileName);
    bool DoIt();

    const int _numPop;
    const int _numKids;
    bool _silent;

    bool _failed;
    int _numGen;
    Indi _best;
    double _avgCost;

private:
    void SelectBest();
    bool ShouldTerminate();
    void SelectForMating();

    Evaluator* const _eval;
    int* const _matingSeq;
    KOpt* _kopt;
    Cross* _cross;
    Indi* _pop;
    int _stagnGen;
};
