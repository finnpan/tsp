/*
License: see tsp.h
 */

#pragma once

#include "tsp.h"

/* TODO:
 * 1. Diversity preservation: Entropy
 * 2. E-sets Type: Block2
 */

class Cross {
public:
    Cross (const Evaluator* e, int nPop);
    ~Cross ();

    void DoIt (Indi& kid, Indi& pa2, int nKids);

private:
    void SetABcycle (const Indi& pa1, const Indi& pa2, int nKids);
    void FormABcycle ();
    void ChangeSol (Indi& kid, int idx, int type);
    EvalType MakeCompleteSol (Indi& kid);
    void MakeUnit ();
    void BackToPa1 (Indi& kid);
    void GoToBest (Indi& kid);

private:
    const Evaluator* const _eval;
    const int _numCity;
    const int _numPop;
    const int _maxNumABcycle;

//////////////////////////////////////////////////////////////
    int       _numABcycle;
    int**     _ABcycle;
    int*      _permuABCycle;
    EvalType* _gainABcycle;

    int       _numABcycleInEset;
    int*      _ABCycleInEset;
    int       _numAppliedCycle;
    int*      _appliedCycle;
    int       _numBestAppliedCycle;
    int*      _bestAppliedCycle;

    int **_nearData;
    int *_koritsu, *_bunki, *_koriInv, *_bunInv;
    int _koritsuMany, _bunkiMany;
    int _stAppear;
    int *_route;
    int _flagSt, _flagCircle, _prType;
    int _posiCurr;
    int *_cycle;
//////////////////////////////////////////////////////////////
    int       _numModiEdge;
    int**     _modiEdge;
    int       _numBestModiEdge;
    int**     _bestModiEdge;

    int*      _path;
    int*      _posi;

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
};

class EAXGA {
public:
    EAXGA ();
    ~EAXGA ();

    void Define (const char* tspFileName);
    void DoIt ();

    int    _numPop;
    int    _numKids;
    bool   _silent;

    int    _numGen;
    Indi   _best;
    double _avgCost;

private:
    void SelectBest ();
    bool ShouldTerminate ();
    void SelectForMating ();

    Evaluator* const _eval;
    KOpt*            _kopt;
    Cross*           _cross;
    Indi*            _pop;
    int*             _matingSeq;
    int              _stagnGen;
};
