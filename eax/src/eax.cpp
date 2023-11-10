/*
License: see tsp.h
 */

#include "eax.h"

EAX::EAX () :
    _eval(new Evaluator()),
    _indexForMating(nullptr),
    _curPop(nullptr),
    _kopt(nullptr),
    _cross(nullptr)
{
    _numPop = 100;
    _numKids = 30;
}

EAX::~EAX ()
{
    delete [] _indexForMating;
    delete [] _curPop;
    delete _kopt;
    delete _cross;
    delete _eval;
}

void EAX::Define ()
{
    _eval->SetInstance(_fileNameTSP);

    const int n = _eval->_numCity;
    _best.Define(n);

    _indexForMating = new int [_numPop + 1];

    _curPop = new Indi [_numPop];
    for (int i = 0; i < _numPop; ++i) {
        _curPop[i].Define(n);
    }

    _kopt = new Kopt(_eval);
    _cross = new Cross(_eval, _numPop);
}

void EAX::DoIt ()
{
    InitPop();
    Init();

    while (1) {
        SetAverageBest();
        if (!_silent) {
            printf("%d: %lld %lf\n",
                   _curNumGen,
                   (long long)_bestValue,
                   _averageValue);
        }

        if (TerminationCondition()) break;

        SelectForMating();

        for (int s =0; s < _numPop; ++s) {
            GenerateKids(s);
        }
        ++_curNumGen;
    }
}

void EAX::Init ()
{
    _accumurateNumCh = 0;
    _curNumGen = 0;
    _stagBest = 0;
}

bool EAX::TerminationCondition ()
{
    if (_averageValue - _bestValue < 0.001) {
        return true;
    }

    if (_stagBest > int(1500 / _numKids)) {
        return true;
    }

    return false;
}

void EAX::SetAverageBest ()
{
    EvalType stockBest = _best._cost;

    _averageValue = 0.0;
    _bestIndex = 0;
    _bestValue = _curPop[0]._cost;

    for (int i = 0; i < _numPop; ++i) {
        _averageValue += _curPop[i]._cost;
        if (_curPop[i]._cost < _bestValue) {
            _bestIndex = i;
            _bestValue = _curPop[i]._cost;
        }
    }

    _best = _curPop[_bestIndex];
    _averageValue /= (double)_numPop;

    if (_best._cost < stockBest) {
        _stagBest = 0;
    } else {
        _stagBest++;
    }
}

void EAX::InitPop ()
{
    for (int i = 0; i < _numPop; ++i) {
        _kopt->DoIt(_curPop[i]);
    }
}

void EAX::SelectForMating ()
{
    for (int i = 0; i < _numPop; i++) {
        _indexForMating[i] = i;
    }
    std::shuffle(_indexForMating, _indexForMating+_numPop, *_eval->_rand);
    _indexForMating[_numPop] = _indexForMating[0];
}

void EAX::GenerateKids (int s)
{
    _cross->SetParents(_curPop[_indexForMating[s]],
                       _curPop[_indexForMating[s+1]], _numKids);

    _cross->DoIt(_curPop[_indexForMating[s]],
                 _curPop[_indexForMating[s+1]], _numKids, 1);

    _accumurateNumCh += _cross->fNumGeneratedCh;
}

Cross::Cross (const Evaluator* e, int nPop) :
    _numCity(e->_numCity),
    _maxNumABcycle(2000),
    _eval(e),
    _numPop(nPop)
{
    /* Set an appropriate value (2000 is usually enough) */

    const int n = _numCity;
    _nearData = new int* [n];
    for (int j = 0; j < n; ++j) {
        _nearData[j] = new int [5];
    }

    _abCycle = new int* [_maxNumABcycle];
    for (int j = 0; j < _maxNumABcycle; ++j) {
        _abCycle[j] = new int [2*n + 4];
    }

    _koritsu = new int [n];
    _bunki = new int [n];
    _koriInv = new int [n];
    _bunInv = new int [n];
    _checkKoritsu = new int [n];
    _route = new int [2*n + 1];
    _permu = new int [_maxNumABcycle];

    _c = new int [2*n+4];

    // Speed Up Start
    _order = new int [n];
    _inv = new int [n];
    _segment = new int* [n];
    for (int j = 0; j < n; ++j) {
        _segment[j] = new int [2];
    }
    _segUnit = new int [n];
    _segPosiList = new int[n];
    _linkAPosi = new int [n];
    _linkBPosi = new int* [n];
    for (int j = 0; j < n; ++j) {
        _linkBPosi[j] = new int [2];
    }
    _posiSeg = new int [n];
    _numElementInUnit = new int [n];
    _centerUnit = new int [n];
    for (int j = 0; j < n; ++j) {
        _centerUnit[j] = 0;
    }
    _listCenterUnit = new int [n+2];
    _gainAB = new EvalType [n];
    _modiEdge = new int* [n];
    for (int j = 0; j < n; ++j) {
        _modiEdge[j] = new int [4];
    }
    _bestModiEdge = new int* [n];
    for (int j = 0; j < n; ++j) {
        _bestModiEdge[j] = new int [4];
    }
    _appliedCylce = new int [n];
    _bestAppliedCylce = new int [n];
    // Speed Up End

    _abCycleInEset = new int [_maxNumABcycle];
}

Cross::~Cross ()
{
    const int n = _numCity;
    delete [] _koritsu;
    delete [] _bunki;
    delete [] _koriInv;
    delete [] _bunInv;
    delete [] _checkKoritsu;
    delete [] _route;
    delete [] _permu;

    for (int j = 0; j < n; ++j) {
        delete[] _nearData[j];
    }
    delete[] _nearData;

    for (int j = 0; j < _maxNumABcycle; ++j) {
        delete[] _abCycle[j];
    }
    delete[] _abCycle;

    delete [] _c;

    // Speed Up Start
    delete [] _order;
    delete [] _inv;

    for (int j = 0; j < n; ++j) {
        delete[] _segment[j];
    }
    delete[] _segment;
    delete[] _segUnit;
    delete [] _segPosiList;
    delete [] _linkAPosi;
    for (int j = 0; j < n; ++j) {
        delete[] _linkBPosi[j];
    }
    delete [] _linkBPosi;
    delete [] _posiSeg;
    delete [] _numElementInUnit;
    delete [] _centerUnit;
    delete [] _listCenterUnit;
    delete [] _gainAB;

    for (int j = 0; j < n; ++j) {
        delete[] _modiEdge[j];
    }
    delete [] _modiEdge;
    for (int j = 0; j < n; ++j) {
        delete[] _bestModiEdge[j];
    }
    delete [] _bestModiEdge;
    delete [] _appliedCylce;
    delete [] _bestAppliedCylce;
    // Speed Up End

    delete [] _abCycleInEset;
}

void Cross::SetParents (const Indi& tPa1, const Indi& tPa2, int numKids)
{
    SetABcycle(tPa1, tPa2, numKids);

    int curr, next, st, pre;
    st = 0;
    curr=-1;
    next = st;
    for (int i = 0; i < _numCity; ++i) {
        pre=curr;
        curr=next;
        if (tPa1._link[curr][0] != pre) {
            next = tPa1._link[curr][0];
        } else {
            next=tPa1._link[curr][1];
        }

        _order[i] = curr;
        _inv[curr] = i;
    }

    assert(next == st);
}

void Cross::DoIt (Indi& tKid, Indi& tPa2, int numKids, int flagP)
{
    int Num;
    int jnum;
    EvalType gain;
    EvalType BestGain;
    double pointMax, point;

    if (numKids <= _numABcycle) {
        Num = numKids;
    } else {
        Num = _numABcycle;
    }

    for (int i = 0; i < _numABcycle; i++) {
        _permu[i] = i;
    }
    std::shuffle(_permu, _permu+_numABcycle, *_eval->_rand);

    fNumGeneratedCh = 0;
    pointMax = 0.0;
    BestGain = 0;
    int flagImp = 0;

    for (int j =0; j < Num; ++j) {
        _numABcycleInEset = 0;
        jnum = _permu[j];
        _abCycleInEset[_numABcycleInEset++] = jnum;

        _numSPL = 0;
        gain = 0;
        _numAppliedCycle = 0;
        _numModiEdge = 0;

        _numAppliedCycle = _numABcycleInEset;
        for (int k = 0; k < _numAppliedCycle; ++k) {
            _appliedCylce[k] = _abCycleInEset[k];
            jnum = _appliedCylce[k];
            ChangeSol(tKid, jnum, flagP);
            gain += _gainAB[jnum];
        }

        MakeUnit();
        MakeCompleteSol(tKid);
        gain += _gainModi;

        ++fNumGeneratedCh;

        point = (double)gain;
        tKid._cost = tKid._cost - gain;

        // if (pointMax < point) {
        if (pointMax < point &&
            tKid._cost != tPa2._cost) {
            pointMax = point;
            BestGain = gain;
            flagImp = 1;

            _numBestAppliedCycle = _numAppliedCycle;
            for (int s = 0; s < _numBestAppliedCycle; ++s) {
                _bestAppliedCylce[s] = _appliedCylce[s];
            }

            _numBestModiEdge = _numModiEdge;
            for (int s = 0; s < _numBestModiEdge; ++s) {
                _bestModiEdge[s][0] = _modiEdge[s][0];
                _bestModiEdge[s][1] = _modiEdge[s][1];
                _bestModiEdge[s][2] = _modiEdge[s][2];
                _bestModiEdge[s][3] = _modiEdge[s][3];
            }
        }

        BackToPa1(tKid);
        tKid._cost = tKid._cost + gain;
    }

    if (flagImp == 1) {
        GoToBest(tKid);
        tKid._cost = tKid._cost - BestGain;
    }
}

void Cross::SetABcycle (const Indi& tPa1, const Indi& tPa2, int numKids)
{
    const int n = _numCity;
    _bunkiMany=0; _koritsuMany=0;
    for (int j = 0; j < n ; ++j) {
        _nearData[j][1]=tPa1._link[j][0];
        _nearData[j][3]=tPa1._link[j][1];

        _nearData[j][0] = 2;

        _koritsu[_koritsuMany]=j;
        _koritsuMany++;

        _nearData[j][2]=tPa2._link[j][0];
        _nearData[j][4]=tPa2._link[j][1];
    }
    for (int j = 0; j < n; ++j) {
        _checkKoritsu[j]=-1;
        _koriInv[_koritsu[j]]=j;
    }

    /**************************************************/

    _numABcycle=0;
    _flagSt=1;
    int r = 0, pr = 0, st = 0, ci = 0;
    while (_koritsuMany!=0) {
        if (_flagSt==1) {
            _posiCurr=0;
            r=(*_eval->_rand)()%_koritsuMany;
            st=_koritsu[r];
            _checkKoritsu[st]=_posiCurr;
            _route[_posiCurr]=st;
            ci=st;
            _prType=2;
        } else if (_flagSt==0) {
            ci=_route[_posiCurr];
        }

        _flagCircle=0;
        while (_flagCircle==0) {
            _posiCurr++;
            pr=ci;

            switch(_prType) {
                case 1:
                    ci=_nearData[pr][_posiCurr%2+1];
                    break;
                case 2:
                    r=(*_eval->_rand)()%2;
                    ci=_nearData[pr][_posiCurr%2+1+2*r];
                    if (r==0) {
                        std::swap(_nearData[pr][_posiCurr%2+1],
                                  _nearData[pr][_posiCurr%2+3]);
                    }
                    break;
                case 3:
                    ci=_nearData[pr][_posiCurr%2+3];
            }

            _route[_posiCurr]=ci;

            if (_nearData[ci][0]==2) {
                if (ci==st) {
                    if (_checkKoritsu[st]==0) {
                        if ((_posiCurr-_checkKoritsu[st])%2==0) {
                            if (_nearData[st][_posiCurr%2+1]==pr) {
                                std::swap(_nearData[ci][_posiCurr%2+1],
                                          _nearData[ci][_posiCurr%2+3]);
                            }
                            _stAppear = 1;
                            FormABcycle();
                            if (_numABcycle == numKids) { goto LLL; }
                            if (_numABcycle == _maxNumABcycle) { goto LLL; }

                            _flagSt=0;
                            _flagCircle=1;
                            _prType=1;
                        } else {
                            std::swap(_nearData[ci][_posiCurr%2+1],
                                      _nearData[ci][_posiCurr%2+3]);
                            _prType=2;
                        }
                        _checkKoritsu[st]=_posiCurr;
                    } else {
                        _stAppear = 2;
                        FormABcycle();
                        if (_numABcycle == numKids) { goto LLL; }
                        if (_numABcycle == _maxNumABcycle) { goto LLL; }

                        _flagSt=1;
                        _flagCircle=1;
                    }
                } else if (_checkKoritsu[ci]==-1) {
                    _checkKoritsu[ci]=_posiCurr;
                    if (_nearData[ci][_posiCurr%2+1]==pr) {
                        std::swap(_nearData[ci][_posiCurr%2+1],
                                  _nearData[ci][_posiCurr%2+3]);
                    }
                    _prType=2;
                } else if (_checkKoritsu[ci]>0) {
                    std::swap(_nearData[ci][_posiCurr%2+1],
                            _nearData[ci][_posiCurr%2+3]);
                    if ((_posiCurr-_checkKoritsu[ci])%2==0) {
                        _stAppear = 1;
                        FormABcycle();
                        if (_numABcycle == numKids) { goto LLL; }
                        if (_numABcycle == _maxNumABcycle) { goto LLL; }

                        _flagSt=0;
                        _flagCircle=1;
                        _prType=1;
                    } else {
                        std::swap(_nearData[ci][(_posiCurr+1)%2+1],
                                  _nearData[ci][(_posiCurr+1)%2+3]);
                        _prType=3;
                    }
                }
            } else if (_nearData[ci][0]==1) {
                if (ci==st) {
                    _stAppear = 1;
                    FormABcycle();
                    if (_numABcycle == numKids) { goto LLL; }
                    if (_numABcycle == _maxNumABcycle) { goto LLL; }

                    _flagSt=1;
                    _flagCircle=1;
                } else {
                    _prType=1;
                }
            }
        }
    }

    while (_bunkiMany!=0) {
        _posiCurr=0;
        r=(*_eval->_rand)()%_bunkiMany;
        st=_bunki[r];
        _route[_posiCurr]=st;
        ci=st;

        _flagCircle=0;
        while (_flagCircle==0) {
            pr=ci;
            _posiCurr++;
            ci=_nearData[pr][_posiCurr%2+1];
            _route[_posiCurr]=ci;
            if (ci==st) {
                _stAppear = 1;
                FormABcycle();
                if (_numABcycle == numKids) { goto LLL; }
                if (_numABcycle == _maxNumABcycle) { goto LLL; }
                _flagCircle=1;
            }
        }
    }

    LLL: ;

    if (_numABcycle == _maxNumABcycle) {
        printf("Error: _maxNumABcycle(%d) must be increased\n", _maxNumABcycle);
        exit(1);
    }
}

void Cross::FormABcycle ()
{
    int j;
    int st_count;
    int edge_type;
    int st, ci, stock;
    int cem;
    EvalType diff;

    if (_posiCurr%2==0) {
        edge_type=1;
    } else {
        edge_type=2;
    }
    st=_route[_posiCurr];
    cem=0;
    _c[cem]=st;

    st_count=0;
    while (1) {
        cem++;
        _posiCurr--;
        ci=_route[_posiCurr];
        if (_nearData[ci][0]==2) {
            _koritsu[_koriInv[ci]]=_koritsu[_koritsuMany-1];
            _koriInv[_koritsu[_koritsuMany-1]]=_koriInv[ci];
            _koritsuMany--;
            _bunki[_bunkiMany]=ci;
            _bunInv[ci]=_bunkiMany;
            _bunkiMany++;
        } else if (_nearData[ci][0]==1) {
            _bunki[_bunInv[ci]]=_bunki[_bunkiMany-1];
            _bunInv[_bunki[_bunkiMany-1]]=_bunInv[ci];
            _bunkiMany--;
        }

        _nearData[ci][0]--;
        if (ci==st) { st_count++; }
        if (st_count==_stAppear) { break; }
        _c[cem]=ci;
    }

    if (cem==2) {
        return;
    }

    _abCycle[_numABcycle][0]=cem;

    if (edge_type==2) {
        stock=_c[0];
        for (int j=0;j<cem-1;j++) { _c[j]=_c[j+1]; }
        _c[cem-1]=stock;
    }

    for (int j=0;j<cem;j++) {
        _abCycle[_numABcycle][j+2]=_c[j];
    }
    _abCycle[_numABcycle][1]=_c[cem-1];
    _abCycle[_numABcycle][cem+2]=_c[0];
    _abCycle[_numABcycle][cem+3]=_c[1];

    _c[cem] = _c[0];
    _c[cem+1] = _c[1];
    diff = 0;
    for (j = 0; j < cem/2; ++j) {
        diff += _eval->_cost[_c[2*j]][_c[1+2*j]] -
                _eval->_cost[_c[1+2*j]][_c[2+2*j]];
    }
    _gainAB[_numABcycle] = diff;
    ++_numABcycle;
}

void Cross::ChangeSol (Indi& tKid, int ABnum, int type)
{
    const int n = _numCity;
    int j;
    int cem,r1,r2,b1,b2;
    int po_r1, po_r2, po_b1, po_b2;

    cem=_abCycle[ABnum][0];
    _c[0]=_abCycle[ABnum][0];

    if (type==2) {
        for (j=0;j<cem+3;j++) {
            _c[cem+3-j]=_abCycle[ABnum][j+1];
        }
    } else {
        for (j=1;j<=cem+3;j++) {
            _c[j]=_abCycle[ABnum][j];
        }
    }

    for (j=0;j<cem/2;j++) {
        r1=_c[2+2*j];r2=_c[3+2*j];
        b1=_c[1+2*j];b2=_c[4+2*j];

        if (tKid._link[r1][0]==r2) {
            tKid._link[r1][0]=b1;
        } else {
            tKid._link[r1][1]=b1;
        }
        if (tKid._link[r2][0]==r1) {
            tKid._link[r2][0]=b2;
        } else {
            tKid._link[r2][1]=b2;
        }

        po_r1 = _inv[r1];
        po_r2 = _inv[r2];
        po_b1 = _inv[b1];
        po_b2 = _inv[b2];

        /* FIXME
         * Using eil101.tsp, sometime heap-buffer-overflow will occur.
         */
        if (po_r1 == 0 && po_r2 == n-1) {
            if (_numSPL >= _numCity) {
                printf("heap-buffer-overflow1: _numSPL = %d\n", _numSPL);
            }
            _segPosiList[_numSPL++] = po_r1;
        } else if (po_r1 == n-1 && po_r2 == 0) {
            if (_numSPL >= _numCity) {
                printf("heap-buffer-overflow2: _numSPL = %d\n", _numSPL);
            }
            _segPosiList[_numSPL++] = po_r2;
        } else if (po_r1 < po_r2) {
            if (_numSPL >= _numCity) {
                printf("heap-buffer-overflow3: _numSPL = %d\n", _numSPL);
            }
            _segPosiList[_numSPL++] = po_r2;
        } else if (po_r2 < po_r1) {
            if (_numSPL >= _numCity) {
                printf("heap-buffer-overflow4: _numSPL = %d\n", _numSPL);
            }
            _segPosiList[_numSPL++] = po_r1;
        } else {
            assert(1 == 2);
        }

        _linkBPosi[po_r1][1] = _linkBPosi[po_r1][0];
        _linkBPosi[po_r2][1] = _linkBPosi[po_r2][0];
        _linkBPosi[po_r1][0] = po_b1;
        _linkBPosi[po_r2][0] = po_b2;
    }
}

void Cross::MakeCompleteSol (Indi& tKid)
{
    const int n = _numCity;
    int j,j1,j2;
    int st,pre,curr,next,a,b,c,d,aa,bb,a1,b1;
    int min_unit_city;
    int near_num;
    int center_un;
    int select_un;
    EvalType diff,max_diff;
    int nearMax;

    _gainModi = 0;

    while (_numUnit != 1) {
        min_unit_city = n + 12345;
        for (int u = 0; u < _numUnit; ++u) {
            if (_numElementInUnit[u] < min_unit_city) {
                center_un = u;
                min_unit_city = _numElementInUnit[u];
            }
        }

        st = -1;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == center_un) {
                int posi = _segment[s][0];
                st = _order[posi];
            }
        }
        assert(st != -1);

        curr = -1;
        next = st;
        _numElementInCU = 0;
        while (1) {
            pre = curr;
            curr = next;
            _centerUnit[curr] = 1;
            _listCenterUnit[_numElementInCU] = curr;
            ++_numElementInCU;

            if (tKid._link[curr][0] != pre) {
                next = tKid._link[curr][0];
            } else {
                next = tKid._link[curr][1];
            }

            if (next == st) { break; }
        }
        _listCenterUnit[_numElementInCU] = _listCenterUnit[0];
        _listCenterUnit[_numElementInCU+1] = _listCenterUnit[1];

        assert(_numElementInCU == _numElementInUnit[center_un]);

        max_diff = std::numeric_limits<EvalType>::min();
        a1 = -1; b1 = -1;
        /*
         *N_near (see Step 5.3 in Section 2.2 of the Online Supplement)
         * nearMax must be smaller than or equal to eva->_maxNumNear (kopt.cpp)
         */
        nearMax = 10;

    RESTART:;
        for (int s = 1; s <= _numElementInCU; ++s) {
            a = _listCenterUnit[s];
            for (near_num = 1; near_num <= nearMax; ++near_num) {
                c = _eval->_near[a][near_num];
                if (_centerUnit[c] == 0) {
                    for (j1 = 0; j1 < 2; ++j1) {
                        b = _listCenterUnit[s-1+2*j1];
                        for (j2 = 0; j2 < 2; ++j2) {
                            d = tKid._link[c][j2];
                            diff = _eval->_cost[a][b] + _eval->_cost[c][d] -
                                   _eval->_cost[a][c] - _eval->_cost[b][d];
                            if (diff > max_diff) {
                                aa = a; bb = b; a1 = c; b1 = d;
                                max_diff = diff;
                            }
                            diff = _eval->_cost[a][b] + _eval->_cost[d][c] -
                                   _eval->_cost[a][d] - _eval->_cost[b][c];
                            if (diff > max_diff) {
                                aa = a; bb = b; a1 = d; b1 = c;
                                max_diff = diff;
                            }
                        }
                    }
                }
            }
        }

        /* This value must also be changed if nearMax is chenged above */
        if (a1 == -1 && nearMax == 10) {
            nearMax = 50;
            goto RESTART;
        } else if (a1 == -1 && nearMax == 50) {
            int r = (*_eval->_rand)() % (_numElementInCU - 1);
            a = _listCenterUnit[r];
            b = _listCenterUnit[r+1];
            for (j = 0; j < n; ++j) {
                if (_centerUnit[j] == 0) {
                    aa = a; bb = b;
                    a1 = j;
                    b1 = tKid._link[j][0];
                    break;
                }
            }
            max_diff = _eval->_cost[aa][bb] + _eval->_cost[a1][b1] -
                       _eval->_cost[a][a1] - _eval->_cost[b][b1];
        }

        if (tKid._link[aa][0] == bb) {
            tKid._link[aa][0]=a1;
        } else {
            tKid._link[aa][1] = a1;
        }
        if (tKid._link[bb][0] == aa) {
            tKid._link[bb][0] = b1;
        } else {
            tKid._link[bb][1] = b1;
        }
        if (tKid._link[a1][0] == b1) {
            tKid._link[a1][0] = aa;
        } else {
            tKid._link[a1][1] = aa;
        }
        if (tKid._link[b1][0] == a1) {
            tKid._link[b1][0] = bb;
        } else {
            tKid._link[b1][1] = bb;
        }

        _modiEdge[_numModiEdge][0] = aa;
        _modiEdge[_numModiEdge][1] = bb;
        _modiEdge[_numModiEdge][2] = a1;
        _modiEdge[_numModiEdge][3] = b1;
        ++_numModiEdge;

        _gainModi += max_diff;

        int posi_a1 = _inv[a1];
        select_un = -1;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segment[s][0] <= posi_a1 && posi_a1 <=  _segment[s][1]) {
                select_un = _segUnit[s];
                break;
            }
        }
        assert(select_un != -1);

        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == select_un) {
                _segUnit[s] = center_un;
            }
        }
        _numElementInUnit[center_un] += _numElementInUnit[select_un];

        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == _numUnit - 1) {
                _segUnit[s] = select_un;
            }
        }
        _numElementInUnit[select_un] = _numElementInUnit[_numUnit - 1];
        --_numUnit;

        for (int s = 0; s < _numElementInCU; ++s) {
            c = _listCenterUnit[s];
            _centerUnit[c] = 0;
        }
    }
}

void Cross::MakeUnit ()
{
    const int n = _numCity;
    int flag = 1;
    for (int s = 0; s < _numSPL; ++s) {
        if (_segPosiList[s] == 0) {
            flag = 0;
            break;
        }
    }
    if (flag == 1) {
        _segPosiList[_numSPL++] = 0;

        _linkBPosi[n-1][1]  = _linkBPosi[n-1][0];
        _linkBPosi[0][1] = _linkBPosi[0][0];
        _linkBPosi[n-1][0] = 0;
        _linkBPosi[0][0] = n-1;
    }

    std::sort(_segPosiList, _segPosiList + _numSPL);

    _numSeg = _numSPL;
    for (int s = 0; s < _numSeg-1; ++s) {
        _segment[s][0] = _segPosiList[s];
        _segment[s][1] = _segPosiList[s+1]-1;
    }

    _segment[_numSeg-1][0] = _segPosiList[_numSeg-1];
    _segment[_numSeg-1][1] = n - 1;

    for (int s = 0; s < _numSeg; ++s) {
        _linkAPosi[_segment[s][0]] = _segment[s][1];
        _linkAPosi[_segment[s][1]] = _segment[s][0];
        _posiSeg[_segment[s][0]] = s;
        _posiSeg[_segment[s][1]] = s;
    }

    for (int s = 0; s < _numSeg; ++s) {
        _segUnit[s] = -1;
    }
    _numUnit = 0;

    int p_st, p1, p2, p_next, p_pre;
    int segNum;

    while (1) {
        flag = 0;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segUnit[s] == -1) {
                p_st = _segment[s][0];
                p_pre = -1;
                p1 = p_st;
                flag = 1;
                break;
            }
        }
        if (flag == 0) {
            break;
        }

        while (1) {
            segNum = _posiSeg[p1];
            _segUnit[segNum] = _numUnit;

            p2 = _linkAPosi[p1];
            p_next = _linkBPosi[p2][0];
            if (p1 == p2) {
                if (p_next == p_pre) {
                    p_next = _linkBPosi[p2][1];
                }
            }

            if (p_next == p_st) {
                ++_numUnit;
                break;
            }

            p_pre = p2;
            p1 = p_next;
        }
    }

    for (int s = 0; s < _numUnit; ++s) {
        _numElementInUnit[s] = 0;
    }

    int unitNum = -1;
    int tmpNumSeg = -1;
    for (int s = 0; s < _numSeg; ++s) {
        if (_segUnit[s] != unitNum) {
            ++tmpNumSeg;
            _segment[tmpNumSeg][0] = _segment[s][0];
            _segment[tmpNumSeg][1] = _segment[s][1];
            unitNum = _segUnit[s];
            _segUnit[tmpNumSeg] = unitNum;
            _numElementInUnit[unitNum] +=
            _segment[s][1] - _segment[s][0] + 1;
        } else {
            _segment[tmpNumSeg][1] = _segment[s][1];
            _numElementInUnit[unitNum] +=
            _segment[s][1] - _segment[s][0] + 1;
        }
    }
    _numSeg = tmpNumSeg + 1;
}

void Cross::BackToPa1 (Indi& tKid)
{
    int aa, bb, a1, b1;
    int jnum;

    for (int s = _numModiEdge -1; s >= 0; --s) {
        aa = _modiEdge[s][0];
        a1 = _modiEdge[s][1];   // $B$3$3$rJQ99$KCm0U(B
        bb = _modiEdge[s][2];   // $B$3$3$rJQ99$KCm0U(B
        b1 = _modiEdge[s][3];

        if (tKid._link[aa][0] == bb) {
            tKid._link[aa][0] = a1;
        } else {
            tKid._link[aa][1] = a1;
        }
        if (tKid._link[b1][0] == a1) {
            tKid._link[b1][0] = bb;
        } else {
            tKid._link[b1][1] = bb;
        }
        if (tKid._link[bb][0] == aa) {
            tKid._link[bb][0] = b1;
        } else {
            tKid._link[bb][1] = b1;
        }
        if (tKid._link[a1][0] == b1) {
            tKid._link[a1][0] = aa;
        } else {
            tKid._link[a1][1] = aa;
        }
    }

    for (int s = 0; s < _numAppliedCycle; ++s) {
        jnum = _appliedCylce[s];
        ChangeSol(tKid, jnum, 2);
    }
}

void Cross::GoToBest (Indi& tKid)
{
    int aa, bb, a1, b1;
    int jnum;

    for (int s = 0; s < _numBestAppliedCycle; ++s) {
        jnum = _bestAppliedCylce[s];
        ChangeSol(tKid, jnum, 1);
    }

    for (int s = 0; s < _numBestModiEdge; ++s) {
        aa = _bestModiEdge[s][0];
        bb = _bestModiEdge[s][1];
        a1 = _bestModiEdge[s][2];
        b1 = _bestModiEdge[s][3];

        if (tKid._link[aa][0] == bb) {
            tKid._link[aa][0]=a1;
        } else {
            tKid._link[aa][1] = a1;
        }
        if (tKid._link[bb][0] == aa) {
            tKid._link[bb][0] = b1;
        } else {
            tKid._link[bb][1] = b1;
        }
        if (tKid._link[a1][0] == b1) {
            tKid._link[a1][0] = aa;
        } else {
            tKid._link[a1][1] = aa;
        }
        if (tKid._link[b1][0] == a1) {
            tKid._link[b1][0] = bb;
        } else {
            tKid._link[b1][1] = bb;
        }
    }
}

void Cross::CheckValid (Indi& indi)
{
    const int n = _numCity;
    int curr, pre, next, st;
    int count;

    st = 0;
    curr = -1;
    next = st;

    count = 0;
    while (1) {
        pre = curr;
        curr = next;
        ++count;
        if (indi._link[curr][0] != pre) {
            next = indi._link[curr][0];
        } else {
            next = indi._link[curr][1];
        }

        if (next == st) { break; }

        if (count > n) {
            printf("Error: Invalid = %d\n", count);
            break;
        }
    }
    if (count != n) {
        printf("Error: Invalid = %d\n", count);
    }
}
