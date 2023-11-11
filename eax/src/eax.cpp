/*
License: see tsp.h
 */

#include "eax.h"

EAXGA::EAXGA () :
    _eval(new Evaluator()),
    _kopt(nullptr),
    _cross(nullptr),
    _pop(nullptr),
    _matingSeq(nullptr)
{
    _numPop = 100;
    _numKids = 30;
    _silent = false;
}

EAXGA::~EAXGA ()
{
    delete _kopt;
    delete _cross;
    delete [] _pop;
    delete [] _matingSeq;
    delete _eval;
}

void EAXGA::Define (const char* tspFileName)
{
    _eval->SetInstance(tspFileName);

    const int n = _eval->_numCity;
    _best.Define(n);

    _kopt = new KOpt(_eval);
    _cross = new Cross(_eval, _numPop);
    _pop = new Indi [_numPop];
    for (int i = 0; i < _numPop; ++i) {
        _pop[i].Define(n);
    }
    _matingSeq = new int [_numPop + 1];

}

void EAXGA::DoIt ()
{
    _numGen = 0;
    _stagnGen = 0;
    _best.MakeRand(_eval);

    for (int i = 0; i < _numPop; ++i) {
        _kopt->DoIt(_pop[i]);
    }

    while (1) {
        SelectBest();
        if (!_silent) {
            printf("%d: %lld %lf\n", _numGen, (long long)_best._cost, _avgCost);
        }

        if (ShouldTerminate()) { break; }

        SelectForMating();
        for (int s = 0; s < _numPop; ++s) {
            _cross->DoIt(_pop[_matingSeq[s]], _pop[_matingSeq[s+1]], _numKids);
        }

        ++_numGen;
    }
}

void EAXGA::SelectBest ()
{
    EvalType stockBest = _best._cost;

    _avgCost = 0.0;
    int bestIndex = 0;
    EvalType bestCost = _pop[0]._cost;

    for (int i = 0; i < _numPop; ++i) {
        _avgCost += _pop[i]._cost;
        if (_pop[i]._cost < bestCost) {
            bestIndex = i;
            bestCost = _pop[i]._cost;
        }
    }

    _best = _pop[bestIndex];
    _avgCost /= _numPop;

    if (_best._cost < stockBest) {
        _stagnGen = 0;
    } else {
        _stagnGen++;
    }
}

bool EAXGA::ShouldTerminate ()
{
    if (_avgCost - _best._cost < 0.001) {
        return true;
    }

    if (_stagnGen > (1500 / _numKids)) {
        return true;
    }

    return false;
}

void EAXGA::SelectForMating ()
{
    for (int i = 0; i < _numPop; i++) {
        _matingSeq[i] = i;
    }
    std::shuffle(_matingSeq, _matingSeq+_numPop, *_eval->_rand);
    _matingSeq[_numPop] = _matingSeq[0];
}

Cross::Cross (const Evaluator* e, int nPop) :
    _eval(e),
    _numCity(e->_numCity),
    _numPop(nPop),
    _maxNumABcycle(2000)
{
    const int n = _numCity;
    _nearData = new int* [n];
    for (int j = 0; j < n; ++j) {
        _nearData[j] = new int [5];
    }

    _ABcycle = new int* [_maxNumABcycle];
    for (int j = 0; j < _maxNumABcycle; ++j) {
        _ABcycle[j] = new int [2*n + 4];
    }

    _koritsu = new int [n];
    _bunki = new int [n];
    _koriInv = new int [n];
    _bunInv = new int [n];
    _route = new int [2*n + 1];
    _permuABCycle = new int [_maxNumABcycle];

    _cycle = new int [2*n+4];

    _path = new int [n];
    _posi = new int [n];
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
    _gainABcycle = new EvalType [_maxNumABcycle];
    _modiEdge = new int* [n];
    for (int j = 0; j < n; ++j) {
        _modiEdge[j] = new int [4];
    }
    _bestModiEdge = new int* [n];
    for (int j = 0; j < n; ++j) {
        _bestModiEdge[j] = new int [4];
    }
    _appliedCycle = new int [n];
    _bestAppliedCycle = new int [n];
    _ABCycleInEset = new int [n];
}

Cross::~Cross ()
{
    const int n = _numCity;
    delete [] _koritsu;
    delete [] _bunki;
    delete [] _koriInv;
    delete [] _bunInv;
    delete [] _route;
    delete [] _permuABCycle;

    for (int j = 0; j < n; ++j) {
        delete[] _nearData[j];
    }
    delete[] _nearData;

    for (int j = 0; j < _maxNumABcycle; ++j) {
        delete[] _ABcycle[j];
    }
    delete[] _ABcycle;

    delete [] _cycle;

    delete [] _path;
    delete [] _posi;

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
    delete [] _gainABcycle;

    for (int j = 0; j < n; ++j) {
        delete[] _modiEdge[j];
    }
    delete [] _modiEdge;
    for (int j = 0; j < n; ++j) {
        delete[] _bestModiEdge[j];
    }
    delete [] _bestModiEdge;
    delete [] _appliedCycle;
    delete [] _bestAppliedCycle;

    delete [] _ABCycleInEset;
}

void Cross::DoIt (Indi& kid, Indi& pa2, int nKids)
{
    EvalType bestGain = 0, gain;
    int idx;

    SetABcycle(kid, pa2, nKids);

    if (nKids > _numABcycle) {
        nKids = _numABcycle;
    }

    for (int i = 0; i < _numABcycle; i++) {
        _permuABCycle[i] = i;
    }
    std::shuffle(_permuABCycle, _permuABCycle+_numABcycle, *_eval->_rand);

    for (int j = 0; j < nKids; ++j) {
        _numABcycleInEset = 0;
        idx = _permuABCycle[j];
        _ABCycleInEset[_numABcycleInEset++] = idx;

        gain = 0;
        _numModiEdge = 0;
        _numSPL = 0;

        _numAppliedCycle = _numABcycleInEset;
        for (int k = 0; k < _numAppliedCycle; ++k) {
            _appliedCycle[k] = _ABCycleInEset[k];
            idx = _appliedCycle[k];
            ChangeSol(kid, idx, 1);
            gain += _gainABcycle[idx];
        }

        MakeUnit();
        gain += MakeCompleteSol(kid);
        kid._cost = kid._cost - gain;

        if (bestGain < gain) {
            bestGain = gain;

            _numBestAppliedCycle = _numAppliedCycle;
            memcpy(_bestAppliedCycle, _appliedCycle, sizeof(int)*_numAppliedCycle);

            _numBestModiEdge = _numModiEdge;
            for (int s = 0; s < _numBestModiEdge; ++s) {
                memcpy(_bestModiEdge[s], _modiEdge[s], sizeof(int)*4);
            }
        }

        BackToPa1(kid);
        kid._cost = kid._cost + gain;
    }

    if (bestGain != 0) {
        GoToBest(kid);
        kid._cost = kid._cost - bestGain;
    }
}

void Cross::SetABcycle (const Indi& pa1, const Indi& pa2, int nKids)
{
    const int n = _numCity;
    int* checkKoritsu = _eval->_buf;
    int curr = -1, next = 0, prev = -1;
    _bunkiMany=0; _koritsuMany=0;

    for (int j = 0; j < n ; ++j) {
        _nearData[j][1]=pa1._link[j][0];
        _nearData[j][3]=pa1._link[j][1];
        _nearData[j][0] = 2;
        _nearData[j][2]=pa2._link[j][0];
        _nearData[j][4]=pa2._link[j][1];

        _koritsu[_koritsuMany++]=j;
        _koriInv[_koritsu[j]]=j;
        checkKoritsu[j]=-1;

        // init _path and _posi for pa1
        prev = curr;
        curr = next;
        if (pa1._link[curr][0] != prev) {
            next = pa1._link[curr][0];
        } else {
            next = pa1._link[curr][1];
        }
        _path[j] = curr;
        _posi[curr] = j;
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
            checkKoritsu[st]=_posiCurr;
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
                    if (checkKoritsu[st]==0) {
                        if ((_posiCurr-checkKoritsu[st])%2==0) {
                            if (_nearData[st][_posiCurr%2+1]==pr) {
                                std::swap(_nearData[ci][_posiCurr%2+1],
                                          _nearData[ci][_posiCurr%2+3]);
                            }
                            _stAppear = 1;
                            FormABcycle();
                            if (_numABcycle == nKids) { goto LLL; }
                            if (_numABcycle == _maxNumABcycle) { goto LLL; }

                            _flagSt=0;
                            _flagCircle=1;
                            _prType=1;
                        } else {
                            std::swap(_nearData[ci][_posiCurr%2+1],
                                      _nearData[ci][_posiCurr%2+3]);
                            _prType=2;
                        }
                        checkKoritsu[st]=_posiCurr;
                    } else {
                        _stAppear = 2;
                        FormABcycle();
                        if (_numABcycle == nKids) { goto LLL; }
                        if (_numABcycle == _maxNumABcycle) { goto LLL; }

                        _flagSt=1;
                        _flagCircle=1;
                    }
                } else if (checkKoritsu[ci]==-1) {
                    checkKoritsu[ci]=_posiCurr;
                    if (_nearData[ci][_posiCurr%2+1]==pr) {
                        std::swap(_nearData[ci][_posiCurr%2+1],
                                  _nearData[ci][_posiCurr%2+3]);
                    }
                    _prType=2;
                } else if (checkKoritsu[ci]>0) {
                    std::swap(_nearData[ci][_posiCurr%2+1],
                            _nearData[ci][_posiCurr%2+3]);
                    if ((_posiCurr-checkKoritsu[ci])%2==0) {
                        _stAppear = 1;
                        FormABcycle();
                        if (_numABcycle == nKids) { goto LLL; }
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
                    if (_numABcycle == nKids) { goto LLL; }
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
                if (_numABcycle == nKids) { goto LLL; }
                if (_numABcycle == _maxNumABcycle) { goto LLL; }
                _flagCircle=1;
            }
        }
    }

    LLL: ;

    if (_numABcycle == _maxNumABcycle) {
        printf("Warning: _maxNumABcycle(%d) must be increased\n", _maxNumABcycle);
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
    _cycle[cem]=st;

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
        _cycle[cem]=ci;
    }

    if (cem==2) {
        return;
    }

    _ABcycle[_numABcycle][0]=cem;

    if (edge_type==2) {
        stock=_cycle[0];
        for (int j=0;j<cem-1;j++) { _cycle[j]=_cycle[j+1]; }
        _cycle[cem-1]=stock;
    }

    for (int j=0;j<cem;j++) {
        _ABcycle[_numABcycle][j+2]=_cycle[j];
    }
    _ABcycle[_numABcycle][1]=_cycle[cem-1];
    _ABcycle[_numABcycle][cem+2]=_cycle[0];
    _ABcycle[_numABcycle][cem+3]=_cycle[1];

    _cycle[cem] = _cycle[0];
    _cycle[cem+1] = _cycle[1];
    diff = 0;
    for (j = 0; j < cem/2; ++j) {
        diff += _eval->_cost[_cycle[2*j]][_cycle[1+2*j]] -
                _eval->_cost[_cycle[1+2*j]][_cycle[2+2*j]];
    }
    _gainABcycle[_numABcycle] = diff;
    ++_numABcycle;
}

void Cross::ChangeSol (Indi& kid, int idx, int type)
{
    const int n = _numCity;
    int j;
    int cem,r1,r2,b1,b2;
    int po_r1, po_r2, po_b1, po_b2;

    cem=_ABcycle[idx][0];
    _cycle[0]=_ABcycle[idx][0];

    if (type==2) {
        for (j=0;j<cem+3;j++) {
            _cycle[cem+3-j]=_ABcycle[idx][j+1];
        }
    } else {
        for (j=1;j<=cem+3;j++) {
            _cycle[j]=_ABcycle[idx][j];
        }
    }

    for (j=0;j<cem/2;j++) {
        r1=_cycle[2+2*j];r2=_cycle[3+2*j];
        b1=_cycle[1+2*j];b2=_cycle[4+2*j];

        if (kid._link[r1][0]==r2) {
            kid._link[r1][0]=b1;
        } else {
            kid._link[r1][1]=b1;
        }
        if (kid._link[r2][0]==r1) {
            kid._link[r2][0]=b2;
        } else {
            kid._link[r2][1]=b2;
        }

        po_r1 = _posi[r1];
        po_r2 = _posi[r2];
        po_b1 = _posi[b1];
        po_b2 = _posi[b2];

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

EvalType Cross::MakeCompleteSol (Indi& kid)
{
    EvalType gainModi = 0;

    const int n = _numCity;
    int j,j1,j2;
    int st,pre,curr,next,a,b,c,d,aa,bb,a1,b1;
    int min_unit_city;
    int near_num;
    int center_un;
    int select_un;
    EvalType diff,max_diff;
    int nearMax;

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
                st = _path[posi];
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

            if (kid._link[curr][0] != pre) {
                next = kid._link[curr][0];
            } else {
                next = kid._link[curr][1];
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
                            d = kid._link[c][j2];
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
                    b1 = kid._link[j][0];
                    break;
                }
            }
            max_diff = _eval->_cost[aa][bb] + _eval->_cost[a1][b1] -
                       _eval->_cost[a][a1] - _eval->_cost[b][b1];
        }

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0]=a1;
        } else {
            kid._link[aa][1] = a1;
        }
        if (kid._link[bb][0] == aa) {
            kid._link[bb][0] = b1;
        } else {
            kid._link[bb][1] = b1;
        }
        if (kid._link[a1][0] == b1) {
            kid._link[a1][0] = aa;
        } else {
            kid._link[a1][1] = aa;
        }
        if (kid._link[b1][0] == a1) {
            kid._link[b1][0] = bb;
        } else {
            kid._link[b1][1] = bb;
        }

        _modiEdge[_numModiEdge][0] = aa;
        _modiEdge[_numModiEdge][1] = bb;
        _modiEdge[_numModiEdge][2] = a1;
        _modiEdge[_numModiEdge][3] = b1;
        ++_numModiEdge;

        gainModi += max_diff;

        int posi_a1 = _posi[a1];
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

    return gainModi;
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

void Cross::BackToPa1 (Indi& kid)
{
    int aa, bb, a1, b1;
    int idx;

    for (int s = _numModiEdge - 1; s >= 0; --s) {
        aa = _modiEdge[s][0];
        a1 = _modiEdge[s][1];
        bb = _modiEdge[s][2];
        b1 = _modiEdge[s][3];

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0] = a1;
        } else {
            kid._link[aa][1] = a1;
        }
        if (kid._link[b1][0] == a1) {
            kid._link[b1][0] = bb;
        } else {
            kid._link[b1][1] = bb;
        }
        if (kid._link[bb][0] == aa) {
            kid._link[bb][0] = b1;
        } else {
            kid._link[bb][1] = b1;
        }
        if (kid._link[a1][0] == b1) {
            kid._link[a1][0] = aa;
        } else {
            kid._link[a1][1] = aa;
        }
    }

    for (int s = 0; s < _numAppliedCycle; ++s) {
        idx = _appliedCycle[s];
        ChangeSol(kid, idx, 2);
    }
}

void Cross::GoToBest (Indi& kid)
{
    int aa, bb, a1, b1;
    int idx;

    for (int s = 0; s < _numBestAppliedCycle; ++s) {
        idx = _bestAppliedCycle[s];
        ChangeSol(kid, idx, 1);
    }

    for (int s = 0; s < _numBestModiEdge; ++s) {
        aa = _bestModiEdge[s][0];
        bb = _bestModiEdge[s][1];
        a1 = _bestModiEdge[s][2];
        b1 = _bestModiEdge[s][3];

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0]=a1;
        } else {
            kid._link[aa][1] = a1;
        }
        if (kid._link[bb][0] == aa) {
            kid._link[bb][0] = b1;
        } else {
            kid._link[bb][1] = b1;
        }
        if (kid._link[a1][0] == b1) {
            kid._link[a1][0] = aa;
        } else {
            kid._link[a1][1] = aa;
        }
        if (kid._link[b1][0] == a1) {
            kid._link[b1][0] = bb;
        } else {
            kid._link[b1][1] = bb;
        }
    }
}
