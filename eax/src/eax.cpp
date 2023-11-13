/*
License: see tsp.h
 */

#include "eax.h"

EAXGA::EAXGA()
    : _eval(new Evaluator()), _kopt(nullptr), _cross(nullptr), _pop(nullptr),
      _matingSeq(nullptr)
{
    _numPop = 100;
    _numKids = 30;
    _silent = false;
}

EAXGA::~EAXGA()
{
    delete _kopt;
    delete _cross;
    delete[] _pop;
    delete[] _matingSeq;
    delete _eval;
}

void EAXGA::Define(const char* tspFileName)
{
    _eval->SetInstance(tspFileName);

    const int n = _eval->_numCity;
    _best.Define(n);

    _kopt = new KOpt(_eval);
    _cross = new Cross(_eval, _numPop);
    _pop = new Indi[_numPop];
    for (int i = 0; i < _numPop; ++i) {
        _pop[i].Define(n);
    }
    _matingSeq = new int[_numPop + 1];
}

void EAXGA::DoIt()
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

        if (ShouldTerminate()) {
            break;
        }

        SelectForMating();
        for (int s = 0; s < _numPop; ++s) {
            _cross->DoIt(_pop[_matingSeq[s]], _pop[_matingSeq[s + 1]],
                         _numKids);
        }

        ++_numGen;
    }
}

void EAXGA::SelectBest()
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

bool EAXGA::ShouldTerminate()
{
    if (_avgCost - _best._cost < 0.001) {
        return true;
    }

    if (_stagnGen > (1500 / _numKids)) {
        return true;
    }

    return false;
}

void EAXGA::SelectForMating()
{
    for (int i = 0; i < _numPop; i++) {
        _matingSeq[i] = i;
    }
    std::shuffle(_matingSeq, _matingSeq + _numPop, *_eval->_rand);
    _matingSeq[_numPop] = _matingSeq[0];
}

Cross::Cross(const Evaluator* e, int nPop)
    : _eval(e), _numCity(e->_numCity), _numPop(nPop), _maxNumABcycle(2000)
{
    const int n = _numCity;
    _nearData = new int*[n];
    for (int j = 0; j < n; ++j) {
        _nearData[j] = new int[5];
    }

    _ABcycle = new int*[_maxNumABcycle];
    for (int j = 0; j < _maxNumABcycle; ++j) {
        _ABcycle[j] = new int[2 * n + 4];
    }

    _koritsu = new int[n];
    _bunki = new int[n];
    _koriInv = new int[n];
    _bunInv = new int[n];
    _route = new int[2 * n + 1];
    _permuABCycle = new int[_maxNumABcycle];

    _cycle = new int[2 * n + 4];

    _path = new int[n];
    _posi = new int[n];
    _segment = new int*[n];
    for (int j = 0; j < n; ++j) {
        _segment[j] = new int[2];
    }
    _segUnit = new int[n];
    _segPosiList = new int[n];
    _linkAPosi = new int[n];
    _linkBPosi = new int*[n];
    for (int j = 0; j < n; ++j) {
        _linkBPosi[j] = new int[2];
    }
    _posiSeg = new int[n];
    _numElementInUnit = new int[n];
    _centerUnit = new int[n];
    for (int j = 0; j < n; ++j) {
        _centerUnit[j] = 0;
    }
    _listCenterUnit = new int[n + 2];
    _gainABcycle = new EvalType[_maxNumABcycle];
    _modiEdge = new int*[n];
    for (int j = 0; j < n; ++j) {
        _modiEdge[j] = new int[4];
    }
    _bestModiEdge = new int*[n];
    for (int j = 0; j < n; ++j) {
        _bestModiEdge[j] = new int[4];
    }
}

Cross::~Cross()
{
    const int n = _numCity;
    delete[] _koritsu;
    delete[] _bunki;
    delete[] _koriInv;
    delete[] _bunInv;
    delete[] _route;
    delete[] _permuABCycle;

    for (int j = 0; j < n; ++j) {
        delete[] _nearData[j];
    }
    delete[] _nearData;

    for (int j = 0; j < _maxNumABcycle; ++j) {
        delete[] _ABcycle[j];
    }
    delete[] _ABcycle;

    delete[] _cycle;

    delete[] _path;
    delete[] _posi;

    for (int j = 0; j < n; ++j) {
        delete[] _segment[j];
    }
    delete[] _segment;
    delete[] _segUnit;
    delete[] _segPosiList;
    delete[] _linkAPosi;
    for (int j = 0; j < n; ++j) {
        delete[] _linkBPosi[j];
    }
    delete[] _linkBPosi;
    delete[] _posiSeg;
    delete[] _numElementInUnit;
    delete[] _centerUnit;
    delete[] _listCenterUnit;
    delete[] _gainABcycle;

    for (int j = 0; j < n; ++j) {
        delete[] _modiEdge[j];
    }
    delete[] _modiEdge;
    for (int j = 0; j < n; ++j) {
        delete[] _bestModiEdge[j];
    }
    delete[] _bestModiEdge;
}

void Cross::DoIt(Indi& kid, Indi& pa2, int nKids)
{
    EvalType bestGain = 0, gain;
    int bestAppliedCycle, appliedCycle;

    SetABcycle(kid, pa2, nKids);

    if (nKids > _numABcycle) {
        nKids = _numABcycle;
    }

    for (int i = 0; i < _numABcycle; i++) {
        _permuABCycle[i] = i;
    }
    std::shuffle(_permuABCycle, _permuABCycle + _numABcycle, *_eval->_rand);

    for (int j = 0; j < nKids; ++j) {
        gain = 0;
        _numModiEdge = 0;
        _numSPL = 0;

        appliedCycle = _permuABCycle[j];
        ChangeSol(kid, appliedCycle, false /*reverse*/);
        gain += _gainABcycle[appliedCycle];

        MakeUnit();
        gain += MakeCompleteSol(kid);
        kid._cost = kid._cost - gain;

        if (bestGain < gain) {
            bestGain = gain;
            bestAppliedCycle = appliedCycle;
            _numBestModiEdge = _numModiEdge;
            for (int s = 0; s < _numBestModiEdge; ++s) {
                memcpy(_bestModiEdge[s], _modiEdge[s], sizeof(int) * 4);
            }
        }

        BackToPa1(kid, appliedCycle);
        kid._cost = kid._cost + gain;
    }

    if (bestGain != 0) {
        GoToBest(kid, bestAppliedCycle);
        kid._cost = kid._cost - bestGain;
    }
}

void Cross::SetABcycle(const Indi& pa1, const Indi& pa2, int nKids)
{
    const int n = _numCity;
    int* checkKoritsu = _eval->_buf;
    int curr = -1, next = 0, prev = -1;
    _bunkiMany = 0;
    _koritsuMany = 0;

    for (int j = 0; j < n; ++j) {
        _nearData[j][1] = pa1._link[j][0];
        _nearData[j][3] = pa1._link[j][1];
        _nearData[j][0] = 2;
        _nearData[j][2] = pa2._link[j][0];
        _nearData[j][4] = pa2._link[j][1];

        _koritsu[_koritsuMany++] = j;
        _koriInv[_koritsu[j]] = j;
        checkKoritsu[j] = -1;

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

    _numABcycle = 0;
    int flagSt = 1;
    int prType = 2;
    int flagCircle = 0;
    int posiCurr = 0;
    int r = 0, pr = 0, st = 0, ci = 0;
    while (_koritsuMany != 0) {
        if (flagSt == 1) {
            posiCurr = 0;
            r = (*_eval->_rand)() % _koritsuMany;
            st = _koritsu[r];
            checkKoritsu[st] = posiCurr;
            _route[posiCurr] = st;
            ci = st;
            prType = 2;
        } else if (flagSt == 0) {
            ci = _route[posiCurr];
        }

        flagCircle = 0;
        while (flagCircle == 0) {
            posiCurr++;
            pr = ci;

            switch (prType) {
                case 1:
                    ci = _nearData[pr][posiCurr % 2 + 1];
                    break;
                case 2:
                    r = (*_eval->_rand)() % 2;
                    ci = _nearData[pr][posiCurr % 2 + 1 + 2 * r];
                    if (r == 0) {
                        std::swap(_nearData[pr][posiCurr % 2 + 1],
                                  _nearData[pr][posiCurr % 2 + 3]);
                    }
                    break;
                case 3:
                    ci = _nearData[pr][posiCurr % 2 + 3];
            }

            _route[posiCurr] = ci;

            if (_nearData[ci][0] == 2) {
                if (ci == st) {
                    if (checkKoritsu[st] == 0) {
                        if ((posiCurr - checkKoritsu[st]) % 2 == 0) {
                            if (_nearData[st][posiCurr % 2 + 1] == pr) {
                                std::swap(_nearData[ci][posiCurr % 2 + 1],
                                          _nearData[ci][posiCurr % 2 + 3]);
                            }
                            FormABcycle(1, posiCurr);
                            if (_numABcycle == nKids) {
                                goto LLL;
                            }
                            if (_numABcycle == _maxNumABcycle) {
                                goto LLL;
                            }

                            flagSt = 0;
                            flagCircle = 1;
                            prType = 1;
                        } else {
                            std::swap(_nearData[ci][posiCurr % 2 + 1],
                                      _nearData[ci][posiCurr % 2 + 3]);
                            prType = 2;
                        }
                        checkKoritsu[st] = posiCurr;
                    } else {
                        FormABcycle(2, posiCurr);
                        if (_numABcycle == nKids) {
                            goto LLL;
                        }
                        if (_numABcycle == _maxNumABcycle) {
                            goto LLL;
                        }

                        flagSt = 1;
                        flagCircle = 1;
                    }
                } else if (checkKoritsu[ci] == -1) {
                    checkKoritsu[ci] = posiCurr;
                    if (_nearData[ci][posiCurr % 2 + 1] == pr) {
                        std::swap(_nearData[ci][posiCurr % 2 + 1],
                                  _nearData[ci][posiCurr % 2 + 3]);
                    }
                    prType = 2;
                } else if (checkKoritsu[ci] > 0) {
                    std::swap(_nearData[ci][posiCurr % 2 + 1],
                              _nearData[ci][posiCurr % 2 + 3]);
                    if ((posiCurr - checkKoritsu[ci]) % 2 == 0) {
                        FormABcycle(1, posiCurr);
                        if (_numABcycle == nKids) {
                            goto LLL;
                        }
                        if (_numABcycle == _maxNumABcycle) {
                            goto LLL;
                        }

                        flagSt = 0;
                        flagCircle = 1;
                        prType = 1;
                    } else {
                        std::swap(_nearData[ci][(posiCurr + 1) % 2 + 1],
                                  _nearData[ci][(posiCurr + 1) % 2 + 3]);
                        prType = 3;
                    }
                }
            } else if (_nearData[ci][0] == 1) {
                if (ci == st) {
                    FormABcycle(1, posiCurr);
                    if (_numABcycle == nKids) {
                        goto LLL;
                    }
                    if (_numABcycle == _maxNumABcycle) {
                        goto LLL;
                    }

                    flagSt = 1;
                    flagCircle = 1;
                } else {
                    prType = 1;
                }
            }
        }
    }

    while (_bunkiMany != 0) {
        posiCurr = 0;
        r = (*_eval->_rand)() % _bunkiMany;
        st = _bunki[r];
        _route[posiCurr] = st;
        ci = st;

        flagCircle = 0;
        while (flagCircle == 0) {
            pr = ci;
            posiCurr++;
            ci = _nearData[pr][posiCurr % 2 + 1];
            _route[posiCurr] = ci;
            if (ci == st) {
                FormABcycle(1, posiCurr);
                if (_numABcycle == nKids) {
                    goto LLL;
                }
                if (_numABcycle == _maxNumABcycle) {
                    goto LLL;
                }
                flagCircle = 1;
            }
        }
    }

LLL:;

    if (_numABcycle == _maxNumABcycle) {
        printf("Warning: _maxNumABcycle(%d) must be increased\n",
               _maxNumABcycle);
    }
}

void Cross::FormABcycle(int stAppear, int& posiCurr)
{
    const int st = _route[posiCurr];
    int st_count = 0;
    int cem = 0;
    _cycle[cem] = st;

    while (1) {
        cem++;
        posiCurr--;
        int ci = _route[posiCurr];
        if (_nearData[ci][0] == 2) {
            _koritsu[_koriInv[ci]] = _koritsu[_koritsuMany - 1];
            _koriInv[_koritsu[_koritsuMany - 1]] = _koriInv[ci];
            _koritsuMany--;
            _bunki[_bunkiMany] = ci;
            _bunInv[ci] = _bunkiMany;
            _bunkiMany++;
        } else if (_nearData[ci][0] == 1) {
            _bunki[_bunInv[ci]] = _bunki[_bunkiMany - 1];
            _bunInv[_bunki[_bunkiMany - 1]] = _bunInv[ci];
            _bunkiMany--;
        }

        _nearData[ci][0]--;
        if (ci == st) {
            st_count++;
        }
        if (st_count == stAppear) {
            break;
        }
        _cycle[cem] = ci;
    }

    if (cem == 2) {
        return;
    }

    _ABcycle[_numABcycle][0] = cem;

    if (posiCurr % 2 != 0) {
        int stock = _cycle[0];
        for (int j = 0; j < cem - 1; j++) {
            _cycle[j] = _cycle[j + 1];
        }
        _cycle[cem - 1] = stock;
    }

    for (int j = 0; j < cem; j++) {
        _ABcycle[_numABcycle][j + 2] = _cycle[j];
    }
    _ABcycle[_numABcycle][1] = _cycle[cem - 1];
    _ABcycle[_numABcycle][cem + 2] = _cycle[0];
    _ABcycle[_numABcycle][cem + 3] = _cycle[1];

    _cycle[cem] = _cycle[0];
    _cycle[cem + 1] = _cycle[1];
    EvalType diff = 0;
    for (int j = 0; j < cem / 2; ++j) {
        diff += _eval->_cost[_cycle[2 * j]][_cycle[1 + 2 * j]]
                - _eval->_cost[_cycle[1 + 2 * j]][_cycle[2 + 2 * j]];
    }
    _gainABcycle[_numABcycle] = diff;
    ++_numABcycle;
}

void Cross::ChangeSol(Indi& kid, int idx, bool reverse, bool updateSeg)
{
    const int n = _numCity;
    int cem, r1, r2, b1, b2;

    cem = _ABcycle[idx][0];
    _cycle[0] = _ABcycle[idx][0];

    if (reverse) {
        for (int j = 0; j < cem + 3; j++) {
            _cycle[cem + 3 - j] = _ABcycle[idx][j + 1];
        }
    } else {
        for (int j = 1; j <= cem + 3; j++) {
            _cycle[j] = _ABcycle[idx][j];
        }
    }

    for (int j = 0; j < cem / 2; j++) {
        r1 = _cycle[2 + 2 * j];
        r2 = _cycle[3 + 2 * j];
        b1 = _cycle[1 + 2 * j];
        b2 = _cycle[4 + 2 * j];

        if (kid._link[r1][0] == r2) {
            kid._link[r1][0] = b1;
        } else {
            kid._link[r1][1] = b1;
        }
        if (kid._link[r2][0] == r1) {
            kid._link[r2][0] = b2;
        } else {
            kid._link[r2][1] = b2;
        }

        if (updateSeg) {
            if (_numSPL >= n) {
                printf("Error: numSPL reach max(%d) in ChangeSol!", n);
                exit(1);
            }
            if (_posi[r1] == 0 && _posi[r2] == n - 1) {
                _segPosiList[_numSPL++] = _posi[r1];
            } else if (_posi[r1] == n - 1 && _posi[r2] == 0) {
                _segPosiList[_numSPL++] = _posi[r2];
            } else if (_posi[r1] < _posi[r2]) {
                _segPosiList[_numSPL++] = _posi[r2];
            } else if (_posi[r2] < _posi[r1]) {
                _segPosiList[_numSPL++] = _posi[r1];
            } else {
                assert(0);
            }

            _linkBPosi[_posi[r1]][1] = _linkBPosi[_posi[r1]][0];
            _linkBPosi[_posi[r2]][1] = _linkBPosi[_posi[r2]][0];
            _linkBPosi[_posi[r1]][0] = _posi[b1];
            _linkBPosi[_posi[r2]][0] = _posi[b2];
        }
    }
}

void Cross::MakeUnit()
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
        if (_numSPL >= n) {
            printf("Error: numSPL reach max(%d) in MakeUnit!", n);
            exit(1);
        }
        _segPosiList[_numSPL++] = 0;

        _linkBPosi[n - 1][1] = _linkBPosi[n - 1][0];
        _linkBPosi[0][1] = _linkBPosi[0][0];
        _linkBPosi[n - 1][0] = 0;
        _linkBPosi[0][0] = n - 1;
    }

    std::sort(_segPosiList, _segPosiList + _numSPL);

    _numSeg = _numSPL;
    for (int s = 0; s < _numSeg - 1; ++s) {
        _segment[s][0] = _segPosiList[s];
        _segment[s][1] = _segPosiList[s + 1] - 1;
    }

    _segment[_numSeg - 1][0] = _segPosiList[_numSeg - 1];
    _segment[_numSeg - 1][1] = n - 1;

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
            _numElementInUnit[unitNum] += _segment[s][1] - _segment[s][0] + 1;
        } else {
            _segment[tmpNumSeg][1] = _segment[s][1];
            _numElementInUnit[unitNum] += _segment[s][1] - _segment[s][0] + 1;
        }
    }
    _numSeg = tmpNumSeg + 1;
}

EvalType Cross::MakeCompleteSol(Indi& kid)
{
    EvalType gainModi = 0;
    constexpr int NearMaxDef = 10;
    int center_un = 0;

    int st, prev, curr, next, a, b, c, d, aa = 0, bb = 0, a1 = 0, b1 = 0;

    while (_numUnit != 1) {
        int min_unit_city = _numCity + 12345;
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
            prev = curr;
            curr = next;
            _centerUnit[curr] = 1;
            _listCenterUnit[_numElementInCU] = curr;
            ++_numElementInCU;

            if (kid._link[curr][0] != prev) {
                next = kid._link[curr][0];
            } else {
                next = kid._link[curr][1];
            }

            if (next == st) {
                break;
            }
        }
        _listCenterUnit[_numElementInCU] = _listCenterUnit[0];
        _listCenterUnit[_numElementInCU + 1] = _listCenterUnit[1];

        assert(_numElementInCU == _numElementInUnit[center_un]);

        EvalType max_diff = std::numeric_limits<EvalType>::min();
        EvalType diff = 0;
        a1 = -1;
        b1 = -1;
        /*
         *N_near (see Step 5.3 in Section 2.2 of the Online Supplement)
         * nearMax must be smaller than or equal to eva->_maxNumNear (kopt.cpp)
         */
        assert(NearMaxDef <= _eval->_maxNumNear);
        int nearMax = NearMaxDef;

    RESTART:;
        for (int s = 1; s <= _numElementInCU; ++s) {
            a = _listCenterUnit[s];
            for (int near_num = 1; near_num <= nearMax; ++near_num) {
                c = _eval->_near[a][near_num];
                if (_centerUnit[c] == 0) {
                    for (int j1 = 0; j1 < 2; ++j1) {
                        b = _listCenterUnit[s - 1 + 2 * j1];
                        for (int j2 = 0; j2 < 2; ++j2) {
                            d = kid._link[c][j2];
                            diff = _eval->_cost[a][b] + _eval->_cost[c][d]
                                   - _eval->_cost[a][c] - _eval->_cost[b][d];
                            if (diff > max_diff) {
                                aa = a;
                                bb = b;
                                a1 = c;
                                b1 = d;
                                max_diff = diff;
                            }
                            diff = _eval->_cost[a][b] + _eval->_cost[d][c]
                                   - _eval->_cost[a][d] - _eval->_cost[b][c];
                            if (diff > max_diff) {
                                aa = a;
                                bb = b;
                                a1 = d;
                                b1 = c;
                                max_diff = diff;
                            }
                        }
                    }
                }
            }
        }

        /* This value must also be changed if nearMax is chenged above */
        if (a1 == -1 && nearMax == NearMaxDef) {
            nearMax = _eval->_maxNumNear;
            goto RESTART;
        } else if (a1 == -1 && nearMax == _eval->_maxNumNear) {
            int r = (*_eval->_rand)() % (_numElementInCU - 1);
            a = _listCenterUnit[r];
            b = _listCenterUnit[r + 1];
            for (int j = 0; j < _numCity; ++j) {
                if (_centerUnit[j] == 0) {
                    aa = a;
                    bb = b;
                    a1 = j;
                    b1 = kid._link[j][0];
                    break;
                }
            }
            max_diff = _eval->_cost[aa][bb] + _eval->_cost[a1][b1]
                       - _eval->_cost[a][a1] - _eval->_cost[b][b1];
        }

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0] = a1;
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
        int select_un = -1;
        for (int s = 0; s < _numSeg; ++s) {
            if (_segment[s][0] <= posi_a1 && posi_a1 <= _segment[s][1]) {
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

void Cross::BackToPa1(Indi& kid, int appliedCycle)
{
    int aa, bb, a1, b1;
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

    ChangeSol(kid, appliedCycle, true /*reverse*/, false /*updateSeg*/);
}

void Cross::GoToBest(Indi& kid, int bestAppliedCycle)
{
    int aa, bb, a1, b1;

    ChangeSol(kid, bestAppliedCycle, false /*reverse*/, false /*updateSeg*/);

    for (int s = 0; s < _numBestModiEdge; ++s) {
        aa = _bestModiEdge[s][0];
        bb = _bestModiEdge[s][1];
        a1 = _bestModiEdge[s][2];
        b1 = _bestModiEdge[s][3];

        if (kid._link[aa][0] == bb) {
            kid._link[aa][0] = a1;
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
