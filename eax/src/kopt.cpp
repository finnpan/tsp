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

#include <iostream>
#include <random>
#include <assert.h>

#include "twoleveltree.h"
#include "kopt.h"

Kopt::Kopt (int N, const Evaluator* e)
{
    _n = N;
    eval = e;

    _link = new int* [_n];
    for (int i = 0; i < _n; ++i) {
        _link[i] = new int [2];
    }
    _ordCity = new int [_n];
    _ordSeg = new int [_n];
    _segCity = new int [_n];
    _orient = new int [_n];
    _linkSeg = new int* [_n];
    for (int i = 0; i < _n; ++i) {
        _linkSeg[i] = new int [2];
    }
    _sizeSeg = new int [_n];
    _citySeg = new int* [_n];
    for (int i = 0; i < _n; ++i) {
        _citySeg[i] = new int [2];
    }

    _t = new int [5];

    _activeV = new int [_n];
    _invNearList = new int* [_n];
    for (int i = 0; i < _n; ++i) {
        _invNearList[i] = new int [500];
    }
    _numOfINL = new int [_n];

    _array = new int [_n+2];
    _checkN = new int [_n];
    _b = new int [_n];
    _gene = new int [_n];

    _tree = new TwoLevelTree(_n);

    SetInvNearList();
}

Kopt::~Kopt ()
{
    for (int i = 0; i < _n; ++i) {
        delete [] _link[i];
    }
    delete [] _link;

    for (int i = 0; i < _n; ++i) {
        delete [] _linkSeg[i];
    }
    delete [] _linkSeg;

    for (int i = 0; i < _n; ++i) {
        delete [] _citySeg[i];
    }
    delete [] _citySeg;

    for (int i = 0; i < _n; ++i) {
        delete [] _invNearList[i];
    }
    delete [] _invNearList;

    delete [] _ordCity;
    delete [] _ordSeg;
    delete [] _segCity;
    delete [] _orient;
    delete [] _sizeSeg;
    delete [] _t;
    delete [] _activeV;
    delete [] _numOfINL;
    delete [] _array;
    delete [] _checkN;
    delete [] _b;
    delete [] _gene;

    delete _tree;
}

void Kopt::SetInvNearList ()
{
    assert(eval->_nearNumMax == 50);

    int c;
    for (int i = 0; i < _n; ++i) {
        _numOfINL[i] = 0;
    }

    for (int i = 0; i < _n; ++i) {
        for (int k = 0; k < 50; ++k) {
            c = eval->_nearCity[i][k];
            if (_numOfINL[c] < 500) {
                _invNearList[c][_numOfINL[c]++] = i;
            } else {
                printf("Check _numOfINL[c] < 500) in kopt.cpp \n");
                fflush(stdout);
            }
        }
    }
}

void Kopt::TransIndiToTree_tlt (const Indi& indi)
{
    int cur = 0;
    for (int t = 0; t < _n; ++t) {
        _checkN[t] = cur;
        cur = indi._link[cur][1];
    }
    _tree->SetTour(_checkN, _checkN + _n - 1);
}

void Kopt::TransTreeToIndi_tlt (Indi& indi)
{
    for (int t = 0; t < _n; ++t) {
        indi._link[t][0] = _tree->GetPrev(t);
        indi._link[t][1] = _tree->GetNext(t);
    }
    eval->DoIt(indi);
}

void Kopt::Sub_tlt ()
{
    EvalType d1, d2;
    int t[5];

    for (int t = 0; t < _n; ++t) { _activeV[t] = 1; }

    LLL1: t[0] = rand()%_n;
    t[1] = t[0];

    while (1) {  // t1's loop
        t[1] = _tree->GetNext(t[1]);
        if (_activeV[t[1]] == 0) { goto EEE; }

        t[2] = _tree->GetPrev(t[1]);
        for (int num1 = 1; num1 < 50; ++num1) {
            t[4] = eval->_nearCity[t[1]][num1];
            t[3] = _tree->GetPrev(t[4]);
            d1 = eval->_edgeDis[t[1]][t[2]] - eval->_edgeDis[t[1]][t[4]];

            if (d1 > 0) {
                d2 = d1 + eval->_edgeDis[t[3]][t[4]] - eval->_edgeDis[t[3]][t[2]];
                if (d2 > 0) {
                    _tree->Flip(t[1], t[2], t[4], t[3]);
                    for (int a = 1; a <= 4; ++a) {
                        for (int k = 0; k < _numOfINL[t[a]]; ++k) {
                            _activeV[_invNearList[t[a]][k]] = 1;
                        }
                    }
                    goto LLL1;
                }
            } else {
                break;
            }
        }

        t[2] = _tree->GetNext(t[1]);
        for (int num1 = 1; num1 < 50; ++num1) {
            t[4] = eval->_nearCity[t[1]][num1];
            t[3] = _tree->GetNext(t[4]);
            d1 = eval->_edgeDis[t[1]][t[2]] - eval->_edgeDis[t[1]][t[4]];

            if (d1 > 0) {
                d2 = d1 + eval->_edgeDis[t[3]][t[4]] - eval->_edgeDis[t[3]][t[2]];
                if (d2 > 0) {
                    _tree->Flip(t[1], t[2], t[4], t[3]);
                    for (int a = 1; a <= 4; ++a) {
                        for (int k = 0; k < _numOfINL[t[a]]; ++k) {
                            _activeV[_invNearList[t[a]][k]] = 1;
                        }
                    }
                    goto LLL1;
                }
            } else {
                break;
            }
        }

        _activeV[t[1]] = 0;
    EEE:;
        if (t[1] == t[0]) { break; }
    }
}

void Kopt::TransIndiToTree (const Indi& indi)
{
    int num;
    int size;
    int orient;

    _array[1] = 0;
    for (int i = 2; i <= _n; ++i) {
        _array[i] = indi._link[_array[i-1]][1];
    }
    _array[0] = _array[_n];
    _array[_n+1] = _array[1];

    num = 1;
    _numSeg = 0;

    while (1) {
        orient = 1;
        size = 0;

        _orient[_numSeg] = orient;
        _ordSeg[_numSeg] = _numSeg;

        _link[_array[num]][0] = -1;
        _link[_array[num]][1] = _array[num+1];
        _ordCity[_array[num]] = size;
        _segCity[_array[num]] = _numSeg;
        _citySeg[_numSeg][this->Turn(orient)] = _array[num];
        ++num;
        ++size;

        for (int i = 0; i < (int)sqrt(_n)-1; ++i) {
            if (num == _n) {
                break;
            }
            _link[_array[num]][0] = _array[num-1];
            _link[_array[num]][1] = _array[num+1];
            _ordCity[_array[num]] = size;
            _segCity[_array[num]] = _numSeg;
            ++num;
            ++size;
        }

        if (num == _n-1) {
            _link[_array[num]][0] = _array[num-1];
            _link[_array[num]][1] = _array[num+1];
            _ordCity[_array[num]] = size;
            _segCity[_array[num]] = _numSeg;
            ++num;
            ++size;
        }

        _link[_array[num]][0] = _array[num-1];
        _link[_array[num]][1] = -1;
        _ordCity[_array[num]] = size;
        _segCity[_array[num]] = _numSeg;
        _citySeg[_numSeg][orient] = _array[num];
        ++num;
        ++size;

        _sizeSeg[_numSeg] = size;
        ++_numSeg;

        if (num == _n+1) {
            break;
        }
    }

    for (int s = 1; s < _numSeg-1; ++s) {
        _linkSeg[s][0] = s-1;
        _linkSeg[s][1] = s+1;
    }
    _linkSeg[0][0] = _numSeg-1;
    _linkSeg[0][1] = 1;
    _linkSeg[_numSeg-1][0] = _numSeg-2;
    _linkSeg[_numSeg-1][1] = 0;

    _fixNumOfSeg = _numSeg;
}

void Kopt::TransTreeToIndi (Indi& indi)
{
    int t_p, t_n;
    for (int t = 0; t < _n; ++t) {
        t_p = this->GetPrev(t);
        t_n = this->GetNext(t);

        indi._link[t][0] = t_p;
        indi._link[t][1] = t_n;
    }
    eval->DoIt(indi);
}

void Kopt::DoIt (Indi& Indi)
{
#if 0
    this->TransIndiToTree(Indi);
    //  this->CheckDetail();           // Check
    //  this->CheckValid();            // Check
    this->Sub();
    this->TransTreeToIndi(Indi);
#else
    this->TransIndiToTree_tlt(Indi);
    this->Sub_tlt();
    this->TransTreeToIndi_tlt(Indi);
#endif
}

void Kopt::Sub ()
{
    int t1_st;
    EvalType dis1, dis2;

    for (int t = 0; t < _n; ++t) {
        _activeV[t] = 1;
    }

    LLL1: t1_st = rand()%_n;
    _t[1] = t1_st;

    while (1) {  // t1's loop
        _t[1] = this->GetNext(_t[1]);
        if (_activeV[_t[1]] == 0) {
            goto EEE;
        }

        ////
        _flagRev = 0;
        _t[2] = this->GetPrev(_t[1]);
        for (int num1 = 1; num1 < 50; ++num1) {
            _t[4] = eval->_nearCity[_t[1]][num1];
            _t[3] = this->GetPrev(_t[4]);
            dis1 = eval->_edgeDis[_t[1]][_t[2]] - eval->_edgeDis[_t[1]][_t[4]];

            if (dis1 > 0) {
                dis2 = dis1 +
                    eval->_edgeDis[_t[3]][_t[4]] - eval->_edgeDis[_t[3]][_t[2]];

                if (dis2 > 0) {
                    this->IncrementImp();
                    for (int a = 1; a <= 4; ++a) {
                        for (int k = 0; k < _numOfINL[_t[a]]; ++k) {
                            _activeV[this->_invNearList[_t[a]][k]] = 1;
                        }
                    }
                    goto LLL1;
                }
            } else {
                break;
            }
        }

        ////
        _flagRev = 1;
        _t[2] = this->GetNext(_t[1]);
        for (int num1 = 1; num1 < 50; ++num1) {
            _t[4] = eval->_nearCity[_t[1]][num1];
            _t[3] = this->GetNext(_t[4]);
            dis1 = eval->_edgeDis[_t[1]][_t[2]] - eval->_edgeDis[_t[1]][_t[4]];

            if (dis1 > 0) {
                dis2 = dis1 +
                        eval->_edgeDis[_t[3]][_t[4]] - eval->_edgeDis[_t[3]][_t[2]];

                if (dis2 > 0) {
                    this->IncrementImp();
                    for (int a = 1; a <= 4; ++a) {
                        for (int k = 0; k < _numOfINL[_t[a]]; ++k) {
                            _activeV[this->_invNearList[_t[a]][k]] = 1;
                        }
                    }
                    goto LLL1;
                }
            } else {
                break;
            }
        }

        _activeV[_t[1]] = 0;
    EEE:;
        if (_t[1] == t1_st) {
            break;
        }
    }
}

int Kopt::GetNext (int t)
{
    int t_n, seg, orient;

    seg = _segCity[t];
    orient = _orient[seg];

    t_n = _link[t][orient];
    if (t_n == -1) {
        seg = _linkSeg[seg][orient];
        orient = Turn(_orient[seg]);
        t_n = _citySeg[seg][orient];
    }
    return t_n;
}

int Kopt::GetPrev (int t)
{
    int t_p, seg, orient;

    seg = _segCity[t];
    orient = _orient[seg];

    t_p = _link[t][this->Turn(orient)];
    if (t_p == -1) {
        seg = _linkSeg[seg][Turn(orient)];
        orient = _orient[seg];
        t_p = _citySeg[seg][orient];
    }
    return t_p;
}

int Kopt::Turn (int &orient)
{
    assert(orient == 0 || orient == 1);
    if (orient == 0) {
        return 1;
    } else if (orient == 1) {
        return 0;
    } else {
        assert(1 == 2);
    }
    return 0;
}

void Kopt::IncrementImp ()
{
    int t1_s, t1_e, t2_s, t2_e;
    int seg_t1_s, seg_t1_e, seg_t2_s, seg_t2_e;
    int ordSeg_t1_s, ordSeg_t1_e, ordSeg_t2_s, ordSeg_t2_e;
    int orient_t1_s, orient_t1_e, orient_t2_s, orient_t2_e;
    int numOfSeg1, numOfSeg2;
    int curr;
    int ord;

    int flag_t2e_t1s;
    int flag_t2s_t1e;
    int length_t1s_seg;
    int length_t1e_seg;
    int seg;

    // Seg1: b->d path
    // Seg2: c->a path

    if (_flagRev == 0) {
        t1_s = _t[1];
        t1_e = _t[3];
        t2_s = _t[4];
        t2_e = _t[2];
    } else if (_flagRev == 1) {
        t1_s = _t[2];
        t1_e = _t[4];
        t2_s = _t[3];
        t2_e = _t[1];
    }

    seg_t1_s = _segCity[t1_s];
    ordSeg_t1_s = _ordSeg[seg_t1_s];
    orient_t1_s = _orient[seg_t1_s];
    seg_t1_e = _segCity[t1_e];
    ordSeg_t1_e = _ordSeg[seg_t1_e];
    orient_t1_e = _orient[seg_t1_e];
    seg_t2_s = _segCity[t2_s];
    ordSeg_t2_s = _ordSeg[seg_t2_s];
    orient_t2_s = _orient[seg_t2_s];
    seg_t2_e = _segCity[t2_e];
    ordSeg_t2_e = _ordSeg[seg_t2_e];
    orient_t2_e = _orient[seg_t2_e];

    //////////////////// Type1 ////////////////////////
    if ((seg_t1_s == seg_t1_e) &&
       (seg_t1_s == seg_t2_s) && (seg_t1_s == seg_t2_e)) {

        if ((_orient[seg_t1_s] == 1 && (_ordCity[t1_s] > _ordCity[t1_e])) ||
            (_orient[seg_t1_s] == 0 && (_ordCity[t1_s] < _ordCity[t1_e]))) {
            std::swap(t1_s, t2_s);
            std::swap(t1_e, t2_e);
            std::swap(seg_t1_s, seg_t2_s);
            std::swap(seg_t1_e, seg_t2_e);
            std::swap(ordSeg_t1_s, ordSeg_t2_s);
            std::swap(ordSeg_t1_e, ordSeg_t2_e);
            std::swap(orient_t1_s, orient_t2_s);
            std::swap(orient_t1_e, orient_t2_e);
        }

        curr = t1_s;
        ord = _ordCity[t1_e];

        while (1) {
            std::swap(_link[curr][0], _link[curr][1]);
            _ordCity[curr] = ord;
            if (curr == t1_e) {
                break;
            }
            curr = _link[curr][Turn(orient_t1_s)];
            if (orient_t1_s == 0) {
                ++ord;
            } else {
                --ord;
            }
        }

        _link[t2_e][orient_t1_s] = t1_e;
        _link[t2_s][Turn(orient_t1_s)] = t1_s;
        _link[t1_s][orient_t1_s] = t2_s;
        _link[t1_e][Turn(orient_t1_s)] = t2_e;

        //    this->CheckDetail();              // Check
        //    this->CheckValid();               // Check
        return;
    }
    //////////////////// Type1 ///////////////////////

    if (ordSeg_t1_e >= ordSeg_t1_s) {
        numOfSeg1 = ordSeg_t1_e - ordSeg_t1_s + 1;
    } else {
        numOfSeg1 = ordSeg_t1_e - ordSeg_t1_s + 1 + _numSeg;
    }
    if (ordSeg_t2_e >= ordSeg_t2_s) {
        numOfSeg2 = ordSeg_t2_e - ordSeg_t2_s + 1;
    } else {
        numOfSeg2 = ordSeg_t2_e - ordSeg_t2_s + 1 + _numSeg;
    }

    if (numOfSeg1 > numOfSeg2) {
        std::swap(numOfSeg1, numOfSeg2);
        std::swap(t1_s, t2_s);
        std::swap(t1_e, t2_e);
        std::swap(seg_t1_s, seg_t2_s);
        std::swap(seg_t1_e, seg_t2_e);
        std::swap(ordSeg_t1_s, ordSeg_t2_s);
        std::swap(ordSeg_t1_e, ordSeg_t2_e);
        std::swap(orient_t1_s, orient_t2_s);
        std::swap(orient_t1_e, orient_t2_e);
    }

    if (_link[t2_e][orient_t2_e] == -1) {
        flag_t2e_t1s = 1;
    } else {
        flag_t2e_t1s = 0;
    }
    if (_link[t2_s][this->Turn(orient_t2_s)] == -1) {
        flag_t2s_t1e = 1;
    } else {
        flag_t2s_t1e = 0;
    }

    length_t1s_seg = abs(_ordCity[t2_e] -
                _ordCity[_citySeg[seg_t2_e][orient_t2_e]]);
    length_t1e_seg = abs(_ordCity[t2_s] -
                _ordCity[_citySeg[seg_t2_s][this->Turn(orient_t2_s)]]);

    ///////////////////// Type2 /////////////////
    if (seg_t1_s == seg_t1_e) {
        if (flag_t2e_t1s == 1 && flag_t2s_t1e == 1) {
            orient_t1_s = Turn(_orient[seg_t1_s]);
            _orient[seg_t1_s] = orient_t1_s;
            _citySeg[seg_t1_s][orient_t1_s] = t1_s;
            _citySeg[seg_t1_s][Turn(orient_t1_s)] = t1_e;
            _linkSeg[seg_t1_s][orient_t1_s] = seg_t2_s;
            _linkSeg[seg_t1_s][Turn(orient_t1_s)] = seg_t2_e;

            //      this->CheckDetail();              // Check
            //      this->CheckValid();               // Check
            return;
        }

        if (flag_t2e_t1s == 0 && flag_t2s_t1e == 1) {
            curr = t1_e;
            ord = _ordCity[t1_s];
            while (1) {
                std::swap(_link[curr][0], _link[curr][1]);
                _ordCity[curr] = ord;
                if (curr == t1_s) {
                    break;
                }
                curr = _link[curr][orient_t2_e];
                if (orient_t2_e == 0) {
                    --ord;
                } else {
                    ++ord;
                }
            }

            _link[t2_e][orient_t2_e] = t1_e;
            _link[t1_s][orient_t2_e] = -1;
            _link[t1_e][Turn(orient_t2_e)] = t2_e;
            _citySeg[seg_t2_e][orient_t2_e] = t1_s;
            //      this->CheckDetail();              // Check
            //      this->CheckValid();               // Check
            return;
        }

        if (flag_t2e_t1s == 1 && flag_t2s_t1e == 0) {
            curr = t1_s;
            ord = _ordCity[t1_e];
            while (1) {
                std::swap(_link[curr][0], _link[curr][1]);
                _ordCity[curr] = ord;
                if (curr == t1_e) {
                    break;
                }
                curr = _link[curr][Turn(orient_t2_s)];
                if (orient_t2_s == 0) {
                    ++ord;
                } else {
                    --ord;
                }
            }

            _link[t2_s][Turn(orient_t2_s)] = t1_s;
            _link[t1_e][Turn(orient_t2_s)] = -1;
            _link[t1_s][orient_t2_s] = t2_s;
            _citySeg[seg_t2_s][Turn(orient_t2_s)] = t1_e;
            //      this->CheckDetail();              // Check
            //      this->CheckValid();               // Check
            return;
        }
    }

    ///////////////////// Type3 /////////////////
    if (flag_t2e_t1s == 1) {
        _linkSeg[seg_t1_s][Turn(orient_t1_s)] = seg_t2_s;
    } else {
        seg_t1_s = _numSeg++;
        orient_t1_s = orient_t2_e;
        _link[t1_s][Turn(orient_t1_s)] = -1;
        _link[_citySeg[seg_t2_e][orient_t2_e]][orient_t1_s] = -1;
        _orient[seg_t1_s] = orient_t1_s;
        _sizeSeg[seg_t1_s] = length_t1s_seg;
        _citySeg[seg_t1_s][Turn(orient_t1_s)] = t1_s;
        _citySeg[seg_t1_s][orient_t1_s] = _citySeg[seg_t2_e][orient_t2_e];
        _linkSeg[seg_t1_s][Turn(orient_t1_s)] = seg_t2_s;
        _linkSeg[seg_t1_s][orient_t1_s] = _linkSeg[seg_t2_e][orient_t2_e];
        seg = _linkSeg[seg_t2_e][orient_t2_e];
        _linkSeg[seg][Turn(_orient[seg])] = seg_t1_s;
    }

    if (flag_t2s_t1e == 1) {
        _linkSeg[seg_t1_e][orient_t1_e] = seg_t2_e;
    } else {
        seg_t1_e = _numSeg++;
        orient_t1_e = orient_t2_s;
        _link[t1_e][orient_t1_e] = -1;
        _link[_citySeg[seg_t2_s][Turn(orient_t2_s)]][Turn(orient_t1_e)] = -1;
        _orient[seg_t1_e] = orient_t1_e;
        _sizeSeg[seg_t1_e] = length_t1e_seg;
        _citySeg[seg_t1_e][orient_t1_e] = t1_e;
        _citySeg[seg_t1_e][Turn(orient_t1_e)] =
                    _citySeg[seg_t2_s][Turn(orient_t2_s)];
        _linkSeg[seg_t1_e][orient_t1_e] = seg_t2_e;
        _linkSeg[seg_t1_e][Turn(orient_t1_e)] =
                    _linkSeg[seg_t2_s][Turn(orient_t2_s)];
        seg = _linkSeg[seg_t2_s][Turn(orient_t2_s)];
        _linkSeg[seg][_orient[seg]] = seg_t1_e;
    }

    _link[t2_e][orient_t2_e] = -1;
    _sizeSeg[seg_t2_e] -= length_t1s_seg;
    _citySeg[seg_t2_e][orient_t2_e] = t2_e;
    _linkSeg[seg_t2_e][orient_t2_e] = seg_t1_e;
    _link[t2_s][Turn(orient_t2_s)] = -1;
    _sizeSeg[seg_t2_s] -= length_t1e_seg;
    _citySeg[seg_t2_s][Turn(orient_t2_s)] = t2_s;
    _linkSeg[seg_t2_s][Turn(orient_t2_s)] = seg_t1_s;

    seg = seg_t1_e;
    while (1) {
        _orient[seg] = Turn(_orient[seg]);
        if (seg == seg_t1_s) {
            break;
        }
        seg = _linkSeg[seg][_orient[seg]];
    }

    if (_sizeSeg[seg_t2_e] < length_t1s_seg) {
        seg = _linkSeg[seg_t2_e][Turn(_orient[seg_t2_e])];
        _linkSeg[seg][_orient[seg]] = seg_t1_s;
        seg = _linkSeg[seg_t2_e][_orient[seg_t2_e]];
        _linkSeg[seg][Turn(_orient[seg])] = seg_t1_s;
        seg = _linkSeg[seg_t1_s][Turn(_orient[seg_t1_s])];
        _linkSeg[seg][_orient[seg]] = seg_t2_e;
        seg = _linkSeg[seg_t1_s][_orient[seg_t1_s]];
        _linkSeg[seg][Turn(_orient[seg])] = seg_t2_e;

        std::swap(_orient[seg_t2_e], _orient[seg_t1_s]);
        std::swap(_sizeSeg[seg_t2_e], _sizeSeg[seg_t1_s]);
        std::swap(_citySeg[seg_t2_e][0], _citySeg[seg_t1_s][0]);
        std::swap(_citySeg[seg_t2_e][1], _citySeg[seg_t1_s][1]);
        std::swap(_linkSeg[seg_t2_e][0], _linkSeg[seg_t1_s][0]);
        std::swap(_linkSeg[seg_t2_e][1], _linkSeg[seg_t1_s][1]);
        std::swap(seg_t2_e, seg_t1_s);
    }

    if (_sizeSeg[seg_t2_s] < length_t1e_seg) {
        seg = _linkSeg[seg_t2_s][Turn(_orient[seg_t2_s])];
        _linkSeg[seg][_orient[seg]] = seg_t1_e;
        seg = _linkSeg[seg_t2_s][_orient[seg_t2_s]];
        _linkSeg[seg][Turn(_orient[seg])] = seg_t1_e;
        seg = _linkSeg[seg_t1_e][Turn(_orient[seg_t1_e])];
        _linkSeg[seg][_orient[seg]] = seg_t2_s;
        seg = _linkSeg[seg_t1_e][_orient[seg_t1_e]];
        _linkSeg[seg][Turn(_orient[seg])] = seg_t2_s;

        std::swap(_orient[seg_t2_s], _orient[seg_t1_e]);
        std::swap(_sizeSeg[seg_t2_s], _sizeSeg[seg_t1_e]);
        std::swap(_citySeg[seg_t2_s][0], _citySeg[seg_t1_e][0]);
        std::swap(_citySeg[seg_t2_s][1], _citySeg[seg_t1_e][1]);
        std::swap(_linkSeg[seg_t2_s][0], _linkSeg[seg_t1_e][0]);
        std::swap(_linkSeg[seg_t2_s][1], _linkSeg[seg_t1_e][1]);
        std::swap(seg_t2_s, seg_t1_e);
    }

    while (_numSeg > _fixNumOfSeg) {
        if (_sizeSeg[_linkSeg[_numSeg-1][0]] <
            _sizeSeg[_linkSeg[_numSeg-1][1]]) {
            this->CombineSeg(_linkSeg[_numSeg-1][0], _numSeg-1);
        } else {
            this->CombineSeg(_linkSeg[_numSeg-1][1], _numSeg-1);
        }
    }

    int ordSeg = 0;
    seg = 0;

    while (1) {
        _ordSeg[seg] = ordSeg;
        ++ordSeg;

        seg = _linkSeg[seg][_orient[seg]];
        if (seg == 0) {
            break;
        }
    }

    // this->CheckDetail();              // Check
    // this->CheckValid();               // Check
    return;
}

void Kopt::CombineSeg (int segL, int segS)
{
    int seg;
    int t_s, t_e, direction; t_s = 0; t_e = 0; direction = 0;
    int ord; ord = 0;
    int increment; increment = 0;
    int curr, next;

    if (_linkSeg[segL][_orient[segL]] == segS) {
        _link[_citySeg[segL][_orient[segL]]][_orient[segL]] =
        _citySeg[segS][Turn(_orient[segS])];

        _link[_citySeg[segS][Turn(_orient[segS])]][Turn(_orient[segS])] =
                    _citySeg[segL][_orient[segL]];
        ord = _ordCity[_citySeg[segL][_orient[segL]]];

        _citySeg[segL][_orient[segL]] = _citySeg[segS][_orient[segS]];
        _linkSeg[segL][_orient[segL]] = _linkSeg[segS][_orient[segS]];
        seg = _linkSeg[segS][_orient[segS]];
        _linkSeg[seg][Turn(_orient[seg])] = segL;

        t_s = _citySeg[segS][Turn(_orient[segS])];
        t_e = _citySeg[segS][_orient[segS]];
        direction = _orient[segS];

        if (_orient[segL] == 1) {
            increment = 1;
        } else {
            increment = -1;
        }
    } else if (_linkSeg[segL][Turn(_orient[segL])] == segS) {
        _link[_citySeg[segL][Turn(_orient[segL])]][Turn(_orient[segL])] =
        _citySeg[segS][_orient[segS]];

        _link[_citySeg[segS][_orient[segS]]][_orient[segS]] =
                    _citySeg[segL][Turn(_orient[segL])];
        ord = _ordCity[_citySeg[segL][Turn(_orient[segL])]];

        _citySeg[segL][Turn(_orient[segL])] =
                    _citySeg[segS][Turn(_orient[segS])];
        _linkSeg[segL][Turn(_orient[segL])] =
                    _linkSeg[segS][Turn(_orient[segS])];
        seg = _linkSeg[segS][Turn(_orient[segS])];
        _linkSeg[seg][_orient[seg]] = segL;

        t_s = _citySeg[segS][_orient[segS]];
        t_e = _citySeg[segS][Turn(_orient[segS])];
        direction = Turn(_orient[segS]);

        if (_orient[segL] == 1) {
            increment = -1;
        } else {
            increment = 1;
        }
    }
    curr = t_s;
    ord = ord + increment;

    while (1) {
        _segCity[curr] = segL;
        _ordCity[curr] = ord;

        next = _link[curr][direction];
        if (_orient[segL] != _orient[segS]) {
            std::swap(_link[curr][0], _link[curr][1]);
        }

        if (curr == t_e) {
            break;
        }
        curr = next;
        ord += increment;
    }
    _sizeSeg[segL] += _sizeSeg[segS];
    --_numSeg;
}

void Kopt::CheckDetail ()
{
    int seg, seg_p, seg_n;
    int ord, ord_p, ord_n;
    int orient;
    int curr;

    seg = 0;

    for (int s = 0; s < _numSeg; ++s) {
        seg = s;
        orient = _orient[seg];
        seg_p = _linkSeg[seg][this->Turn(orient)];
        seg_n = _linkSeg[seg][orient];

        ord = _ordSeg[seg];
        ord_p = ord - 1 ;
        if (ord_p < 0) {
            ord_p = _numSeg - 1;
        }
        ord_n = ord + 1;
        if (ord_n >= _numSeg) {
            ord_n = 0;
        }

        assert(ord_p == _ordSeg[seg_p]);
        assert(ord_n == _ordSeg[seg_n]);

        curr = _citySeg[s][0];
        int count = 0;

        while (1) {
            ++count;
            if (curr == _citySeg[s][1]) {
                break;
            }
            curr = _link[curr][1];
            assert(curr != -1);
        }
        assert(count == _sizeSeg[s]);
    }

    int t;
    int t_n, t_p, t_s, t_e;

    for (t = 0; t < _n; ++t) {
        seg = _segCity[t];
        orient = _orient[seg];
        t_s = _citySeg[seg][0];
        t_e = _citySeg[seg][1];

        t_p = _link[t][0];
        t_n = _link[t][1];

        if (t == t_s) {
            assert(t_p == -1);
        } else {
            assert(t_p != -1);
            assert(t == _link[t_p][1]);
            assert(seg == _segCity[t_p]);
            assert(_ordCity[t] == _ordCity[t_p] + 1);
        }

        if (t == t_e) {
            assert(t_n == -1);
        } else{
            assert(t_n != -1);
            assert(t == _link[t_n][0]);
            assert(seg == _segCity[t_n]);
            assert(_ordCity[t] == _ordCity[t_n] - 1);
        }
    }
}

void Kopt::CheckValid ()
{
    int t_st, t_c, t_n;
    int count;
    int Invalid = 0;

    for (int i = 0; i < _n; ++i) {
        _checkN[i] = 0;
    }

    t_st = rand() % _n;
    t_n = t_st;

    count = 0;
    while (1) {
        t_c = t_n;
        _checkN[t_c] = 1;
        ++count;

        t_n = this->GetNext(t_c);

        if (t_n == t_st) {
            break;
        }

        if (count == _n+1) {
            Invalid = 1;
            break;
        }
    }

    for (int i = 0; i < _n; ++i) {
        if (_checkN[i] == 0) {
           Invalid = 1;
        }
    }

    if (Invalid == 1) {
        printf("Invalid \n"); fflush(stdout);
        assert(1 == 2);
    }
}

void Kopt::MakeRandSol (Indi& indi)
{
    int r;

    for (int j = 0; j < _n; ++j) {
        _b[j] = j;
    }

    for (int i = 0; i < _n; ++i) {
        r = rand() % (_n-i);
        _gene[i] = _b[r];
        _b[r] = _b[_n-i-1];
    }

    for (int j2 = 1 ; j2 < _n-1; ++j2) {
        indi._link[_gene[j2]][0] = _gene[j2-1];
        indi._link[_gene[j2]][1] = _gene[j2+1];
    }
    indi._link[_gene[0]][0] = _gene[_n-1];
    indi._link[_gene[0]][1] = _gene[1];
    indi._link[_gene[_n-1]][0] = _gene[_n-2];
    indi._link[_gene[_n-1]][1] = _gene[0];

    eval->DoIt(indi);
}
