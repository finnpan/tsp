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

#include <cstdio>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <random>

#include "flipper.h"

using EvalType = int32_t;

class Indi {
public:
    Indi();
    ~Indi();

    Indi(const Indi&) = delete;
    Indi& operator=(const Indi&);
    Indi(Indi&&) = delete;
    Indi& operator=(Indi&&) = delete;

    void Define(int n);
    void ToArr(int* arr, int* arrInv = nullptr) const;

    int _numCity;
    int** _link;
    EvalType _cost;
};

class Evaluator {
public:
    Evaluator();
    ~Evaluator();

    void DoIt(Indi& indi) const;
    bool SetInstance(const char filename[]);
    void MakeRand(Indi& indi) const;

    std::random_device* const _rDev;
    std::mt19937* const _rand;
    const int _maxNumNear;
    int _numCity;
    EvalType** _cost;
    int** _near;
    int* _routeBuf;

private:
    void Allocate(int n);
    void Deallocate(int n);

private:
    double* _x;
    double* _y;
};

class KOpt {
public:
    KOpt(const Evaluator* e);
    ~KOpt();

    void DoIt(Indi& indi);

private:
    void SetInvNearList();
    void TransIndiToTree(const Indi& indi);
    void TransTreeToIndi(Indi& indi);
    void Local_search_2_opt_neighborhood();

private:
    const Evaluator* const _eval;
    const int _numCity;
    const int _maxNumINL;
    int* _numINL;
    int** _invNearList;
    CClk_flipper _flipper;
};
