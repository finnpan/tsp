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

class TwoLevelTree;
class Kopt {
public:
    Kopt (int N, const Evaluator* e);
    ~Kopt ();

    /* Apply a local search with the 2-opt neighborhood */
    void DoIt (Indi& Indi);
    /* Set a random tour */
    void MakeRandSol (Indi& indi);

private:
    void TransIndiToTree (const Indi& indi);
    void TransTreeToIndi (Indi& indi);
    void Sub ();

    void SetInvNearList ();

    void TransIndiToTree_tlt (const Indi& indi);
    void TransTreeToIndi_tlt (Indi& indi);
    void Sub_tlt ();

    int GetNext (int t);
    int GetPrev (int t);
    void IncrementImp ();
    void CombineSeg (int segL, int segS);

    void CheckDetail ();
    void CheckValid ();

    int Turn (int &orient);

private:
    int _n;

    int _fixNumOfSeg;
    int _numSeg;
    int **_link;
    int *_segCity;
    int *_ordCity;

    int *_ordSeg;
    int *_orient;
    int **_linkSeg;
    int *_sizeSeg;
    int **_citySeg;

    int *_t;
    int _flagRev;

    int *_activeV;
    int **_invNearList;
    int *_numOfINL;

    int *_array;
    int *_checkN;
    int *_gene;
    int *_b;

    const Evaluator* eval;
    TwoLevelTree* _tree;
};
