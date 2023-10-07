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

using EvalType = long long;

namespace Utils {
	void Permutation (int* a, int n);
};

class Indi {
public:
    Indi ();
    ~Indi ();
    void Define (int N);
    Indi& operator= (const Indi& src);
    bool operator== (const Indi& indi2);

    int _n;                    /* Number of cities */
    int** _link;               /* _link[i][]: two vertices adjacent to i */
    EvalType _cost; /* Tour length of */
};

class Evaluator {
public:
    Evaluator ();
    ~Evaluator ();
    void SetInstance (const char filename[]);
    void DoIt (Indi& indi) const;
    bool CheckValid (int* array, EvalType value);

    int _nearNumMax;       /* Maximum number of k (see below) */
    int **_nearCity;       /* NearCity[i][k]: k-th nearest city from */
    EvalType  **_edgeDis;  /* EdgeDis[i][j]: distance between i and j */
    int _nCity;             /* Number of cities */
    double *_x;             /* x[i]: x-coordinate of */
    double *_y;             /* y[i]: x-coordinate of */
    int* _checkedN;
};
