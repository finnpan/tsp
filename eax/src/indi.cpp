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
#include <assert.h>
#include <random>
#include <string.h>

#include "indi.h"

void Utils::Permutation (int* a, int n)
{
    int r = 0;
    for (int i = 0; i < n; ++i) {
        a[i] = i;
    }
    for (int i = 0; i < n - 1; ++i) {
        r = rand() % (n - (i + 1)) + i + 1;
        std::swap(a[i], a[r]);
    }
}

Indi::Indi ()
{
    _n = 0;
    _link = NULL;
    _cost = 0;
}

Indi::~Indi ()
{
    for (int i = 0; i < _n; ++i) {
        delete[] _link[i];
    }
    delete[] _link;
}

void Indi::Define (int N)
{
    _n = N;

    _link = new int* [_n];
    for (int i = 0; i < _n; ++i) {
        _link[i] = new int [2];
    }
}

Indi& Indi::operator= (const Indi& src)
{
    _n = src._n;

    for (int i = 0; i < _n; ++i) {
        for (int j = 0; j < 2; ++j) {
            _link[i][j] = src._link[i][j];
        }
    }
    _cost = src._cost;

    return *this;
}

bool Indi::operator== (const Indi& src)
{
    int curr, next, pre;

    if (_n != src._n) {
        return false;
    }
    if (_cost != src._cost) {
        return false;
    }

    curr = 0;
    pre = -1;
    for (int i = 0; i < _n; ++i) {
        if (_link[curr][0] == pre) {
            next = _link[curr][1];
        } else {
            next = _link[curr][0];
        }

        if (src._link[curr][0] != next && src._link[curr][1] != next) {
            return false;
        }

        pre = curr;
        curr = next;
    }

    return true;
}

Evaluator::Evaluator ()
{
    _edgeDis = NULL;
    _nearCity = NULL;
    _nCity = 0;
    _nearNumMax = 50;
    _checkedN = nullptr;
}

Evaluator::~Evaluator ()
{
    for (int i = 0; i < _nCity; ++i) {
        delete [] _edgeDis[i];
    }
    delete [] _edgeDis;
    for (int i = 0; i < _nCity; ++i) {
        delete [] _nearCity[i];
    }
    delete [] _nearCity;

    delete [] _x;
    delete [] _y;

    if (_checkedN) {
        delete [] _checkedN;
    }
}

void Evaluator::SetInstance (const char filename[])
{
    FILE* fp;
    int n;
    char word[80], type[80];

    fp = fopen(filename, "r");
    if (!fp) {
        printf("%s open failed\n", filename);
        assert(0);
    }

    ////// read instance //////
    while (1) {
        if (fscanf(fp, "%s", word) == EOF) {
            break;
        }

        if (strcmp(word, "DIMENSION") == 0) {
            fscanf(fp, "%s", word);
            assert(strcmp(word, ":") == 0);
            fscanf(fp, "%d", &_nCity);
        }

        if (strcmp(word, "EDGE_WEIGHT_TYPE") == 0) {
            fscanf(fp, "%s", word);
            assert(strcmp(word, ":") == 0);
            fscanf(fp, "%s", type);
        }

        if(strcmp(word, "NODE_COORD_SECTION") == 0) {
            break;
        }
    }
    if (strcmp(word, "NODE_COORD_SECTION") != 0) {
        printf("Error in reading the instance\n");
        exit(0);
    }

    _x = new double [_nCity];
    _y = new double [_nCity];
    _checkedN = new int[_nCity];

    for (int i = 0; i < _nCity; ++i) {
        fscanf(fp, "%d", &n);
        assert(i+1 == n);
        fscanf(fp, "%s", word);
        _x[i] = atof(word);
        fscanf(fp, "%s", word);
        _y[i] = atof(word);
    }

    fclose(fp);
    //////////////////////////

    _edgeDis = new EvalType* [_nCity];
    for (int i = 0; i < _nCity; ++i) { _edgeDis[i] = new EvalType[_nCity]; }
    _nearCity = new int* [_nCity];
    for (int i = 0; i < _nCity; ++i) { _nearCity[i] = new int [_nearNumMax+1]; }

    double r = 0;
    if (strcmp(type, "EUC_2D") == 0) {
        for (int i = 0; i < _nCity ; ++i) {
            for (int j = 0; j < _nCity ; ++j) {
                _edgeDis[i][j] = (EvalType)(
                    sqrt((_x[i]-_x[j])*(_x[i]-_x[j])+(_y[i]-_y[j])*(_y[i]-_y[j]))+0.5
                );
            }
        }
    } else if (strcmp(type, "ATT") == 0) {
        for (int i = 0; i < _nCity; ++i) {
            for (int j = 0; j < _nCity; ++j) {
                r = sqrt(((_x[i]-_x[j])*(_x[i]-_x[j])+(_y[i]-_y[j])*(_y[i]-_y[j]))/10.0);
                EvalType t = (EvalType)r;
                if ((double)t < r) {
                    _edgeDis[i][j] = t+1;
                } else {
                    _edgeDis[i][j] = t;
                }
            }
        }
    } else if (strcmp(type, "CEIL_2D") == 0) {
        for(int i = 0; i < _nCity ; ++i) {
            for (int j = 0; j < _nCity ; ++j) {
                _edgeDis[i][j] = (EvalType)ceil(
                        sqrt((_x[i]-_x[j])*(_x[i]-_x[j])+(_y[i]-_y[j])*(_y[i]-_y[j]))
                );
            }
        }
    } else {
        printf("EDGE_WEIGHT_TYPE is not supported\n");
        exit(1);
    }

    int ci;
    int j1 ,j2 ,j3;
    int city_num = 0;
    EvalType min_dis;

    for (ci = 0; ci < _nCity; ++ci) {
        for (j3 = 0; j3 < _nCity; ++j3) { _checkedN[j3] = 0; }
        _checkedN[ci] = 1;
        _nearCity[ci][0] = ci;
        for (j1 = 1; j1 <= _nearNumMax; ++j1) {
            min_dis = std::numeric_limits<EvalType>::max();
            for (j2 = 0; j2 < _nCity; ++j2) {
                if (_edgeDis[ci][j2] <= min_dis && _checkedN[j2] == 0) {
                    city_num = j2;
                    min_dis = _edgeDis[ci][j2];
                }
            }
            _nearCity[ci][j1] = city_num;
            _checkedN[city_num] = 1;
        }
    }
}

void Evaluator::DoIt (Indi& indi) const
{
    EvalType d = 0;
    for (int i = 0; i < _nCity; ++i) {
        d = d + _edgeDis[i][indi._link[i][0]];
        d = d + _edgeDis[i][indi._link[i][1]];
    }
    indi._cost = d/2;
}

bool Evaluator::CheckValid (int* array, EvalType value)
{
    int* check = _checkedN;
    EvalType distance;

    for (int i = 0; i < _nCity; ++i) {
        check[i] = 0;
    }

    for (int i = 0; i < _nCity; ++i) {
        ++check[array[i]-1];
    }

    for (int i = 0; i < _nCity; ++i) {
        if (check[i] != 1) {
            return false;
        }
    }

    distance = 0;
    for (int i = 0; i < _nCity-1; ++i) {
        distance += _edgeDis[array[i]-1][array[i+1]-1];
    }
    distance += _edgeDis[array[_nCity-1]-1][array[0]-1];
    if (distance != value) {
        return false;
    }
    return true;
}
