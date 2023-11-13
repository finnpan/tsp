/*
License: see tsp.h
 */

#include "tsp.h"

Indi::Indi()
    : _numCity(0), _link(nullptr), _cost(std::numeric_limits<EvalType>::max())
{
}

Indi::~Indi()
{
    for (int i = 0; i < _numCity; ++i) {
        delete[] _link[i];
    }
    delete[] _link;
}

void Indi::Define(int n)
{
    _numCity = n;
    _link = new int*[n];
    for (int i = 0; i < n; ++i) {
        _link[i] = new int[2];
    }
}

Indi& Indi::operator=(const Indi& rhs)
{
    if (this != &rhs) {
        _numCity = rhs._numCity;
        for (int i = 0; i < _numCity; ++i) {
            for (int j = 0; j < 2; ++j) {
                _link[i][j] = rhs._link[i][j];
            }
        }
        _cost = rhs._cost;
    }
    return *this;
}

void Indi::ToArr(int* arr, int* arrInv) const
{
    int curr = -1, next = 0, prev = -1;
    for (int i = 0; i < _numCity; ++i) {
        prev = curr;
        curr = next;
        if (_link[curr][0] != prev) {
            next = _link[curr][0];
        } else {
            next = _link[curr][1];
        }
        arr[i] = curr;
        if (arrInv) {
            arrInv[curr] = i;
        }
    }
}

Evaluator::Evaluator()
    : _rDev(new std::random_device()), _rand(new std::mt19937((*_rDev)())),
      _maxNumNear(50), _numCity(0), _cost(nullptr), _near(nullptr),
      _routeBuf(nullptr), _x(nullptr), _y(nullptr)
{
}

Evaluator::~Evaluator()
{
    Deallocate(_numCity);
    delete _rand;
    delete _rDev;
}

void Evaluator::Allocate(int n)
{
    _cost = new EvalType*[n];
    for (int i = 0; i < n; ++i) {
        _cost[i] = new EvalType[n];
    }

    _near = new int*[n];
    for (int i = 0; i < n; ++i) {
        _near[i] = new int[_maxNumNear + 1];
    }

    _routeBuf = new int[n];
    _x = new double[n];
    _y = new double[n];
}

void Evaluator::Deallocate(int n)
{
    for (int i = 0; i < n; ++i) {
        delete[] _cost[i];
    }
    delete[] _cost;

    for (int i = 0; i < n; ++i) {
        delete[] _near[i];
    }
    delete[] _near;

    delete[] _routeBuf;
    delete[] _x;
    delete[] _y;
}

void Evaluator::DoIt(Indi& indi) const
{
    indi._cost = 0;
    for (int i = 0; i < _numCity; ++i) {
        indi._cost += _cost[i][indi._link[i][1]];
    }
}

bool Evaluator::SetInstance(const char filename[])
{
    const int oldCityNum = _numCity;
    char word[80], type[80];

    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: failed to open %s!\n", filename);
        return false;
    }

    /* read instance */
    while (1) {
        if (fscanf(fp, "%s", word) == EOF) {
            break;
        }
        if (strcmp(word, "DIMENSION") == 0) {
            fscanf(fp, "%s", word);
            assert(strcmp(word, ":") == 0);
            fscanf(fp, "%d", &_numCity);
        } else if (strcmp(word, "EDGE_WEIGHT_TYPE") == 0) {
            fscanf(fp, "%s", word);
            assert(strcmp(word, ":") == 0);
            fscanf(fp, "%s", type);
        } else if (strcmp(word, "NODE_COORD_SECTION") == 0) {
            break;
        }
    }
    if (strcmp(word, "NODE_COORD_SECTION") != 0) {
        printf("Error: in reading the instance\n");
        return false;
    }

    if (_numCity <= 0) {
        printf("Error: invalid city number in tsp file!\n");
        return false;
    }

    if (_numCity != oldCityNum) {
        Deallocate(oldCityNum);
        Allocate(_numCity);
    }

    int n;
    for (int i = 0; i < _numCity; ++i) {
        fscanf(fp, "%d", &n);
        assert(i + 1 == n);
        fscanf(fp, "%s", word);
        _x[i] = atof(word);
        fscanf(fp, "%s", word);
        _y[i] = atof(word);
    }

    fclose(fp);

    double r = 0;
    if (strcmp(type, "EUC_2D") == 0) {
        for (int i = 0; i < _numCity; ++i) {
            for (int j = 0; j < _numCity; ++j) {
                r = sqrt((_x[i] - _x[j]) * (_x[i] - _x[j])
                         + (_y[i] - _y[j]) * (_y[i] - _y[j]))
                    + 0.5;
                _cost[i][j] = (EvalType)(r);
            }
        }
    } else if (strcmp(type, "ATT") == 0) {
        for (int i = 0; i < _numCity; ++i) {
            for (int j = 0; j < _numCity; ++j) {
                r = sqrt(((_x[i] - _x[j]) * (_x[i] - _x[j])
                          + (_y[i] - _y[j]) * (_y[i] - _y[j]))
                         / 10.0);
                EvalType t = (EvalType)r;
                if ((double)t < r) {
                    _cost[i][j] = t + 1;
                } else {
                    _cost[i][j] = t;
                }
            }
        }
    } else if (strcmp(type, "CEIL_2D") == 0) {
        for (int i = 0; i < _numCity; ++i) {
            for (int j = 0; j < _numCity; ++j) {
                r = sqrt((_x[i] - _x[j]) * (_x[i] - _x[j])
                         + (_y[i] - _y[j]) * (_y[i] - _y[j]));
                _cost[i][j] = (EvalType)ceil(r);
            }
        }
    } else {
        printf("Error: EDGE_WEIGHT_TYPE is not supported\n");
        return false;
    }

    int nearCity = 0;
    EvalType minCost;

    for (int ci = 0; ci < _numCity; ++ci) {
        for (int j = 0; j < _numCity; ++j) {
            _routeBuf[j] = 0;
        }
        _routeBuf[ci] = 1;
        _near[ci][0] = ci;
        for (int ni = 1; ni <= _maxNumNear; ++ni) {
            minCost = std::numeric_limits<EvalType>::max();
            for (int j = 0; j < _numCity; ++j) {
                if (_cost[ci][j] <= minCost && _routeBuf[j] == 0) {
                    nearCity = j;
                    minCost = _cost[ci][j];
                }
            }
            _near[ci][ni] = nearCity;
            _routeBuf[nearCity] = 1;
        }
    }

    return true;
}

void Evaluator::MakeRand(Indi& indi) const
{
    const int n = _numCity;
    int* route = _routeBuf;

    for (int i = 0; i < n; ++i) {
        route[i] = i;
    }
    std::shuffle(route, route + n, *_rand);

    for (int i = 1; i < n - 1; ++i) {
        indi._link[route[i]][0] = route[i - 1];
        indi._link[route[i]][1] = route[i + 1];
    }
    indi._link[route[0]][0] = route[n - 1];
    indi._link[route[0]][1] = route[1];
    indi._link[route[n - 1]][0] = route[n - 2];
    indi._link[route[n - 1]][1] = route[0];

    DoIt(indi);
}

KOpt::KOpt(const Evaluator* e)
    : _eval(e), _numCity(e->_numCity), _tree(new TwoLevelTree(e->_numCity)),
      _maxNumINL(500)
{
    const int n = _numCity;
    _numINL = new int[n];
    _invNearList = new int*[n];
    for (int i = 0; i < n; ++i) {
        _invNearList[i] = new int[_maxNumINL];
    }

    SetInvNearList();
}

KOpt::~KOpt()
{
    delete _tree;
    delete[] _numINL;
    for (int i = 0; i < _numCity; ++i) {
        delete[] _invNearList[i];
    }
    delete[] _invNearList;
}

void KOpt::DoIt(Indi& indi)
{
    _eval->MakeRand(indi);
    TransIndiToTree(indi);
    Local_search_2_opt_neighborhood();
    TransTreeToIndi(indi);
}

void KOpt::SetInvNearList()
{
    const int maxNumNear = _eval->_maxNumNear;
    const int n = _numCity;
    int c;

    for (int i = 0; i < n; ++i) {
        _numINL[i] = 0;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < maxNumNear; ++j) {
            c = _eval->_near[i][j];
            if (_numINL[c] < _maxNumINL) {
                _invNearList[c][_numINL[c]++] = i;
            } else {
                printf("Warning: check _numINL[c] < %d in kopt.cpp\n",
                       _maxNumINL);
                fflush(stdout);
            }
        }
    }
}

void KOpt::TransIndiToTree(const Indi& indi)
{
    int* route = _eval->_routeBuf;
    indi.ToArr(route);
    _tree->SetTour(route, route + _numCity - 1);
}

void KOpt::TransTreeToIndi(Indi& indi) const
{
    for (int i = 0; i < _numCity; ++i) {
        indi._link[i][0] = _tree->GetPrev(i);
        indi._link[i][1] = _tree->GetNext(i);
    }
    _eval->DoIt(indi);
}

void KOpt::Local_search_2_opt_neighborhood()
{
    const int n = _numCity;
    const int maxNumNear = _eval->_maxNumNear;
    EvalType d1, d2;
    int t[5];

    int* active = _eval->_routeBuf;
    for (int i = 0; i < n; ++i) {
        active[i] = 1;
    }

BBB:
    t[0] = (*_eval->_rand)() % n;
    t[1] = t[0];

    while (1) {  // t1's loop
        t[1] = _tree->GetNext(t[1]);
        if (active[t[1]] == 0) {
            goto EEE;
        }

        for (int rev = 0; rev < 2; rev++) {
            t[2] = (rev == 0) ? _tree->GetPrev(t[1]) : _tree->GetNext(t[1]);
            for (int k = 1; k < maxNumNear; ++k) {
                t[3] = _eval->_near[t[1]][k];
                t[4] = (rev == 0) ? _tree->GetPrev(t[3]) : _tree->GetNext(t[3]);
                d1 = _eval->_cost[t[1]][t[2]] - _eval->_cost[t[1]][t[3]];
                if (d1 > 0) {
                    d2 = d1 + _eval->_cost[t[3]][t[4]]
                         - _eval->_cost[t[2]][t[4]];
                    if (d2 > 0) {
                        _tree->Flip(t[1], t[2], t[3], t[4]);
                        for (int i = 1; i <= 4; ++i) {
                            for (int j = 0; j < _numINL[t[i]]; ++j) {
                                active[_invNearList[t[i]][j]] = 1;
                            }
                        }
                        goto BBB;
                    }
                } else {
                    break;
                }
            }
        }
        active[t[1]] = 0;

    EEE:
        if (t[1] == t[0]) {
            break;
        };
    }
}
