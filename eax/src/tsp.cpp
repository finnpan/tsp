/*
License: see tsp.h
 */

#include "tsp.h"

Indi::Indi () :
    _numCity(0),
    _link(nullptr),
    _cost(0)
{}

Indi::~Indi ()
{
    for (int i = 0; i < _numCity; ++i) {
        delete[] _link[i];
    }
    delete[] _link;
}

void Indi::Define (int n)
{
    _numCity = n;
    _link = new int* [n];
    for (int i = 0; i < n; ++i) {
        _link[i] = new int [2];
    }
}

Indi& Indi::operator= (const Indi& rhs)
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

bool Indi::operator== (const Indi& rhs) const
{
    int curr = 0, next = -1, prev = -1;

    if (_numCity != rhs._numCity) {
        return false;
    }

    for (int i = 0; i < _numCity; ++i) {
        if (_link[curr][0] == prev) {
            next = _link[curr][1];
        } else {
            next = _link[curr][0];
        }
        if (rhs._link[curr][0] != next && rhs._link[curr][1] != next) {
            return false;
        }
        prev = curr;
        curr = next;
    }

    return true;
}

Evaluator::Evaluator () :
    _rDev(new std::random_device()),
    _rand(new std::mt19937((*_rDev)())),
    _maxNumNear(50),
    _numCity(0),
    _cost(nullptr),
    _near(nullptr),
    _buf(nullptr),
    _x(nullptr),
    _y(nullptr)
{}

Evaluator::~Evaluator ()
{
    for (int i = 0; i < _numCity; ++i) {
        delete [] _cost[i];
    }
    delete [] _cost;

    for (int i = 0; i < _numCity; ++i) {
        delete [] _near[i];
    }
    delete [] _near;

    delete [] _buf;
    delete [] _x;
    delete [] _y;
    delete _rand;
    delete _rDev;
}

void Evaluator::DoIt (Indi& indi) const
{
    indi._cost = 0;
    for (int i = 0; i < _numCity; ++i) {
        indi._cost += _cost[i][indi._link[i][0]];
        indi._cost += _cost[i][indi._link[i][1]];
    }
    indi._cost /= 2;
}

void Evaluator::SetInstance (const char filename[])
{
    FILE* fp;
    int n;
    char word[80], type[80];

    fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: failed to open %s!\n", filename);
        assert(0);
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
        exit(0);
    }

    _buf = new int[_numCity];
    _x = new double [_numCity];
    _y = new double [_numCity];

    for (int i = 0; i < _numCity; ++i) {
        fscanf(fp, "%d", &n);
        assert(i+1 == n);
        fscanf(fp, "%s", word);
        _x[i] = atof(word);
        fscanf(fp, "%s", word);
        _y[i] = atof(word);
    }

    fclose(fp);

    _cost = new EvalType* [_numCity];
    for (int i = 0; i < _numCity; ++i) { _cost[i] = new EvalType[_numCity]; }
    _near = new int* [_numCity];
    for (int i = 0; i < _numCity; ++i) { _near[i] = new int [_maxNumNear+1]; }

    double r = 0;
    if (strcmp(type, "EUC_2D") == 0) {
        for (int i = 0; i < _numCity ; ++i) {
            for (int j = 0; j < _numCity ; ++j) {
                r = sqrt((_x[i]-_x[j])*(_x[i]-_x[j])+(_y[i]-_y[j])*(_y[i]-_y[j]))+0.5;
                _cost[i][j] = (EvalType)(r);
            }
        }
    } else if (strcmp(type, "ATT") == 0) {
        for (int i = 0; i < _numCity; ++i) {
            for (int j = 0; j < _numCity; ++j) {
                r = sqrt(((_x[i]-_x[j])*(_x[i]-_x[j])+(_y[i]-_y[j])*(_y[i]-_y[j]))/10.0);
                EvalType t = (EvalType)r;
                if ((double)t < r) {
                    _cost[i][j] = t+1;
                } else {
                    _cost[i][j] = t;
                }
            }
        }
    } else if (strcmp(type, "CEIL_2D") == 0) {
        for(int i = 0; i < _numCity ; ++i) {
            for (int j = 0; j < _numCity ; ++j) {
                r = sqrt((_x[i]-_x[j])*(_x[i]-_x[j])+(_y[i]-_y[j])*(_y[i]-_y[j]));
                _cost[i][j] = (EvalType)ceil(r);
            }
        }
    } else {
        printf("Error: EDGE_WEIGHT_TYPE is not supported\n");
        exit(1);
    }

    int nearCity = 0;
    EvalType minCost;

    for (int ci = 0; ci < _numCity; ++ci) {
        for (int j = 0; j < _numCity; ++j) { _buf[j] = 0; }
        _buf[ci] = 1;
        _near[ci][0] = ci;
        for (int ni = 1; ni <= _maxNumNear; ++ni) {
            minCost = std::numeric_limits<EvalType>::max();
            for (int j = 0; j < _numCity; ++j) {
                if (_cost[ci][j] <= minCost && _buf[j] == 0) {
                    nearCity = j;
                    minCost = _cost[ci][j];
                }
            }
            _near[ci][ni] = nearCity;
            _buf[nearCity] = 1;
        }
    }
}

Kopt::Kopt (const Evaluator* e) :
    _numCity(e->_numCity),
    _eval(e),
    _tree(new TwoLevelTree(e->_numCity)),
    _maxNumINL(500)
{
    const int n = _numCity;
    _numINL = new int [n];
    _invNearList = new int* [n];
    for (int i = 0; i < n; ++i) {
        _invNearList[i] = new int [_maxNumINL];
    }

    SetInvNearList();
}

Kopt::~Kopt ()
{
    delete _tree;
    delete [] _numINL;
    for (int i = 0; i < _numCity; ++i) {
        delete [] _invNearList[i];
    }
    delete [] _invNearList;
}

void Kopt::DoIt (Indi& indi)
{
    MakeRandSol(indi);
    TransIndiToTree(indi);
    Local_search_2_opt_neighborhood();
    TransTreeToIndi(indi);
}

void Kopt::MakeRandSol (Indi& indi) const
{
    const int n = _numCity;
    int* route = _eval->_buf;

    for (int i = 0; i < n; ++i) {
        route[i] = i;
    }
    std::shuffle(route, route+n, *_eval->_rand);

    for (int i = 1 ; i < n-1; ++i) {
        indi._link[route[i]][0] = route[i-1];
        indi._link[route[i]][1] = route[i+1];
    }
    indi._link[route[0]][0] = route[n-1];
    indi._link[route[0]][1] = route[1];
    indi._link[route[n-1]][0] = route[n-2];
    indi._link[route[n-1]][1] = route[0];

    _eval->DoIt(indi);
}

void Kopt::SetInvNearList ()
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
                printf("Warning: check _numINL[c] < %d in kopt.cpp\n", _maxNumINL);
                fflush(stdout);
            }
        }
    }
}

void Kopt::TransIndiToTree (const Indi& indi)
{
    const int n = _numCity;
    int* route = _eval->_buf;
    int c = 0;
    for (int i = 0; i < n; ++i) {
        route[i] = c;
        c = indi._link[c][1];
    }
    _tree->SetTour(route, route + n - 1);
}

void Kopt::TransTreeToIndi (Indi& indi) const
{
    for (int i = 0; i < _numCity; ++i) {
        indi._link[i][0] = _tree->GetPrev(i);
        indi._link[i][1] = _tree->GetNext(i);
    }
    _eval->DoIt(indi);
}

void Kopt::Local_search_2_opt_neighborhood ()
{
    const int n = _numCity;
    const int maxNumNear = _eval->_maxNumNear;
    EvalType d1, d2;
    int t[5];

    int* active = _eval->_buf;
    for (int i = 0; i < n; ++i) { active[i] = 1; }

    LLL1: t[0] = (*_eval->_rand)()%n;
    t[1] = t[0];

    while (1) {  // t1's loop
        t[1] = _tree->GetNext(t[1]);
        if (active[t[1]] == 0) { goto EEE; }

        t[2] = _tree->GetPrev(t[1]);
        for (int k = 1; k < maxNumNear; ++k) {
            t[4] = _eval->_near[t[1]][k];
            t[3] = _tree->GetPrev(t[4]);
            d1 = _eval->_cost[t[1]][t[2]] - _eval->_cost[t[1]][t[4]];

            if (d1 > 0) {
                d2 = d1 + _eval->_cost[t[3]][t[4]] - _eval->_cost[t[3]][t[2]];
                if (d2 > 0) {
                    _tree->Flip(t[1], t[2], t[4], t[3]);
                    for (int i = 1; i <= 4; ++i) {
                        for (int j = 0; j < _numINL[t[i]]; ++j) {
                            active[_invNearList[t[i]][j]] = 1;
                        }
                    }
                    goto LLL1;
                }
            } else {
                break;
            }
        }

        t[2] = _tree->GetNext(t[1]);
        for (int k = 1; k < maxNumNear; ++k) {
            t[4] = _eval->_near[t[1]][k];
            t[3] = _tree->GetNext(t[4]);
            d1 = _eval->_cost[t[1]][t[2]] - _eval->_cost[t[1]][t[4]];

            if (d1 > 0) {
                d2 = d1 + _eval->_cost[t[3]][t[4]] - _eval->_cost[t[3]][t[2]];
                if (d2 > 0) {
                    _tree->Flip(t[1], t[2], t[4], t[3]);
                    for (int i = 1; i <= 4; ++i) {
                        for (int j = 0; j < _numINL[t[i]]; ++j) {
                            active[_invNearList[t[i]][j]] = 1;
                        }
                    }
                    goto LLL1;
                }
            } else {
                break;
            }
        }

        active[t[1]] = 0;
    EEE:;
        if (t[1] == t[0]) { break; }
    }
}
