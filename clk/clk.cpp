/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*               CODE FOR TESTING ITERATED LIN-KERNIGHAN                    */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 17, 1995                                                    */
/*                                                                          */
/*  For a short describtion see usage ()                                    */
/*                                                                          */
/****************************************************************************/

#include "linkern.h"

#define BIGDOUBLE (1e30)

static int run_silently = 1;
static int number_runs = 0;
static char *nodefile = (char *) NULL;


int main (int ac, char **av)
{
    int rval = 0;
    int ncount;
    double val, best;
    double startzeit, kzeit;
    int edgeNum, *edgeList = (int *) NULL;
    int *incycle = (int *) NULL, *outcycle = (int *) NULL;
    int* *near = (int* *) NULL;
    CCdatagroup dat;
    int in_repeater;
    const int kNum = 50;
    double minCost;
    int nearCity;
    const int eNum = 10;
    int eCnt, listIdx;
    int iter;

    unsigned int seed = (unsigned int) time (0);
    srand(seed);

    printf ("Chained Lin-Kernighan with seed %d\n", seed);
    fflush (stdout);

    if (ac != 2) {
        printf ("please input tsp file\n");
        return 1;
    }
    nodefile = av[1];

    startzeit = CCutil_zeit ();

    if (CCutil_gettsplib (nodefile, &ncount, &dat)) {
        fprintf (stderr, "could not read the TSPLIB file\n");
        rval = 1;
        goto CLEANUP;
    }

    in_repeater = ncount;

    incycle = CC_SAFE_MALLOC (ncount, int);
    if (!incycle) {
        rval = 1;
        goto CLEANUP;
    }

    /* find k-nearest */
    kzeit = CCutil_zeit ();
    near = CC_SAFE_MALLOC (ncount, int*);
    for (int ci = 0; ci < ncount; ++ci) {
        near[ci] = CC_SAFE_MALLOC (kNum, int);
    }
    for (int ci = 0; ci < ncount; ++ci) {
        for (int j = 0; j < ncount; ++j) {
            incycle[j] = 0;
        }
        incycle[ci] = 1;
        for (int ni = 0; ni < kNum; ++ni) {
            minCost = BIGDOUBLE;
            nearCity = -1;
            for (int cj = 0; cj < ncount; ++cj) {
                int cost = CCutil_dat_edgelen(ci, cj, &dat);
                if (cost <= minCost && incycle[cj] == 0) {
                    nearCity = cj;
                    minCost = cost;
                }
            }
            assert(nearCity != -1);
            near[ci][ni] = nearCity;
            incycle[nearCity] = 1;
        }
    }
    printf ("Time to find k-nearest: %.2f\n", CCutil_zeit () - kzeit);
    fflush (stdout);

    /* generate edge by k-nearest */
    kzeit = CCutil_zeit ();
    edgeNum = 0;
    for (int ci = 0; ci < ncount; ++ci) {
        eCnt = 0;
        for (int ni = 0; ni < kNum && eCnt < eNum; ++ni) {
            if (near[ci][ni] > ci) {
                eCnt++;
                edgeNum++;
            }
        }
    }
    assert(edgeNum > 0);
    edgeList = CC_SAFE_MALLOC (2*edgeNum, int);
    listIdx = 0;
    for (int ci = 0; ci < ncount; ++ci) {
        eCnt = 0;
        for (int ni = 0; ni < kNum && eCnt < eNum; ++ni) {
            if (near[ci][ni] > ci) {
                edgeList[listIdx] = ci;
                edgeList[listIdx+1] = near[ci][ni];
                eCnt++;
                listIdx += 2;
            }
        }
    }
    assert(listIdx == 2*edgeNum);
    printf ("Time to generate edge: %.2f\n", CCutil_zeit () - kzeit);
    fflush (stdout);

    CCutil_randcycle (ncount, incycle);

    outcycle = CC_SAFE_MALLOC (ncount, int);
    if (!outcycle) {
        rval = 1;
        goto CLEANUP;
    }

    iter = 0;
    best = BIGDOUBLE;
    do {
        printf ("\nStarting Run %d\n", iter);
        if (CClinkern_tour (ncount, &dat, edgeNum, edgeList, 100000000,
                in_repeater, incycle, outcycle, &val, run_silently)) {
            fprintf (stderr, "CClinkern_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
        if (val < best) {
            best = val;
        }
    } while (++iter < number_runs);
    printf ("Overall Best Cycle: %.0f\n", val);
    printf ("Total Running Time: %.2f\n", CCutil_zeit () - startzeit);
    fflush (stdout);

CLEANUP:

    for (int ci = 0; ci < ncount; ++ci) {
        CC_IFFREE(near[ci], int);
    }
    CC_IFFREE(near, int*);

    CC_IFFREE (edgeList, int);
    CC_IFFREE (incycle, int);
    CC_IFFREE (outcycle, int);
    CCutil_freedatagroup (&dat);

    return rval;
}
