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

#ifndef __UTIL_H
#define __UTIL_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <alloc.h>


#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

double CCutil_zeit (void);

void CCutil_randcycle (int ncount, int *cyc);

typedef struct CCdatagroup {
    int    (*edgelen) (int i, int j, struct CCdatagroup *dat);
    double  *x;
    double  *y;
    int      norm;
} CCdatagroup;

int
    CCutil_dat_edgelen (int i, int j, CCdatagroup *dat),
    CCutil_gettsplib (char *datname, int *ncount, CCdatagroup *dat);

void
    CCutil_init_datagroup (CCdatagroup *dat),
    CCutil_freedatagroup (CCdatagroup *dat);

#define CC_KD_NORM_TYPE    128          /* Kdtrees work      */
#define CC_X_NORM_TYPE     256          /* Old nearest works */
#define CC_NORM_BITS      (CC_KD_NORM_TYPE | CC_X_NORM_TYPE)
#define CC_EUCLIDEAN_CEIL (0 |   CC_KD_NORM_TYPE)
#define CC_EUCLIDEAN      (1 |   CC_KD_NORM_TYPE)
#define CC_ATT            (2 |    CC_X_NORM_TYPE)

#endif /* __UTIL_H */
