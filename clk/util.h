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

#define CC_EUCLIDEAN_CEIL (0)
#define CC_EUCLIDEAN      (1)
#define CC_ATT            (2)

#endif /* __UTIL_H */
