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

#ifndef __KDTREE_H
#define __KDTREE_H

#include "mpool.h"
#include "util.h"

typedef struct CCkdnode {
    double cutval;
    struct CCkdnode *loson;
    struct CCkdnode *hison;
    struct CCkdnode *father;
    struct CCkdnode *next;
    struct CCkdbnds *bnds;
    int              lopt;
    int              hipt;
    char             bucket;
    char             empty;
    char             cutdim;
} CCkdnode;

typedef struct CCkdtree {
    CCkdnode        *root;
    CCkdnode       **bucketptr;
    int             *perm;
    mpool           *pool;
} CCkdtree;

typedef struct CCkdbnds {
    double           x[2];
    double           y[2];
    struct CCkdbnds *next;
} CCkdbnds;


void
    CCkdtree_free (CCkdtree *kt);

int
    CCkdtree_build (CCkdtree *kt, int ncount, CCdatagroup *dat, double *wcoord),
    CCkdtree_quadrant_k_nearest (CCkdtree *kt, int ncount, int k,
        CCdatagroup *dat, double *wcoord, int wantlist, int *ocount,
        int **olist, int silent),
    CCkdtree_node_k_nearest (CCkdtree *kt, int ncount, int n, int k,
        CCdatagroup *dat, double *wcoord, int *list),
    CCkdtree_node_quadrant_k_nearest (CCkdtree *kt, int ncount, int n, int k,
        CCdatagroup *dat, double *wcoord, int *list);

#endif  /* __KDTREE_H */
