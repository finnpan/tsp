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
/*  int CCkdtree_build (CCkdtree *kt, int ncount, CCdatagroup *dat,         */
/*      double *wcoord)                                */
/*    -When called, intree should point to a CCkdtree struct that the       */
/*     funtion will load with the tree it builds. The wcoord array          */
/*     is used for node weights (like in Held-Karp), it can be NULL.        */
/*     The node weights must be nonegative (for cutoffs).                   */
/*                                                                          */
/*  void CCkdtree_free (CCkdtree *kt)                                       */
/*    -Frees the space (including the ptrs) used by kt.                     */
/*                                                                          */
/*    NOTES:                                                                */
/*       On a 32 bit machine, a CCkdtree on n nodes needs about 52n         */
/*     bytes of memory. CCkdtree_build will return 1 if an error            */
/*     occurs (most likely running out of memory).                          */
/*                                                                          */
/****************************************************************************/

#include "kdtree.h"

#define CUTOFF 5
#define BNDS_DEPTH 5   /* When bnds info is recorded */
#define BIGDOUBLE (1e30)


static void
    kdtree_free_work (CCkdnode *p, CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world),
    kdtree_free_world (CCptrworld *kdnode_world, CCptrworld *kdbnds_world);
static unsigned char
    findmaxspread (int l, int u, CCkdtree *thetree, double *datx,
           double *daty, double *datw);
static CCkdnode
    *build (int l, int u, int *depth, double *current_bnds_x,
           double *current_bnds_y, CCkdtree *thetree, double *datx,
           double *daty, double *datw);

static void
    CCutil_rselect (int *arr, int l, int r, int m, double *coord);


CC_PTRWORLD_ROUTINES (CCkdnode, kdnodealloc, kdnode_bulk_alloc, kdnodefree)
CC_PTRWORLD_LEAKS_ROUTINE (CCkdnode, kdnode_check_leaks, empty, char)

CC_PTRWORLD_ROUTINES (CCkdbnds, kdbndsalloc, kdbnds_bulk_alloc, kdbndsfree)
CC_PTRWORLD_LEAKS_ROUTINE (CCkdbnds, kdbnds_check_leaks, x[0], double)

int CCkdtree_build (CCkdtree *intree, int ncount, CCdatagroup *dat,
        double *wcoord)
{
    int i;
    int depth;
    double current_bnds_x[2];
    double current_bnds_y[2];
    CCkdtree *thetree;

    CCptrworld_init (&intree->kdnode_world);
    CCptrworld_init (&intree->kdbnds_world);

    if (wcoord != (double *) NULL) {
        for (i = 0; i < ncount; i++) {
            if (wcoord[i] < -0.00000001) {
                fprintf (stderr, "Cannot build with negative node weights\n");
                return 1;
            }
        }
    }

    thetree = intree;
    thetree->perm = CC_SAFE_MALLOC (ncount, int);
    if (!thetree->perm)
        return 1;
    for (i = 0; i < ncount; i++)
        thetree->perm[i] = i;

    thetree->bucketptr = CC_SAFE_MALLOC (ncount, CCkdnode *);
    if (!thetree->bucketptr) {
        CC_FREE (thetree->perm, int);
        return 1;
    }

    depth = 0;
    current_bnds_x[0] = -BIGDOUBLE;
    current_bnds_x[1] =  BIGDOUBLE;
    current_bnds_y[0] = -BIGDOUBLE;
    current_bnds_y[1] =  BIGDOUBLE;

    thetree->root = build (0, ncount - 1, &depth, current_bnds_x,
                     current_bnds_y, thetree, dat->x, dat->y, wcoord);
    if (!(thetree->root)) {
        fprintf (stderr, "Unable to build CCkdtree\n");
        CC_FREE (thetree->perm, int);
        CC_FREE (thetree->bucketptr, CCkdnode *);
        return 1;
    } else {
        thetree->root->father = (CCkdnode *) NULL;
        return 0;
    }
}

void CCkdtree_free (CCkdtree *kt)
{
    if (kt->perm)
        CC_FREE (kt->perm, int);
    if (kt->bucketptr)
        CC_FREE (kt->bucketptr, CCkdnode *);
    kdtree_free_work (kt->root, &kt->kdnode_world, &kt->kdbnds_world);
    kt->root = (CCkdnode *) NULL;

    kdtree_free_world (&kt->kdnode_world, &kt->kdbnds_world);
}

static void kdtree_free_world (CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world)
{
    int total, onlist;

    if (kdnode_check_leaks (kdnode_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding kdnodes\n",
                 total - onlist);
    }
    if (kdbnds_check_leaks (kdbnds_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding kdbnds\n", total - onlist);
    }
    CCptrworld_delete (kdnode_world);
    CCptrworld_delete (kdbnds_world);
}

static void kdtree_free_work (CCkdnode *p, CCptrworld *kdnode_world,
        CCptrworld *kdbnds_world)
{
    if (p->bucket) {
        if (p->bnds)
            kdbndsfree (kdbnds_world, p->bnds);
        kdnodefree (kdnode_world, p);
    } else {
        kdtree_free_work (p->loson, kdnode_world, kdbnds_world);
        kdtree_free_work (p->hison, kdnode_world, kdbnds_world);
        if (p->bnds)
            kdbndsfree (kdbnds_world, p->bnds);
        kdnodefree (kdnode_world, p);
    }
}

static CCkdnode *build (int l, int u, int *depth, double *current_bnds_x,
        double *current_bnds_y, CCkdtree *thetree, double *datx, double *daty,
        double *datw)
{
    CCkdnode *p;
    int i, m;
    double savebnd;

    (*depth)++;
    p = kdnodealloc (&thetree->kdnode_world);
    if (!p) {
        (*depth)--;
        return (CCkdnode *) NULL;
    }
    p->empty = 0;

    if (u - l + 1 < CUTOFF) {
        p->bucket = 1;
        p->lopt = l;
        p->hipt = u;
        for (i = l; i <= u; i++)
            thetree->bucketptr[thetree->perm[i]] = p;
        p->bnds = (CCkdbnds *) NULL;
    } else {
        p->bucket = 0;
        if (!((*depth) % BNDS_DEPTH)) {
            p->bnds = kdbndsalloc (&thetree->kdbnds_world);
            if (!p->bnds) {
                (*depth)--;
                kdnodefree (&thetree->kdbnds_world, p);
                return (CCkdnode *) NULL;
            }
            p->bnds->x[0] = current_bnds_x[0];
            p->bnds->x[1] = current_bnds_x[1];
            p->bnds->y[0] = current_bnds_y[0];
            p->bnds->y[1] = current_bnds_y[1];
        } else {
            p->bnds = (CCkdbnds *) NULL;
        }

        p->cutdim = findmaxspread (l, u, thetree, datx, daty, datw);
        m = (l + u) / 2;
        switch (p->cutdim) {
        case 0:
            CCutil_rselect (thetree->perm, l, u, m, datx);
            p->cutval = datx[thetree->perm[m]];

            savebnd = current_bnds_x[1];
            current_bnds_x[1] = p->cutval;
            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_x[1] = savebnd;

            savebnd = current_bnds_x[0];
            current_bnds_x[0] = p->cutval;
            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_x[0] = savebnd;

            break;
        case 1:
            CCutil_rselect (thetree->perm, l, u, m, daty);
            p->cutval = daty[thetree->perm[m]];

            savebnd = current_bnds_y[1];
            current_bnds_y[1] = p->cutval;
            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_y[1] = savebnd;

            savebnd = current_bnds_y[0];
            current_bnds_y[0] = p->cutval;
            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }
            current_bnds_y[0] = savebnd;

            break;
        case 2:
            CCutil_rselect (thetree->perm, l, u, m, datw);
            p->cutval = datw[thetree->perm[m]];

            p->loson = build (l, m, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw);
            if (!p->loson) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }

            p->hison = build (m + 1, u, depth, current_bnds_x, current_bnds_y,
                              thetree, datx, daty, datw);
            if (!p->hison) {
                (*depth)--;
                kdnodefree (&thetree->kdnode_world, p);
                return (CCkdnode *) NULL;
            }

            break;
        }
        p->loson->father = p;
        p->hison->father = p;
    }
    (*depth)--;
    return p;
}

static unsigned char findmaxspread (int l, int u, CCkdtree *thetree,
        double *datx, double *daty, double *datw)
{
    int i;
    double xmax, xmin, xval, xspread;
    double ymax, ymin, yval, yspread;
    double wmax, wmin, wval, wspread;

    wmax = (double) 0.0;
    wmin = (double) 0.0;

    if (datw != (double *) NULL) {
        wmin = datw[thetree->perm[l]];
        wmax = wmin;
    }
    xmin = datx[thetree->perm[l]];
    xmax = xmin;
    ymin = daty[thetree->perm[l]];
    ymax = ymin;
    for (i = l + 1; i <= u; i++) {
        xval = datx[thetree->perm[i]];
        if (xval < xmin)
            xmin = xval;
        else if (xval > xmax)
            xmax = xval;
        yval = daty[thetree->perm[i]];
        if (yval < ymin)
            ymin = yval;
        else if (yval > ymax)
            ymax = yval;
        if (datw != (double *) NULL) {
            wval = datw[thetree->perm[i]];
            if (wval < wmin)
                wmin = wval;
            else if (wval > wmax)
                wmax = wval;
        }
    }

    xspread = xmax - xmin;
    yspread = ymax - ymin;

    if (datw != (double *) NULL) {
        wspread = (wmax - wmin);
        if (xspread >= yspread && xspread >= wspread)
            return (unsigned char) 0;
        else if (yspread >= xspread && yspread >= wspread)
            return (unsigned char) 1;
        else {
            return (unsigned char) 2;
        }
    } else {
        if (xspread >= yspread)
            return (unsigned char) 0;
        else
            return (unsigned char) 1;
    }
}




/****************************************************************************/
/*  int CCkdtree_quadrant_k_nearest (CCkdtree *kt, int ncount, int k,       */
/*      CCdatagroup *dat, double *wcoord,                                   */
/*      int wantlist, int *ocount, int **olist, int silent)                 */
/*    RETURNS the quadrant k-nearest neighbor graph.                        */
/*                                                                          */
/*  int CCkdtree_node_k_nearest (CCkdtree *kt, int ncount, int n, int k,    */
/*      CCdatagroup *dat, double *wcoord, int *list)                        */
/*    RETURNS the k nearest points to point n.                              */
/*      -The k points are return in list (and list must be allocated by     */
/*       calling routine.                                                   */
/*      -kt is a pointer to a CCkdtree previously built by                  */
/*       CCkdtree_build.                                                    */
/*                                                                          */
/*  int CCkdtree_node_quadrant_k_nearest (CCkdtree *kt, int ncount,         */
/*      int n, int k, CCdatagroup *dat, double *wcoord, int *list,)         */
/*    RETURNS the quadrant k nearest point to point n.                      */
/*      -see CCkdtree_node_k_nearest.                                       */
/*                                                                          */
/*    NOTES:                                                                */
/*       If memory is tight, use CCkdtree_node_k_nearest to get the         */
/*    edges one node at a time. (CCkdtree_k_nearest () builds a hash        */
/*    table to avoid duplicate edges, and it will use 8 * nedges            */
/*    bytes.); all other routines return 0 if successful and 1 otherwise.   */
/*                                                                          */
/****************************************************************************/


#define BIGDOUBLE (1e30)

typedef struct shortedge {
    int end;
    double length;
} shortedge;

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

#define Fedgelen(n1, n2)                                                     \
    (datw != (double *) NULL ?                                               \
      CCutil_dat_edgelen ((n1), (n2), dat)                                   \
            + datw[(n1)] + datw[(n2)] :                                      \
      CCutil_dat_edgelen ((n1), (n2), dat))

#define dtrunc(x) (((x)>0.0)?floor(x):ceil(x))

static void
    node_k_nearest_work (CCkdtree *thetree, CCdatagroup *dat, double *datw,
        CCkdnode *p, int target, int num, shortedge *nearlist,
        double *worst_on_list, CCkdbnds *box);
static int
    run_kdtree_k_nearest (CCkdtree *kt, int ncount, int k, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ocount, int **olist, int doquad,
        int silent),
    put_in_table (int i, int j, int *added, intptr **table,
        CCptrworld *intptr_world),
    q_run_it (CCkdtree *thetree, CCdatagroup *dat, double *datw, int *llist,
         int *lcount, int *list, int target, int num, CCkdbnds *box),
    run_kdtree_node_k_nearest (CCkdtree *thetree, CCdatagroup *dat,
         double *datw, int *list, int target, int num, CCkdbnds *box),
    ball_in_bounds (CCdatagroup *dat, CCkdbnds *bnds, int n, double dist);


CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulk_alloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

int CCkdtree_quadrant_k_nearest (CCkdtree *kt, int ncount, int k,
        CCdatagroup *dat, double *wcoord, int wantlist, int *ocount,
        int **olist, int silent)
{
    return run_kdtree_k_nearest (kt, ncount, k, dat, wcoord,
                                 wantlist, ocount, olist, 1, silent);
}


static int run_kdtree_k_nearest (CCkdtree *kt, int ncount, int k,
        CCdatagroup *dat, double *wcoord, int wantlist, int *ocount,
        int **olist, int doquad, int silent)
{
    int i, n;
    intptr *ip, *ipnext;
    int total, onlist;
    CCkdtree localkt, *mykt;
    int added, ntotal = 0;
    int rval = 0;
    int *list = (int *) NULL;
    int goal = (doquad ? (4 * k) : k);
    int newtree = 0;
    intptr **table = (intptr **) NULL;
    CCptrworld intptr_world;

    CCptrworld_init (&intptr_world);

    if (wcoord != (double *) NULL) {
        for (i = 0; i < ncount; i++) {
            if (wcoord[i] < -0.00000001) {
                fprintf (stderr, "Cannot CCkdtree with negative node weights\n");
                return 1;
            }
        }
    }

    if (wantlist) {
        *ocount = 0;
        *olist = (int *) NULL;
    }

    if (kt == (CCkdtree *) NULL) {
        if (CCkdtree_build (&localkt, ncount, dat, wcoord)) {
            fprintf (stderr, "Unable to build CCkdtree\n");
            return 1;
        }
        mykt = &localkt;
        newtree = 1;
    } else {
        mykt = kt;
    }


    table = CC_SAFE_MALLOC (ncount, intptr *);
    if (!table) {
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < ncount; i++)
        table[i] = (intptr *) NULL;
    list = CC_SAFE_MALLOC (goal, int);
    if (!list) {
        rval = 1;
        goto CLEANUP;
    }

    for (n = 0; n < ncount; n++) {
        if (doquad) {
            if (CCkdtree_node_quadrant_k_nearest (mykt, ncount, n, k, dat,
                       wcoord, list)) {
                rval = 1;
                goto CLEANUP;
            }
        } else {
            if (CCkdtree_node_k_nearest (mykt, ncount, n, k, dat, wcoord,
                                         list)) {
                rval = 1;
                goto CLEANUP;
            }
        }
        for (i = 0; i < goal; i++) {
            if (list[i] != -1) {
                if (put_in_table (n, list[i], &added, table, &intptr_world))  {
                    rval = 1;
                    goto CLEANUP;
                } else {
                    ntotal += added;
                }
            }
        }
/*
        if (n == 0) {
            printf ("Neighbors of Node %d (%d, %d) :\n", n,
                                      (int) dat->x[n], (int) dat->y[n]);
            for (i = 0; i < goal; i++) {
                if (list[i] != -1) {
                    printf ("%d  %d (%d, %d)\n", list[i],
                      CCutil_dat_edgelen (n, list[i], dat),
                      (int) dat->x[list[i]], (int) dat->y[list[i]]);
                }
            }
        }
*/
        if (!silent) {
            if (n % 1000 == 999) {
                printf ("."); fflush (stdout);
            }
        }
    }

    if (!silent) {
        printf (" %d edges\n", ntotal); fflush (stdout);
    }

    if (wantlist) {
        int j = 0;
        *olist = CC_SAFE_MALLOC (2 * ntotal, int);
        if (!(*olist)) {
            rval = 1;
            goto CLEANUP;
        }
        *ocount = ntotal;
        for (i = 0; i < ncount; i++) {
            for (ip = table[i]; ip; ip = ipnext) {
                ipnext =  ip->next;
                (*olist)[j++] = i;
                (*olist)[j++] = ip->this;
                intptrfree (&intptr_world, ip);
            }
            table[i] = (intptr *) NULL;
        }
    } else {
        for (i = 0; i < ncount; i++) {
            intptr_listfree (&intptr_world, table[i]);
            table[i] = (intptr *) NULL;
        }
    }
    if (intptr_check_leaks (&intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs in kdnear\n",
                 total - onlist);
    }

CLEANUP:

    CCptrworld_delete (&intptr_world);
    if (table)
        CC_FREE(table, intptr *);
    if (list)
        CC_FREE (list, int);
    if (newtree)
        CCkdtree_free (&localkt);

    return rval;
}

static int put_in_table (int i, int j, int *added, intptr **table,
        CCptrworld *intptr_world)
{
    intptr *ip;

    if (j < i) {
        int temp;
        CC_SWAP(i, j, temp);
    }

    for (ip = table[i]; ip; ip = ip->next) {
        if (ip->this == j) {
            *added = 0;
            return 0;
        }
    }
    if (intptr_listadd (&table[i], j, intptr_world)) {
        *added = 0;
        return 1;
    }
    *added = 1;
    return 0;
}

int CCkdtree_node_quadrant_k_nearest (CCkdtree *kt, int ncount, int n, int k,
            CCdatagroup *dat, double *wcoord, int *list)
{
    CCkdbnds localbnds;
    int i, lcount = 0;
    int *llist = (int *) NULL;
    int rval = 0;
    CCkdtree localkt;
    CCkdtree *thetree;
    int newtree = 0;

    if (kt == (CCkdtree *) NULL) {
        if (CCkdtree_build (&localkt, ncount, dat, wcoord)) {
            fprintf (stderr, "Unable to build CCkdtree\n");
            return 1;
        }
        thetree = &localkt;
        newtree = 1;
    } else {
        thetree = kt;
    }

    llist = CC_SAFE_MALLOC (k, int);
    if (!llist) {
        rval = 1;
        goto CLEANUP;
    }

    localbnds.x[0] = dat->x[n];
    localbnds.x[1] = BIGDOUBLE;
    localbnds.y[0] = dat->y[n];
    localbnds.y[1] = BIGDOUBLE;
    if (q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                  &localbnds)) {
        fprintf (stderr, "run_kdtree_node_k_nearest failed\n");
        rval = 1;
        goto CLEANUP;
    }

    localbnds.x[0] = dat->x[n];
    localbnds.x[1] = BIGDOUBLE;
    localbnds.y[0] = -BIGDOUBLE;
    localbnds.y[1] = dat->y[n];
    if (q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                  &localbnds)) {
        fprintf (stderr, "run_kdtree_node_k_nearest failed\n");
        rval = 1;
        goto CLEANUP;
    }

    localbnds.x[0] = -BIGDOUBLE;
    localbnds.x[1] = dat->x[n];
    localbnds.y[0] = -BIGDOUBLE;
    localbnds.y[1] = dat->y[n];
    if (q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                  &localbnds)) {
        fprintf (stderr, "run_kdtree_node_k_nearest failed\n");
        rval = 1;
        goto CLEANUP;
    }

    localbnds.x[0] = -BIGDOUBLE;
    localbnds.x[1] = dat->x[n];
    localbnds.y[0] = dat->y[n];
    localbnds.y[1] = BIGDOUBLE;
    if (q_run_it (thetree, dat, wcoord, llist, &lcount, list, n, k,
                  &localbnds)) {
        fprintf (stderr, "run_kdtree_node_k_nearest failed\n");
        rval = 1;
        goto CLEANUP;
    }

    for (i = lcount; i < (4 * k); i++)
        list[i] = -1;

CLEANUP:

    CC_FREE (llist, int);
    if (newtree)
        CCkdtree_free (&localkt);

    return rval;
}

static int q_run_it (CCkdtree *thetree, CCdatagroup *dat, double *datw,
        int *llist, int *lcount, int *list, int target, int num, CCkdbnds *box)
{
    int i, j;

    if (run_kdtree_node_k_nearest (thetree, dat, datw, llist, target, num,
                                   box))
        return 1;
    for (i = 0; i < num; i++) {
        if (llist[i] != -1) {
            for (j = 0; j < *lcount; j++)
                if (list[j] == llist[i])
                    break;
            if (j == *lcount)
                list[(*lcount)++] = llist[i];
        }
    }
    return 0;
}

int CCkdtree_node_k_nearest (CCkdtree *kt, int ncount, int n, int k,
        CCdatagroup *dat, double *wcoord, int *list)
{
    CCkdtree localkt;
    CCkdtree *thetree;
    int newtree = 0;
    int rval = 0;

    if (kt == (CCkdtree *) NULL) {
        if (CCkdtree_build (&localkt, ncount, dat, wcoord)) {
            fprintf (stderr, "Unable to build CCkdtree\n");
            return 1;
        }
        thetree = &localkt;
        newtree = 1;
    } else {
        thetree = kt;
    }

    rval = run_kdtree_node_k_nearest (thetree, dat, wcoord, list, n, k,
                                      (CCkdbnds *) NULL);
    if (newtree)
        CCkdtree_free (&localkt);
    return rval;
}

static int run_kdtree_node_k_nearest (CCkdtree *thetree, CCdatagroup *dat,
        double *datw, int *list, int target, int num, CCkdbnds *box)
{
    int i;
    CCkdnode *p, *lastp;
    double diff;
    shortedge *nearlist = (shortedge *) NULL;
    double worst_on_list = BIGDOUBLE;

    nearlist = CC_SAFE_MALLOC (num + 1, shortedge);
    if (!nearlist) {
        return 1;
    }
    for (i = 0; i < num; i++)
        nearlist[i].length = BIGDOUBLE;
    nearlist[num].length = -BIGDOUBLE;

/*
    To do top down search just use:

        node_k_nearest_work (thetree->root);
*/

    p = thetree->bucketptr[target];
    node_k_nearest_work (thetree, dat, datw, p, target, num, nearlist,
                         &worst_on_list, box);
    while (1) {
        lastp = p;
        p = p->father;
        if (p == (CCkdnode *) NULL)
            break;
        switch (p->cutdim) {
        case 0:
            diff = p->cutval - dat->x[target];
            if (lastp == p->loson) {    /* So target is on low side */
               if (worst_on_list > dtrunc(diff))
                   if (box == (CCkdbnds *) NULL || p->cutval <= box->x[1])
                       node_k_nearest_work (thetree, dat, datw, p->hison,
                              target, num, nearlist, &worst_on_list, box);
            } else {
               if (worst_on_list > dtrunc(-diff))
                   if (box == (CCkdbnds *) NULL || p->cutval >= box->x[0])
                       node_k_nearest_work (thetree, dat, datw, p->loson,
                              target, num, nearlist, &worst_on_list, box);
            }
            break;
        case 1:
            diff = p->cutval - dat->y[target];
            if (lastp == p->loson) {
               if (worst_on_list > dtrunc(diff))
                   if (box == (CCkdbnds *) NULL || p->cutval <= box->y[1])
                       node_k_nearest_work (thetree, dat, datw, p->hison,
                              target, num, nearlist, &worst_on_list, box);
            } else {
               if (worst_on_list > dtrunc(-diff))
                   if (box == (CCkdbnds *) NULL || p->cutval >= box->y[0])
                       node_k_nearest_work (thetree, dat, datw, p->loson,
                              target, num, nearlist, &worst_on_list, box);
            }
            break;
        case 2:
            if (lastp == p->loson) {
                if (worst_on_list > p->cutval + datw[target])
                    node_k_nearest_work (thetree, dat, datw, p->hison,
                              target, num, nearlist, &worst_on_list, box);
            } else {
                node_k_nearest_work (thetree, dat, datw, p->loson,
                              target, num, nearlist, &worst_on_list, box);
            }
            break;
        }
        if (datw == (double *) NULL && p->bnds &&
               ball_in_bounds (dat, p->bnds, target, worst_on_list))
              /* Doing extra check for box with quad-nearest appears to slow */
              /* things down.                                                */
            break;
    }

    int ntot = 0;
    for (i = 0; i < num; i++) {
        if (nearlist[i].length < BIGDOUBLE)
            list[ntot++] = nearlist[i].end;
    }
    if (ntot < num) {
        if (box == (CCkdbnds *) NULL) {
            fprintf (stderr, "WARNING: There do not exist %d neighbors\n",
                        num);
        }
        for (i = ntot; i < num; i++)
            list[i] = -1;
    }


    CC_FREE (nearlist, shortedge);
    return 0;
}

static void node_k_nearest_work (CCkdtree *thetree, CCdatagroup *dat,
        double *datw, CCkdnode *p, int target, int num, shortedge *nearlist,
        double *worst_on_list, CCkdbnds *box)
{
    int i, k;
    double val, thisx, thisdist;

    if (p->bucket) {
        for (i = p->lopt; i <= p->hipt; i++) {
            if (thetree->perm[i] != target) {
                if (box == (CCkdbnds *) NULL ||
                    (dat->x[thetree->perm[i]] >= box->x[0] &&
                    dat->x[thetree->perm[i]] <= box->x[1] &&
                    dat->y[thetree->perm[i]] >= box->y[0] &&
                    dat->y[thetree->perm[i]] <= box->y[1])) {
                    thisdist = Fedgelen (thetree->perm[i], target);
                    if (*worst_on_list > thisdist) {
                        for (k = 0; nearlist[k+1].length > thisdist; k++) {
                            nearlist[k].end = nearlist[k + 1].end;
                            nearlist[k].length = nearlist[k + 1].length;
                        }
                        nearlist[k].length = thisdist;
                        nearlist[k].end = thetree->perm[i];
                        *worst_on_list = nearlist[0].length;
                    }
                }
            }
        }
    } else {
        val = p->cutval;
        switch (p->cutdim) {
        case 0:
            thisx = dat->x[target];
            if (thisx < val) {
                node_k_nearest_work (thetree, dat, datw, p->loson, target, num,
                        nearlist, worst_on_list, box);
                /* Truncation for floating point coords */
                if (*worst_on_list > dtrunc(val - thisx))
                    if (box == (CCkdbnds *) NULL || val >= box->x[0])
                        node_k_nearest_work (thetree, dat, datw, p->hison,
                               target, num, nearlist, worst_on_list, box);
            } else {
                node_k_nearest_work (thetree, dat, datw, p->hison,
                               target, num, nearlist, worst_on_list, box);
                if (*worst_on_list > dtrunc(thisx - val))
                    if (box == (CCkdbnds *) NULL || val <= box->x[1])
                        node_k_nearest_work (thetree, dat, datw, p->loson,
                               target, num, nearlist, worst_on_list, box);
            }
            break;
        case 1:
            thisx = dat->y[target];
            if (thisx < val) {
                node_k_nearest_work (thetree, dat, datw, p->loson, target,
                               num, nearlist, worst_on_list, box);
                if (*worst_on_list > dtrunc(val - thisx))
                    if (box == (CCkdbnds *) NULL || val >= box->y[0])
                        node_k_nearest_work (thetree, dat, datw, p->hison,
                               target, num, nearlist, worst_on_list, box);
            } else {
                node_k_nearest_work (thetree, dat, datw, p->hison, target,
                               num, nearlist, worst_on_list, box);
                if (*worst_on_list > dtrunc(thisx - val))
                    if (box == (CCkdbnds *) NULL || val <= box->y[1])
                        node_k_nearest_work (thetree, dat, datw, p->loson,
                               target, num, nearlist, worst_on_list, box);
            }
            break;
        case 2:
            thisx = datw[target];
            node_k_nearest_work (thetree, dat, datw, p->loson, target,
                               num, nearlist, worst_on_list, box);
            if (*worst_on_list > val + thisx)
                node_k_nearest_work (thetree, dat, datw, p->hison, target,
                               num, nearlist, worst_on_list, box);
            break;
        }
    }
}

static int ball_in_bounds (CCdatagroup *dat, CCkdbnds *bnds, int n,
        double dist)
{
    if (dtrunc(dat->x[n] - bnds->x[0]) < dist ||
        dtrunc(bnds->x[1] - dat->x[n]) < dist ||
        dtrunc(dat->y[n] - bnds->y[0]) < dist ||
        dtrunc(bnds->y[1] - dat->y[n]) < dist)
        return 0;
    return 1;
}



/****************************************************************************/
/*  void CCutil_rselect (int *arr, int l, int r, int m,                     */
/*      double *coord)                                                      */
/*    arr - permutation that will be rearranged                             */
/*    l,r - specify the range of arr that we are interested in              */
/*    m - is the index into l,r that is the break point for the perm        */
/*    coord - gives the keys that determine the ordering                    */
/*                                                                          */
/****************************************************************************/

static void
    select_split (int *arr, int n, double v, int *start, int *end,
           double *coord),
    select_sort (int *arr, int n, double *coord),
    select_sort_dsample (double *samp, int n);

/* NSAMPLES should be odd */
#define NSAMPLES 3
#define SORTSIZE 20

void CCutil_rselect (int *arr, int l, int r, int m, double *coord)
{
    double samplevals[NSAMPLES];
    int i;
    int st, en;
    int n;

    arr += l;
    n = r - l + 1;
    m -= l;

    while (n > SORTSIZE) {
        for (i = 0; i < NSAMPLES; i++) {
            samplevals[i] = coord[arr[rand() % n]];
        }
        select_sort_dsample (samplevals, NSAMPLES);
        select_split (arr, n, samplevals[(NSAMPLES - 1) / 2], &st, &en, coord);
        if (st > m) {
            n = st;
        } else if (en <= m) {
            arr += en;
            n -= en;
            m -= en;
        } else {
            return;
        }
    }

    select_sort (arr, n, coord);
    return;
}

static void select_split (int *arr, int n, double v, int *start, int *end,
                          double *coord)
{
    int i, j, k;
    int t;

    i = 0;
    j = k = n;

    while (i < j) {
        if (coord[arr[i]] < v) {
            i++;
        } else if (coord[arr[i]] == v) {
            j--;
            CC_SWAP (arr[i], arr[j], t);
        } else {
            j--;
            k--;
            t = arr[i];
            arr[i] = arr[j];
            arr[j] = arr[k];
            arr[k] = t;
        }
    }
    *start = j;
    *end = k;
    return;
}

static void select_sort (int *arr, int n, double *coord)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++) {
        t = arr[i];
        for (j = i; j > 0 && coord[arr[j - 1]] > coord[t]; j--) {
            arr[j] = arr[j - 1];
        }
        arr[j] = t;
    }
}

static void select_sort_dsample (double *samp, int n)
{
    int i, j;
    double t;

    for (i = 1; i < n; i++) {
        t = samp[i];
        for (j = i; j > 0 && samp[j - 1] > t; j--) {
            samp[j] = samp[j - 1];
        }
        samp[j] = t;
    }
}
