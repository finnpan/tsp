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

#include "util.h"


double CCutil_zeit (void)
{
    return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}

void CCutil_randcycle (int ncount, int *cyc)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++) cyc[i] = i;
    for (i = ncount; i > 1; i--) {
        k = rand() % i;
        CC_SWAP (cyc[i - 1], cyc[k], temp);
    }
}



static int
    edgelen_nonorm (int i, int j, CCdatagroup *dat),
    euclid_edgelen (int i, int j, CCdatagroup *dat),
    euclid_ceiling_edgelen (int i, int j, CCdatagroup *dat),
    att_edgelen (int i, int j, CCdatagroup *dat),
    dat_setnorm (CCdatagroup *dat, int norm);

int CCutil_dat_edgelen (int i, int j, CCdatagroup *dat)
{
    return (dat->edgelen)(i, j, dat);
}

int dat_setnorm (CCdatagroup *dat, int norm)
{
    switch (norm) {
    case CC_EUCLIDEAN_CEIL:
        dat->edgelen = euclid_ceiling_edgelen;
        break;
    case CC_EUCLIDEAN:
        dat->edgelen = euclid_edgelen;
        break;
    case CC_ATT:
        dat->edgelen = att_edgelen;
        break;
    default:
        fprintf (stderr, "ERROR:  Unknown NORM %d.\n", norm);
        return 1;
    }
    dat->norm = norm;

    return 0;
}

int CCutil_gettsplib(char *datname, int *ncount, CCdatagroup *dat)
{
    char buf[256], key[256], field[256];
    char *p;
    FILE *in;
    int norm = -1;

    CCutil_init_datagroup (dat);
    *ncount = -1;

    if ((in = fopen (datname, "r")) == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        return 1;
    }

    while (fgets (buf, 254, in) != (char *) NULL) {
        p = buf;
        while (*p != '\0') {
            if (*p == ':')
                *p = ' ';
            p++;
        }
        p = buf;
        if (sscanf (p, "%s", key) != EOF) {
            p += strlen (key);
            while (*p == ' ')
                p++;
            if (!strcmp (key, "NAME")) {
                printf ("Problem Name: %s", p);
            } else if (!strcmp (key, "TYPE")) {
                printf ("Problem Type: %s", p);
                if (sscanf (p, "%s", field) == EOF || strcmp (field, "TSP")) {
                    fprintf (stderr, "Not a TSP problem\n");
                    return 1;
                }
            } else if (!strcmp (key, "COMMENT")) {
                printf ("%s", p);
            } else if (!strcmp (key, "DIMENSION")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in DIMENSION line\n");
                    return 1;
                }
                *ncount = atoi (field);
                printf ("Number of Nodes: %d\n", *ncount);
            } else if (!strcmp (key, "EDGE_WEIGHT_TYPE")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in EDGE_WEIGHT_TYPE line\n");
                    return 1;
                }
                if (!strcmp (field, "EUC_2D")) {
                    norm = CC_EUCLIDEAN;
                    printf ("Rounded Euclidean Norm (CC_EUCLIDEAN)\n");
                } else if (!strcmp (field, "ATT")) {
                    norm = CC_ATT;
                    printf ("ATT Norm (CC_ATT)\n");
                } else if (!strcmp (field, "CEIL_2D")) {
                    norm = CC_EUCLIDEAN_CEIL;
                    printf ("Rounded Up Euclidean Norm (CC_EUCLIDEAN_CEIL)\n");
                } else {
                    fprintf (stderr, "ERROR: Not set up for norm %s\n", field);
                    return 1;
                }
                if (dat_setnorm (dat, norm)) {
                    fprintf (stderr, "ERROR: Couldn't set norm %d\n", norm);
                    return 1;
                }
            } else if (!strcmp (key, "NODE_COORD_SECTION")) {
                int i;
                if (*ncount <= 0) {
                    fprintf (stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (dat->x != (double *) NULL) {
                    fprintf (stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    CCutil_freedatagroup (dat);
                    return 1;
                }
                if (1) {
                    dat->x = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->x) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    dat->y = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->y) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0; i < *ncount; i++) {
                        fscanf (in, "%*d %lf %lf", &(dat->x[i]), &(dat->y[i]));
                    }
                } else {
                    fprintf (stderr, "ERROR: Node coordinates with norm %d?\n",
                                 norm);
                    return 1;
                }
            }
        }
    }
    fclose (in);

    if (dat->x == (double *) NULL) {
        fprintf (stderr, "ERROR: Didn't find the data\n");
        return 1;
    } else {
        return 0;
    }
}

static int edgelen_nonorm (int i, int j, CCdatagroup *dat)
{
    fprintf (stderr, "CCutil_dat_edgelen has been called with no norm set\n");
    fprintf (stderr, "This is a FATAL ERROR\n");
    if (i != 0 || j != 0 || dat != (CCdatagroup *) NULL) {
        /* so the compiler won't complain about unused variables */
        fprintf (stderr, "This is a FATAL ERROR\n");
        exit (1);
    }
    return -1;
}

static int euclid_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    int temp;

    temp = (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
    return temp;
}

static int euclid_ceiling_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    return (int) (ceil (sqrt (t1 * t1 + t2 * t2)));
}

static int att_edgelen (int i, int j, CCdatagroup *dat)
{
    double xd = dat->x[i] - dat->x[j];
    double yd = dat->y[i] - dat->y[j];
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);

    int k = (int) rij;
    double tij = (double) k;

    int dij;
    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}

void CCutil_init_datagroup (CCdatagroup *dat)
{
    dat->x = (double *) NULL;
    dat->y = (double *) NULL;
    dat->norm = 0;
    dat->edgelen = edgelen_nonorm;
}

void CCutil_freedatagroup (CCdatagroup *dat)
{
    CC_IFFREE (dat->x, double);
    CC_IFFREE (dat->y, double);
}



/****************************************************************************/
/*  void *CCutil_allocrus (size_t size)                                     */
/*    RETURNS a pointer to an allocated block of "size" memory.             */
/*                                                                          */
/*  void CCutil_freerus (void *ptr)                                         */
/*    FREES ptr.                                                            */
/*                                                                          */
/*  void *CCutil_reallocrus (void *ptr, size_t size)                        */
/*    REALLOCS ptr to size bytes.                                           */
/*                                                                          */
/*  int CCutil_reallocrus_scale (void **pptr, int *pnnum, int count,        */
/*      double scale, size_t size)                                          */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int *pnnum (a reference to the number of objects in the               */
/*                allocated space)                                          */
/*    int count (a minimum value for the new nnum)                          */
/*    double scale (a scale factor to apply to nnum)                        */
/*    int size (the size of objects to be realloced)                        */
/*    RETURNS 0 if *pptr was successfully changed to point to at            */
/*            least max(*pnnum*scale, *pnnum+1000, count) objects.          */
/*            *pnnum is changed to the new object count.                    */
/*            Otherwise, prints an error message, leaves *pptr and          */
/*            *pnnum alone, and returns nonzero.                            */
/*                                                                          */
/*  int CCutil_reallocrus_count (void **pptr, int count,                    */
/*      size_t size)                                                        */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int count (number of objects to be realloced)                         */
/*    int size (the size of the objects to be realloced)                    */
/*    RETURNS 0 is successful, and 1 if the realloc failed.                 */
/*                                                                          */
/*  CCbigchunkptr *CCutil_bigchunkalloc (void)                              */
/*         RETURNS a CCbigchunkptr with the "this_one" field loaded with a  */
/*                 a pointer to a bigchunk of memory.                       */
/*    NOTES:                                                                */
/*       The idea is to use bigchunks (the size of a bigchunk is defined    */
/*       by CC_BIGCHUNK in util.h) to supply local routines with memory     */
/*       for ptrs, so the memory can be shared with other                   */
/*       local routines.                                                    */
/*                                                                          */
/*  CCutil_bigchunkfree (CCbigchunkptr *bp)                                 */
/*    ACTION: Frees a CCbigchunkptr.                                        */
/*                                                                          */
/*  void CCptrworld_init (CCptrworld *world)                                */
/*     initialize a CCptrworld with 1 reference                             */
/*                                                                          */
/*  void CCptrworld_add (CCptrworld *world)                                 */
/*     add a reference to a CCptrworld                                      */
/*                                                                          */
/*  void CCptrworld_delete (CCptrworld *world)                              */
/*     delete a reference to a ptrworld, and free if no more references     */
/*                                                                          */
/****************************************************************************/

typedef struct CCbigchunk {
    char space[CC_BIGCHUNK];
    CCbigchunkptr ptr;
} CCbigchunk;

void *CCutil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf (stderr, "Warning: 0 bytes allocated\n");
    }

    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void CCutil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }

    free (p);
}

void *CCutil_reallocrus (void *ptr, size_t size)
{
    void *newptr;

    if (!ptr) {
        return CCutil_allocrus (size);
    } else {
        newptr = (void *) realloc (ptr, size);
        if (!newptr) {
            fprintf (stderr, "Out of memory.  Tried to grow to %d bytes\n",
                     (int) size);
        }
        return newptr;
    }
}

int CCutil_reallocrus_scale (void **pptr, int *pnnum, int count, double scale,
        size_t size)
{
    int newsize = (int) (((double) *pnnum) * scale);
    void *p;

    if (newsize < *pnnum+1000) newsize = *pnnum+1000;
    if (newsize < count) newsize = count;
    p = CCutil_reallocrus (*pptr, newsize * size);
    if (!p) {
        return 1;
    } else {
        *pptr = p;
        *pnnum = newsize;
        return 0;
    }
}

int CCutil_reallocrus_count (void **pptr, int count, size_t size)
{
    void *p = CCutil_reallocrus (*pptr, count * size);

    if (!p) {
        return 1;
    } else {
        *pptr = p;
        return 0;
    }
}

CCbigchunkptr *CCutil_bigchunkalloc (void)
{
    CCbigchunk *p = CC_SAFE_MALLOC (1, CCbigchunk);

    if (p == (CCbigchunk *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_bigchunkalloc\n");
        return (CCbigchunkptr *) NULL;
    }
    p->ptr.this_chunk = p;
    p->ptr.this_one = (void *) p->space;
    return &(p->ptr);
}

void CCutil_bigchunkfree (CCbigchunkptr *bp)
{
    /* This copy is necessary since CC_FREE zeros its first argument */
    CCbigchunk *p = bp->this_chunk;

    CC_FREE (p, CCbigchunk);
}

void CCptrworld_init (CCptrworld *world)
{
    world->refcount = 1;
    world->freelist = (void *) NULL;
    world->chunklist = (CCbigchunkptr *) NULL;
}

void CCptrworld_add (CCptrworld *world)
{
    world->refcount++;
}

void CCptrworld_delete (CCptrworld *world)
{
    world->refcount--;
    if (world->refcount <= 0) {
        CCbigchunkptr *bp, *bpnext;

        for (bp = world->chunklist ; bp; bp = bpnext) {
            bpnext = bp->next;
            CCutil_bigchunkfree (bp);
        }
        world->chunklist = (CCbigchunkptr *) NULL;
        world->freelist = (void *) NULL;
        world->refcount = 0;
    }
}
