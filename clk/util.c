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


/****************************************************************************/
/*                                                                          */
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 2, 1995 (cofeb16)                                        */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
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




/****************************************************************************/
/*                                                                          */
/*                         SORTING ROUTINES                                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*   Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/*   DATE:  February 24, 1994                                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCutil_rselect (int *arr, int l, int r, int m,                     */
/*      double *coord, CCrandstate *rstate)                                 */
/*    arr - permutation that will be rearranged                             */
/*    l,r - specify the range of arr that we are interested in              */
/*    m - is the index into l,r that is the break point for the perm        */
/*    coord - gives the keys that determine the ordering                    */
/*                                                                          */
/****************************************************************************/

#define BITS_PER_PASS (8)


static void
    select_split (int *arr, int n, double v, int *start, int *end,
           double *coord),
    select_sort (int *arr, int n, double *coord),
    select_sort_dsample (double *samp, int n);


/**********  Median - Select Routines **********/

/* NSAMPLES should be odd */
#define NSAMPLES 3
#define SORTSIZE 20


void CCutil_rselect (int *arr, int l, int r, int m, double *coord,
        CCrandstate *rstate)
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
            samplevals[i] = coord[arr[CCutil_lprand (rstate) % n]];
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



/****************************************************************************/
/*                                                                          */
/*              MACHINE INDEPENDENT RANDOM NUMBER GENERATOR                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  DIMACS  (modified for TSP)                                 */
/*  Date: February 7, 1995  (cofeb16)                                       */
/*        September 18, 2001  (billenium fix)                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCutil_sprand (int seed, CCrandstate *r)                           */
/*    - Call once to initialize the generator.                              */
/*                                                                          */
/*  int CCutil_lprand (CCrandstate *r)                                      */
/*    - Returns an integer in the range 0 to CC_PRANDMAX - 1.               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*    NOTES (from DIMACS):                                                  */
/*        This file contains a set of c-language functions for generating   */
/*    uniform integers.   This is a COMPLETELY PORTABLE generator. It will  */
/*    give IDENTICAL sequences of random numbers for any architecture with  */
/*    at least 30-bit integers, regardless of the integer representation,   */
/*    INT_MAX value, or roundoff/truncation method, etc.                    */
/*        This Truly Remarkable RNG is described more fully in              */
/*    J. Bentley's column, ``The Software Exploratorium ''. It is based on  */
/*    one in Knuth, Vol 2, Section 3.2.2 (Algorithm A).                     */
/*                                                                          */
/****************************************************************************/



void CCutil_sprand (int seed, CCrandstate *r)
{
    int i, ii;
    int last, next;
    int *arr = r->arr;

    seed %= CC_PRANDMAX;
    if (seed < 0) seed += CC_PRANDMAX;

    arr[0] = last = seed;
    next = 1;
    for (i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        arr[ii] = next;
        next = last - next;
        if (next < 0)
            next += CC_PRANDMAX;
        last = arr[ii];
    }
    r->a = 0;
    r->b = 24;
    for (i = 0; i < 165; i++)
        last = CCutil_lprand (r);
}


int CCutil_lprand (CCrandstate *r)
{
    int t;

    if (r->a-- == 0)
        r->a = 54;
    if (r->b-- == 0)
        r->b = 54;

    t = r->arr[r->a] - r->arr[r->b];

    if (t < 0)
        t += CC_PRANDMAX;

    r->arr[r->a] = t;

    return t;
}


double CCutil_zeit (void)
{
    return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}




/****************************************************************************/
/*                                                                          */
/*                        PORTABLE GETOPT                                   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: 1993 (?) (fmfeb02)                                                */
/*  Modified: 15 February 1995 (bico)  - added warning                      */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_bix_getopt (int argc, char **argv, const char *def,          */
/*      int *p_optind, char **p_optarg)                                     */
/*     parse an argument list                                               */
/*                                                                          */
/****************************************************************************/


int CCutil_bix_getopt (int ac, char **av, const char *def, int *p_optind,
        char **p_optarg)
{
    int c;
    char *sp = av[*p_optind];
    char bwarn[2];

    if (*p_optind < 1 || *p_optind >= ac) {
        *p_optind = ac;
        return (EOF);
    }
    if ((int) *sp != (int) '-')
        return (EOF);
    if ((int) *(sp + 1) == (int) '-') {
        (*p_optind)++;
        return (EOF);
    }
    (av[*p_optind])++;
    sp++;
    while ((int) *sp != (int) *def && (int) *def != (int) '\0')
            def++;
    if ((int) *def == (int) '\0') {
        *p_optind = ac;
        bwarn[0] = *sp;                          /* Bico: February 8, 1995 */
        bwarn[1] = '\0';
        printf ("Illegal option: -%s\n", bwarn);
        return CC_BIX_GETOPT_UNKNOWN;
    }
    if ((int) *(def + 1) != (int) ':') {
        c = *sp;
        if ((int) *(sp + 1) != (int) '\0')
            *sp = '-';
        else
            (*p_optind)++;
        return (c);
    } else {
        if ((int) *(sp + 1) != (int) '\0') {
            *p_optarg = sp + 1;
            c = *sp;
            (*p_optind)++;
            return (c);
        } else if (*p_optind >= ac - 1) {
            *p_optind = ac;
            return (EOF);
        } else {
            *p_optarg = av[*p_optind + 1];
            c = *sp;
            *p_optind += 2;
            return (c);
        }
    }
}




/****************************************************************************/
/*                                                                          */
/*                   SOME DATA READING ROUTINES                             */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 2, 1995                                                     */
/*  Changes: 17.7.96 (Bico)                                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_gettsplib (char *datname, int *ncount, CCdatagroup *dat)     */
/*    READS an xxx.tsp TSPLIB file, and returns the dat structure to        */
/*            generate edge lengths.                                        */
/*     -datname should be the name of a TSPLIB xxx.tsp file.                */
/*     -ncount returns the number of nodes.                                 */
/*     -dat returns the data.                                               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*     Functions prototyped in util.h. Functions return 0 when they         */
/*    succeed and nonzero when they fail (usually do to bad filenames or    */
/*    not enough memory).                                                   */
/*     The TSPLIB reader works for all problems in TSPLIB_1.2, but does     */
/*    not include all of the options listed in Reinelt's orginal TSPLIB     */
/*    paper. It returns a failure on linhp318.tsp, since there is no        */
/*    place for fixed edges in our edge length dat structure.               */
/*                                                                          */
/****************************************************************************/


#define MATRIX_LOWER_DIAG_ROW  0
#define MATRIX_UPPER_ROW       1
#define MATRIX_UPPER_DIAG_ROW  2
#define MATRIX_FULL_MATRIX     3

int CCutil_gettsplib(char *datname, int *ncount, CCdatagroup *dat)
{
    char buf[256], key[256], field[256];
    char *p;
    FILE *in;
    int matrixform = MATRIX_LOWER_DIAG_ROW;
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
                if (!strcmp (field, "EXPLICIT")) {
                    norm = CC_MATRIXNORM;
                    printf ("Explicit Lengths (CC_MATRIXNORM)\n");
                } else if (!strcmp (field, "EUC_2D")) {
                    norm = CC_EUCLIDEAN;
                    printf ("Rounded Euclidean Norm (CC_EUCLIDEAN)\n");
                } else if (!strcmp (field, "EUC_3D")) {
                    norm = CC_EUCLIDEAN_3D;
                    printf ("Rounded Euclidean 3D Norm (CC_EUCLIDEAN_3D)\n");
                } else if (!strcmp (field, "MAX_2D")) {
                    norm = CC_MAXNORM;
                    printf ("Max Norm (CC_MAXNORM)\n");
                } else if (!strcmp (field, "MAN_2D")) {
                    norm = CC_MANNORM;
                    printf ("Max Norm (CC_MAXNORM)\n");
                } else if (!strcmp (field, "GEO")) {
                    norm = CC_GEOGRAPHIC;
                    printf ("Geographical Norm (CC_GEOGRAPHIC)\n");
                } else if (!strcmp (field, "GEOM")) {
                    norm = CC_GEOM;
                    printf ("Geographical Norm in Meters (CC_GEOM)\n");
                } else if (!strcmp (field, "ATT")) {
                    norm = CC_ATT;
                    printf ("ATT Norm (CC_ATT)\n");
                } else if (!strcmp (field, "CEIL_2D")) {
                    norm = CC_EUCLIDEAN_CEIL;
                    printf ("Rounded Up Euclidean Norm (CC_EUCLIDEAN_CEIL)\n");
                } else if (!strcmp (field, "DSJRAND")) {
                    norm = CC_DSJRANDNORM;
                    printf ("David Johnson Random Norm (CC_DSJRANDNORM)\n");
                } else {
                    fprintf (stderr, "ERROR: Not set up for norm %s\n", field);
                    return 1;
                }
                if (CCutil_dat_setnorm (dat, norm)) {
                    fprintf (stderr, "ERROR: Couldn't set norm %d\n", norm);
                    return 1;
                }
            } else if (!strcmp (key, "EDGE_WEIGHT_FORMAT")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in EDGE_WEIGHT_FORMAT line\n");
                    return 1;
                }
                if (!strcmp (field, "LOWER_DIAG_ROW")) {
                    matrixform = MATRIX_LOWER_DIAG_ROW;
                } else if (!strcmp (field, "UPPER_ROW")) {
                    matrixform = MATRIX_UPPER_ROW;
                } else if (!strcmp (field, "UPPER_DIAG_ROW")) {
                    matrixform = MATRIX_UPPER_DIAG_ROW;
                } else if (!strcmp (field, "FULL_MATRIX")) {
                    matrixform = MATRIX_FULL_MATRIX;
                } else if (strcmp (field, "FUNCTION")) {
                    fprintf (stderr, "Cannot handle format: %s\n", field);
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
                if ((norm & CC_NORM_SIZE_BITS) == CC_D2_NORM_SIZE) {
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
                } else if ((norm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
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
                    dat->z = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->z) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0; i < *ncount; i++) {
                        fscanf (in, "%*d %lf %lf %lf",
                               &(dat->x[i]), &(dat->y[i]), &(dat->z[i]));
                    }
                } else {
                    fprintf (stderr, "ERROR: Node coordinates with norm %d?\n",
                                 norm);
                    return 1;
                }
            } else if (!strcmp (key, "EDGE_WEIGHT_SECTION")) {
                int i, j;
                if (*ncount <= 0) {
                    fprintf (stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (dat->adj != (int **) NULL) {
                    fprintf (stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    CCutil_freedatagroup (dat);
                    return 1;
                }
                if ((norm & CC_NORM_SIZE_BITS) == CC_MATRIX_NORM_SIZE) {
                    dat->adj = CC_SAFE_MALLOC (*ncount, int *);
                    dat->adjspace = CC_SAFE_MALLOC ((*ncount)*(*ncount+1)/2,
                                                    int);
                    if (dat->adj == (int **) NULL ||
                        dat->adjspace == (int *) NULL) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0, j = 0; i < *ncount; i++) {
                        dat->adj[i] = dat->adjspace + j;
                        j += (i+1);
                    }
                    if (matrixform == MATRIX_LOWER_DIAG_ROW) {
                        for (i = 0; i < *ncount; i++) {
                            for (j = 0; j <= i; j++)
                                fscanf (in, "%d", &(dat->adj[i][j]));
                        }
                    } else if (matrixform == MATRIX_UPPER_ROW ||
                               matrixform == MATRIX_UPPER_DIAG_ROW ||
                               matrixform == MATRIX_FULL_MATRIX) {
                        int **tempadj = (int **) NULL;
                        int *tempadjspace = (int *) NULL;
                        tempadj = CC_SAFE_MALLOC (*ncount, int *);
                        tempadjspace = CC_SAFE_MALLOC ((*ncount) * (*ncount),
                                                       int);
                        if (tempadj == (int **) NULL ||
                            tempadjspace == (int *) NULL) {
                            CC_IFFREE (tempadj, int *);
                            CC_IFFREE (tempadjspace, int);
                            CCutil_freedatagroup (dat);
                            return 1;
                        }
                        for (i = 0; i < *ncount; i++) {
                            tempadj[i] = tempadjspace + i * (*ncount);
                            if (matrixform == MATRIX_UPPER_ROW) {
                                tempadj[i][i] = 0;
                                for (j = i + 1; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            } else if (matrixform == MATRIX_UPPER_DIAG_ROW) {
                                for (j = i; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            } else {
                                for (j = 0; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            }
                        }
                        for (i = 0; i < *ncount; i++) {
                            for (j = 0; j <= i; j++)
                                dat->adj[i][j] = tempadj[j][i];
                        }
                        CC_FREE (tempadjspace, int);
                        CC_FREE (tempadj, int *);
                    }
                } else {
                    fprintf (stderr, "ERROR: Matrix with norm %d?\n",
                             norm);
                    return 1;
                }
            } else if (!strcmp (key, "FIXED_EDGES_SECTION")) {
                fprintf (stderr, "ERROR: Not set up for fixed edges\n");
                return 1;
            }
        }
    }
    fclose (in);

    if (dat->x == (double *) NULL && dat->adj == (int **) NULL) {
        fprintf (stderr, "ERROR: Didn't find the data\n");
        return 1;
    } else {
        return 0;
    }
}




/****************************************************************************/
/*                                                                          */
/*                           DHEAP ROUTINES                                 */
/*                                                                          */
/*                                                                          */
/*                              TSP CODE                                    */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 9, 1995                                                  */
/*  Reference: R.E. Tarjan, Data Structures and Network Algorithms          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_dheap_init (CCdheap *h, int k)                               */
/*        -h should point to a CCdheap struct.                              */
/*        -k the max number of elements in the dheap.                       */
/*                                                                          */
/*  void CCutil_dheap_free (CCdheap *h)                                     */
/*    -frees the spaces allocated by CCutil_dheap_init                      */
/*                                                                          */
/*  int CCutil_dheap_resize (CCdheap *h, int newsize)                       */
/*    -REALLOCs h so it can contain newsize elements.                       */
/*    -returns -1 if it can't resize the heap.                              */
/*                                                                          */
/*  int CCutil_dheap_findmin (CCdheap *h)                                   */
/*    -returns the index of the element with min value h->key[i]            */
/*    -returns -1 if no elements in heap.                                   */
/*                                                                          */
/*  int CCutil_dheap_insert (CCdheap *h, int i)                             */
/*    -inserts the element with index i (so its key should be loaded        */
/*     beforehand in h->key[i]).                                            */
/*                                                                          */
/*  void CCutil_dheap_delete (CCdheap *h, int i)                            */
/*    -deletes the element with index i.                                    */
/*                                                                          */
/*  int CCutil_dheap_deletemin (CCdheap *h)                                 */
/*    -returns the min element in the heap, and deletes the min element     */
/*    -returns -1 if no elements in heap.                                   */
/*                                                                          */
/*  void CCutil_dheap_changekey (CCdheap *h, int i, double newkey)          */
/*    -changes the key of the element with index i to newkey.               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*      A k-element heap will malloc 16k bytes of memory. If memory is      */
/*  tight, using integer keys (instead of doubles), brings it down to       */
/*  12k bytes, and if arbitrary deletions are not required, with a little   */
/*  rewriting, the h->loc field can be eliminated, bring the space down     */
/*  to 8k bytes.                                                            */
/*      These routines work with indices into the h->key array, so in       */
/*  some cases, you will need to maintain a separate names array to know    */
/*  what element belongs to index i. For an example, see the k_nearest      */
/*  code in kdnear.c.                                                       */
/*                                                                          */
/****************************************************************************/

#define HEAP_D 3
#define HEAP_UP(x) (((x)-1)/HEAP_D)
#define HEAP_DOWN(x) (((x)*HEAP_D)+1)


static void
    dheap_siftup (CCdheap *h, int i, int x),
    dheap_siftdown (CCdheap *h, int i, int x);

static int
    dheap_minchild (int x, CCdheap *h);



int CCutil_dheap_init (CCdheap *h, int k)
{
    h->loc = (int *) NULL;
    h->key = (double *) NULL;
    h->entry = CC_SAFE_MALLOC (k, int);
    if (!h->entry)
        return 1;
    h->loc = CC_SAFE_MALLOC (k, int);
    if (!h->loc) {
        CC_FREE (h->entry, int);
        return 1;
    }
    h->key = CC_SAFE_MALLOC (k, double);
    if (!h->key) {
        CC_FREE (h->entry, int);
        CC_FREE (h->loc, int);
        return 1;
    }
    h->total_space = k;
    h->size = 0;
    return 0;
}

void CCutil_dheap_free (CCdheap *h)
{
    CC_IFFREE (h->entry, int);
    CC_IFFREE (h->loc, int);
    CC_IFFREE (h->key, double);
}

int CCutil_dheap_resize (CCdheap *h, int newsize)
{
    if (newsize < h->size || newsize < h->total_space) return 0;
    if (CCutil_reallocrus_count ((void **) &(h->key), newsize,
                                 sizeof (double))) {
        return -1;
    }
    if (CCutil_reallocrus_count ((void **) &(h->entry), newsize,
                                 sizeof (int))) {
        return -1;
    }
    if (CCutil_reallocrus_count ((void **) &(h->loc), newsize, sizeof (int))) {
        return -1;
    }
    h->total_space = newsize;

    return 0;
}

int CCutil_dheap_findmin (CCdheap *h)
{
    if (h->size == 0)
        return -1;
    else
        return h->entry[0];
}

int CCutil_dheap_insert (CCdheap *h, int i)
{
    if (h->size >= h->total_space) {
        fprintf (stderr, "Error - heap already full\n");
        return 1;
    }
    h->size++;
    dheap_siftup (h, i, h->size - 1);
    return 0;
}

void CCutil_dheap_delete (CCdheap *h, int i)
{
    int j;

    h->size--;
    j = h->entry[h->size];
    h->entry[h->size] = -1;

    if (j != i) {
        if (h->key[j] <= h->key[i]) {
            dheap_siftup (h, j, h->loc[i]);
        } else {
            dheap_siftdown (h, j, h->loc[i]);
        }
    }
}

int  CCutil_dheap_deletemin (CCdheap *h)
{
    int i;

    if (h->size == 0)
        return -1;
    else {
        i = h->entry[0];
        CCutil_dheap_delete (h, i);
        return i;
    }
}

void CCutil_dheap_changekey (CCdheap *h, int i, double newkey)
{
    if (newkey < h->key[i]) {
        h->key[i] = newkey;
        dheap_siftup (h, i, h->loc[i]);
    } else if (newkey > h->key[i]) {
        h->key[i] = newkey;
        dheap_siftdown (h, i, h->loc[i]);
    }
}

static void dheap_siftup (CCdheap *h, int i, int x)
{
    int p;

    p = HEAP_UP (x);
    while (x && h->key[h->entry[p]] > h->key[i]) {
        h->entry[x] = h->entry[p];
        h->loc[h->entry[p]] = x;
        x = p;
        p = HEAP_UP (p);
    }
    h->entry[x] = i;
    h->loc[i] = x;
}

static void dheap_siftdown (CCdheap *h, int i, int x)
{
    int c;

    c = dheap_minchild (x, h);

    while (c >= 0 && h->key[h->entry[c]] < h->key[i]) {
        h->entry[x] = h->entry[c];
        h->loc[h->entry[c]] = x;
        x = c;
        c = dheap_minchild (c, h);
    }
    h->entry[x] = i;
    h->loc[i] = x;
}

static int dheap_minchild (int x, CCdheap *h)
{
    int c = HEAP_DOWN (x);
    int cend;
    double minval;
    int minloc;

    if (c >= h->size)
        return -1;
    minval = h->key[h->entry[c]];
    minloc = c;
    cend = c + HEAP_D;
    if (h->size < cend)
        cend = h->size;
    for (c++; c < cend; c++) {
        if (h->key[h->entry[c]] < minval) {
            minval = h->key[h->entry[c]];
            minloc = c;
        }
    }
    return minloc;
}




/****************************************************************************/
/*                                                                          */
/*       FUNCTIONS FOR COMPUTING EDGE LENGTHS FOR GEOMETRIC PROBLEMS        */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: Summer 1994                                                       */
/*        Modified - March 2, 1995                                          */
/*                 - October 5, 1995 (Bico)                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_dat_setnorm (CCdatagroup *dat, int norm)                     */
/*     NOTES:                                                               */
/*         Supported norms (with defs in edgelen.h) are:                    */
/*             CC_MAXNORM  -  the L-infinity norm                           */
/*             CC_EUCLIDEAN_CEIL - the norm for the plaXXXX problems        */
/*             CC_EUCLIDEAN - rounded L-2 norm                              */
/*             CC_EUCLIDEAN_3D - rounded L-2 norm in 3 space                */
/*             CC_USER - a norm specified by the user                       */
/*             CC_GEOGRAPHIC - distances on a sphere (Groetshel and         */
/*                             Holland)                                     */
/*             CC_GEOM - sphere in meters, coords in decimal degrees,       */
/*                             slight modification of CC_GEOGRAPHIC         */
/*             CC_ATT - pseudo-Euclidean norm for att532                    */
/*             CC_MATRIXNORM - complete graph (lower + diagonal matrix)     */
/*             CC_DSJRANDNORM - random edgelengths                          */
/*             CC_CRYSTAL - Bland-Shallcross xray norm                      */
/*                     - The coordinates generated for CC_CRYSTAL problems  */
/*                (in CCutil_getdata.c) have been diveded by the motor      */
/*                speeds (this makes the edgelen function faster) and       */
/*                scaled by CRYSTAL_SCALE (currently 10000) and rounded to  */
/*                the nearest integer (this lets the edgelen function       */
/*                produce integer lengths without further rounding). The    */
/*                result is a closer approximation to the Bland -           */
/*                Shallcross floating point length function than that       */
/*                given in TSPLIB_1.2.                                      */
/*             CC_SPARSE - a sparse graph                                   */
/*             CC_RHMAPx - where x = 1, 2, 3, 4, 5 one of 5 RH mapping      */
/*                norms.                                                    */
/*                                                                          */
/*         If CCUTIL_EDGELEN_FUNCTIONPTR has been defined in util.h,        */
/*         then CCutil_dat_edgelen is a pointer to a function instead of    */
/*         a function.  This saves a function call and results in           */
/*         improved performance on some machines for edgelen-intensive      */
/*         routines like linkern.  The function pointer is set by           */
/*         CCutil_dat_setnorm.                                              */
/*                                                                          */
/*         IMPORTANT: This means that if CCUTIL_EDGELEN_FUNCTIONPTR is set  */
/*         and you have more than one CCdatagroup, you must call            */
/*         CCutil_dat_setnorm whenever you switch from using one            */
/*         CCdatagroup to the other.  IF YOU DON'T DO THIS, EDGELEN WILL    */
/*         RETURN INCORRECT RESULTS.  For this reason,                      */
/*         CCUTIL_EDGELEN_FUNCTIONPTR should only be set with extreme       */
/*         caution.                                                         */
/*                                                                          */
/*         NOTE: CCUTIL_EDGELEN_FUNCTIONPTR does not work with the          */
/*         subdivision code for parallel TSP.                               */
/*                                                                          */
/*    To define a user norm, you must perform the following steps:          */
/*    1.  In util.h, define the struct CCdata_user to contain the data      */
/*        necessary for the computation of edge lengths.                    */
/*    2.  In edgelen.c, write the init_userdat and free_userdat functions   */
/*        which initialize and free a CCdata_user structure.                */
/*    3.  In edgelen.c, write the user_edgelen function which               */
/*        computes the length of the edge for node i to node j, using the   */
/*        userdat field of the CCdatagroup argument (userdat is of type     */
/*        CCdata_user).                                                     */
/*    4.  In getdata.c, write the build_user, read_user_text,               */
/*        read_user_binary, readmaster_user, and writemaster_user           */
/*        routines.  read_user_text reads the data file which provides      */
/*        the data for computing the edge lengths.  build_user and          */
/*        read_user_binary are optional routines which build random         */
/*        datasets and read binary datafiles.  writemaster_user writes a    */
/*        binary version of that data to the master file, and               */
/*        readmaster_user reads that same data.  See the comments before    */
/*        those routines in getdata for more details on what they should    */
/*        do.                                                               */
/*    5.  In getdata.c, write permute_user, which permutes the data to      */
/*        reflect a permutation of the nodes.                               */
/*                                                                          */
/*  int CCutil_dat_edgelen (int i, int j, CCdatagroup *dat)                 */
/*     compute the length of an edge                                        */
/*                                                                          */
/*  void CCutil_dat_getnorm (CCdatagroup *dat, int *norm)                   */
/*     get the norm of a CCdatagroup                                        */
/*                                                                          */
/*  void CCutil_init_datagroup (CCdatagroup *dat)                           */
/*     initialize a CCdatagroup                                             */
/*                                                                          */
/*  void CCutil_freedatagroup (CCdatagroup *dat)                            */
/*     free a CCdatagroup                                                   */
/*                                                                          */
/****************************************************************************/

static double
    dtrunc (double);

static int
    edgelen_nonorm (int i, int j, CCdatagroup *dat),
    max_edgelen (int i, int j, CCdatagroup *dat),
    man_edgelen (int i, int j, CCdatagroup *dat),
    euclid_edgelen (int i, int j, CCdatagroup *dat),
    euclid_ceiling_edgelen (int i, int j, CCdatagroup *dat),
    euclid3d_edgelen (int i, int j, CCdatagroup *dat),
    geographic_edgelen (int i, int j, CCdatagroup *dat),
    geom_edgelen (int i, int j, CCdatagroup *dat),
    att_edgelen (int i, int j, CCdatagroup *dat),
    dsjrand_edgelen (int i, int j, CCdatagroup *dat),
    crystal_edgelen (int i, int j, CCdatagroup *dat),
    matrix_edgelen (int i, int j, CCdatagroup *dat),
    sparse_edgelen (int i, int j, CCdatagroup *dat),
    user_edgelen (int i, int j, CCdatagroup *dat),
    rhmap1_edgelen (int i, int j, CCdatagroup *dat),
    rhmap2_edgelen (int i, int j, CCdatagroup *dat),
    rhmap3_edgelen (int i, int j, CCdatagroup *dat),
    rhmap4_edgelen (int i, int j, CCdatagroup *dat),
    rhmap5_edgelen (int i, int j, CCdatagroup *dat),
    toroidal_edgelen (int i, int j, CCdatagroup *dat);

static void
    init_userdat (CCdata_user *userdat),
    free_userdat (CCdata_user *userdat),
    init_rhdata (CCdata_rhvector *rhdat),
    free_rhdata (CCdata_rhvector *rhdat);


#ifndef M_PI
#define M_PI 3.14159265358979323846264
#endif

#ifdef  CCUTIL_EDGELEN_FUNCTIONPTR

int (*CCutil_dat_edgelen) (int i, int j, CCdatagroup *dat) = edgelen_nonorm;

#else /* CCUTIL_EDGELEN_FUNCTIONPTR */

int CCutil_dat_edgelen (int i, int j, CCdatagroup *dat)
{
    if (dat->ndepot) {
        if (i >= dat->orig_ncount) {
            return dat->depotcost[j];
        } else if (j >= dat->orig_ncount) {
            return dat->depotcost[i];
        }
    }
    return (dat->edgelen)(i, j, dat);
}

#endif /* CCUTIL_EDGELEN_FUNCTIONPTR */


int CCutil_dat_setnorm (CCdatagroup *dat, int norm)
{
    switch (norm) {
    case CC_EUCLIDEAN_CEIL:
        dat->edgelen = euclid_ceiling_edgelen;
        break;
    case CC_EUCLIDEAN:
        dat->edgelen = euclid_edgelen;
        break;
    case CC_MAXNORM:
        dat->edgelen = max_edgelen;
        break;
    case CC_MANNORM:
        dat->edgelen = man_edgelen;
        break;
    case CC_EUCLIDEAN_3D:
        dat->edgelen = euclid3d_edgelen;
        break;
    case CC_USER:
        dat->edgelen = user_edgelen;
        break;
    case CC_GEOGRAPHIC:
        dat->edgelen = geographic_edgelen;
        break;
    case CC_GEOM:
        dat->edgelen = geom_edgelen;
        break;
    case CC_ATT:
        dat->edgelen = att_edgelen;
        break;
    case CC_MATRIXNORM:
        dat->edgelen = matrix_edgelen;
        break;
    case CC_DSJRANDNORM:
        dat->edgelen = dsjrand_edgelen;
        break;
    case CC_CRYSTAL:
        dat->edgelen = crystal_edgelen;
        break;
    case CC_SPARSE:
        dat->edgelen = sparse_edgelen;
        break;
    case CC_RHMAP1:
        dat->edgelen = rhmap1_edgelen;
        break;
    case CC_RHMAP2:
        dat->edgelen = rhmap2_edgelen;
        break;
    case CC_RHMAP3:
        dat->edgelen = rhmap3_edgelen;
        break;
    case CC_RHMAP4:
        dat->edgelen = rhmap4_edgelen;
        break;
    case CC_RHMAP5:
        dat->edgelen = rhmap5_edgelen;
        break;
    case CC_EUCTOROIDAL:
        dat->edgelen = toroidal_edgelen;
        break;
    default:
        fprintf (stderr, "ERROR:  Unknown NORM %d.\n", norm);
        return 1;
    }
    dat->norm = norm;

#ifdef CCUTIL_EDGELEN_FUNCTIONPTR
    CCutil_dat_edgelen = dat->edgelen;
#endif /* CCUTIL_EDGELEN_FUNCTIONPTR */

    return 0;
}

void CCutil_dat_getnorm (CCdatagroup *dat, int *norm)
{
    (*norm) = dat->norm;
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

/* Several variables that would normally be called y1 and y2 are called
   yy1 and yy2 to avoid conflict with the bessel functions */

static int max_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];

    if (t1 < 0)
        t1 *= -1;
    if (t2 < 0)
        t2 *= -1;
    t1 += 0.5;
    t2 += 0.5;

    return (int) (t1 < t2 ? t2 : t1);
}

static int man_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];

    if (t1 < 0)
        t1 *= -1;
    if (t2 < 0)
        t2 *= -1;

    return (int) (t1 + t2 + 0.5);
}


static int euclid_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    int temp;

    temp = (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
    return temp;
}

static int toroidal_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j];
    double t2 = dat->y[i] - dat->y[j];
    int temp;

    if (t1 < 0) t1 = -t1;
    if (t2 < 0) t2 = -t2;
    if (dat->gridsize - t1 < t1) t1 = dat->gridsize - t1;
    if (dat->gridsize - t2 < t2) t2 = dat->gridsize - t2;
    temp = (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
    return temp;
}

static int euclid3d_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
    double t3 = dat->z[i] - dat->z[j];
    int temp;

    temp = (int) (sqrt (t1 * t1 + t2 * t2 + t3 * t3) + 0.5);
    return temp;
}

static int euclid_ceiling_edgelen (int i, int j, CCdatagroup *dat)
{
    double t1 = dat->x[i] - dat->x[j], t2 = dat->y[i] - dat->y[j];
/*
    int rd;
    double max;

    max = sqrt (t1 * t1 + t2 * t2);
    rd = (int) max;
    return (((max - rd) > .000000001) ? rd + 1 : rd);
*/
    return (int) (ceil (sqrt (t1 * t1 + t2 * t2)));
}

#define GH_PI (3.141592)

static int geographic_edgelen (int i, int j, CCdatagroup *dat)
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;
    double x1 = dat->x[i], x2 = dat->x[j], yy1 = dat->y[i], yy2 = dat->y[j];

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (yy1);
    min = yy1 - deg;
    longi = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (yy2);
    min = yy2 - deg;
    longj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
                + 1.0);
    return dd;
}

static int geom_edgelen (int i, int j, CCdatagroup *dat)
{
    double lati, latj, longi, longj;
    double q1, q2, q3, q4, q5;

    lati = M_PI * dat->x[i] / 180.0;
    latj = M_PI * dat->x[j] / 180.0;

    longi = M_PI * dat->y[i] / 180.0;
    longj = M_PI * dat->y[j] / 180.0;

    q1 = cos (latj) * sin(longi - longj);
    q3 = sin((longi - longj)/2.0);
    q4 = cos((longi - longj)/2.0);
    q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
    q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
    return (int) (6378388.0 * atan2(sqrt(q1*q1 + q2*q2), q5) + 1.0);
}

#if 0
static int geom_edgelen (int i, int j, CCdatagroup *dat)
{
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;

    lati = M_PI * (dat->x[i] / 180.0);
    latj = M_PI * (dat->x[j] / 180.0);

    longi = M_PI * (dat->y[i] / 180.0);
    longj = M_PI * (dat->y[j] / 180.0);

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378388.0 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
                + 1.0);
    return dd;
}
#endif

static int att_edgelen (int i, int j, CCdatagroup *dat)
{
    double xd = dat->x[i] - dat->x[j];
    double yd = dat->y[i] - dat->y[j];
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}

static double dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

static int dsjrand_edgelen (int i, int j, CCdatagroup *dat)
{
    int di = (int) dat->x[i];
    int dj = (int) dat->x[j];
    int x, y, z;

    x = di&dj;
    y = di|dj;
    z = dat->dsjrand_param;

    x *= z;
    y *= x;
    z *= y;

    z ^= dat->dsjrand_param;

    x *= z;
    y *= x;
    z *= y;

    x = ((di+dj)^z)&0x7fffffff;
    return (int)(x * dat->dsjrand_factor);
}

#define CRYSTAL_SCALE 10000

#define CRYSTAL_FLIP_TOL ((180 * CRYSTAL_SCALE * 4) / 5)
#define CRYSTAL_NEEDS_FLIP(x) ((x) > (CRYSTAL_FLIP_TOL))
#define CRYSTAL_FLIP(x) ((2 * (CRYSTAL_FLIP_TOL)) - (x))

static int crystal_edgelen (int i, int j, CCdatagroup *dat)
{
    double w, w1;

    w = dat->x[i] - dat->x[j];
    if (w < 0)
        w = -w;
    w1 = dat->y[i] - dat->y[j];
    if (w1 < 0)
        w1 = -w1;
    if (CRYSTAL_NEEDS_FLIP (w1))
        w1 = CRYSTAL_FLIP (w1);
    if (w < w1)
        w = w1;
    w1 = dat->z[i] - dat->z[j];
    if (w1 < 0)
        w1 = -w1;
    if (w < w1)
        w = w1;

    return (int) w;
}

static int matrix_edgelen (int i, int j, CCdatagroup *dat)
{
    if (i > j)
        return (dat->adj[i])[j];
    else
        return (dat->adj[j])[i];
}

static int sparse_edgelen (int i, int j, CCdatagroup *dat)
{
    int *adj;
    int k, deg;

    if (i > j) {
        CC_SWAP (i, j, k);
    }
    adj = dat->adj[i];
    deg = dat->degree[i];

    for (k = 0; k < deg; k++) {
        if (adj[k] == j) {
            return dat->len[i][k];
        }
    }
    return dat->default_len;
}

void CCutil_init_datagroup (CCdatagroup *dat)
{
    dat->x = (double *) NULL;
    dat->y = (double *) NULL;
    dat->z = (double *) NULL;
    dat->adj = (int **) NULL;
    dat->adjspace = (int *) NULL;
    dat->len      = (int **) NULL;
    dat->lenspace = (int *) NULL;
    dat->degree   = (int *) NULL;
    dat->norm = 0;
    dat->dsjrand_param = 1;
    dat->dsjrand_factor = 1.0;
    dat->default_len = 100000;
    dat->sparse_ecount = 0;
    dat->edgelen = edgelen_nonorm;
    init_userdat (&dat->userdat);
    init_rhdata (&dat->rhdat);
    dat->ndepot = 0;
    dat->orig_ncount = 0;
    dat->depotcost = (int *) NULL;
    dat->orig_names = (int *) NULL;
}

void CCutil_freedatagroup (CCdatagroup *dat)
{
    CC_IFFREE (dat->x, double);
    CC_IFFREE (dat->y, double);
    CC_IFFREE (dat->z, double);
    CC_IFFREE (dat->adj, int *);
    CC_IFFREE (dat->adjspace, int);
    CC_IFFREE (dat->len, int *);
    CC_IFFREE (dat->lenspace, int);
    CC_IFFREE (dat->degree, int);
    free_userdat (&dat->userdat);
    free_rhdata (&dat->rhdat);
    CC_IFFREE (dat->depotcost, int);
    CC_IFFREE (dat->orig_names, int);
}

static void init_userdat (CCdata_user *userdat)
{
    userdat->x = (double *) NULL;
    userdat->y = (double *) NULL;
}

static void free_userdat (CCdata_user *userdat)
{
    CC_IFFREE (userdat->x, double);
    CC_IFFREE (userdat->y, double);
}

static int user_edgelen (int i, int j, CCdatagroup *dat)
{
    double dw = dat->userdat.x[i] - dat->userdat.x[j];
    double dw1 = dat->userdat.y[i] - dat->userdat.y[j];
    static const double ibm_xmult[7] = {1062.5,
        300.0,
        300.0,
        250.0,
        300.0,
        1000.0,
        154.6};
    static const double ibm_xadd[7] = {155.0 - 0.01 * 1062.5,
        197.5 - 0.05 * 300.0,
        212.5 - 0.10 * 300.0,
        227.5 - 0.15 * 250.0,
        240.5 - 0.20 * 300.0,
        255.0 - 0.25 * 1000.0,
        305.0 - 0.30 * 154.6};
    static const double ibm_ymult[7] = {1062.5,
        450.0,
        350.0,
        250.0,
        300.0,
        900.0,
        157.7};
    static const double ibm_yadd[7] = {150.0 - 0.01 * 1062.5,
        192.5 - 0.05 * 450.0,
        215.0 - 0.10 * 350.0,
        232.5 - 0.15 * 250.0,
        245.5 - 0.20 * 300.0,
        250.0 - 0.25 * 900.0,
        295.0 - 0.30 * 157.7};

    if (dw < 0.0)
        dw = -dw;
    dw /= 25400.0;
    if (dw <= 0.01) {
        dw *= 15500.0;
    } else if (dw >= 0.30) {
        dw = dw * 154.6 + (305.0 - 0.3 * 154.6);
    } else {
        dw = dw * ibm_xmult[(int) (dw / 0.05)] +
            ibm_xadd[(int) (dw / 0.05)];
    }
    if (dw1 < 0.0)
        dw1 = -dw1;
    dw1 /= 25400.0;
    if (dw1 <= 0.01) {
        dw1 *= 15000.0;
    } else if (dw1 >= 0.30) {
        dw1 = dw1 * 157.7 + (295.0 - 0.3 * 157.7);
    } else {
        dw1 = dw1 * ibm_ymult[(int) (dw1 / 0.05)] +
            ibm_yadd[(int) (dw1 / 0.05)];
    }
    if (dw < dw1)
        dw = dw1;
    return (int) dw;
}

static void init_rhdata (CCdata_rhvector *rhdat)
{
    rhdat->space = (char *) NULL;
    rhdat->vectors = (char **) NULL;
    rhdat->rhlength = 0;
    rhdat->dist_00 = 0;
    rhdat->dist_01 = 0;
    rhdat->dist_02 = 0;
    rhdat->dist_22 = 0;
    rhdat->p       = 0.0;
}

static void free_rhdata (CCdata_rhvector *rhdat)
{
    CC_IFFREE (rhdat->space, char);
    CC_IFFREE (rhdat->vectors, char *);
    rhdat->rhlength = 0;
}

static int rhmap1_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int rhlength = dat->rhdat.rhlength;
    char *v1 = vectors[i];
    char *v2 = vectors[j];
    int n;
    int sum = 0;
    int dist_00 = dat->rhdat.dist_00;
    int dist_01 = dat->rhdat.dist_01;
    int dist_02 = dat->rhdat.dist_02;
    int dist_12 = dat->rhdat.dist_12;
    int dist_22 = dat->rhdat.dist_22;

    if (v1 == (char *) NULL || v2 == (char *) NULL) return 0;

    for (n=0; n<rhlength; n++) {
        if (v1[n] == 2) {
            if (v2[n] == 0) sum += dist_02;
            else if (v2[n] == 1) sum += dist_12;
            else sum += dist_22;
        } else {
            if (v1[n] == v2[n]) sum += dist_00;
            else if (v2[n] != 2) sum += dist_01;
            else if (v1[n] == 0) sum += dist_02;
            else sum += dist_12;
        }
    }
    return sum;
}

static int rhmap2_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int rhlength = dat->rhdat.rhlength;
    char *v1 = vectors[i];
    char *v2 = vectors[j];
    int n;
    double sum = 0;
    double p = dat->rhdat.p;

    if (v1 == (char *) NULL || v2 == (char *) NULL) return 0;

    for (n=0; n<rhlength; n++) {
        if (v1[n] == 0) {
            if (v2[n] == 1) sum += 1;
            else if (v2[n] == 2) sum += p;
        } else if (v1[n] == 1) {
            if (v2[n] == 0) sum += 1;
            else if (v2[n] == 2) sum += (1-p);
        } else {
            if (v2[n] == 0) sum += p;
            else if (v2[n] == 1) sum += (1-p);
            else sum += 2*p*(1-p);
        }
    }
    return (int) (sum * 100.0 + 0.5);
}

#define MAX_DIST 1000

static int rhmap3_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int len = dat->rhdat.rhlength;
    char *first = vectors[i];
    char *second = vectors[j];
    int xindex;
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    int n = 0;
    double P = dat->rhdat.p;
    double Q = 1.0 - P;
    double trans;
    double temp;
    double theta;
    double term;

    if (first == (char *) NULL) {
        if (second == (char *) NULL) return 0;
        first = second;
        second = (char *) NULL;
    }

    if (second == (char *) NULL) {
        for (xindex = 0; xindex < len; xindex++) {
            if (first[xindex] == 1) a++;
            else if (first[xindex] == 0) b++;
        }
        trans = pow(sqrt(P), (double) a) * pow(sqrt(Q), (double) b);
    } else {
        for (xindex = 0; xindex < len; xindex++) {
            if ((first[xindex] != 2) || (second[xindex] != 2)) {
                n++;
                if ((first[xindex] == 1) && (second[xindex] == 1)) a++;
                if ((first[xindex] == 1) && (second[xindex] == 0)) b++;
                if ((first[xindex] == 0) && (second[xindex] == 1)) c++;
                if ((first[xindex] == 0) && (second[xindex] == 0)) d++;
            }
        }
        if (n == 0) return MAX_DIST;
        temp = (n - (a*P) - (d*Q));
        term = (4.0 * n * P * Q * (b+c));

        if (term >= (temp * temp)) return MAX_DIST;

        theta = (temp - sqrt((temp * temp) - term)) / (2.0 * n * P * Q);

        if (theta >= 1.0) return MAX_DIST;

        trans = pow((1.0 - (theta*P)),  (double) d) *
                pow((1.0 - (theta*Q)),  (double) a) *
                pow((theta * sqrt(P*Q)), (double)  (b+c));

    }
    return ((int) (-10.0 * log10(trans)));
}

static int rhmap4_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int len = dat->rhdat.rhlength;
    char *first = vectors[i];
    char *second = vectors[j];
    int xindex;
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    int n = 0;
    double P = dat->rhdat.p;
    double Q = 1.0 - P;
    double trans;
    double temp;
    double theta;
    double term;

    if (first == (char *) NULL) {
        if (second == (char *) NULL) return 0;
        first = second;
        second = (char *) NULL;
    }

    if (second == (char *) NULL) {
        for (xindex = 0; xindex < len; xindex++) {
            if (first[xindex] == 1) a++;
            else if (first[xindex] == 0) b++;
        }
        trans = pow(sqrt(P), (double) a) * pow(sqrt(Q), (double) b);
    } else {
        for (xindex = 0; xindex < len; xindex++) {
            if ((first[xindex] != 2) && (second[xindex] != 2)) {
                n++;
                if ((first[xindex] == 1) && (second[xindex] == 1)) a++;
                if ((first[xindex] == 1) && (second[xindex] == 0)) b++;
                if ((first[xindex] == 0) && (second[xindex] == 1)) c++;
                if ((first[xindex] == 0) && (second[xindex] == 0)) d++;
            }
        }
        if (n == 0) return MAX_DIST;
        temp = (n - (a*P) - (d*Q));
        term = (4.0 * n * P * Q * (b+c));

        if (term >= (temp * temp)) return MAX_DIST;

        theta = (temp - sqrt((temp * temp) - term)) / (2.0 * n * P * Q);

        if (theta >= 1.0) return MAX_DIST;

        trans = pow((1.0 - (theta*P)),  (double) d) *
                pow((1.0 - (theta*Q)),  (double) a) *
                pow((theta * sqrt(P*Q)),  (double) (b+c));

    }
    return ((int) (-10.0 * log10(trans)));
}

static int rhmap5_edgelen (int i, int j, CCdatagroup *dat)
{
    char **vectors = dat->rhdat.vectors;
    int rhlength = dat->rhdat.rhlength;
    char *v1 = vectors[i];
    char *v2 = vectors[j];
    int n;
    int mis = 0;
    int cnt = 0;

    if (v1 == (char *) NULL || v2 == (char *) NULL) return 0;

    for (n=0; n<rhlength; n++) {
        if (v1[n] != 2 && v2[n] != 2) {
            cnt++;
            if (v1[n] != v2[n]) mis++;
        }
    }
    if (cnt == 0) return 0;
    else return (int) (10.0 * rhlength * mis) / cnt;
}
