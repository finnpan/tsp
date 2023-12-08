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
/*                         SORTING ROUTINES                                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*   Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/*   DATE:  February 24, 1994                                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void CCutil_int_perm_quicksort (int *perm, int *len, int n)             */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void CCutil_double_perm_quicksort (int *perm, double *len, int n)       */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void CCutil_rselect (int *arr, int l, int r, int m,                     */
/*      double *coord, CCrandstate *rstate)                                 */
/*    arr - permutation that will be rearranged                             */
/*    l,r - specify the range of arr that we are interested in              */
/*    m - is the index into l,r that is the break point for the perm        */
/*    coord - gives the keys that determine the ordering                    */
/*                                                                          */
/****************************************************************************/

#include "config.h"
#include "util.h"

#define BITS_PER_PASS (8)


static void
    select_split (int *arr, int n, double v, int *start, int *end,
           double *coord),
    select_sort (int *arr, int n, double *coord),
    select_sort_dsample (double *samp, int n);


void CCutil_int_perm_quicksort (int *perm, int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1)
        return;

    CC_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        CC_SWAP (perm[i], perm[j], temp);
    }
    CC_SWAP (perm[0], perm[j], temp);

    CCutil_int_perm_quicksort (perm, len, j);
    CCutil_int_perm_quicksort (perm + i, len, n - i);
}


void CCutil_double_perm_quicksort (int *perm, double *len, int n)
{
    int i, j, temp;
    double t;

    if (n <= 1)
        return;

    CC_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        CC_SWAP (perm[i], perm[j], temp);
    }
    CC_SWAP (perm[0], perm[j], temp);

    CCutil_double_perm_quicksort (perm, len, j);
    CCutil_double_perm_quicksort (perm + i, len, n - i);
}


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
/*                    EDGE LIST UTILITY ROUTINES                            */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by: Applegate, Bixby, Chvatal, and Cook                         */
/*  Date: February 8, 1995                                                  */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_edge_to_cycle (int ncount, int *elist, int *yesno,           */
/*      int *cyc)                                                           */
/*    CONVERTS an edgelist to a cycle.                                      */
/*     -ncount is the number of nodes.                                      */
/*     -elist is an edgelist in end1 end2 format.                           */
/*     -yesno returns 1 if elist describes a tour and 0 otherwise.          */
/*     -cyc returns the cycle in permutation format if it is not NULL       */
/*      (if cyc is not NULL, then it should point to an array of            */
/*      length at least ncount).                                            */
/*     Returns a nonzero value if there was an error.                       */
/*                                                                          */
/****************************************************************************/


int CCutil_edge_to_cycle (int ncount, int *elist, int *yesno, int *cyc)
{
    int *Lside, *Rside;
    int i, k, end1, end2, prev, this, next, start, okfirst, first = 0;
    int rval = 0;

    *yesno = 0;

    Lside = CC_SAFE_MALLOC (ncount, int);
    Rside = CC_SAFE_MALLOC (ncount, int);
    if (!Lside || !Rside) {
        fprintf (stderr, "out of memory in CCutil_edge_to_cycle\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        Lside[i] = Rside[i] = -1;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        end1 = elist[k++];
        end2 = elist[k++];
        if (Lside[end1] == -1)
            Lside[end1] = end2;
        else
            Rside[end1] = end2;
        if (Lside[end2] == -1)
            Lside[end2] = end1;
        else
            Rside[end2] = end1;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        end1 = elist[k++];
        end2 = elist[k++];
        if (Lside[end1] == -1 || Rside[end1] == -1 ||
            Lside[end2] == -1 || Rside[end2] == -1) {
            *yesno = 0;  goto CLEANUP;
        }
    }
    start = elist[0];
    prev = -2;
    this = start;
    k = 0;
    okfirst = 0;
    do {
        if (this == first)
           okfirst = 1;
        if (Lside[this] != prev)
            next = Lside[this];
        else
            next = Rside[this];
        prev = this;
        this = next;
        k++;
    } while (next != start && k < ncount);

    if (k != ncount || !okfirst) {
        *yesno = 0;  goto CLEANUP;
    }

    *yesno = 1;

    if (cyc) {
        start = first;
        prev = -2;
        this = start;
        k = 0;
        do {
            cyc[k++] = this;
            if (Lside[this] != prev)
                next = Lside[this];
            else
                next = Rside[this];
            prev = this;
            this = next;
        } while (next != start && k < ncount);
    }


CLEANUP:

    CC_IFFREE (Lside, int);
    CC_IFFREE (Rside, int);

    return rval;
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
