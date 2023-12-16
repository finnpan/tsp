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
