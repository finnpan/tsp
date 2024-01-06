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
#define CC_BIX_GETOPT_UNKNOWN -3038


static int run_silently = 1;
static int kick_type = CC_LK_WALK_KICK;
static int number_runs = 0;
static char *nodefile = (char *) NULL;


int main (int, char **);
static void usage (char *f);
static int
   print_command (int ac, char **av),
   parseargs (int, char **),
   CCutil_bix_getopt (int ac, char **av, const char *def,
                       int *p_optind, char **p_optarg);


int main (int ac, char **av)
{
    int ncount;
    double val, best;
    double startzeit, kzeit;
    int edgeNum, *edgeList = (int *) NULL;
    int *incycle = (int *) NULL, *outcycle = (int *) NULL;
    int* *near = (int* *) NULL;
    CCdatagroup dat;

    int rval = print_command (ac, av);
    if (rval) {
        fprintf (stderr, "%s\n", "print_command failed");
        goto CLEANUP;
    }

    unsigned int seed = (unsigned int) time (0);
    srand(seed);

    if (parseargs (ac, av)) return 1;

    printf ("Chained Lin-Kernighan with seed %d\n", seed);
    fflush (stdout);

    if (!nodefile) {
        usage (av[0]);
        return 1;
    }

    startzeit = CCutil_zeit ();

    if (CCutil_gettsplib (nodefile, &ncount, &dat)) {
        fprintf (stderr, "could not read the TSPLIB file\n");
        rval = 1;
        goto CLEANUP;
    }

    int in_repeater = ncount;

    incycle = CC_SAFE_MALLOC (ncount, int);
    if (!incycle) {
        rval = 1;
        goto CLEANUP;
    }

    /* find k-nearest */
    kzeit = CCutil_zeit ();
    const int kNum = 50;
    double minCost;
    int nearCity;
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
    const int eNum = 10;
    int eCnt, listIdx;
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

    int iter = 0;
    best = BIGDOUBLE;
    do {
        printf ("\nStarting Run %d\n", iter);
        if (CClinkern_tour (ncount, &dat, edgeNum, edgeList, 100000000,
                in_repeater, incycle, outcycle, &val, run_silently, kick_type)) {
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


static int parseargs (int ac, char **av)
{
    int c, k;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "a:bBD:E:g:G:h:k:lI:K:N:o:q:Qr:R:s:S:t:y:Y:", &boptind, &boptarg)) != EOF)
        switch (c) {
        case 'K':
            k = atoi (boptarg);
            if (k == CC_LK_RANDOM_KICK)         kick_type = CC_LK_RANDOM_KICK;
            else if (k == CC_LK_CLOSE_KICK)     kick_type = CC_LK_CLOSE_KICK;
            else if (k == CC_LK_WALK_KICK)      kick_type = CC_LK_WALK_KICK;
            else fprintf (stderr, "unknown kick type, using default\n");
            break;
        case 'Q':
            run_silently++;
            break;
        case 'r':
            number_runs = atoi (boptarg);
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    if (boptind < ac)
        nodefile = av[boptind++];

    if (boptind > ac) {
        usage (av[0]);
        return 1;
    }
    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "usage: %s [- see below -] [tsplib_file or dat_file]\n", f);
    fprintf (stderr, "   -K #  kick (%d-Random, %d-Close, %d-Random_Walk [default])\n",
           CC_LK_RANDOM_KICK, CC_LK_CLOSE_KICK, CC_LK_WALK_KICK);
    fprintf (stderr, "   -r #  number of runs\n");
    fprintf (stderr, "   -Q    run silently\n");
}

static int print_command (int ac, char **av)
{
    int rval = 0;
    int i, cmdlen = 0;
    char *cmdout = (char *) NULL;

    for (i=0; i<ac; i++) {
        cmdlen += strlen(av[i]) + 1;
    }
    cmdout = CC_SAFE_MALLOC (cmdlen, char);
    if (!cmdout) {
        fprintf (stderr, "%s\n", "out of memory in print_command");
        rval = 1;
        goto CLEANUP;
    }

    cmdlen = 0;
    for (i=0; i<ac; i++) {
        strcpy (cmdout + cmdlen, av[i]);
        cmdlen += strlen(av[i]);
        cmdout[cmdlen] = ' ';
        cmdlen++;
    }
    cmdout[cmdlen-1] = '\0';
    printf ("%s\n", cmdout); fflush (stdout);

CLEANUP:

    CC_IFFREE (cmdout, char);
    return rval;
}

static int CCutil_bix_getopt (int ac, char **av, const char *def,
        int *p_optind, char **p_optarg)
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
