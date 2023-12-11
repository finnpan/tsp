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
#include "kdtree.h"

#define BIGDOUBLE (1e30)
#define CC_BIX_GETOPT_UNKNOWN -3038


static int norm = CC_EUCLIDEAN;
static int run_silently = 1;
static int kick_type = CC_LK_WALK_KICK;
static int number_runs = 0;
static char *nodefile = (char *) NULL;


int
   main (int, char **);
static int
   print_command (int ac, char **av),
   parseargs (int, char **);
static void
   usage (char *f);
static int
    CCutil_bix_getopt (int ac, char **av, const char *def,
                       int *p_optind, char **p_optarg);



int main (int ac, char **av)
{
    int k, ncount;
    double val, best;
    double startzeit;
    int tempcount, *templist;
    int *incycle = (int *) NULL, *outcycle = (int *) NULL;
    CCdatagroup dat;
    int rval = 0, in_repeater;
    int quadtry = 3;

    rval = print_command (ac, av);
    if (rval) {
        fprintf (stderr, "%s\n", "print_command failed");
        goto CLEANUP;
    }

    unsigned int seed = (unsigned int) time (0);
    srand(seed);

    if (parseargs (ac, av))
        return 1;

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
    norm = dat.norm;

    in_repeater = ncount;

    incycle = CC_SAFE_MALLOC (ncount, int);
    if (!incycle) {
        rval = 1;
        goto CLEANUP;
    }

    if ((norm & CC_NORM_BITS) == CC_KD_NORM_TYPE) {
        CCkdtree localkt;
        double kzeit = CCutil_zeit ();

        if (CCkdtree_build (&localkt, ncount, &dat, (double *) NULL)) {
            fprintf (stderr, "CCkdtree_build failed\n");
            rval = 1;
            goto CLEANUP;
        }
        printf ("Time to build kdtree: %.2f\n", CCutil_zeit () - kzeit);
        fflush (stdout);

        kzeit = CCutil_zeit ();
        if (CCkdtree_quadrant_k_nearest (&localkt, ncount, quadtry,
                &dat, (double *) NULL, 1, &tempcount, &templist,run_silently)) {
            fprintf (stderr, "CCkdtree-quad nearest code failed\n");
            rval = 1;
            goto CLEANUP;
        }
        if (!run_silently) {
            printf ("Time to find quad %d-nearest: %.2f\n",
                    quadtry, CCutil_zeit () - kzeit);
            fflush (stdout);
        }

        kzeit = CCutil_zeit ();
        CCutil_randcycle (ncount, incycle);
        if (!run_silently) {
            printf ("Time to grow tour: %.2f\n",
                    CCutil_zeit () - kzeit);
            fflush (stdout);
        }
        CCkdtree_free (&localkt);
    } else if ((norm & CC_NORM_BITS) == CC_X_NORM_TYPE) {
        // do nothing
    } else {
        assert(0);
    }

    outcycle = CC_SAFE_MALLOC (ncount, int);
    if (!outcycle) {
        rval = 1;
        goto CLEANUP;
    }

    k = 0;
    best = BIGDOUBLE;
    do {
        printf ("\nStarting Run %d\n", k);
        if (CClinkern_tour (ncount, &dat, tempcount, templist, 100000000,
                in_repeater, incycle, outcycle, &val, run_silently, kick_type)) {
            fprintf (stderr, "CClinkern_tour failed\n");
            rval = 1;
            goto CLEANUP;
        }
        if (val < best) {
            best = val;
        }
    } while (++k < number_runs);
    printf ("Overall Best Cycle: %.0f\n", val);
    printf ("Total Running Time: %.2f\n", CCutil_zeit () - startzeit);
    fflush (stdout);

CLEANUP:

#ifndef BIG_PROBLEM
    CC_IFFREE (templist, int);
#endif
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
            else if (k == CC_LK_GEOMETRIC_KICK) kick_type = CC_LK_GEOMETRIC_KICK;
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
    fprintf (stderr, "   -K #  kick (%d-Random, %d-Geometric, %d-Close, %d-Random_Walk [default])\n",
           CC_LK_RANDOM_KICK, CC_LK_GEOMETRIC_KICK,
           CC_LK_CLOSE_KICK, CC_LK_WALK_KICK);
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
