#include "LKH.h"

/*
 * The PrintParameters function prints the problem parameters to
 * standard output.
*/

void PrintParameters()
{
    printff("ASCENT_CANDIDATES = %d\n", AscentCandidates);
    printff("BACKBONE_TRIALS = %d\n", BackboneTrials);
    printff("BACKTRACKING = %s\n", Backtracking ? "YES" : "NO");
    if (Excess >= 0)
        printff("EXCESS = %g\n", Excess);
    else
        printff("# EXCESS =\n");
    printff("EXTRA_CANDIDATES = %d %s\n", ExtraCandidates,
            ExtraCandidateSetSymmetric ? "SYMMETRIC" : "");
    printff("GAIN23 = %s\n", Gain23Used ? "YES" : "NO");
    printff("GAIN_CRITERION = %s\n", GainCriterionUsed ? "YES" : "NO");
    if (InitialPeriod >= 0)
        printff("INITIAL_PERIOD = %d\n", InitialPeriod);
    else
        printff("# INITIAL_PERIOD =\n");
    printff("INITIAL_STEP_SIZE = %d\n", InitialStepSize);
    printff("KICK_TYPE = %d\n", KickType);
    printff("KICKS = %d\n", Kicks);
    if (MaxBreadth == INT_MAX)
        printff("# MAX_BREADTH =\n");
    else
        printff("MAX_BREADTH = %d\n", MaxBreadth);
    printff("MAX_CANDIDATES = %d %s\n", MaxCandidates,
            CandidateSetSymmetric ? "SYMMETRIC" : "");
    if (MaxSwaps >= 0)
        printff("MAX_SWAPS = %d\n", MaxSwaps);
    else
        printff("# MAX_SWAPS =\n");
    if (MaxTrials >= 0)
        printff("MAX_TRIALS = %d\n", MaxTrials);
    else
        printff("# MAX_TRIALS =\n");
    printff("MOVE_TYPE = %d\n", MoveType);
    printff("%sNONSEQUENTIAL_MOVE_TYPE = %d\n", PatchingA > 1 ? "" : "# ",
            NonsequentialMoveType);
    if (Optimum == MINUS_INFINITY)
        printff("# OPTIMUM =\n");
    else
        printff("OPTIMUM = " GainFormat "\n", Optimum);
    printff("PATCHING_A = %d %s\n", PatchingA,
            PatchingARestricted ? "RESTRICTED"
            : PatchingAExtended ? "EXTENDED"
                                : "");
    printff("PATCHING_C = %d %s\n", PatchingC,
            PatchingCRestricted ? "RESTRICTED"
            : PatchingCExtended ? "EXTENDED"
                                : "");
    printff("PRECISION = %d\n", Precision);
    printff("%sPROBLEM_FILE = %s\n", ProblemFileName ? "" : "# ",
            ProblemFileName ? ProblemFileName : "");
    printff("RESTRICTED_SEARCH = %s\n", RestrictedSearch ? "YES" : "NO");
    printff("RUNS = %d\n", Runs);
    printff("SEED = %u\n", Seed);
    printff("STOP_AT_OPTIMUM = %s\n", StopAtOptimum ? "YES" : "NO");
    printff("SUBGRADIENT = %s\n", Subgradient ? "YES" : "NO");
    printff("SUBSEQUENT_MOVE_TYPE = %d\n",
            SubsequentMoveType == 0 ? MoveType : SubsequentMoveType);
    printff("SUBSEQUENT_PATCHING = %s\n", SubsequentPatching ? "YES" : "NO");
    if (TimeLimit == DBL_MAX)
        printff("# TIME_LIMIT =\n");
    else
        printff("TIME_LIMIT = %0.1f\n", TimeLimit);
    if (TotalTimeLimit == DBL_MAX)
        printff("# TOTAL_TIME_LIMIT =\n");
    else
        printff("TOTAL_TIME_LIMIT = %0.1f\n", TotalTimeLimit);
    printff("TRACE_LEVEL = %d\n\n", TraceLevel);
}
