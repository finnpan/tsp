#include "LKH.h"

/*
 * This file contains the main function of the program.
 */

int main(int argc, char* argv[])
{
    GainType Cost;
    double Time, LastTime;

    /* Read the specification of the problem */
    if (argc >= 2) ParameterFileName = argv[1];
    ReadParameters();
    StartTime = LastTime = GetTime();
    MaxMatrixDimension = 20000;
    MergeWithTour = MergeWithTourIPT;
    ReadProblem();

    AllocateStructures();
    CreateCandidateSet();
    InitializeStatistics();

    if (Norm != 0)
        BestCost = PLUS_INFINITY;
    else {
        /* The ascent has solved the problem! */
        Optimum = BestCost = (GainType)LowerBound;
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        Runs = 0;
    }

    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        if (LastTime - StartTime >= TotalTimeLimit) {
            if (TraceLevel >= 1) printff("*** Time limit exceeded ***\n");
            Run--;
            break;
        }
        Cost = FindTour(); /* using the Lin-Kernighan heuristic */
        if (Run > 1) Cost = MergeTourWithBestTour();
        if (Cost < BestCost) {
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
        }
        if (Cost < Optimum) {
            Optimum = Cost;
            printff("*** New optimum = " GainFormat " ***\n", Optimum);
        }
        Time = fabs(GetTime() - LastTime);
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
            printff("Run %d: Cost = " GainFormat, Run, Cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff(", Gap = %0.4f%%", 100.0 * (Cost - Optimum) / Optimum);
            printff(", Time = %0.2f sec. %s\n\n", Time,
                    Cost < Optimum    ? "<"
                    : Cost == Optimum ? "="
                                      : "");
        }
        SRandom(++Seed);
    }
    PrintStatistics();
    return EXIT_SUCCESS;
}
