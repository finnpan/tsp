#include <ctime>

#define private public
#include "tsp.h"
#undef private

int main (int argc, char* argv[])
{
    if (argc < 2) {
        printf("Error: please input instance file!\n");
        return 0;
    }

    Evaluator eval;
    eval.SetInstance(argv[1]);

    KOpt kopt(&eval);

    EvalType e = 0;
    clock_t t = 0;
    Indi indi;
    indi.Define(eval._numCity);

    constexpr int times = 20;
    for (int i = 0; i < times; i++) {
        clock_t start = clock();
        kopt.DoIt(indi);
        t += clock() - start;
        e += indi._cost;
    }

    printf("\e[1;32m[ ==== PASSED ==== ] KOpt tested: ");
    printf("iter = %d, avg = %5.5lf, time = %5.5f\n",
           times, (double)e/times, (double)(t)/CLOCKS_PER_SEC);
    printf("\e[0m");

    return 0;
}
