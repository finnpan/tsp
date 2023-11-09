#include <ctime>

#define private public
#include "tsp.h"
#undef private

int main (int argc, char* argv[])
{
    if (argc < 2) {
        printf("please input instance file!\n");
        return 0;
    }

    Evaluator eval;
    eval.SetInstance(argv[1]);

    Kopt kopt(&eval);

    EvalType e = 0;
    clock_t t = 0;
    Indi indi;
    indi.Define(eval._nCity);

    constexpr int times = 100;
    for (int i = 0; i < times; i++) {
        clock_t start = clock();
        kopt.DoIt(indi);
        t += clock() - start;
        e += indi._cost;
    }

    printf("[ ==== PASSED ==== ] kopt tested: iter = %d, avg = %5.5lf, time = %5.5f\n",
           times,
           (double)e/times,
           (double)(t)/CLOCKS_PER_SEC);

    return 0;
}
