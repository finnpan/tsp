#include "eax.h"

int main (int argc, char* argv[])
{
    if (argc < 2) {
        printf("Error: please input instance file!\n");
        return 0;
    }

    EAXGA ga;
    ga._silent = true;
    ga.Define(argv[1]);
    ga.DoIt();

    printf("\e[1;32m[ ==== PASSED ==== ] EAXGA tested: ");
    printf("iter = %d, best = %lld, avg = %lf\n",
           ga._numGen, (long long)ga._best._cost, ga._avgCost);
    printf("\e[0m");

    return 0;
}
