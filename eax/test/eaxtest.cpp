#include "eax.h"

int main (int argc, char* argv[])
{
    if (argc < 2) {
        printf("Error: please input instance file!\n");
        return 0;
    }

    EAX eax;
    eax._silent = true;
    eax._fileNameTSP = argv[1];
    eax.Define();
    eax.DoIt();

    printf("\e[1;32m[ ==== PASSED ==== ] eax tested: ");
    printf("iter = %d, best = %lld, avg = %lf\n",
           eax._curNumGen, (long long)eax._bestValue, eax._averageValue);
    printf("\e[0m");

    return 0;
}
