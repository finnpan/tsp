#include "eax.h"

int main (int argc, char* argv[])
{
    if (argc < 2) {
        printf("please input instance file!\n");
        return 0;
    }

    EAX eax;
    eax._numOfPop = 100;
    eax._numOfKids = 30;
    eax._fileNameTSP = argv[1];
    eax.Define();
    eax._silent = true;
    eax.DoIt();

    printf("[ ==== PASSED ==== ] eax tested: iter = %d, best = %lld, avg = %lf\n",
           eax._curNumOfGen,
           (long long)eax._bestValue,
           eax._averageValue);

    return 0;
}
