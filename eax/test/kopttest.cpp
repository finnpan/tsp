#include <stdio.h>
#include <ctime>
#include <cstring>

#define private public
#include "kopt.h"
#undef private

#include "eax.h"

int main (int argc, char* argv[])
{
    if (argc < 2) {
        printf("please input instance file!\n");
        return 0;
    }

    EAX* eax = new EAX();
    eax->_numOfPop = 100;
    eax->_numOfKids = 30;
    eax->_fileNameTSP = argv[1];
    eax->Define();

    // test kopt
    Indi indi[2];
    indi[0].Define(eax->_eval->_nCity);
    indi[1].Define(eax->_eval->_nCity);

    EvalType e[2];
    clock_t t[2];
    memset(e, 0, sizeof(e));
    memset(t, 0, sizeof(t));

    printf("start testing kopt...\n");
    for (int i = 0; i < 1000; i++) {
        clock_t start;
        eax->_kopt->MakeRandSol(indi[0]);

        if (i % 100 == 0) {
            printf("%3d, old ...\n", i);
        }
        start = clock();
        eax->_kopt->TransIndiToTree(indi[0]);
        eax->_kopt->Sub();
        eax->_kopt->TransTreeToIndi(indi[1]);
        t[0] += clock() - start;
        e[0] += indi[1]._cost;

        if (i % 100 == 0) {
            printf("%3d, new ...\n", i);
        }
        start = clock();
        eax->_kopt->TransIndiToTree_tlt(indi[0]);
        eax->_kopt->Sub_tlt();
        eax->_kopt->TransTreeToIndi_tlt(indi[1]);
        t[1] += clock() - start;
        e[1] += indi[1]._cost;
    }

    printf("old kopt: cost = %lld, time = %5.5f\n", e[0], (double)(t[0])/CLOCKS_PER_SEC);
    printf("new kopt: cost = %lld, time = %5.5f\n", e[1], (double)(t[1])/CLOCKS_PER_SEC);

    delete eax;

    return 0;
}
