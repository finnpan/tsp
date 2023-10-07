#include <stdio.h>
#include <ctime>

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

#if 1
    eax->DoIt();
#else
    // test kopt
    Indi indi1, indi2;
    indi1.Define(eax->_eval->_nCity);
    indi2.Define(eax->_eval->_nCity);

    EvalType e1 = 0, e2 = 0;
    clock_t start, end;
    clock_t t1 = 0, t2 = 0;

    printf("start testing kopt...\n");
    for (int i = 0; i < 10; i++) {
        eax->_kopt->MakeRandSol(indi1);

        start = clock();
        eax->_kopt->TransIndiToTree(indi1);
        eax->_kopt->Sub();
        eax->_kopt->TransTreeToIndi(indi2);
        end = clock();
        e1 += indi2._cost;
        t1 += end - start;

        start = clock();
        eax->_kopt->TransIndiToTree2(indi1);
        eax->_kopt->Sub2();
        eax->_kopt->TransTreeToIndi2(indi2);
        end = clock();
        e2 += indi2._cost;
        t2 += end - start;
    }

    printf("cost1 = %lld, time1 = %f\n", e1, (double)(t1)/CLOCKS_PER_SEC);
    printf("cost2 = %lld, time2 = %f\n", e2, (double)(t2)/CLOCKS_PER_SEC);
#endif

    delete eax;

    return 0;
}
