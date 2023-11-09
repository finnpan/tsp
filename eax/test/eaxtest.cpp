#include <stdio.h>

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

    eax->DoIt();

    delete eax;

    return 0;
}
