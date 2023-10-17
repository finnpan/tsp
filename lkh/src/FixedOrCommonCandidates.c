#include "LKH.h"

/* 
 * The FixedCandidates function returns the number of fixed or
 * common candidate edges emanating from a given node, N.
 */

int FixedCandidates(Node * N)
{
    int Count = 0;

    Count = N->FixedTo2 ? 2 : N->FixedTo1 ? 1 : 0;
    if (Count > 2)
        eprintf("Node %d has more than two required candidate edges",
                N->Id);
    return Count;
}
