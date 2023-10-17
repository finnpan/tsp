#include "Segment.h"
#include "LKH.h"

/*
 * The SegmentSize function returns the number of nodes in the 
 * tour segment between two given nodes in the current direction. 
 * Note, however, that if the two-level tree is used,
 * the number of nodes is only approximate (for efficiency reasons).
 * 
 * Time complexity: O(1).
 */

#ifdef TWO_LEVEL_TREE

int SegmentSize(Node * ta, Node * tb)
{
    Segment *Pa, *Pb;
    int nLeft, nMid, nRight;

    Pa = ta->Parent;
    Pb = tb->Parent;
    if (Pa == Pb) {
        int n = Reversed == Pa->Reversed ? tb->Rank - ta->Rank :
            ta->Rank - tb->Rank;
        return (n < 0 ? n + Dimension : n) + 1;
    }
    nLeft =
        Reversed ==
        Pa->Reversed ? Pa->Last->Rank - ta->Rank : ta->Rank -
        Pa->First->Rank;
    if (nLeft < 0)
        nLeft += Pa->Size;
    nMid = !Reversed ? Pb->Rank - Pa->Rank : Pa->Rank - Pb->Rank;
    if (nMid < 0)
        nMid += Groups;
    nMid = nMid == 2 ? (!Reversed ? Pa->Suc : Pa->Pred)->Size
        : (nMid - 1) * GroupSize;
    nRight =
        Reversed ==
        Pb->Reversed ? tb->Rank -
        Pb->First->Rank : Pb->Last->Rank - tb->Rank;
    if (nRight < 0)
        nRight += Pb->Size;
    return nLeft + nMid + nRight + 2;
}

#endif
