#include "LKH.h"

/*
 * The AddTourCandidates function extends the candidate set with tour
 * edges given in the tour files.
 *   
 * The function is called from GenerateCandidateSet and OrderCandidateSet.  
*/

void AddTourCandidates()
{
    Node *Na;

    /* Add fixed edges */
    Na = FirstNode;
    do {
        if (Na->FixedTo1)
            AddCandidate(Na, Na->FixedTo1, D(Na, Na->FixedTo1), 0);
        if (Na->FixedTo2)
            AddCandidate(Na, Na->FixedTo2, D(Na, Na->FixedTo2), 0);
    }
    while ((Na = Na->Suc) != FirstNode);
}
