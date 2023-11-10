/*
License: see tree.h
 */

#include "tree.h"

TwoLevelTree::TwoLevelTree (int cityNum, int originCity) :
    _segNum(static_cast<int>(std::sqrt(cityNum)) + 1),
    _cityNum(cityNum), _originCity(originCity)
{
    assert(_cityNum > 0);
    assert(_originCity >= 0);
    assert(ParentNum() > 1);
    _parentNodes = new ParentNode[ParentNum()];
    _childNodes = new ChildNode[ChildNum()];
    _nominalSegLen = cityNum / ParentNum();
}

TwoLevelTree::~TwoLevelTree ()
{
    if (_parentNodes) {
        delete[] _parentNodes;
        _parentNodes = nullptr;
    }
    if (_childNodes) {
        delete[] _childNodes;
        _childNodes = nullptr;
    }
}

TwoLevelTree& TwoLevelTree::operator= (const TwoLevelTree& other) noexcept
{
    if (this == &other) {
        return *this;
    }
    delete[] _parentNodes;
    _parentNodes = new ParentNode[other.ParentNum()];
    delete[] _childNodes;
    _childNodes = new ChildNode[other.ChildNum()];
    _segNum = other._segNum;
    _cityNum = other._cityNum;
    _originCity = other._originCity;
    _nominalSegLen = other._nominalSegLen;

    auto* origNode = other.OriginCityNode();
    ChildNode::Iter begin(origNode, 0, _cityNum),
                    last(origNode->prev, _cityNum-1, _cityNum);
    SetTour(begin, last);
    return *this;
}

TwoLevelTree& TwoLevelTree::operator= (TwoLevelTree&& other)
{
    if (this == &other) {
        return *this;
    }
    _parentNodes = other._parentNodes;
    _childNodes = other._childNodes;
    _segNum = other._segNum;
    _cityNum = other._cityNum;
    _originCity = other._originCity;
    _nominalSegLen = other._nominalSegLen;
    other._parentNodes = nullptr;
    other._childNodes = nullptr;
    other._segNum = 0;
    other._cityNum = 0;
    other._originCity = 0;
    other._nominalSegLen = 0;
    return *this;
}

bool TwoLevelTree::Between (const ChildNode* a, const ChildNode* b,
                            const ChildNode* c) const
{
    assert(a != b && a != c && b != c);
    auto *pa = a->parent, *pb = b->parent, *pc = c->parent;

    // all same parents: in a single segment.
    if (pa == pb && pb == pc) {
        if (pa->reverse) {
            if (c->id < a->id) {
                return b->id < a->id && b->id > c->id;
            } else {
                return b->id < a->id || b->id > c->id;
            }
        } else {
            if (c->id > a->id) {
                return b->id > a->id && b->id < c->id;
            } else {
                return b->id > a->id || b->id < c->id;
            }
        }
    }

    /*
     * all three parents are distinct.
     * note that the parents are in a cyclcal list
     */
    if (pa != pb && pa != pc && pb != pc) {
        if (pc->id > pa->id) {
            return pb->id > pa->id && pb->id < pc->id;
        } else {
            return pb->id > pa->id || pb->id < pc->id;
        }
    }

    /*
     * two nodes share the parent, one different.
     * only true if we can reach v from u before leaving a the segment
     */
    if (pa == pb) {
        return IsPathInSingleSegment(a, b);
    } else if (pb == pc) {
        return IsPathInSingleSegment(b, c);
    } else if (pa == pc) {
        return !IsPathInSingleSegment(a, c);
    }

    return false;
}

/*
 * Remove two arcs (a, b) and (c, d), and add two others (a, c) and (b, d).
 * The two arcs should both be in forward or backward orientation.
 * Either (a, d) or (b, c) is reversed. The smaller one is preferred.
 */
void TwoLevelTree::Flip (ChildNode* a, ChildNode* b, ChildNode* c, ChildNode* d)
{
    bool isForward = GetNext(a) == b;
    assert((GetNext(c) == d) == isForward);
    assert(!((a == c) && (b == d)));
    if (b == c || d == a) { // in this case, even after flip, still the same
        return;
    }

    // reverse the old subpath (b, c) or (d, a)
    // and reconnection of (a, c)  and (b, d) is also performed automatically
    // thus, no need to delete the old arcs explicitly
    // we tend to reverse the shorter path for possibly reduced computation cost
    if (IsApproximatelyShorter(b, c, d, a)) {
        if (isForward) {
            Reverse(b, c);
        } else {
            Reverse(c, b);
        }
    } else {
        if (isForward) {
            Reverse(d, a);
        } else {
            Reverse(a, d);
        }
    }
}

/*
 * Perform a double-bridge move. Suppose the next node in a forward tour
 * of the four arguments are an, bn, cn and dn respectively.
 * Then, (a, an), (b, bn), (c, cn) and (d, dn) are removed.
 * New arcs (a, cn), (b, dn), (c, an), and (d, bn) are inserted. Note:
 * (1) The arguments a, b, c, d should be given in a forward tour order,
 * and there should be at least one other node between each pairs of them.
 * (2) Any two of a, b, c, d should lie in diffrent segments.
 */
void TwoLevelTree::DoubleBridgeMove (ChildNode* a, ChildNode* b,
                                     ChildNode* c, ChildNode* d)
{
    assert(Between(a, b, c));
    assert(Between(b, c, d));
    assert(Between(c, d, a));
    assert(Between(d, a, b));
    assert(a->parent != b->parent);
    assert(a->parent != c->parent);
    assert(a->parent != d->parent);
    assert(b->parent != c->parent);
    assert(b->parent != d->parent);
    assert(c->parent != d->parent);
    auto *an = GetNext(a), *bn = GetNext(b),
         *cn = GetNext(c), *dn = GetNext(d);

    // (1) split and merge to make all the above segment boundaries
    ChildNode* arr[4] = {a, b, c, d};
    for (auto* p : arr) {
        // special case: p and pn are already boundary nodes,
        // i.e., [.....p] -> [pn....]
        if (p->parent == GetNext(p)->parent) {
            SplitAndMerge(p, false, Direction::forward);
        }

#ifndef NDEBUG
        assert((p == p->parent->begin ||
                p == p->parent->end));
        auto* q = GetNext(p);
        assert((q == q->parent->begin ||
                q == q->parent->end));
        assert(p->parent->next == q->parent);
#endif
    }

    // (2) reconnect. Note that p and q are both segment boundary nodes.
    // must be connected in the right order
    ConnectForward(a, cn);
    ConnectForward(d, bn);
    ConnectForward(c, an);
    ConnectForward(b, dn);

    // (3) each segment itself is not changed due to reconnection
    // However, the order of the segments is changed and re-id is needed
    auto* p = HeadParentNode();
    int id = 0;
    do {
        p->id = id;
        id++;
        p = p->next;
    } while (p != HeadParentNode());
}

/*
 * To facilitate implementation, we use implicit rebalance here, because
 * in practice the complicated full rebalance is empricially unnecessary.
 */
void TwoLevelTree::Reverse (ChildNode* a, ChildNode* b)
{
    if (a == b || GetNext(b) == a) {
        return;
    }
    // (1) the path is contained in a single segment
    if (IsPathInSingleSegment(a, b)) {
        ReverseSegment(a, b);
        return;
    }

    // (2) multiple segments are involved,
    // we simply split and merge to make complete segments
    // split and merge the smaller half heuristically
    SplitAndMerge_0(a);
    // whether a and b are in a single segment now
    if (IsPathInSingleSegment(a, b)) {
        ReverseSegment(a, b);
        return;
    }

    // split and merge b similarly
    SplitAndMerge_1(a, b);
    // whether a and b are in a single segment now
    if (IsPathInSingleSegment(a, b)) {
        ReverseSegment(a, b);
        return;
    }

    // now the forward path a ----> b contains multiple complete segments
    // even after merge suppose s1 [a...] [....] [....] [....b] s2
    // note: in a forward path, we always go to .next for the parent node
    assert(a == a->parent->ForwardBeginNode());
    assert(b == b->parent->BackwardBeginNode());

    ParentNode *q = nullptr, *p = nullptr, *tmp = nullptr;
    auto* s1 = a->parent->prev;
    auto* s2 = b->parent->next;

    // (a) each segment between a and b should be reversed
    // (b) reverse the positions of the segments and reconnect them
    // between s1 and s2. We also need to update the ID
    // and adjust the connections in the ends of each segment
    p = s1;
    q = b->parent;
    while (q != s1) {
        tmp = q->prev;
        q->reverse = !q->reverse;
        ReverseSeg_0(p, q);
        p = q;
        q = tmp;
    }
    ReverseSeg_0(p, s2);

    // now p is s2
    assert((p->id + 1) % ParentNum() == p->next->id);
}

void TwoLevelTree::ReverseSegment (ChildNode* a, ChildNode* b)
{
    assert(a->parent == b->parent);
    auto* parent = a->parent;
    // if exactly a complete segment
    if ((a == parent->begin && b == parent->end) ||
        (b == parent->begin && a == parent->end)) {
        ReverseCompleteSegment(a, b);
    } else { // only a part of the segment
        auto pathLen = std::abs(a->id - b->id) + 1;  // IDs are consecutive
        if (pathLen <= _nominalSegLen * 3 / 4) {
            ReversePartialSegment(a, b);
        } else { // split at a and b and merge with their neighbors
            // leave a and b in the original segment
            // to make a complete segment for reversion
            SplitAndMerge(a, false, Direction::backward);
            SplitAndMerge(b, false, Direction::forward);
            ReverseCompleteSegment(a, b);
        }
    }
}

void TwoLevelTree::ReverseCompleteSegment (ChildNode* a, ChildNode* b)
{
    assert(a->parent == b->parent);
    auto* const parent = a->parent;
    assert(a == parent->ForwardBeginNode() && b == parent->BackwardBeginNode());

    // the following line cannot work directly in the flip method when (b, c)
    // is a single segment. auto prevA = GetPrev(a), nextB = GetNext(b);
    auto* prevA = parent->prev->ForwardEndNode();
    auto* nextB = parent->next->ForwardBeginNode();
    parent->reverse = !parent->reverse;
    // repair the 4 connections to the neighbor segments
    // prevA now should go to b
    if (prevA->parent->reverse) {
        prevA->prev = b;
    } else {
        prevA->next = b;
    }
    // a should now go to nextB
    if (parent->reverse) {
        a->prev = nextB;
    } else {
        a->next = nextB;
    }
    // nextB should go back to a
    if (nextB->parent->reverse) {
        nextB->next = a;
    } else {
        nextB->prev = a;
    }
    // b should go back to prevA
    if (parent->reverse) {
        b->next = prevA;
    } else {
        b->prev = prevA;
    }
}

void TwoLevelTree::ReversePartialSegment (ChildNode* a, ChildNode* b)
{
    // we need change the connections and the IDs,
    // and possibly the segment endpoints
    assert(a->parent == b->parent);
    auto* const parent = a->parent;
    auto *prevA = GetPrev(a), *nextB = GetNext(b);
    auto partialSegLen = std::abs(a->id - b->id) + 1;
    ChildNode *q = nullptr, *p = nullptr, *tmp = nullptr;
    int id = 0;

    // now we reconstruct the connections from
    // prevA -> b .. -> a -> nextB along the forward direction
    p = prevA;
    q = b;
    while (q != prevA) {
        tmp = GetPrev(q);
        // connection p with q, p -> q on forward tour
        ConnectArcForward(p, q);
        p = q;
        q = tmp;
    }
    ConnectArcForward(p, nextB);

    // if one of them is originally an endpoint (at most one can be)
    if (a == parent->begin) {
        parent->begin = b;
    } else if (a == parent->end) {
        parent->end = b;
    } else if (b == parent->begin) {
        parent->begin = a;
    } else if (b == parent->end) {
        parent->end = a;
    }

    // relabel the IDs for the forward path b --> a.
    // Note ID is numbered according to node.next.
    if (parent->reverse) { // a --next-- --next-- b
        if (a == parent->begin) {
            id = b->next->id - partialSegLen;
        } else {
            id = a->prev->id + 1;
        }
        RelabelId(a, b, id);
    } else { // b --next-- --next-- a
        if (b == parent->begin) {
            id = a->next->id - partialSegLen;
        } else {
            id = b->prev->id + 1;
        }
        RelabelId(b, a, id);
    }
}

/*
 * Split a segment at s, and merge one half to its neighbor segment
 * specified by the direction. If includeSelf is true, then the node s
 * is merged to its neighbor; otherwise, it stays in its own segment.
 */
void TwoLevelTree::SplitAndMerge (ChildNode* s, bool includeSelf, Direction dir)
{
    auto* const parent = s->parent;
    ParentNode* neighborParent = (dir == Direction::forward)?
                                  parent->next : parent->prev;

    // the new boundary of the parent segment after being splitted
    ChildNode* boundary = nullptr;
    // the node of neighbor segment to be connected to splited nodes
    ChildNode *q = nullptr, *p = nullptr, *tmp = nullptr;
    int deltaId = 0, cnt = 0;

    if (dir == Direction::forward) {
        boundary = includeSelf? GetPrev(s) : s;
        p = neighborParent->ForwardBeginNode();
        deltaId = neighborParent->reverse? 1 : -1;

        // find last one
        q = parent->ForwardEndNode();

        // merge these nodes to the neighbor
        while (q != boundary) {
            tmp = GetPrev(q);
            q->parent = neighborParent;
            ConnectArcForward(q, p);
            // relabel the newly merged part in the neighbor segment
            q->id = p->id + deltaId;
            p = q;
            ++cnt;
            q = tmp;
        }
        if (neighborParent->reverse) {
            neighborParent->end = p;
        } else {
            neighborParent->begin = p;
        }

        // repair the boundary of the old segment
        ConnectArcForward(boundary, p);
        if (parent->reverse) {
            parent->begin = boundary;
        } else {
            parent->end = boundary;
        }
    } else {
        boundary = includeSelf? GetNext(s) : s;
        p = neighborParent->BackwardBeginNode();
        deltaId = neighborParent->reverse? -1 : 1;

        // find last one
        q = parent->BackwardEndNode();

        // merge these nodes to the neighbor
        while (q != boundary) {
            tmp = GetNext(q);
            q->parent = neighborParent;
            ConnectArcForward(p, q);
            // relabel the newly merged part in the neighbor segment
            q->id = p->id + deltaId;
            p = q;
            ++cnt;
            q = tmp;
        }
        if (neighborParent->reverse) {
            neighborParent->begin = p;
        } else {
            neighborParent->end = p;
        }

        // repair the boundary of the old segment
        ConnectArcForward(p, boundary);
        if (parent->reverse) {
            parent->end = boundary;
        } else {
            parent->begin = boundary;
        }
    }

    neighborParent->size += cnt;
    parent->size -= cnt;
    assert(parent->size > 0); // we cannot leave an empty segment
}
