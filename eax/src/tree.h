/*
Original License:
  MIT License
  Copyright (c) 2019 Shuhua Gao
  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:
  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

Repo:
  https://github.com/ShuhuaGao/two-level-tree
Commit:
  805f00d6c16e525e7aabb66c064eb950e891fb1b
 */

#pragma once

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <utility>

class TwoLevelTree {
public:
    enum class Direction { forward, backward };

    class ChildNode;
    class ParentNode {
    public:
        bool reverse = false;
        int size = 0;
        int id = 0;
        ParentNode* prev = nullptr;
        ParentNode* next = nullptr;
        ChildNode* begin = nullptr;
        ChildNode* end = nullptr;

        ChildNode* ForwardEndNode() const { return reverse ? begin : end; }
        ChildNode* ForwardBeginNode() const { return reverse ? end : begin; }
        ChildNode* BackwardBeginNode() const { return reverse ? begin : end; }
        ChildNode* BackwardEndNode() const { return reverse ? end : begin; }
    };

    class ChildNode {
    public:
        int id = 0;
        int city = -1;
        ChildNode* prev = nullptr;
        ChildNode* next = nullptr;
        ParentNode* parent = nullptr;

        class Iter {
            const ChildNode* _node;
            int _idx;
            const int _len;

        public:
            Iter(const ChildNode* node, int idx, int len)
                : _node(node), _idx(idx), _len(len)
            {
            }
            int operator*() { return _node->city; }
            int operator-(const Iter& oth) const { return _idx - oth._idx; }
            Iter& operator++()
            {
                if (_idx < _len - 1) {
                    _node = _node->next;
                    ++_idx;
                }
                return *this;
            }
            Iter operator++(int)
            {
                Iter cpy = *this;
                ++(*this);
                return cpy;
            }
            Iter& operator--()
            {
                if (_idx > 0) {
                    _node = _node->prev;
                    --_idx;
                }
                return *this;
            }
            Iter operator--(int)
            {
                Iter cpy = *this;
                --(*this);
                return cpy;
            }
        };
    };

private:
    // store the double linked list in an array for fast access
    ParentNode* _parentNodes = nullptr;
    ChildNode* _childNodes = nullptr;
    int _segNum = 0;
    int _cityNum = 0;
    int _originCity = 0;
    int _nominalSegLen = 0;

public:
    TwoLevelTree(int cityNum, int originCity = 0);
    ~TwoLevelTree();
    TwoLevelTree(const TwoLevelTree& other) noexcept { *this = other; }
    TwoLevelTree& operator=(const TwoLevelTree& other) noexcept;
    TwoLevelTree(TwoLevelTree&& other) { *this = std::move(other); }
    TwoLevelTree& operator=(TwoLevelTree&& other);

    int GetNext(int current) const { return GetNext(GetNode(current))->city; }
    int GetPrev(int current) const { return GetPrev(GetNode(current))->city; }

    bool Between(int a, int b, int c) const
    {
        return Between(GetNode(a), GetNode(b), GetNode(c));
    }

    void Flip(int a, int b, int c, int d)
    {
        Flip(GetNode(a), GetNode(b), GetNode(c), GetNode(d));
    }

    void DoubleBridgeMove(int a, int b, int c, int d)
    {
        DoubleBridgeMove(GetNode(a), GetNode(b), GetNode(c), GetNode(d));
    }

    int ParentNum() const { return _segNum; }
    int ChildNum() const { return _cityNum; }
    int OriginCity() const { return _originCity; }

    bool IsCityValid(int city) const
    {
        return city >= _originCity && city < _originCity + _cityNum;
    }

    template <typename iterator>
    void SetTour(iterator begin, iterator last);

private:
    bool Between(const ChildNode* a, const ChildNode* b,
                 const ChildNode* c) const;

    void Flip(ChildNode* a, ChildNode* b, ChildNode* c, ChildNode* d);

    void DoubleBridgeMove(ChildNode* a, ChildNode* b, ChildNode* c,
                          ChildNode* d);

    void Reverse(ChildNode* a, ChildNode* b);

    void ReverseSegment(ChildNode* a, ChildNode* b);

    void ReverseCompleteSegment(ChildNode* a, ChildNode* b);

    void ReversePartialSegment(ChildNode* a, ChildNode* b);

    void SplitAndMerge(ChildNode* s, bool includeSelf, Direction dir);

    const ChildNode* GetNode(int city) const
    {
        assert(IsCityValid(city));
        return &_childNodes[city - _originCity];
    }

    ChildNode* GetNode(int city)
    {
        auto* ct = const_cast<const TwoLevelTree*>(this);
        return const_cast<ChildNode*>(ct->GetNode(city));
    }

    const ParentNode* GetParentNode(int city) const
    {
        return GetNode(city)->parent;
    }

    ParentNode* GetParentNode(int city)
    {
        auto* ct = const_cast<const TwoLevelTree*>(this);
        return const_cast<ParentNode*>(ct->GetParentNode(city));
    }

    const ParentNode* HeadParentNode() const { return &_parentNodes[0]; }

    ParentNode* HeadParentNode()
    {
        auto* ct = const_cast<const TwoLevelTree*>(this);
        return const_cast<ParentNode*>(ct->HeadParentNode());
    }

    const ParentNode* TailParentNode() const
    {
        return &_parentNodes[_segNum - 1];
    }

    ParentNode* TailParentNode()
    {
        auto* ct = const_cast<const TwoLevelTree*>(this);
        return const_cast<ParentNode*>(ct->TailParentNode());
    }

    const ChildNode* OriginCityNode() const { return GetNode(_originCity); }

    ChildNode* OriginCityNode()
    {
        auto* ct = const_cast<const TwoLevelTree*>(this);
        return const_cast<ChildNode*>(ct->OriginCityNode());
    }

    const ChildNode* GetNext(const ChildNode* current) const
    {
        if (current->parent->reverse) {
            return current->prev;
        }
        return current->next;
    }

    ChildNode* GetNext(const ChildNode* current)
    {
        auto* ct = const_cast<const TwoLevelTree*>(this);
        return const_cast<ChildNode*>(ct->GetNext(current));
    }

    const ChildNode* GetPrev(const ChildNode* current) const
    {
        if (current->parent->reverse) {
            return current->next;
        }
        return current->prev;
    }

    ChildNode* GetPrev(const ChildNode* current)
    {
        auto* ct = const_cast<const TwoLevelTree*>(this);
        return const_cast<ChildNode*>(ct->GetPrev(current));
    }

    bool IsApproximatelyShorter(ChildNode* a, ChildNode* b, ChildNode* c,
                                ChildNode* d) const
    {
        int abSegNum = CountSegmants(a, b);
        int cdSegNum = CountSegmants(c, d);
        // consider number of segments firstly due to imbalance
        if (abSegNum != cdSegNum) {
            return abSegNum < cdSegNum;
        }
        int lenA = std::abs(a->id - a->parent->ForwardBeginNode()->id);
        int lenB = std::abs(b->id - b->parent->ForwardEndNode()->id);
        int lenC = std::abs(c->id - c->parent->ForwardBeginNode()->id);
        int lenD = std::abs(d->id - d->parent->ForwardEndNode()->id);
        return lenA + lenB > lenC + lenD;
    }

    // forward path a --> b
    int CountSegmants(ChildNode* a, ChildNode* b) const
    {
        int n = ParentNum();
        auto apid = a->parent->id, bpid = b->parent->id;
        // how many segments are involved in the forward path a --> b
        if (apid == bpid) {
            // whether the forward a --> b is in the single segment
            if (IsPathInSingleSegment(a, b)) {
                return 1;
            }
            return n;
        } else if (bpid > apid) {
            return bpid - apid + 1;
        }
        return bpid + n - apid + 1;
    }

    bool IsPathInSingleSegment(const ChildNode* a, const ChildNode* b) const
    {
        if (a->parent == b->parent) {
            if (a->parent->reverse) {
                return a->id > b->id;
            }
            return a->id < b->id;
        }
        return false;
    }

    void ConnectArcForward(ChildNode* p, ChildNode* q) const
    {
        if (p->parent->reverse) {
            p->prev = q;
        } else {
            p->next = q;
        }
        if (q->parent->reverse) {
            q->next = p;
        } else {
            q->prev = p;
        }
    }

    void ConnectForward(ChildNode* p, ChildNode* q) const
    {
        ConnectArcForward(p, q);
        p->parent->next = q->parent;
        q->parent->prev = p->parent;
    }

    void RelabelId(ChildNode* a, ChildNode* b, int idA) const
    {
        assert(a->parent == b->parent);
        a->id = idA;
        while (a != b) {
            a->next->id = a->id + 1;
            a = a->next;
        }
    }

    void SplitAndMerge_0(ChildNode* a)
    {
        if (a == a->parent->ForwardBeginNode()) {
            return;
        }
        auto* forwEndA = a->parent->ForwardEndNode();
        int forwHalfLenA = std::abs(forwEndA->id - a->id) + 1;
        if (forwHalfLenA <= a->parent->size / 2) {
            SplitAndMerge(a, true, Direction::forward);
        } else {
            SplitAndMerge(a, false, Direction::backward);
        }
    }

    void SplitAndMerge_1(ChildNode* a, ChildNode* b)
    {
        if (b == b->parent->BackwardBeginNode()) {
            return;
        }
        // to handle the special cases: [......b..] -> [a......]
        // (i.e., reverse almost a full circle)
        if (b->parent->next == a->parent) {
            SplitAndMerge(b, true, Direction::backward);
            return;
        }
        auto* backEndB = b->parent->BackwardEndNode();
        auto backHalfLenB = std::abs(backEndB->id - b->id) + 1;
        if (backHalfLenB <= b->parent->size / 2) {
            SplitAndMerge(b, true, Direction::backward);
        } else {
            SplitAndMerge(b, false, Direction::forward);
        }
    }

    void ReverseSeg_0(ParentNode* p, ParentNode* q) const
    {
        // p -> q forward
        p->next = q;
        q->prev = p;
        // all parents nodes are placed in a cyclic list
        q->id = (p->id + 1) % ParentNum();
        // the neighbor nodes of p and q segments should be connected properly
        ConnectArcForward(p->ForwardEndNode(), q->ForwardBeginNode());
    }
};

template <typename iterator>
void TwoLevelTree::SetTour(iterator begin, iterator last)
{
    const int tourLen = last - begin + 1;
    (void)tourLen;
    assert(tourLen == _cityNum);
    const int segNum = ParentNum();
    const int segLen = _cityNum / segNum;

    ParentNode* parent;
    int s, e, c, city;
    ChildNode* node;
    iterator iter = begin;

    for (int segId = 0; segId < segNum; segId++) {
        // first build the parent for this segment
        parent = &_parentNodes[segId];
        parent->id = segId;
        parent->prev =
            (segId > 0) ? &_parentNodes[segId - 1] : TailParentNode();
        parent->next =
            (segId + 1 < segNum) ? &_parentNodes[segId + 1] : HeadParentNode();
        parent->reverse = false;
        // this segment range in the given order tour (the end excluded)
        s = segId * segLen;
        e = s + segLen;
        // the last segment takes all the remaining cities
        if (segId == segNum - 1) {
            e = _cityNum;
        }

        // build the segment node one by one
        parent->size = e - s;
        parent->begin = GetNode(*iter);
        c = s;
        while (c < e) {
            city = *iter;
            node = GetNode(city);
            node->city = city;
            node->parent = parent;
            // cycle tour
            if (c == 0) {
                node->prev = GetNode(*last);
            } else {
                node->prev = GetNode(*(--iter));
                ++iter;  // restore
            }
            if (c + 1 == _cityNum) {
                node->next = GetNode(*begin);
            } else {
                node->next = GetNode(*(++iter));
                --iter;  // restore
            }
            node->id = c - s;
            ++c;
            ++iter;
        }
        if (c == _cityNum) {
            parent->end = GetNode(*last);
        } else {
            parent->end = GetNode(*(--iter));
            ++iter;  // restore
        }
    }
}
