/*
License: see tree.h
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>
#include <cassert>

#define private public
#include "tree.h"
#undef private

// whether a and b are neighbors and a is before b on a forward tour
static bool IsOrderedNeighbor(const TwoLevelTree& tree, int a, int b)
{
    return tree.GetNext(tree.GetNode(a)) == tree.GetNode(b)
           && tree.GetNode(a) == tree.GetPrev(tree.GetNode(b));
}

static bool Between(const TwoLevelTree& tree, int a, int b, int c)
{
    return tree.Between(tree.GetNode(a), tree.GetNode(b), tree.GetNode(c));
}

static std::vector<int> GetTourViaParent(TwoLevelTree& tree, int start_city)
{
    std::vector<int> ans;
    auto p = tree.GetParentNode(start_city);
    do {
        auto q = p->ForwardBeginNode();
        while (q != p->ForwardEndNode()) {
            ans.push_back(q->city);
            q = tree.GetNext(q);
        }
        ans.push_back(q->city);
        p = p->next;
    } while (p != tree.GetParentNode(start_city));
    return ans;
}

static std::vector<int>
GetTour(const TwoLevelTree& tree, int start = -1,
        TwoLevelTree::Direction dir = TwoLevelTree::Direction::forward)
{
    std::vector<int> tour(tree._cityNum);
    if (start < 0) {
        start = tree._originCity;
    }
    auto* node = tree.GetNode(start);
    for (int i = 0; i < tree._cityNum; ++i) {
        tour[i] = node->city;
        if (dir == TwoLevelTree::Direction::forward) {
            node = tree.GetNext(node);
        } else {
            node = tree.GetPrev(node);
        }
    }
    return tour;
}

void Move2Opt(TwoLevelTree& tree, int t1, int t2, int t3, int t4)
{
    tree.Flip(t1, t2, t4, t3);
}

void Undo2OptMove(TwoLevelTree& tree, int t1, int t2, int t3, int t4)
{
    //Move2Opt(tree, t2, t3, t4, t1);
    tree.Flip(t2, t3, t1, t4);
}

static std::vector<int> ActualSegmentSizes(const TwoLevelTree& tree,
                                           int start = -1)
{
    if (tree.IsCityValid(start)) {
        std::vector<int> ans;
        ans.reserve(tree.ParentNum());
        auto* parent = tree.GetParentNode(start);
        auto* p = parent;
        do {
            ans.push_back(p->size);
            p = p->next;
        } while (p != parent);
        return ans;
    } else {
        std::vector<int> ans(tree.ParentNum());
        for (int i = 0; i < tree.ParentNum(); i++) {
            ans[i] = tree._parentNodes[i].size;
        }
        return ans;
    }
}

TEST(Build_tree_from_an_ordered_list_of_cities, 1)
{
    int cityNum = 67, start_city = 2;
    std::vector<int> order(cityNum);
    std::iota(order.begin(), order.end(), start_city);
    std::shuffle(order.begin(), order.end(), std::mt19937{123});
    TwoLevelTree tree{cityNum, start_city};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);

    EXPECT_EQ(tree._cityNum, cityNum);
    EXPECT_EQ(tree.ParentNum(), (int)std::sqrt(cityNum) + 1);

    // traverse the tour
    for (int i = 0; i < cityNum; i++) {
        auto city = order[i];
        auto node = tree.GetNode(city);
        EXPECT_EQ(node->city, city);
        if (i + 1 < cityNum) {
            EXPECT_EQ(node->next->city, order[i + 1]);
        } else {
            EXPECT_EQ(node->next->city, order[0]);
        }
        if (i - 1 >= 0) {
            EXPECT_EQ(node->prev->city, order[i - 1]);
        } else {
            EXPECT_EQ(node->prev->city, order.back());
        }
    }

    // a tour is a cycle
    // (1) next orientation
    int count_city = 0;
    auto node = tree.OriginCityNode();
    do {
        node = node->next;
        count_city++;
    } while (node != tree.OriginCityNode());
    EXPECT_EQ(count_city, cityNum);
    // (2) prev orientation
    node = tree.OriginCityNode();
    count_city = 0;
    do {
        node = node->prev;
        count_city++;
    } while (node != tree.OriginCityNode());
    EXPECT_EQ(count_city, cityNum);

    // check each segment
    count_city = 0;
    auto parent_node = tree.HeadParentNode();
    do {
        EXPECT_EQ(parent_node->end->next, parent_node->next->begin);
        EXPECT_EQ(parent_node->begin->prev, parent_node->prev->end);
        count_city += parent_node->size;
        parent_node = parent_node->next;

    } while (parent_node != tree.HeadParentNode());
    EXPECT_EQ(count_city, cityNum);

    // // in the initial tour, all segments are splitted in order
    auto first_city = order.front(), last_city = order.back();
    EXPECT_EQ(tree.GetNode(first_city)->parent, tree.HeadParentNode());
    EXPECT_EQ(tree.GetNode(last_city)->parent, tree.TailParentNode());
    EXPECT_EQ(tree.TailParentNode()->next, tree.HeadParentNode());
    EXPECT_EQ(tree.HeadParentNode()->prev, tree.TailParentNode());
}

TEST(Prev_Next_andBetween, 1)
{
    // no reversal yet in this test

    int cityNum = 10, origin = 1;
    std::vector<int> order = {3, 6, 8, 4, 1, 2, 5, 9, 10, 7};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);

    // 1. prev and next
    for (int i : {4, 8, 2, 5}) {
        auto city = order[i], prev_city = order[i - 1],
             next_city = order[i + 1];
        EXPECT_EQ(tree.GetNext(tree.GetNode(city)), tree.GetNode(next_city));
        EXPECT_EQ(tree.GetPrev(tree.GetNode(city)), tree.GetNode(prev_city));
    }

    auto first_city_node = tree.GetNode(order.front());
    auto last_city_node = tree.GetNode(order.back());
    EXPECT_EQ(tree.GetNext(last_city_node), first_city_node);
    EXPECT_EQ(tree.GetPrev(first_city_node), last_city_node);

    // 2. between
    auto Between = [&tree](int a, int b, int c) {
        return tree.Between(tree.GetNode(a), tree.GetNode(b), tree.GetNode(c));
    };
    EXPECT_TRUE(Between(3, 6, 8));
    EXPECT_TRUE(Between(8, 4, 1));
    EXPECT_TRUE(Between(3, 8, 10));
    EXPECT_TRUE(Between(3, 5, 7));
    EXPECT_TRUE(Between(9, 7, 3));
    EXPECT_TRUE(Between(6, 1, 3));
    EXPECT_TRUE(Between(10, 7, 5));
    EXPECT_TRUE(Between(6, 8, 3));
    EXPECT_TRUE(Between(7, 3, 6));
    EXPECT_TRUE(Between(7, 3, 10));
    EXPECT_TRUE(Between(5, 10, 1));
    EXPECT_TRUE(Between(4, 1, 2));
    EXPECT_TRUE(Between(3, 1, 7));
    EXPECT_TRUE(Between(2, 10, 1));
    EXPECT_TRUE(Between(10, 4, 1));
    EXPECT_FALSE(Between(6, 4, 8));
    EXPECT_FALSE(Between(6, 4, 8));
    EXPECT_FALSE(Between(10, 3, 7));
    EXPECT_FALSE(Between(10, 1, 8));
    EXPECT_FALSE(Between(3, 7, 9));
    EXPECT_FALSE(Between(1, 4, 2));
    EXPECT_FALSE(Between(6, 3, 10));
}

TEST(Reverse_exactly_a_complete_segment, 1)
{
    int cityNum = 14, origin = 1;
    std::vector<int> order = {11, 13, 6, 8, 4, 1, 2, 5, 9, 10, 7, 12, 14, 3};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);

    EXPECT_EQ(tree.GetNode(11)->parent->id, 0);
    EXPECT_EQ(tree.GetNode(13)->parent->id, 0);
    EXPECT_EQ(tree.GetNode(6)->parent->id, 0);
    EXPECT_EQ(tree.GetNode(8)->parent->id, 1);
    EXPECT_EQ(tree.GetNode(14)->parent->id, 3);
    EXPECT_EQ(tree.GetNode(7)->parent->id, 3);
    EXPECT_EQ(tree.GetNode(3)->parent->id, 3);

    tree.Reverse(tree.GetNode(8), tree.GetNode(1));  //[8, 4, 1]
    EXPECT_TRUE(tree.GetParentNode(8)->reverse);
    EXPECT_TRUE(tree.GetParentNode(4)->reverse);
    EXPECT_TRUE(tree.GetParentNode(1)->reverse);
    EXPECT_EQ(tree.GetNext(tree.GetNode(6)), tree.GetNode(1));
    EXPECT_EQ(tree.GetNext(tree.GetNode(8)), tree.GetNode(2));
    EXPECT_EQ(tree.GetNext(tree.GetNode(4)), tree.GetNode(8));
    EXPECT_EQ(tree.GetNext(tree.GetNode(1)), tree.GetNode(4));
    EXPECT_TRUE(
        tree.Between(tree.GetNode(1), tree.GetNode(4), tree.GetNode(8)));
    // only change the reversal bit:
    // the beginning and ending nodes remain unchanged
    EXPECT_EQ(tree.GetParentNode(4)->begin, tree.GetNode(8));
    EXPECT_EQ(tree.GetParentNode(4)->end, tree.GetNode(1));
    std::vector<int> ans = {11, 13, 6, 1, 4, 8, 2, 5, 9, 10, 7, 12, 14, 3};
    EXPECT_EQ(GetTour(tree, 11), ans);
    auto node = tree.OriginCityNode();
    do {
        EXPECT_EQ(tree.GetNext(tree.GetPrev(node)), node);
        EXPECT_EQ(tree.GetPrev(tree.GetNext(node)), node);
        node = tree.GetNext(node);
    } while (node != tree.OriginCityNode());
    EXPECT_FALSE(
        tree.Between(tree.GetNode(6), tree.GetNode(13), tree.GetNode(1)));

    tree.Reverse(tree.GetNode(11), tree.GetNode(6));  //[11, 13, 6]
    EXPECT_TRUE(tree.GetParentNode(11)->reverse);
    EXPECT_EQ(tree.GetNext(tree.GetNode(11)), tree.GetNode(1));
    EXPECT_EQ(tree.GetPrev(tree.GetNode(13)), tree.GetNode(6));
    EXPECT_EQ(tree.GetNext(tree.GetNode(13)), tree.GetNode(11));
    EXPECT_EQ(tree.GetPrev(tree.GetNode(6)), tree.GetNode(3));
    EXPECT_TRUE(
        tree.Between(tree.GetNode(6), tree.GetNode(13), tree.GetNode(1)));
    ans = {6, 13, 11, 1, 4, 8, 2, 5, 9, 10, 7, 12, 14, 3};
    EXPECT_EQ(GetTour(tree, 6), ans);
    node = tree.OriginCityNode();
    do {
        EXPECT_EQ(tree.GetNext(tree.GetPrev(node)), node);
        EXPECT_EQ(tree.GetPrev(tree.GetNext(node)), node);
        node = tree.GetPrev(node);
    } while (node != tree.OriginCityNode());

    tree.Reverse(tree.GetNode(10), tree.GetNode(3));  // [10, 7, 12, 14, 3]
    EXPECT_EQ(tree.GetPrev(tree.GetNode(3)), tree.GetNode(9));
    EXPECT_EQ(tree.GetPrev(tree.GetNode(10)), tree.GetNode(7));
    EXPECT_EQ(tree.GetNext(tree.GetNode(14)), tree.GetNode(12));
    ans = {6, 13, 11, 1, 4, 8, 2, 5, 9, 3, 14, 12, 7, 10};
    EXPECT_EQ(GetTour(tree, 6), ans);
    node = tree.OriginCityNode();
    do {
        EXPECT_EQ(tree.GetNext(tree.GetPrev(node)), node);
        EXPECT_EQ(tree.GetPrev(tree.GetNext(node)), node);
        node = tree.GetPrev(node);
    } while (node != tree.OriginCityNode());

    tree.Reverse(tree.GetNode(6), tree.GetNode(11));  // [6, 13, 11]
    EXPECT_EQ(tree.GetParentNode(11)->reverse, false);
    EXPECT_EQ(tree.GetPrev(tree.GetNode(11)), tree.GetNode(10));
    EXPECT_EQ(tree.GetPrev(tree.GetNode(13)), tree.GetNode(11));
    EXPECT_EQ(tree.GetNext(tree.GetNode(6)), tree.GetNode(1));
    ans = {11, 13, 6, 1, 4, 8, 2, 5, 9, 3, 14, 12, 7, 10};
    EXPECT_EQ(GetTour(tree, 11), ans);
    node = tree.OriginCityNode();
    do {
        EXPECT_EQ(tree.GetNext(tree.GetPrev(node)), node);
        EXPECT_EQ(tree.GetPrev(tree.GetNext(node)), node);
        node = tree.GetPrev(node);
    } while (node != tree.OriginCityNode());
    EXPECT_EQ(GetTourViaParent(tree, 1),
              GetTour(tree, tree.GetParentNode(1)->ForwardBeginNode()->city));
}

TEST(Reverse_a_partial_segment_with_no_split_and_merge, 1)
{
    int cityNum = 23, origin = 1;
    std::vector<int> order = {11, 13, 6,  8,  4,  1,  2,  5,  9,  10, 7, 12,
                              14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);
    // assert segment
    std::vector<int> expected_segment_sizes{4, 4, 4, 4, 7};
    EXPECT_EQ(ActualSegmentSizes(tree), expected_segment_sizes);
    // partial reverse. The nominal segment length is 4,
    // thus if the partial segment has a
    // length <= 3, it is reversed in the segment with not split-merge.
    tree.Reverse(tree.GetNode(4), tree.GetNode(2));  // reverse [4, 1, 2]
    EXPECT_TRUE(IsOrderedNeighbor(tree, 8, 2));
    EXPECT_TRUE(IsOrderedNeighbor(tree, 4, 5));
    EXPECT_TRUE(IsOrderedNeighbor(tree, 2, 1));
    EXPECT_TRUE(IsOrderedNeighbor(tree, 1, 4));
    EXPECT_FALSE(IsOrderedNeighbor(tree, 4, 1));
    EXPECT_EQ(tree.GetParentNode(1)->begin, tree.GetNode(2));
    EXPECT_EQ(tree.GetParentNode(1)->end, tree.GetNode(5));
    std::vector<int> expected_tour = {11, 13, 6,  8,  2,  1,  4,  5,
                                      9,  10, 7,  12, 14, 3,  15, 16,
                                      17, 18, 20, 19, 23, 22, 21};
    EXPECT_EQ(GetTour(tree, 11), expected_tour);
    auto node = tree.GetParentNode(1)->begin;
    auto end = tree.GetParentNode(1)->end;
    while (node != end) {
        EXPECT_EQ(node->next->id - node->id, 1);
        node = node->next;
    }

    tree.Reverse(tree.GetNode(20), tree.GetNode(23));  // reverse [20, 19, 23]
    EXPECT_TRUE(IsOrderedNeighbor(tree, 20, 22));
    EXPECT_TRUE(IsOrderedNeighbor(tree, 23, 19));
    EXPECT_TRUE(IsOrderedNeighbor(tree, 18, 23));
    EXPECT_EQ(tree.GetParentNode(17)->begin, tree.GetNode(17));
    EXPECT_EQ(tree.GetParentNode(20)->end, tree.GetNode(21));
    // check IDs
    node = tree.GetParentNode(17)->begin;
    end = tree.GetParentNode(17)->end;
    while (node != end) {
        EXPECT_EQ(node->next->id - node->id, 1);
        node = node->next;
    }
    expected_tour = {11, 13, 6,  8,  2,  1,  4,  5,  9,  10, 7, 12,
                     14, 3,  15, 16, 17, 18, 23, 19, 20, 22, 21};
    EXPECT_EQ(GetTour(tree, 11), expected_tour);

    // let's check first reverse a whole segment
    // reverse [17, 18, 23, 19, 20, 22, 21]
    tree.Reverse(tree.GetNode(17), tree.GetNode(21));
    expected_tour = {11, 13, 6,  8,  2,  1,  4,  5,  9,  10, 7, 12,
                     14, 3,  15, 16, 21, 22, 20, 19, 23, 18, 17};
    EXPECT_EQ(GetTour(tree, 11), expected_tour);
    EXPECT_TRUE(tree.GetParentNode(23)->reverse);
    tree.Reverse(tree.GetNode(23), tree.GetNode(17));  // reverse [23, 18, 17]
    expected_tour = {11, 13, 6,  8,  2,  1,  4,  5,  9,  10, 7, 12,
                     14, 3,  15, 16, 21, 22, 20, 19, 17, 18, 23};
    EXPECT_EQ(GetTour(tree, 11), expected_tour);
    EXPECT_TRUE(IsOrderedNeighbor(tree, 17, 18));
    EXPECT_TRUE(IsOrderedNeighbor(tree, 23, 11));
    EXPECT_TRUE(IsOrderedNeighbor(tree, 19, 17));
    EXPECT_EQ(tree.GetNext(tree.GetNode(21)), tree.GetNode(22));
    EXPECT_EQ(tree.GetPrev(tree.GetNode(11)), tree.GetNode(23));
    EXPECT_TRUE(Between(tree, 11, 22, 23));
    EXPECT_TRUE(Between(tree, 18, 23, 1));
    EXPECT_TRUE(Between(tree, 5, 7, 3));
    EXPECT_FALSE(Between(tree, 15, 18, 22));
    // check IDs
    node = tree.GetParentNode(22)->begin;
    end = tree.GetParentNode(22)->end;
    while (node != end) {
        EXPECT_EQ(node->next->id - node->id, 1);
        node = node->next;
    }
    EXPECT_EQ(GetTourViaParent(tree, 1),
              GetTour(tree, tree.GetParentNode(1)->ForwardBeginNode()->city));
}

TEST(Split_and_merge, 1)
{
    std::vector<int> v;
    int cityNum = 23, origin = 1;
    std::vector<int> order = {11, 13, 6,  8,  4,  1,  2,  5,  9,  10, 7, 12,
                              14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);
    // assert segment
    std::vector<int> expected_segment_sizes{4, 4, 4, 4, 7};
    EXPECT_EQ(ActualSegmentSizes(tree), expected_segment_sizes);

    tree.SplitAndMerge(tree.GetNode(6), true, TwoLevelTree::Direction::forward);
    EXPECT_EQ(tree.GetParentNode(6), tree.GetParentNode(4));
    v = std::vector<int>{2, 6, 4, 4, 7};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    v = std::vector<int>{6,  8,  4,  1,  2,  5,  9,  10, 7,  12, 14, 3,
                         15, 16, 17, 18, 20, 19, 23, 22, 21, 11, 13};
    EXPECT_TRUE(GetTour(tree, 6) == v);
    v = std::vector<int>{11, 13, 6,  8,  4,  1,  2,  5,  9,  10, 7, 12,
                         14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    EXPECT_TRUE(GetTour(tree, 11) == v);

    // reverse the 2nd segment [6, 8, 4, 1, 2, 5]
    tree.Reverse(tree.GetNode(6), tree.GetNode(5));
    v = std::vector<int>{10, 9,  6,  8,  4,  1,  2,  5, 13, 11, 21, 22,
                         23, 19, 20, 18, 17, 16, 15, 3, 14, 12, 7};
    EXPECT_TRUE(GetTour(tree, 10, TwoLevelTree::Direction::backward) == v);
    v = std::vector<int>{11, 13, 5,  2,  1,  4,  8,  6,  9,  10, 7, 12,
                         14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    EXPECT_TRUE(GetTour(tree, 11) == v);
    EXPECT_TRUE(tree.GetParentNode(4)->reverse);
    // split and merge will not change the tour
    tree.SplitAndMerge(tree.GetNode(4), true, TwoLevelTree::Direction::forward);
    v = std::vector<int>{1,  4,  8,  6,  9,  10, 7,  12, 14, 3, 15, 16,
                         17, 18, 20, 19, 23, 22, 21, 11, 13, 5, 2};
    EXPECT_TRUE(GetTour(tree, 1) == v);
    v = std::vector<int>{2, 3, 7, 4, 7};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    v = std::vector<int>{11, 13, 5,  2,  1,  4,  8,  6,  9,  10, 7, 12,
                         14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    EXPECT_TRUE(GetTour(tree, 11) == v);
    EXPECT_TRUE(tree.GetParentNode(2)->reverse);
    EXPECT_FALSE(tree.GetParentNode(4)->reverse);

    // try backward merge
    //doesn't include 19
    tree.SplitAndMerge(tree.GetNode(19), false,
                       TwoLevelTree::Direction::backward);
    v = std::vector<int>{11, 13, 5,  2,  1,  4,  8,  6,  9,  10, 7, 12,
                         14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    EXPECT_TRUE(GetTour(tree, 11) == v);
    v = std::vector<int>{2, 3, 7, 7, 4};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    EXPECT_EQ(tree.GetParentNode(19)->begin, tree.GetNode(19));
    EXPECT_EQ(tree.GetParentNode(16)->end, tree.GetNode(20));
    EXPECT_TRUE(tree.GetParentNode(2)->reverse);

    // try another backward merge,
    // note that the segment for [5, 2, 1] has the reverse bit
    tree.SplitAndMerge(tree.GetNode(10), true,
                       TwoLevelTree::Direction::backward);
    v = std::vector<int>{11, 13, 5,  2,  1,  4,  8,  6,  9,  10, 7, 12,
                         14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    EXPECT_TRUE(GetTour(tree, 11) == v);
    v = std::vector<int>{2, 8, 2, 7, 4};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    EXPECT_TRUE(tree.GetParentNode(9)->reverse);
    EXPECT_EQ(tree.GetParentNode(9)->end, tree.GetNode(5));
    EXPECT_EQ(tree.GetParentNode(7)->end, tree.GetNode(12));
    EXPECT_EQ(tree.GetParentNode(9)->begin, tree.GetNode(10));
    EXPECT_EQ(tree.GetParentNode(12)->begin, tree.GetNode(7));

    // another one, here the segment containing 5 has the reverse bit
    tree.SplitAndMerge(tree.GetNode(2), true, TwoLevelTree::Direction::forward);
    v = std::vector<int>{11, 13, 5,  2,  1,  4,  8,  6,  9,  10, 7, 12,
                         14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    EXPECT_TRUE(GetTour(tree, 11) == v);
    v = std::vector<int>{2, 1, 9, 7, 4};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    EXPECT_TRUE(tree.GetParentNode(5)->reverse);
    EXPECT_FALSE(tree.GetParentNode(1)->reverse);
    EXPECT_EQ(tree.GetParentNode(5)->begin, tree.GetNode(5));
    EXPECT_EQ(tree.GetParentNode(5)->end, tree.GetNode(5));
    EXPECT_EQ(tree.GetParentNode(12)->begin, tree.GetNode(2));
    EXPECT_EQ(tree.GetParentNode(2)->end, tree.GetNode(12));
    v = std::vector<int>{2,  5, 13, 11, 21, 22, 23, 19, 20, 18, 17, 16,
                         15, 3, 14, 12, 7,  10, 9,  6,  8,  4,  1};
    EXPECT_TRUE(GetTour(tree, 2, TwoLevelTree::Direction::backward) == v);
}

TEST(Reverse_a_partial_segment_with_split_and_merge, 1)
{
    std::vector<int> v;
    int cityNum = 23, origin = 1;
    std::vector<int> order = {11, 13, 6,  8,  4,  1,  2,  5,  9,  10, 7, 12,
                              14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);
    // assert segment
    std::vector<int> expected_segment_sizes{4, 4, 4, 4, 7};
    EXPECT_EQ(ActualSegmentSizes(tree), expected_segment_sizes);

    // the nominal length is 4, if a part to be reversed is > 3,
    // then split and merge is used.
    tree.Reverse(tree.GetNode(18), tree.GetNode(23));
    v = std::vector<int>{22, 21, 11, 13, 6,  8,  4,  1,  2,  5,  9, 10,
                         7,  12, 14, 3,  15, 16, 17, 23, 19, 20, 18};
    EXPECT_TRUE(GetTour(tree, 22) == v);
    v = std::vector<int>{6, 4, 4, 5, 4};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    EXPECT_TRUE(tree.GetParentNode(18)->reverse);
    EXPECT_FALSE(tree.GetParentNode(22)->reverse);

    // how about to reverse [11, 13, 6, 8],
    // note that no forward merging is actually needed
    tree.Reverse(tree.GetNode(11), tree.GetNode(8));
    v = std::vector<int>{8,  6, 13, 11, 4,  1,  2,  5,  9,  10, 7, 12,
                         14, 3, 15, 16, 17, 23, 19, 20, 18, 22, 21};
    EXPECT_TRUE(GetTour(tree, 8) == v);
    EXPECT_TRUE(tree.GetParentNode(22)->reverse);
    EXPECT_EQ(tree.GetParentNode(21)->begin, tree.GetNode(21));
    EXPECT_EQ(tree.GetParentNode(21)->end, tree.GetNode(23));
    EXPECT_TRUE(tree.GetParentNode(8)->reverse);
    v = std::vector<int>{12, 7,  10, 9,  5,  2,  1,  4,  11, 13, 6, 8,
                         21, 22, 18, 20, 19, 23, 17, 16, 15, 3,  14};
    EXPECT_TRUE(GetTour(tree, 12, TwoLevelTree::Direction::backward) == v);

    // reverse [19, 20, 18, 22], whose reverse bit is set
    tree.Reverse(tree.GetNode(19), tree.GetNode(22));
    v = std::vector<int>{21, 8,  6, 13, 11, 4,  1,  2,  5,  9,  10, 7,
                         12, 14, 3, 15, 16, 17, 23, 22, 18, 20, 19};
    EXPECT_TRUE(GetTour(tree, 21) == v);
    v = std::vector<int>{5, 4, 4, 6, 4};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    EXPECT_FALSE(tree.GetParentNode(19)->reverse);
    EXPECT_EQ(tree.GetParentNode(19)->begin, tree.GetNode(22));
    EXPECT_TRUE(tree.GetParentNode(21)->reverse);
    EXPECT_EQ(tree.GetParentNode(21)->begin, tree.GetNode(11));
    EXPECT_EQ(GetTourViaParent(tree, 1),
              GetTour(tree, tree.GetParentNode(1)->ForwardBeginNode()->city));
}

TEST(Reverse_multiple_segments_with_split_and_merge, 1)
{
    std::vector<int> v;
    int cityNum = 23, origin = 1;
    std::vector<int> order = {11, 13, 6,  8,  4,  1,  2,  5,  9,  10, 7, 12,
                              14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);
    // assert segment
    std::vector<int> expected_segment_sizes{4, 4, 4, 4, 7};
    EXPECT_EQ(ActualSegmentSizes(tree), expected_segment_sizes);

    // though a and b are not in the same segment,
    // after split-and-merge, they are.
    tree.Reverse(tree.GetNode(6), tree.GetNode(4));
    v = std::vector<int>{11, 13, 4,  8,  6,  1,  2,  5,  9,  10, 7, 12,
                         14, 3,  15, 16, 17, 18, 20, 19, 23, 22, 21};
    EXPECT_TRUE(GetTour(tree, 11) == v);
    v = std::vector<int>{2, 6, 4, 4, 7};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);
    EXPECT_EQ(tree.GetParentNode(4)->begin, tree.GetNode(4));

    tree.Reverse(tree.GetNode(22), tree.GetNode(8));
    v = std::vector<int>{8, 4,  13, 11, 21, 22, 6,  1,  2,  5,  9, 10,
                         7, 12, 14, 3,  15, 16, 17, 18, 20, 19, 23};
    EXPECT_TRUE(GetTour(tree, 8) == v);
    EXPECT_TRUE(tree.GetParentNode(8)->reverse);
    EXPECT_TRUE(tree.GetParentNode(22)->reverse);
    EXPECT_EQ(tree.GetParentNode(8)->end, tree.GetNode(8));
    EXPECT_FALSE(tree.GetParentNode(23)->reverse);
    v = std::vector<int>{6, 4, 4, 4, 5};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);

    // reverse multiple segments
    // now [8, 4] are deprived
    tree.Reverse(tree.GetNode(13), tree.GetNode(5));
    v = std::vector<int>{5,  2, 1,  6,  22, 21, 11, 13, 9,  10, 7, 12,
                         14, 3, 15, 16, 17, 18, 20, 19, 23, 8,  4};
    EXPECT_TRUE(GetTour(tree, 5) == v);
    EXPECT_FALSE(tree.GetParentNode(22)->reverse);
    EXPECT_TRUE(tree.GetParentNode(2)->reverse);
    v = std::vector<int>{4, 4, 4, 4, 7};
    EXPECT_TRUE(ActualSegmentSizes(tree) == v);

    tree.Reverse(tree.GetNode(6), tree.GetNode(14));
    v = std::vector<int>{13, 11, 21, 22, 6, 3, 15, 16, 17, 18, 20, 19,
                         23, 8,  4,  5,  2, 1, 14, 12, 7,  10, 9};
    EXPECT_TRUE(GetTour(tree, 13) == v);
    v = std::vector<int>{5, 3, 7, 3, 5};
    EXPECT_TRUE(ActualSegmentSizes(tree, 13) == v);

    // let's traverse via the parents
    auto p_parent = tree.HeadParentNode();
    do {
        EXPECT_EQ(p_parent->next->prev, p_parent);
        EXPECT_EQ(p_parent->prev->next, p_parent);
        EXPECT_EQ((p_parent->id + 1) % tree.ParentNum(), p_parent->next->id);
        p_parent = p_parent->next;
    } while (p_parent != tree.HeadParentNode());

    std::vector<int> segment_size_res, node_res;
    p_parent = tree.GetParentNode(13);
    do {
        segment_size_res.push_back(p_parent->size);
        auto p = p_parent->ForwardBeginNode();
        while (p != p_parent->ForwardEndNode()) {
            node_res.push_back(p->city);
            p = tree.GetNext(p);
        }
        node_res.push_back(p->city);
        p_parent = p_parent->next;
    } while (p_parent != tree.GetParentNode(13));
    v = std::vector<int>{5, 3, 7, 3, 5};
    EXPECT_TRUE(segment_size_res == v);
    EXPECT_EQ(node_res, GetTour(tree, 13));

    EXPECT_EQ(GetTourViaParent(tree, 1),
              GetTour(tree, tree.GetParentNode(1)->ForwardBeginNode()->city));
}

TEST(flip, 1)
{
    std::vector<int> v;
    int cityNum = 12, origin = 1;
    std::vector<int> order = {3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);

    tree.Flip(tree.GetNode(3), tree.GetNode(6), tree.GetNode(10),
              tree.GetNode(7));
    v = std::vector<int>{6, 8, 4, 1, 12, 2, 5, 9, 10, 3, 11, 7};
    EXPECT_TRUE(GetTour(tree, 6) == v);
    EXPECT_EQ(GetTourViaParent(tree, 1),
              GetTour(tree, tree.GetParentNode(1)->ForwardBeginNode()->city));

    tree.Reverse(tree.GetNode(4), tree.GetNode(10));
    v = std::vector<int>{6, 8, 10, 9, 5, 2, 12, 1, 4, 3, 11, 7};
    EXPECT_TRUE(GetTour(tree, 6) == v);
    EXPECT_EQ(GetTourViaParent(tree, 1),
              GetTour(tree, tree.GetParentNode(1)->ForwardBeginNode()->city));

    tree.Flip(tree.GetNode(8), tree.GetNode(10), tree.GetNode(7),
              tree.GetNode(6));
    v = std::vector<int>{10, 9, 5, 2, 12, 1, 4, 3, 11, 7, 8, 6};
    EXPECT_TRUE(GetTour(tree, 10) == v);
    EXPECT_EQ(GetTourViaParent(tree, 1),
              GetTour(tree, tree.GetParentNode(1)->ForwardBeginNode()->city));

    auto node = tree.OriginCityNode();
    do {
        EXPECT_EQ(tree.GetNext(tree.GetPrev(node)), node);
        EXPECT_EQ(tree.GetPrev(tree.GetNext(node)), node);
        node = tree.GetPrev(node);
    } while (node != tree.OriginCityNode());

    // check ID of the 2 segment
    node = tree.GetParentNode(2)->begin;
    auto end = tree.GetParentNode(2)->end;
    while (node != end) {
        EXPECT_EQ(node->next->id - node->id, 1);
        node = node->next;
    }

    // backward
    v = std::vector<int>{10, 9, 5, 2, 12, 1, 4, 3, 11, 7, 8, 6};
    EXPECT_TRUE(GetTour(tree, 10) == v);
    tree.Flip(1, 12, 9, 10);
    v = std::vector<int>{1, 9, 5, 2, 12, 10, 6, 8, 7, 11, 3, 4};
    EXPECT_TRUE(GetTour(tree, 1) == v);

    tree.Flip(10, 6, 8, 7);
    v = std::vector<int>{10, 8, 6, 7, 11, 3, 4, 1, 9, 5, 2, 12};
    EXPECT_TRUE(GetTour(tree, 10) == v);
}

TEST(opt_move_and_undo, 1)
{
    std::vector<int> v, v2;
    int cityNum = 12, origin = 1;
    std::vector<int> order = {3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);

    Move2Opt(tree, 5, 9, 3, 11);
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 11, 7, 10, 9};
    EXPECT_TRUE(GetTour(tree, 3) == v);
    // undo
    Undo2OptMove(tree, 5, 9, 3, 11);
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    EXPECT_TRUE(GetTour(tree, 3) == v);
    v = std::vector<int>{12, 2, 5, 9, 10, 7, 11, 3, 6, 8, 4, 1};
    EXPECT_TRUE(GetTour(tree, 12) == v);

    Move2Opt(tree, 12, 2, 7, 10);
    v = std::vector<int>{3, 11, 7, 2, 5, 9, 10, 12, 1, 4, 8, 6};
    v2 = std::vector<int>{3, 6, 8, 4, 1, 12, 10, 9, 5, 2, 7, 11};
    EXPECT_TRUE((GetTour(tree, 3) == v || GetTour(tree, 3) == v2));
    // undo
    Undo2OptMove(tree, 12, 2, 7, 10);
    v = std::vector<int>{12, 2, 5, 9, 10, 7, 11, 3, 6, 8, 4, 1};
    EXPECT_TRUE(GetTour(tree, 12) == v);
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    EXPECT_TRUE(GetTour(tree, 3) == v);
}

TEST(Tree_deep_copy_move_and_independency, 1)
{
    std::vector<int> v;
    int cityNum = 12, origin = 1;
    std::vector<int> order = {3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);

    auto tree2 = tree;
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    EXPECT_TRUE(GetTour(tree2, 3) == v);

    Move2Opt(tree, 5, 9, 3, 11);
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 11, 7, 10, 9};
    EXPECT_TRUE(GetTour(tree, 3) == v);
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    EXPECT_TRUE(GetTour(tree2, 3) == v);
    v = std::vector<int>{5, 9, 10, 7, 11, 3, 6, 8, 4, 1, 12, 2};
    EXPECT_TRUE(GetTour(tree2, 5) == v);

    TwoLevelTree tree3 = std::move(tree);
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 11, 7, 10, 9};
    EXPECT_TRUE(GetTour(tree3, 3) == v);
    Undo2OptMove(tree3, 5, 9, 3, 11);
    v = std::vector<int>{3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    EXPECT_TRUE(GetTour(tree3, 3) == v);
}

TEST(Double_bridge_move, 1)
{
    std::vector<int> v;
    int cityNum = 12, origin = 1;
    std::vector<int> order = {3, 6, 8, 4, 1, 12, 2, 5, 9, 10, 7, 11};
    TwoLevelTree tree{cityNum, origin};
    tree.SetTour(order.begin(), order.begin() + order.size() - 1);

    tree.DoubleBridgeMove(12, 5, 11, 8);
    v = std::vector<int>{2, 5, 4, 1, 12, 3, 6, 8, 9, 10, 7, 11};
    EXPECT_TRUE(GetTour(tree, 2) == v);
    // let's traverse via the parents
    auto p_parent = tree.HeadParentNode();
    int id = 0;
    do {
        EXPECT_EQ(p_parent->id, id);
        EXPECT_EQ(p_parent->next->prev, p_parent);
        EXPECT_EQ(p_parent->prev->next, p_parent);
        EXPECT_EQ((p_parent->id + 1) % tree.ParentNum(), p_parent->next->id);
        p_parent = p_parent->next;
        id++;
    } while (p_parent != tree.HeadParentNode());

    tree.DoubleBridgeMove(3, 9, 2, 4);
    v = std::vector<int>{2, 6, 8, 9, 1, 12, 3, 5, 4, 10, 7, 11};
    EXPECT_TRUE(GetTour(tree, 2) == v);
    auto node = tree.OriginCityNode();
    do {
        EXPECT_EQ(tree.GetNext(tree.GetPrev(node)), node);
        EXPECT_EQ(tree.GetPrev(tree.GetNext(node)), node);
        node = tree.GetNext(node);
    } while (node != tree.OriginCityNode());

    p_parent = tree.HeadParentNode();
    id = 0;
    do {
        EXPECT_EQ(p_parent->id, id);
        EXPECT_EQ(p_parent->next->prev, p_parent);
        EXPECT_EQ(p_parent->prev->next, p_parent);
        EXPECT_EQ((p_parent->id + 1) % tree.ParentNum(), p_parent->next->id);
        p_parent = p_parent->next;
        id++;
    } while (p_parent != tree.HeadParentNode());

    tree.DoubleBridgeMove(5, 11, 6, 1);
    v = std::vector<int>{4, 10, 7, 11, 12, 3, 5, 8, 9, 1, 2, 6};
    EXPECT_TRUE(GetTour(tree, 4) == v);
    p_parent = tree.HeadParentNode();
    int size = 0;
    do {
        EXPECT_EQ(p_parent->next->prev, p_parent);
        EXPECT_EQ(p_parent->prev->next, p_parent);
        EXPECT_EQ((p_parent->id + 1) % tree.ParentNum(), p_parent->next->id);
        if (!p_parent->reverse) {
            if (p_parent->next->reverse) {
                EXPECT_EQ(tree.GetNext(p_parent->end), p_parent->next->end);
            } else {
                EXPECT_EQ(tree.GetNext(p_parent->end), p_parent->next->begin);
            }
        }
        size += p_parent->size;
        p_parent = p_parent->next;
    } while (p_parent != tree.HeadParentNode());
    assert(size == 12);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
