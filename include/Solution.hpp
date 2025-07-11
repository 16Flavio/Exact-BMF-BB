#pragma once

#include <vector>
#include <climits>
#include <ostream>
#include <limits>

using namespace std;

struct Solution {
    int cost;                
    vector<bool> W;             
    vector<bool> H;             

    Solution() : cost(INT_MAX) {}

    Solution(int c, const vector<bool> &w, const vector<bool> &h)
        : cost(c), W(w), H(h) {}
};

ostream& operator<<(ostream &out, const Solution &sol);