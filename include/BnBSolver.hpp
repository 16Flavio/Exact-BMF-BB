#pragma once

#include "BMFInstance.hpp"
#include "Solution.hpp"
#include "BnBNode.hpp"
#include <queue>
#include <chrono>
#include <vector>
#include <functional>
#include <gurobi_c++.h>

using namespace std;

struct NodeCompare {
    bool operator()(const BnBNode& a, const BnBNode& b) const {
        if (a.depth() != b.depth())
            return a.depth() < b.depth(); 
        return a.lb() > b.lb();
    }
};

class BnBSolver {
    const BMFInstance &instance_;
    vector<Solution> bestSolutions_;  
    int bestCost_, exploredNodes_ = 0;
    double elapsedMs_;
    int selectMostPromisingVariable(const BnBNode &node) const;
    int selectFirstFree(const BnBNode &node) const;
    vector<pair<uint64_t, uint64_t>> generate_binary_pairs(const int &r) const ;

    static GRBEnv sharedHalfLeafEnv_;
    static GRBModel sharedHalfLeafModel_;
    static bool halfLeafModelInitialized_;
    void initializeHalfLeafModel() const;
    Solution solveHalfLeaf_Gurobi(const BnBNode& node, bool solveW) const;
    Solution solveHalfLeafW_Gurobi(const BnBNode& node) const;
    Solution solveHalfLeafH_Gurobi(const BnBNode& node) const;

public:
    
    BnBSolver(const BMFInstance &instance, const Solution &initialSol);

    void run();

    const vector<Solution>& getBestSolutions() const;

    double getElapsedTime() const;

    int getExploredNodes() const { return exploredNodes_;}
};