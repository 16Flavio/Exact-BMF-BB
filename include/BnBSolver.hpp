#pragma once

#include "BMFInstance.hpp"
#include "Solution.hpp"
#include "BnBNode.hpp"
#include <queue>
#include <chrono>
#include <vector>
#include <functional>

using namespace std;

struct NodeCompare {
    bool operator()(const BnBNode& a, const BnBNode& b) const {
        if (a.depth() != b.depth())
            return a.depth() < b.depth(); 
        return a.cost() > b.cost();

        /*if (a.lb() != b.lb()){
            return a.lb() > b.lb();
        }
        return a.depth() < b.depth();*/
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

public:
    
    BnBSolver(const BMFInstance &instance, const Solution &initialSol);

    void run();

    const vector<Solution>& getBestSolutions() const;

    double getElapsedTime() const;

    int getExploredNodes() const { return exploredNodes_;}
    
};