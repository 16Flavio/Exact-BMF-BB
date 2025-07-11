#include "BnBSolver.hpp"
#include <cassert>

using namespace std;

BnBSolver::BnBSolver(const BMFInstance &instance, const Solution &initialSol)
    : instance_(instance), bestSolutions_(), bestCost_(initialSol.cost), elapsedMs_(0.0) {
    //bestSolutions_.push_back(initialSol);
}

void BnBSolver::run() {

    cout << "Start!\n";

    int m = instance_.getMatrix().numRows();
    int n = instance_.getMatrix().numCols();
    int r = instance_.getRank();

    using Clock = chrono::high_resolution_clock;
    auto start = Clock::now();

    priority_queue<BnBNode, vector<BnBNode>, NodeCompare> pq;
    
    pq.push(BnBNode(instance_)); 

    auto push_if_valid = [&](const BnBNode& child) {
        if (!(child.lb() > bestCost_ || (child.cost() > bestCost_ && !child.isLeaf())) && child.checkSymmetry()){  
            pq.push(child);
        }/*else{
            cout << "[#Explored: " << exploredNodes
                    << "] Depth: " << child.depth()
                    << ", LB: " << child.lb()
                    << ", Cost: " << child.cost()
                    << ", BestCost: " << bestCost_
                    << ", QueueSize: " << pq.size()
                    << "\n";
            cout << "lb is too high for child\n";
        }*/
    };

    vector<pair<uint64_t, uint64_t>> possibilities = generate_binary_pairs(instance_.getRank());

    cout << "Finished to generate binary pairs!\n";

    while (!pq.empty()) {

        BnBNode node = pq.top(); pq.pop();

        /*int lb = node.lb();
        if (lb > bestCost_ || (node.cost() >= bestCost_  && !node.isLeaf())){
            //cout << "lb is too high\n";
            continue;
        }*/
        if(node.lb() > bestCost_ || (node.cost() > bestCost_ && !node.isLeaf())){
            //cout << "lb is too high for node\n";
            continue;
        }

        if (node.isLeaf()) {

            int cost = node.cost();
            if (cost < bestCost_) {
                bestCost_ = cost;
                bestSolutions_.clear();
            }
            if (cost == bestCost_) {
                Solution sol;
                sol.cost = cost;
                sol.W.resize(m * r);
                sol.H.resize(r * n);
                for (int i = 0; i < m; ++i)
                    for (int k = 0; k < r; ++k)
                        sol.W[i*r + k] = (node.getW(i,k) == 1);
                for (int k = 0; k < r; ++k)
                    for (int j = 0; j < n; ++j)
                        sol.H[k*n + j] = (node.getH(k,j) == 1);
                bestSolutions_.push_back(sol);
                cout << "Leaf found with cost of " << cost << "\n";
                cout << "Actual best cost:" << bestCost_ << "\n";
            }
            continue;
        }else if (node.isHalfLeaf()) {
            if (node.remainingVariablesW() == 0) {
                /*for (int j = 0; j < n; j++) {
                    if (!node.isFixedHColumn(j)) {
                        for (const auto& possibility : possibilities) {
                            uint64_t b = possibility.second;
                            vector<pair<int, int>> assignments;
                            for (int k = 0; k < r; k++) {
                                if (!node.isFixedH(k, j)) {
                                    int idxH = m * r + k * n + j;
                                    int bit = (b >> k) & 1;
                                    assignments.push_back({idxH, bit});
                                }
                            }
                            if (!assignments.empty()) {
                                BnBNode child = node.branchMultiple(assignments);
                                push_if_valid(child);
                            }
                        }
                        break;
                    }
                }*/
                Solution sol = node.computeHalfWFixed();
                int cost = sol.cost;
                if (cost < bestCost_) {
                    bestCost_ = cost;
                    bestSolutions_.clear();
                }
                if (cost == bestCost_) {
                    bestSolutions_.push_back(sol);
                    cout << "Leaf found with cost of " << cost << "\n";
                    cout << "Actual best cost:" << bestCost_ << "\n";
                }
            }else if(node.remainingVariablesH() == 0){
                /*for (int i = 0; i < m; i++) {
                    if (!node.isFixedWRow(i)) {
                        for (const auto& possibility : possibilities) {
                            uint64_t a = possibility.first;
                            vector<pair<int, int>> assignments;
                            for (int k = 0; k < r; k++) {
                                if (!node.isFixedW(i, k)) {
                                    int idxW = i * r + k;
                                    int bit = (a >> k) & 1;
                                    assignments.push_back({idxW, bit});
                                }
                            }
                            if (!assignments.empty()) {
                                BnBNode child = node.branchMultiple(assignments);
                                push_if_valid(child);
                            }
                        }
                        break;
                    }
                }*/
                Solution sol = node.computeHalfHFixed();
                int cost = sol.cost;
                if (cost < bestCost_) {
                    bestCost_ = cost;
                    bestSolutions_.clear();
                }
                if (cost == bestCost_) {
                    bestSolutions_.push_back(sol);
                    cout << "Leaf found with cost of " << cost << "\n";
                    cout << "Actual best cost:" << bestCost_ << "\n";
                }
            }
        }
        else {
            int i = node.depth();
            int j = node.depth();
            
            for (const auto& possibility : possibilities) {
                uint64_t a = possibility.first;
                uint64_t b = possibility.second;
                vector<pair<int, int>> assignments;
                
                for (int k = 0; k < r; k++) {
                    if (!node.isFixedW(i, k)) {
                        int idxW = i * r + k;
                        int bit = (a >> k) & 1;
                        assignments.push_back({idxW, bit});
                    }
                }
                
                for (int k = 0; k < r; k++) {
                    if (!node.isFixedH(k, j)) {
                        int idxH = m * r + k * n + j;
                        int bit = (b >> k) & 1;
                        assignments.push_back({idxH, bit});
                    }
                }
                
                if (!assignments.empty()) {
                    BnBNode child = node.branchMultiple(assignments);
                    cout << "[#Explored: " << exploredNodes_
                    << "] Depth: " << child.depth()
                    << ", LB: " << child.lb()
                    << ", Cost: " << child.cost()
                    << ", BestCost: " << bestCost_
                    << ", QueueSize: " << pq.size()
                    << "\n";
                    child.printWH();
                    push_if_valid(child);
                }
            }
        }

        if (exploredNodes_ % 10 == 0) {
            cout << "[#Explored: " << exploredNodes_
                    << "] Depth: " << node.depth()
                    << ", LB: " << node.lb()
                    << ", Cost: " << node.cost()
                    << ", BestCost: " << bestCost_
                    << ", QueueSize: " << pq.size()
                    << "\n";
            //node.printWH();
        }

        exploredNodes_ ++;

    }

    auto end = Clock::now();
    elapsedMs_ = chrono::duration<double, milli>(end - start).count();
}

int BnBSolver::selectMostPromisingVariable(const BnBNode &node) const {
    int m = instance_.getMatrix().numRows();
    int n = instance_.getMatrix().numCols();
    int r = instance_.getRank();
    const Matrix &X = instance_.getMatrix();

    int bestScore = -1;
    int bestIdx = -1;

    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < r; ++k) {
            if (node.isFixedW(i, k)) continue;  
            int score = 0;
            for (int j = 0; j < n; ++j) {
                if (X.get(i, j)) {
                    if (!node.isFixedH(k, j) || node.getH(k, j) == 1) {
                        score++;
                    }
                }
            }

            int idx = i * r + k;  
            if (score > bestScore) {
                bestScore = score;
                bestIdx = idx;
            }
        }
    }

    for (int k = 0; k < r; ++k) {
        for (int j = 0; j < n; ++j) {
            if (node.isFixedH(k, j)) continue;  
            int score = 0;
            for (int i = 0; i < m; ++i) {
                if (X.get(i, j)) {
                    if (!node.isFixedW(i, k) || node.getW(i, k) == 1) {
                        score++;
                    }
                }
            }

            int idx = m * r + k * n + j; 
            if (score > bestScore) {
                bestScore = score;
                bestIdx = idx;
            }
        }
    }

    return bestIdx;
}

int BnBSolver::selectFirstFree(const BnBNode &node) const {
    int total = node.totalVariables();
    for (int idx = 0; idx < total; ++idx) {
        int m = instance_.getMatrix().numRows();
        int r = instance_.getRank();
        if (idx < m*r) {
            int i = idx / r, k = idx % r;
            if (!node.isFixedW(i,k)) return idx;
        } else {
            int off = idx - m*r;
            int kk = off / instance_.getMatrix().numCols();
            int j  = off % instance_.getMatrix().numCols();
            if (!node.isFixedH(kk,j)) return idx;
        }
    }
    return -1; // leaf
}

vector<pair<uint64_t, uint64_t>> BnBSolver::generate_binary_pairs(const int &r) const {
    vector<pair<uint64_t, uint64_t>> pairs;
    
    if (r < 0 || r > 64) return pairs;  
    
    const uint64_t max_val = (r == 64) ? ~0ULL : (1ULL << r) - 1;
    
    for (uint64_t a = 0; ; a++) {
        for (uint64_t b = 0; ; b++) {
            pairs.emplace_back(a, b);
            if (b == max_val) break;
        }
        if (a == max_val) break;
    }

    return pairs;
}

const vector<Solution>& BnBSolver::getBestSolutions() const {
    return bestSolutions_;
}

double BnBSolver::getElapsedTime() const {
    return elapsedMs_;
}