#pragma once

#include <vector>
#include <limits>
#include <cassert>
#include <string>
#include <unordered_map>
#include <gurobi_c++.h>
#include "BMFInstance.hpp"
#include "Solution.hpp"

using namespace std;

using AssignVec = vector<int8_t>;  

class BnBNode {
    const BMFInstance *instance_;
    int r_, m_, n_, depth_, variableFreeW_, variableFreeH_;
    double lb_;
    int cost_;
    //AssignVec assignments_;
    vector<uint64_t> W_values_, W_mask_;
    vector<uint64_t> H_values_, H_mask_;

    vector<vector<int>> Xtilde_;

    inline int totalW() const { return m_ * r_; }
    inline int totalH() const { return r_ * n_; }
    inline int wordCountW() const { return (totalW() + 63) / 64; }
    inline int wordCountH() const { return (totalH() + 63) / 64; }

    int computeLowerBound() const;
    double computeGurobiLB() const;
    int computePartialCost() const;

    static GRBEnv* sharedEnv;
    static GRBModel* relaxedModel;
    static bool modelInitialized;
    static vector<vector<GRBVar>> W_vars;
    static vector<vector<GRBVar>> H_vars;
    static vector<vector<GRBVar>> Z_vars;
    static vector<vector<GRBVar>> D_vars;
    static vector<vector<vector<GRBVar>>> T_vars;

    static void initializeRelaxedModel(const BMFInstance* instance);
    void updateModelBounds() const;

    static inline double timeLimit_ = 0.0;

public:

    explicit BnBNode(const BMFInstance &instance)
        : instance_(&instance),
          r_(instance.getRank()),
          m_(instance.getMatrix().numRows()),
          n_(instance.getMatrix().numCols()),
          depth_(0),
          cost_(0),
          variableFreeW_(m_*r_),
          variableFreeH_(n_*r_),
          Xtilde_(m_, vector<int>(n_, -1))
{
    W_values_.assign(wordCountW(), 0ULL);
    W_mask_.assign(wordCountW(), 0ULL);
    H_values_.assign(wordCountH(), 0ULL);
    H_mask_.assign(wordCountH(), 0ULL);
    
    if (!modelInitialized) {
        initializeRelaxedModel(instance_);
        modelInitialized = true;
    }

    lb_ = computeGurobiLB();
    //lb_ = computeLowerBound();
}

    BnBNode(const BnBNode &other, const vector<pair<int, int>>& assignments)
        : instance_(other.instance_),
          r_(other.r_), m_(other.m_), n_(other.n_),
          depth_(other.depth_),  
          variableFreeW_(other.variableFreeW_),
          variableFreeH_(other.variableFreeH_),
          cost_(other.cost_),
          W_values_(other.W_values_),
          W_mask_(other.W_mask_),
          H_values_(other.H_values_),
          H_mask_(other.H_mask_),
          Xtilde_(other.Xtilde_)
    {
        int countW = 0, countH = 0;
        bool possible = true;
        for (const auto& assignment : assignments) {
            int idx = assignment.first;
            int value = assignment.second;
            bool isW = (idx < totalW());

            int word, bit;
            uint64_t mask;

            
            if (isW) {
                word = idx / 64;
                bit = idx % 64;
                mask = 1ULL << bit;

                assert(word >= 0 && word < (int)W_mask_.size());
                assert(((W_mask_[word] >> bit) & 1) == 0 && "Variable W déjà fixée");

                if (value)
                    W_values_[word] |= mask;
                else
                    W_values_[word] &= ~mask;

                W_mask_[word] |= mask;
                countW++;

            } else {
                int idxH = idx - totalW();
                word = idxH / 64;
                bit = idxH % 64;
                mask = 1ULL << bit;

                assert(word >= 0 && word < (int)H_mask_.size());
                assert(((H_mask_[word] >> bit) & 1) == 0 && "Variable H déjà fixée");

                if (value)
                    H_values_[word] |= mask;
                else
                    H_values_[word] &= ~mask;

                H_mask_[word] |= mask;
                countH++;
            }
        }
        variableFreeW_ -= countW;
        variableFreeH_ -= countH;

        depth_++;

        if(other.isHalfLeaf()){
            if(other.remainingVariablesW() == 0){
                int j = depth_-1;
                for(int i = 0; i < m_; i++){
                    bool result = false;
                    for(int k = 0; k < r_; k++){
                        if(getW(i,k) == 1 && getH(k,j) == 1){
                            result = true;
                            break;
                        }
                    }
                    Xtilde_[i][j] = result ? 1 : 0;
                    int x_val = instance_->getMatrix().get(i,j) ? 1 : 0;
                    if(x_val != Xtilde_[i][j]){
                        cost_ ++;
                    }
                }
            }else if(other.remainingVariablesH() == 0){
                int i = depth_-1;
                for(int j = 0; j < n_; j++){
                    bool result = false;
                    for(int k = 0; k < r_; k++){
                        if(getW(i,k) == 1 && getH(k,j) == 1){
                            result = true;
                            break;
                        }
                    }
                    Xtilde_[i][j] = result ? 1 : 0;
                    int x_val = instance_->getMatrix().get(i,j) ? 1 : 0;
                    if(x_val != Xtilde_[i][j]){
                        cost_ ++;
                    }
                }
            }
        }else{
            int i = depth_-1;
            int j = depth_-1;
            int x_val;
            for(int l = 0; l < depth_; l++){

                bool result = false;
                for(int k = 0; k < r_; k++){
                    if(getW(i,k) == 1 && getH(k,l) == 1){
                        result = true;
                        break;
                    }
                }
                if(Xtilde_[i][l] == -1){
                    Xtilde_[i][l] = result ? 1 : 0;
                    x_val = instance_->getMatrix().get(i,l) ? 1 : 0;
                    if(x_val != Xtilde_[i][l]){
                        cost_ ++;
                    }
                }

                result = false;
                for(int k = 0; k < r_; k++){
                    if(getW(l,k) == 1 && getH(k,j) == 1){
                        result = true;
                        break;
                    }
                }
                if(Xtilde_[l][j] == -1){
                    Xtilde_[l][j] = result ? 1 : 0;
                    x_val = instance_->getMatrix().get(l,j) ? 1 : 0;
                    if(x_val != Xtilde_[l][j]){
                        cost_ ++;
                    }
                }
            }
        }
        
        lb_ = computeGurobiLB();
        //lb_ = computeLowerBound();
    }

    BnBNode branchMultiple(const vector<pair<int, int>>& assignments) const {
        return BnBNode(*this, assignments);
    }

    bool checkSymmetry() const;

    bool isFixedW(int i, int k) const {
        int idx = i * r_ + k;
        int word = idx / 64;
        int bit  = idx % 64;
        return (W_mask_[word] >> bit) & 1;
    }

    bool isFixedWRow(int i) const {
        int idx = i * r_;
        int word = idx / 64;
        int bit  = idx % 64;
        return (W_mask_[word] >> bit) & 1;
    }

    bool isFixedH(int k, int j) const {
        int idx = k * n_ + j;
        int word = idx / 64;
        int bit  = idx % 64;
        return (H_mask_[word] >> bit) & 1;
    }
    bool isFixedHColumn(int j) const{
        int idx = j;
        int word = idx / 64;
        int bit  = idx % 64;
        return (H_mask_[word] >> bit) & 1;
    }

    int getW(int i, int k) const {
        int idx = i * r_ + k;
        int word = idx / 64;
        int bit  = idx % 64;
        return (W_values_[word] >> bit) & 1;
    }

    int getH(int k, int j) const {
        int idx = k * n_ + j;
        int word = idx / 64;
        int bit  = idx % 64;
        return (H_values_[word] >> bit) & 1;
    }

    int depth() const { return depth_; }
    int totalVariables() const { return m_*r_ + r_*n_; }
    int remainingVariablesW() const {return variableFreeW_;}
    int remainingVariablesH() const {return variableFreeH_;}

    bool isLeaf() const {return (variableFreeH_ == 0 && variableFreeW_ == 0);}
    bool isHalfLeaf() const {return (variableFreeH_ == 0 || variableFreeW_ == 0);}

    double lb() const{return lb_;};

    int cost() const{return cost_;};

    void printWH() const {
        cout << "W matrix (" << m_ << " x " << r_ << "):\n";
        for (int i = 0; i < m_; ++i) {
            for (int k = 0; k < r_; ++k) {
                int idx = i * r_ + k;
                int word = idx / 64;
                int bit  = idx % 64;

                bool is_fixed = (W_mask_[word] >> bit) & 1;
                if (is_fixed) {
                    bool val = (W_values_[word] >> bit) & 1;
                    cout << val << " ";
                } else {
                    cout << ". ";
                }
            }
            cout << "\n";
        }

        cout << "H matrix (" << r_ << " x " << n_ << "):\n";
        for (int k = 0; k < r_; ++k) {
            for (int j = 0; j < n_; ++j) {
                int idx =  k * n_ + j;
                int word = idx / 64;
                int bit  = idx % 64;

                bool is_fixed = (H_mask_[word] >> bit) & 1;
                if (is_fixed) {
                    bool val = (H_values_[word] >> bit) & 1;
                    cout << val << " ";
                } else {
                    cout << ". ";
                }
            }
            cout << "\n";
        }

        cout << "X approximation matrix (" << m_ << " x " << n_ << "):\n";
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j < n_; ++j) {
                if (Xtilde_[i][j] != -1) {
                    cout << Xtilde_[i][j] << " ";
                } else {
                    cout << ". ";
                }
            }
            cout << "\n";
        }
    }

    static void setTimeLimit(double limit) { timeLimit_ = limit; }

    void storeLastSolution();
    void resetStartValues();

    Solution computeHalfHFixed() const;
    Solution computeHalfWFixed() const;

};