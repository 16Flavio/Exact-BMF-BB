#include "BnBNode.hpp"
#include "gurobi_c++.h"
#include <vector>
#include <stdexcept>
#include <math.h>

using namespace std;

GRBEnv* BnBNode::sharedEnv = nullptr;
GRBModel* BnBNode::relaxedModel = nullptr;
bool BnBNode::modelInitialized = false;
vector<vector<GRBVar>> BnBNode::W_vars;
vector<vector<GRBVar>> BnBNode::H_vars;
vector<vector<GRBVar>> BnBNode::Z_vars;
vector<vector<GRBVar>> BnBNode::D_vars;
vector<vector<vector<GRBVar>>> BnBNode::T_vars;

void BnBNode::initializeRelaxedModel(const BMFInstance* instance){
    if (modelInitialized) return;
    
    const Matrix &X = instance->getMatrix();
    int m = X.numRows();
    int n = X.numCols();
    int r = instance->getRank();

    sharedEnv = new GRBEnv();
    sharedEnv->set(GRB_IntParam_OutputFlag, 0);
    relaxedModel = new GRBModel(*sharedEnv);
    
    relaxedModel->set(GRB_IntParam_Presolve, 0);
    //relaxedModel->set(GRB_IntParam_Presolve, 1);
    //relaxedModel->set(GRB_IntParam_Cuts, 0);
    relaxedModel->set(GRB_IntParam_Cuts, 1);
    relaxedModel->set(GRB_DoubleParam_Heuristics, 0);
    relaxedModel->set(GRB_IntParam_Threads, 1);
    relaxedModel->set(GRB_IntParam_MIPFocus, 1);

    W_vars.clear();
    H_vars.clear();
    Z_vars.clear();    
    D_vars.clear();
    T_vars.clear();

    W_vars.resize(m, vector<GRBVar>(r));
    H_vars.resize(r, vector<GRBVar>(n));
    Z_vars.resize(m, vector<GRBVar>(n));
    D_vars.resize(m, vector<GRBVar>(n));
    T_vars.resize(m, vector<vector<GRBVar>>(r, vector<GRBVar>(n)));
    
    for (int i = 0; i < m; i++) {
        for (int k = 0; k < r; k++) {
            W_vars[i][k] = relaxedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        }
    }
    
    for (int k = 0; k < r; k++) {
        for (int j = 0; j < n; j++) {
            H_vars[k][j] = relaxedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        }
    }
    
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            Z_vars[i][j] = relaxedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            D_vars[i][j] = relaxedModel->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            for (int k = 0; k < r; ++k) {
                    T_vars[i][k][j] = relaxedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                }
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < r; ++k) {
            for (int j = 0; j < n; ++j) {
                relaxedModel->addConstr(T_vars[i][k][j] <= W_vars[i][k]);
                relaxedModel->addConstr(T_vars[i][k][j] <= H_vars[k][j]);
                relaxedModel->addConstr(T_vars[i][k][j] >= W_vars[i][k] + H_vars[k][j] - 1);
            }
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            GRBLinExpr sumT = 0;
            for (int k = 0; k < r; ++k) {
                relaxedModel->addConstr(Z_vars[i][j] >= T_vars[i][k][j]);
                sumT += T_vars[i][k][j];
            }
            relaxedModel->addConstr(Z_vars[i][j] <= sumT);
        }
    }

    
    GRBLinExpr obj = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double x_val = X.get(i, j) ? 1.0 : 0.0;
        
            relaxedModel->addConstr(D_vars[i][j] >= (x_val - Z_vars[i][j]));
            relaxedModel->addConstr(D_vars[i][j] >= (Z_vars[i][j] - x_val));
            
            obj += D_vars[i][j];
        }
    }

    relaxedModel->setObjective(obj, GRB_MINIMIZE);
    
    relaxedModel->update();
    relaxedModel->set(GRB_IntParam_DualReductions, 0);
    modelInitialized = true;
}

void BnBNode::updateModelBounds() const {

    if(timeLimit_ > 0.0){
        double depthFactor = exp(-depth_*0.5);
        relaxedModel->set(GRB_DoubleParam_TimeLimit, max(timeLimit_*depthFactor,0.01));
    }

    for (int i = 0; i < m_; i++) {
        for (int k = 0; k < r_; k++) {
            GRBVar var = W_vars[i][k];
            if (isFixedW(i, k)) {
                var.set(GRB_DoubleAttr_LB, getW(i, k));
                var.set(GRB_DoubleAttr_UB, getW(i, k));
            } else {
                var.set(GRB_DoubleAttr_LB, 0.0);
                var.set(GRB_DoubleAttr_UB, 1.0);
            }
        }
    }

    for (int k = 0; k < r_; k++) {
        for (int j = 0; j < n_; j++) {
            GRBVar var = H_vars[k][j];
            if (isFixedH(k, j)) {
                var.set(GRB_DoubleAttr_LB, getH(k, j));
                var.set(GRB_DoubleAttr_UB, getH(k, j));
            } else {
                var.set(GRB_DoubleAttr_LB, 0.0);
                var.set(GRB_DoubleAttr_UB, 1.0);
            }
        }
    }

    relaxedModel->update();
}

double BnBNode::computeGurobiLB() const {
    updateModelBounds();

    relaxedModel->optimize();
    
    int status = relaxedModel->get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL) {
        return relaxedModel->get(GRB_DoubleAttr_ObjVal);
    }
    return numeric_limits<double>::max();
}

int BnBNode::computePartialCost() const {
    const Matrix& X = instance_->getMatrix();
    int error = 0;
    
    if (isLeaf()) {
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j < n_; ++j) {
                int x_val = X.get(i,j) ? 1 : 0;
                if (Xtilde_[i][j] != x_val) {
                    ++error;
                }
            }
        }
    } else {
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j < n_; ++j) {
                int v = Xtilde_[i][j];
                int x_val = X.get(i,j) ? 1 : 0;
                if (v != -1 && v != x_val) {
                    ++error;
                }
            }
        }
    }
    return error;
}

bool BnBNode::checkSymmetry() const {
    vector<pair<int, int>> column_ranges(r_);
    
    for (int k = 0; k < r_; k++) {
        int fixed_sum = 0;
        int free_count = 0;
        for (int i = 0; i < m_; i++) {
            if (isFixedW(i, k)) {
                fixed_sum += getW(i, k);
            } else {
                free_count++;
            }
        }
        column_ranges[k] = {fixed_sum, fixed_sum + free_count};
    }

    for (int k = 0; k < r_ - 1; k++) {
        if (column_ranges[k].second < column_ranges[k+1].first) {
            return false;  
        }
    }
    return true;
}