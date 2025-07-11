#include "BnBNode.hpp"
#include "gurobi_c++.h"
#include <vector>
#include <stdexcept>
#include <math.h>

using namespace std;

Solution BnBNode::computeHalfHFixed() const{
    const Matrix &X = instance_->getMatrix();
    int m = X.numRows();
    int n = X.numCols();
    int r = instance_->getRank();

    for (int i = 0; i < m; i++) {
        for (int k = 0; k < r; k++) {
            if (isFixedW(i, k)) {
                Whalf[i][k].set(GRB_DoubleAttr_LB, getW(i, k));
                Whalf[i][k].set(GRB_DoubleAttr_UB, getW(i, k));
            } else {
                Whalf[i][k].set(GRB_DoubleAttr_LB, 0.0);
                Whalf[i][k].set(GRB_DoubleAttr_UB, 1.0);
            }
        }
    }

    sharedModel_halfleafW->update();
    sharedModel_halfleafW->optimize();

    Solution sol;
    sol.cost = sharedModel_halfleafW->get(GRB_DoubleAttr_ObjVal);
    sol.W.resize(m * r);
    sol.H.resize(r * n);
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < r; ++k)
            sol.W[i*r + k] = (Whalf[i][k].get(GRB_DoubleAttr_X) == 1.0);
    for (int k = 0; k < r; ++k)
        for (int j = 0; j < n; ++j)
            sol.H[k*n + j] = (getH(k,j) == 1);
    
    return sol;

}
Solution BnBNode::computeHalfWFixed() const{
    const Matrix &X = instance_->getMatrix();
    int m = X.numRows();
    int n = X.numCols();
    int r = instance_->getRank();

    for (int k = 0; k < r; k++) {
        for (int j = 0; j < n; ++j) {
            if (isFixedH(k, j)) {
                Hhalf[k][j].set(GRB_DoubleAttr_LB, getH(k, j));
                Hhalf[k][j].set(GRB_DoubleAttr_UB, getW(k, j));
            } else {
                Hhalf[k][j].set(GRB_DoubleAttr_LB, 0.0);
                Hhalf[k][j].set(GRB_DoubleAttr_UB, 1.0);
            }
        }
    }

    sharedModel_halfleafH->update();
    sharedModel_halfleafH->optimize();

    Solution sol;
    sol.cost = sharedModel_halfleafH->get(GRB_DoubleAttr_ObjVal);
    sol.W.resize(m * r);
    sol.H.resize(r * n);
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < r; ++k)
            sol.W[i*r + k] = (getW(i,k) == 1);
    for (int k = 0; k < r; ++k)
        for (int j = 0; j < n; ++j)
            sol.H[k*n + j] = (Hhalf[k][j].get(GRB_DoubleAttr_X) == 1.0);
    
    return sol;
}


void BnBNode::initializeGurobiModel() const {
    //Initialize model for LB
    const Matrix &X = instance_->getMatrix();
    int m = X.numRows();
    int n = X.numCols();
    int r = instance_->getRank();

    if (sharedEnv == nullptr) {
        sharedEnv = new GRBEnv();
        sharedEnv->set(GRB_IntParam_OutputFlag, 0);
        sharedModel = new GRBModel(*sharedEnv);
    }

    sharedEnv->set(GRB_IntParam_OutputFlag, 0);
    sharedEnv->set(GRB_IntParam_LogToConsole, 0);

    if (timeLimit_ > 0.0) {
        double depthFactor = exp(-depth_*0.5);
        sharedModel->set(GRB_DoubleParam_TimeLimit, timeLimit_*depthFactor);
    }

    sharedModel->set(GRB_IntParam_Presolve, 0); 

    sharedModel->set(GRB_IntParam_Cuts, 0);         
    sharedModel->set(GRB_IntParam_MIPFocus, 1);    

    sharedModel->set(GRB_DoubleParam_Heuristics, 0.0); 

    sharedEnv->start();

    Wc.resize(m, vector<GRBVar>(r));
    Hc.resize(r, vector<GRBVar>(n));
    Z.resize(m, vector<GRBVar>(n));
    D.resize(m, vector<GRBVar>(n));

    vector<vector<vector<GRBVar>>> T(m, vector<vector<GRBVar>>(n, vector<GRBVar>(r)));

    for (int i = 0; i < m; i++) {
        for (int k = 0; k < r; k++) {
            Wc[i][k] = sharedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        }
    }

    for (int k = 0; k < r; k++) {
        for (int j = 0; j < n; j++) {
            Hc[k][j] = sharedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Z[i][j] = sharedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            D[i][j] = sharedModel->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < r; k++) {
                T[i][j][k] = sharedModel->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                sharedModel->addConstr(T[i][j][k] <= Wc[i][k]);
                sharedModel->addConstr(T[i][j][k] <= Hc[k][j]);
                sharedModel->addConstr(T[i][j][k] >= Wc[i][k] + Hc[k][j] - 1.0);
            }
            GRBLinExpr sumT = 0;
            for (int k = 0; k < r; k++) {
                sharedModel->addConstr(Z[i][j] >= T[i][j][k]);
                sumT += T[i][j][k];
            }
            sharedModel->addConstr(Z[i][j] <= sumT);

            double x_val = X.get(i,j) ? 1.0 : 0.0;
            sharedModel->addConstr(D[i][j] >= Z[i][j] - x_val);
            sharedModel->addConstr(D[i][j] >= x_val - Z[i][j]);
        }
    }

    GRBLinExpr obj = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            obj += D[i][j];
        }
    }
    sharedModel->setObjective(obj, GRB_MINIMIZE);

    sharedModel->update();
    modelInitialized = true;
    cached_m = m; cached_n = n; cached_r = r;
    //Initialize Model for HalfH
    sharedModel_halfleafH = new GRBModel(*sharedEnv);
    sharedModel_halfleafH->set(GRB_IntParam_OutputFlag, 0);

    Hhalf.resize(r, vector<GRBVar>(n));
    Zhalf.resize(m, vector<GRBVar>(n));
    Dhalf.resize(m, vector<GRBVar>(n));

    for (int k = 0; k < r; k++) {
        for (int j = 0; j < n; j++) {
            Hhalf[k][j] = sharedModel_halfleafH->addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }
    }

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n ; j++){
            Zhalf[i][j] = sharedModel_halfleafH->addVar(0.0, 1.0, 0.0, GRB_BINARY);
            GRBLinExpr z_expr = 0;
            for(int k = 0; k < r; k++){
                z_expr += getW(i,k)*Hhalf[k][j];
            }
            sharedModel_halfleafH->addConstr(Zhalf[i][j] == z_expr);
        }
    }

    for(int i = 0; i < m;i++){
        for(int j = 0; j < n;j++){
            Dhalf[i][j] = sharedModel_halfleafH->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            double xval = X.get(i, j) ? 1.0 : 0.0;
            sharedModel_halfleafH->addConstr(Dhalf[i][j] >= Zhalf[i][j] - xval);
            sharedModel_halfleafH->addConstr(Dhalf[i][j] >= xval - Zhalf[i][j]);
        }
    }

    obj = 0;
    int x_val;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            obj += Dhalf[i][j];
        }
    }
    sharedModel_halfleafH->setObjective(obj, GRB_MINIMIZE);
    sharedModel_halfleafH->update();

    //Initialize Model for HalfW
    sharedModel_halfleafW = new GRBModel(*sharedEnv);
    sharedModel_halfleafW->set(GRB_IntParam_OutputFlag, 0);

    Whalf.resize(m, vector<GRBVar>(r));
    Zhalf.resize(m, vector<GRBVar>(n));
    Dhalf.resize(m, vector<GRBVar>(n));

    for (int i = 0; i < m; i++) {
        for (int k = 0; k < r; k++) {
            Whalf[i][k] = sharedModel_halfleafW->addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }
    }

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n ; j++){
            Zhalf[i][j] = sharedModel_halfleafW->addVar(0.0, 1.0, 0.0, GRB_BINARY);
            GRBLinExpr z_expr = 0;
            for(int k = 0; k < r; k++){
                z_expr += Whalf[i][k] * getH(k,j);
            }
            sharedModel_halfleafW->addConstr(Zhalf[i][j] == z_expr);
        }
    }

    for(int i = 0; i < m;i++){
        for(int j = 0; j < n;j++){
            Dhalf[i][j] = sharedModel_halfleafW->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            double xval = X.get(i, j) ? 1.0 : 0.0;
            sharedModel_halfleafW->addConstr(Dhalf[i][j] >= Zhalf[i][j] - xval);
            sharedModel_halfleafW->addConstr(Dhalf[i][j] >= xval - Zhalf[i][j]);
        }
    }

    obj = 0;
    x_val = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            obj += Dhalf[i][j];
        }
    }
    sharedModel_halfleafW->setObjective(obj, GRB_MINIMIZE);
    sharedModel_halfleafW->update();
}

double BnBNode::computeGurobiLB() const {
    const Matrix &X = instance_->getMatrix();
    int m = X.numRows();
    int n = X.numCols();
    int r = instance_->getRank();

    try {

        for (int i = 0; i < m; i++) {
            for (int k = 0; k < r; k++) {
                if (isFixedW(i, k)) {
                    Wc[i][k].set(GRB_DoubleAttr_LB, getW(i, k));
                    Wc[i][k].set(GRB_DoubleAttr_UB, getW(i, k));
                } else {
                    Wc[i][k].set(GRB_DoubleAttr_LB, 0.0);
                    Wc[i][k].set(GRB_DoubleAttr_UB, 1.0);
                }
            }
        }

        for (int k = 0; k < r; k++) {
            for (int j = 0; j < n; j++) {
                if (isFixedH(k, j)) {
                    Hc[k][j].set(GRB_DoubleAttr_LB, getH(k, j));
                    Hc[k][j].set(GRB_DoubleAttr_UB, getH(k, j));
                } else {
                    Hc[k][j].set(GRB_DoubleAttr_LB, 0.0);
                    Hc[k][j].set(GRB_DoubleAttr_UB, 1.0);
                }
            }
        }

        if (!lastW_.empty()) {
            for (int i = 0; i < m_; i++)
                for (int k = 0; k < r_; k++)
                    Wc[i][k].set(GRB_DoubleAttr_Start, lastW_[i][k]);
            for (int k = 0; k < r_; k++)
                for (int j = 0; j < n_; j++)
                    Hc[k][j].set(GRB_DoubleAttr_Start, lastH_[k][j]);
        }

        sharedModel->update();
        try{
            sharedModel->optimize();
        }catch(GRBException& e){
            const_cast<BnBNode*>(this)->resetStartValues();
            sharedModel->optimize();
        }
        
        const_cast<BnBNode*>(this)->storeLastSolution();

        double objVal = GRB_INFINITY;
        int status = sharedModel->get(GRB_IntAttr_Status);
        
        if (status == GRB_OPTIMAL) {
            objVal = sharedModel->get(GRB_DoubleAttr_ObjVal);
        } 
        else if (status == GRB_SUBOPTIMAL) {
            try {
                objVal = sharedModel->get(GRB_DoubleAttr_ObjVal);
            } 
            catch (...) {
                try {
                    objVal = sharedModel->get(GRB_DoubleAttr_ObjBound);
                } 
                catch (...) {
                    objVal = GRB_INFINITY;
                }
            }
        }
        else if (status == GRB_INFEASIBLE || status == GRB_UNBOUNDED) {
            objVal = GRB_INFINITY;
        }
        else {
            try {
                objVal = sharedModel->get(GRB_DoubleAttr_ObjBound);
            } 
            catch (...) {
                objVal = GRB_INFINITY;
            }
        }

        for (int i = 0; i < m; i++) {
            for (int k = 0; k < r; k++) {
                Wc[i][k].set(GRB_DoubleAttr_LB, 0.0);
                Wc[i][k].set(GRB_DoubleAttr_UB, 1.0);
            }
        }

        for (int k = 0; k < r; k++) {
            for (int j = 0; j < n; j++) {
                Hc[k][j].set(GRB_DoubleAttr_LB, 0.0);
                Hc[k][j].set(GRB_DoubleAttr_UB, 1.0);
            }
        }

        sharedModel->update();
        return objVal;

    } catch (GRBException &e) {
        if (e.getErrorCode() == 10003) {
            return GRB_INFINITY;
        }
        cerr << "Erreur Gurobi (" << e.getErrorCode() << "): " << e.getMessage() << "\n";
        return GRB_INFINITY;
    } catch (...) {
        cerr << "Exception inconnue dans computeGurobiLB\n";
        return GRB_INFINITY;
    }
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
    vector<int> sums(r_, 0);
    vector<bool> complete(r_, true);
    
    for (int k = 0; k < r_; k++) {
        for (int i = 0; i < m_; i++) {
            if (isFixedW(i, k)) {
                sums[k] += getW(i, k);
            } else {
                complete[k] = false;
            }
        }
    }
    
    for (int k = 0; k < r_ - 1; k++) {
        if (complete[k] && complete[k+1] && sums[k] < sums[k+1]) {
            return false;
        }
    }
    
    return true;
}

void BnBNode::storeLastSolution(){
    lastW_.resize(m_, vector<double>(r_));
    lastH_.resize(r_, vector<double>(n_));
    for (int i = 0; i < m_; i++)
        for (int k = 0; k < r_; k++)
            lastW_[i][k] = Wc[i][k].get(GRB_DoubleAttr_X);
    for (int k = 0; k < r_; k++)
        for (int j = 0; j < n_; j++)
            lastH_[k][j] = Hc[k][j].get(GRB_DoubleAttr_X);
}

void BnBNode::resetStartValues() {
    for (int i = 0; i < m_; i++) {
        for (int k = 0; k < r_; k++) {
            if (!isFixedW(i, k)) {  
                Wc[i][k].set(GRB_DoubleAttr_Start, GRB_UNDEFINED);
            }
        }
    }

    for (int k = 0; k < r_; k++) {
        for (int j = 0; j < n_; j++) {
            if (!isFixedH(k, j)) {  
                Hc[k][j].set(GRB_DoubleAttr_Start, GRB_UNDEFINED);
            }
        }
    }

    lastW_.clear();
    lastH_.clear();
}