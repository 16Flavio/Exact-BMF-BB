#pragma once

#include <string>
#include "Matrix.hpp"

using namespace std;

class BMFInstance{

    Matrix X_;
    int r_;

public:

    BMFInstance() : X_(0,0), r_(0) {}

    //Load binary matrix X and fixe the rank r
    void loadFromFile(const string &filename, int r);

    //get the matrix and the rank
    const Matrix& getMatrix() const {return X_;}
    int getRank() const {return r_;}

};