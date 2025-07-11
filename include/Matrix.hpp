#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdint>

using namespace std;

class Matrix{
    
    int m_, n_, wordsPerRow_;
    vector<uint64_t> data_;

    /*Calculate index in data_
    inline int index(int row, int col) const {
        return row * wordsPerRow_ + (col >> 6);
    }
    
    inline uint64_t mask(int col) const{
        return uint64_t(1) << (col & 63);
    }*/

public:

    //Constructor : to get dimension
    Matrix(int rows, int cols);
    //Load Matrix from txt file
    void loadFromFile(const string &filename);

    //Writing and reading into the matrix
    void set(int i, int j, bool val);
    bool get(int i, int j) const;

    const uint64_t* getRowPtr(int i) const;

    //Print the matrix
    void print(ostream &out = cout) const;


    int numRows() const{return m_;}
    int numCols() const{return n_;}

};