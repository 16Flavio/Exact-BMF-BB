#include "Matrix.hpp"
#include <fstream>
#include <sstream>   
#include <stdexcept> 
#include <cassert>
#include <vector>

using namespace std;

Matrix::Matrix(int rows, int cols)
    : m_(rows), n_(cols), wordsPerRow_((cols + 63) / 64),
      data_(rows * wordsPerRow_, 0ull) {}

void Matrix::loadFromFile(const string &filename) {
    ifstream in(filename);
    if (!in) throw runtime_error("Cannot open file: " + filename);

    string line;
    vector<vector<int>> rawRows;

    while (getline(in, line)) {
        istringstream iss(line);
        vector<int> row;
        int bit;
        while (iss >> bit) {
            if (bit != 0 && bit != 1)
                throw runtime_error("Invalid entry in matrix (not 0/1)");
            row.push_back(bit);
        }
        if (!row.empty())
            rawRows.push_back(move(row));
    }

    if (rawRows.empty()) throw runtime_error("Empty matrix");

    m_ = rawRows.size();
    n_ = rawRows[0].size();

    for (const auto &row : rawRows) {
        if (row.size() != static_cast<size_t>(n_))
            throw runtime_error("Inconsistent row size in matrix");
    }

    wordsPerRow_ = (n_ + 63) / 64;
    data_.assign(m_ * wordsPerRow_, 0ULL);

    for (int i = 0; i < m_; ++i) {
        for (int j = 0; j < n_; ++j) {
            set(i, j, rawRows[i][j]);
        }
    }
}

void Matrix::set(int i, int j, bool value) {
    assert(i >= 0 && i < m_ && j >= 0 && j < n_);
    int wordIdx = i * wordsPerRow_ + (j / 64);
    uint64_t mask = uint64_t(1) << (j % 64);
    if (value)
        data_[wordIdx] |= mask;
    else
        data_[wordIdx] &= ~mask;
}

bool Matrix::get(int i, int j) const {
    assert(i >= 0 && i < m_ && j >= 0 && j < n_);
    int wordIdx = i * wordsPerRow_ + (j / 64);
    uint64_t mask = uint64_t(1) << (j % 64);
    return (data_[wordIdx] & mask) != 0;
}

const uint64_t* Matrix::getRowPtr(int i) const {
    assert(i >= 0 && i < m_);
    return &data_[i * wordsPerRow_];
}

void Matrix::print(ostream &out) const {
    for (int i = 0; i < m_; ++i) {
        for (int j = 0; j < n_; ++j) {
            out << (get(i, j) ? '1' : '0') << ' ';
            if (j < n_ - 1) out << ' ';
        }
        out << '\n';
    }
}

