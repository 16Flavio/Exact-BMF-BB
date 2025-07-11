#include "BMFInstance.hpp"
#include <fstream>
#include <stdexcept>

using namespace std;

void BMFInstance::loadFromFile(const string &filename, int r){

    X_.loadFromFile(filename);

    cout << "Specified Matrix :\n";
    X_.print();
    cout << "\n";
    cout << "Specified rank : " << r << "\n";

    if(r <= 0 || r > min(X_.numRows(), X_.numCols())){
        throw invalid_argument("Rank must be between 1 and min(rows, cols)");
    }
    r_ = r;
}