#include <iostream>
#include <string>

#include "BMFInstance.hpp"
#include "Heuristic.hpp"
#include "BnBSolver.hpp"

using namespace std;

int main(int argc, char** argv) {
    if (argc < 5 || argc > 6) {
        cerr << "Usage: " << argv[0] << " <matrix_file> <rank> <threshold> <timelimitLB> [initial_upper_bound]\n";
        return 1;
    }

    const string filename = argv[1];
    const int rank = stoi(argv[2]);
    const int threshold = stoi(argv[3]);
    const double timelimit = stod(argv[4]);
    
    int initialUB = -1;
    if (argc == 6) {
        initialUB = stoi(argv[5]);
    }

    BMFInstance instance;
    try {
        instance.loadFromFile(filename, rank);
    } catch (const exception &e) {
        cerr << "Loading error: " << e.what() << "\n";
        return 1;
    }

    BnBNode::setTimeLimit(timelimit);

    Solution initialSol;
    if(initialUB > 0){
        initialSol.cost = initialUB;
        cout << "Upper bound defined: " << initialUB << "\n";
    }else {
        initialSol = computeInitialHeuristic(instance, threshold);
        cout << "Heuristic - cost: " << initialSol.cost << "\n";
    } 

    BnBSolver solver(instance, initialSol);
    solver.run();

    const auto &bestSols = solver.getBestSolutions();
    cout << "Time to solve (ms): " << solver.getElapsedTime() << "\n";
    cout << "Best cost find: " << bestSols.front().cost << "\n";
    cout << "Node explored: " << solver.getExploredNodes() << "\n";

    for (size_t idx = 0; idx < bestSols.size(); ++idx) {
        cout << "--- Solution #" << idx+1 << " ---\n";
        
        int m = instance.getMatrix().numRows();
        int r = instance.getRank();
        int n = instance.getMatrix().numCols();
        
        cout << "W (" << m << "x" << r << "):\n";
        for (int i = 0; i < m; ++i) {
            for (int k = 0; k < r; ++k)
                cout << (bestSols[idx].W[i*r + k] ? '1' : '0') << ' ';
            cout << '\n';
        }
        
        cout << "H (" << r << "x" << n << "):\n";
        for (int k = 0; k < r; ++k) {
            for (int j = 0; j < n; ++j)
                cout << (bestSols[idx].H[k*n + j] ? '1' : '0') << ' ';
            cout << '\n';
        }
    }

    return 0;
}
