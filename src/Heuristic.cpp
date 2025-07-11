#include "Heuristic.hpp"
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

using vvi = vector<std::vector<int>>;

Solution computeInitialHeuristic(const BMFInstance &instance, int threshold) {
    const Matrix &Xmat = instance.getMatrix();
    int m = Xmat.numRows();
    int n = Xmat.numCols();
    int r = instance.getRank();

    vvi R(m, vector<int>(n));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            R[i][j] = Xmat.get(i, j) ? 1 : 0;

    vector<vector<int>> W2(m, vector<int>(r, 0));
    vector<vector<int>> H2(r, vector<int>(n, 0));

    // Précompute count of ones per row
    vector<int> count(m, 0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            if (R[i][j] == 1) count[i]++;

    for (int k = 0; k < r; ++k) {
        
        int maxi = numeric_limits<int>::min();
        int idx = -1;
        for (int i = 0; i < m; ++i) {
            if (count[i] > maxi) {
                maxi = count[i];
                idx = i;
            }
        }
        if (idx < 0) break;
        count[idx] = 0;
        W2[idx][k] = 1;

        for (int j = 0; j < n; ++j)
            if (R[idx][j] == 1)
                H2[k][j] = 1;

        for (int i = 0; i < m; ++i) {
            if (i == idx) continue;
            int gain = 0, cost = 0;
            for (int j = 0; j < n; ++j) {
                if (R[i][j] == 1 && H2[k][j] == 1) gain++;
                else if (R[i][j] == 0 && H2[k][j] == 1) cost++;
            }
            if (gain > cost && gain >= threshold)
                W2[i][k] = 1;
        }

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                if (W2[i][k] == 1 && H2[k][j] == 1)
                    R[i][j] = 0;

        for (int i = 0; i < m; ++i) {
            if (W2[i][k] == 1) {
                count[i] = 0;
                for (int j = 0; j < n; ++j)
                    if (R[i][j] == 1) count[i]++;
            }
        }
    }

    Solution sol;
    sol.W.resize(m * r);
    sol.H.resize(r * n);
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < r; ++k)
            sol.W[i * r + k] = (W2[i][k] == 1);
    for (int k = 0; k < r; ++k)
        for (int j = 0; j < n; ++j)
            sol.H[k * n + j] = (H2[k][j] == 1);

    int errors = 0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            bool cover = false;
            for (int k = 0; k < r; ++k) {
                if (sol.W[i * r + k] && sol.H[k * n + j]) { cover = true; break; }
            }
            if (cover != static_cast<bool>(Xmat.get(i, j))) errors++;
        }
    }
    sol.cost = errors;
    return sol;
}