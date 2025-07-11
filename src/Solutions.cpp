#include "Solution.hpp"
#include <climits>
#include <iostream>

using namespace std;

ostream& operator<<(ostream &out, const Solution &sol) {
    out << "Cost: " << sol.cost << "\n";
    out << "W: ";
    for (bool b : sol.W) out << (b ? '1' : '0');
    out << "\nH: ";
    for (bool b : sol.H) out << (b ? '1' : '0');
    out << '\n';
    return out;
}