#pragma once

#include "BMFInstance.hpp"
#include "Solution.hpp"

Solution computeInitialHeuristic(const BMFInstance &instance, int threshold = 1);