#include <iostream>
#include <vector>
#include <random>
#include <time.h>

#include "Problem.hpp"
#include "KnapsackSolver.hpp"

using namespace knapsack_solver;
using std::vector;
using std::cout;
using std::endl;

int main() {
    srand(time(NULL) * 1000);
    int instance_size = 17;
    vector<int> ks;
    ks.push_back(25);
    Problem p(instance_size, 1, ks, 10, 10, 0.4);

    cout << p;

    Requirements rq;

    cout << KnapsackSolver::Greedy(p, rq) << endl;

    rq.structureToFind = Requirements::StructureToFind::IGNORE_CONNECTIONS;
    cout << KnapsackSolver::Greedy(p, rq) << endl;
    rq.structureToFind = Requirements::StructureToFind::PATH;

    cout << KnapsackSolver::BranchAndBound(p, rq) << endl;
    cout << KnapsackSolver::BruteForce(p, rq) << endl;

    return 0;
}