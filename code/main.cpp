#include <iostream> // cout, endl
#include <vector>

#include <cstdlib> // std::srand()
#include <ctime> // std::time()

#include "Problem.hpp"
#include "KnapsackSolver.hpp"

using namespace knapsack_solver;
using std::vector;
using std::cout;
using std::endl;

int main() {
    std::srand(std::time(0) * 1000);
    cout << "Random value on [0, " << RAND_MAX << "]: " << std::rand() << '\n';

    Problem::GenerationSettings gs;
    gs.instance_size = 200;
    gs.sub_knapsacks = 1;
    gs.randomize_knapsack_sizes = false;
    gs.knapsack_size_limit_exclusive = 200;
    gs.fixed_knapsack_sizes.push_back(225);
    gs.value_limit_exclusive = 10;
    gs.weight_limit_exclusive = 10;
    gs.connection_density = 0.01;

    Problem::BatchGenerateProblemsJSON(gs, 1, "../tests", "batch1", "problem");

    Problem p("../tests/batch1/problems/problem0.json");
    //cout << p;

    Requirements rq;

    cout << Floyd::Solve(p, rq) << endl;
    
    cout << KnapsackSolver::Greedy(p, rq) << endl;

    rq.structureToFind = Requirements::StructureToFind::IGNORE_CONNECTIONS;
    cout << KnapsackSolver::Greedy(p, rq) << endl;
    rq.structureToFind = Requirements::StructureToFind::PATH;

    /*cout << KnapsackSolver::BranchAndBound(p, rq) << endl;
    PackagedSolution ps = KnapsackSolver::BruteForce(p, rq);
    cout << ps << endl;
    ps.ExportJSON("../tests/json/res.json");*/

    return 0;
}