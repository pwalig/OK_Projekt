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
    gs.instance_size = 18;
    gs.sub_knapsacks = 1;
    gs.randomize_knapsack_sizes = false;
    gs.knapsack_size_limit_exclusive = 25;
    gs.fixed_knapsack_sizes.push_back(30);
    gs.value_limit_exclusive = 10;
    gs.weight_limit_exclusive = 10;
    gs.connection_density = 0.1;

    Problem::Requirements rq;
    rq.structureToFind = Problem::Requirements::StructureToFind::PATH;

    PackagedProblem pp (gs, rq);

    cout << GreedySolver::Solve(pp, {Problem::SortMode::WEIGHT_VALUE_RATIO, 1}) << endl;
    cout << GreedySolver::Solve(pp, {Problem::SortMode::WEIGHT, 1}) << endl;
    cout << GreedySolver::Solve(pp, {Problem::SortMode::VALUE, 1}) << endl;

    BruteForceSolver::Options bop;
    bop.iterative = true;
    pp.requirements.structureToFind = Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS;
    bop.search_order = BruteForceSolver::Options::SearchOrder::ZERO_FIRST;
    cout << BruteForceSolver::Solve(pp, bop) << endl;
    bop.search_order = BruteForceSolver::Options::SearchOrder::ONE_FIRST;
    cout << BruteForceSolver::Solve(pp, bop) << endl;
    bop.search_order = BruteForceSolver::Options::SearchOrder::UNCONSTRAINED;
    cout << BruteForceSolver::Solve(pp, bop) << endl;

    return 0;
}