#pragma once

#include <iostream> // std::ostream, std::ifstream
#include <vector>
#include <string>

#include "Problem.hpp"

namespace knapsack_solver {

class Solution {
    public:
    int max_value;
    std::vector<bool> selected;
    Solution();
    Solution(int InstanceSize);

    // WARNING - destructive method
    // method modifies this solution and remainingSpace vector
    void AddItem(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace);
    
    // WARNING - destructive method
    // method modifies this solution and remainingSpace vector
    void RemoveItem(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace);
    
    // WARNING - destructive method
    // method modifies this solution and remainingSpace vector
    /// @returns true if addition was succesfull, false if item did not fit
    bool AddItemIfFits(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace);
};

class PackagedSolution {
    public:
        std::string algorithm;
        Solution solution;
        double to_optimum_ratio;
        double solve_time;
        std::vector<int> remainingSpaces;
    
    void ExportJSON(const std::string file_name);
};

class KnapsackSolver{
    public:
    static bool Fits(const std::vector<int> & weights, const std::vector<int> & remainingSpace);
    static std::vector<int> CalculateRemainingSpaces(const Solution & solution, const Problem & problem);
};

class BruteForceSolver{
    public:
    struct Options{
        enum class SearchOrder { ZERO_FIRST, ONE_FIRST, RANDOM };
        SearchOrder search_order = SearchOrder::ZERO_FIRST;

        bool iterative = false;
    };

    static PackagedSolution Solve(const Problem & problem, const Requirements & requirements, const Options options);
};

class BranchAndBoundSolver{
    private:
    static int GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode, Solution currentSolution, std::vector<int> remainingSpace);

    static Solution DFSLateFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, std::vector<int> remainingSpace);
    static Solution DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, std::vector<int> remainingSpace, const int & lower_bound);

    public:
    struct Options{
        enum class BoundingFunction { NONE, CONTINOUS, ACYCLIC, BASE_DYNAMIC };
        BoundingFunction bounding_function = BoundingFunction::NONE;
        bool late_fit = false;
    };

    static Solution BnBLateFitPath(const Problem & problem);
    static Solution BnBEarlyFitPath(const Problem & problem);

    static PackagedSolution Solve(const Problem & problem, const Requirements & requirements, const Options & options);
};

class GreedySolver{
    public:
    struct Options{
        Problem::SortMode sort_mode = Problem::SortMode::WEIGHT_VALUE_RATIO;
        int buffor = 1;
    };

    static PackagedSolution Solve(const Problem & problem, const Requirements & requirements, const Options & options);
    static Solution GreedyUniversal(const Problem & problem, const Requirements & requirements, const Options & options);
    static Solution GreedyIgnoreConnections(const Problem & problem, const Options & options);
    static Solution GreedyPath(const Problem & problem, const Options & options);
};

class FloydSolver{
    public:
    static PackagedSolution Solve(const Problem & problem, const Requirements & requirements);
    static Solution Connected(const Problem & problem);
};


} // namespace knapsack_solver

std::ostream& operator<<(std::ostream& os, const knapsack_solver::Solution& s);

std::ostream& operator<<(std::ostream& os, const knapsack_solver::PackagedSolution& ps);