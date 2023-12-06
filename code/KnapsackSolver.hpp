#pragma once

#include <iostream> // std::ostream, std::ifstream
#include <vector>
#include <string>

#include "Problem.hpp"

namespace knapsack_solver {



class Solution {
    bool IsPathDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & length) const;
    bool IsCycleDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & start, const int & length) const;

    public:
    int max_value;
    std::vector<bool> selected;
    std::vector<int> remainingSpace;
    bool valid;

    Solution();
    Solution(const int & InstanceSize, const std::vector<int> & available_space);

    void AddItem(const Problem & problem, const int & selected_item_id);
    void AddItemForce(const Problem & problem, const int & selected_item_id);
    void RemoveItem(const Problem & problem, const int & selected_item_id);
    bool Fits(const Problem & problem, const int & selected_item_id);
    
    /// @returns true if addition was succesfull, false if item did not fit or was already in solution
    bool AddItemIfFits(const Problem & problem, const int & selected_item_id);

    bool IsPath(const Problem & problem) const;
    bool IsCycle(const Problem & problem) const;
    bool IsTree(const Problem & problem) const;
    bool IsConnected(const Problem & problem) const;
    bool IsCyclePossible(const Problem & problem) const;
};



class Validation{
    public:
    struct ValidationStatus{
        bool undergone = false;
        bool valid = true;
        bool value = true;
        bool remaining_space = true;
        bool fit = true;
        bool structure = true;
        bool self_valid = true;
        //bool quality;
     };
    static std::vector<int> CalculateRemainingSpaces(const Solution & solution, const Problem & problem);
    static int CalculateMaxValue(const Solution & solution, const Problem & problem);
    static ValidationStatus Validate(const Solution & solution, const PackagedProblem & problem);
    static int GoalFunction(const Solution & solution, const PackagedProblem & problem);
};



class PackagedSolution {
    public:
        std::string algorithm;
        Solution solution;
        double to_optimum_ratio;
        double solve_time;
        Validation::ValidationStatus validation_status;
    
    void ExportJSON(const std::string file_name);
};



class BruteForceSolver{
    public:
    struct Options{
        enum class SearchOrder { ZERO_FIRST, ONE_FIRST, RANDOM, GRAY_CODE, UNCONSTRAINED };
        SearchOrder search_order = SearchOrder::UNCONSTRAINED;

        bool iterative = true;
        bool late_fit = true;
    };

    private:
    static Solution Max(const Solution & a, const Solution & b);
    static Solution SolutionFromNumber(int num, const Problem & problem);
    static Solution DFS(const Problem & problem, Solution currentSolution, const Options::SearchOrder & search_order, const int & depth);
    static Solution DFS(const Problem & problem, Solution currentSolution, const Options::SearchOrder & search_order, const bool & add, const int & depth);

    public:
    static PackagedSolution Solve(const PackagedProblem & problem, const Options options);
    static Solution Iterative(const PackagedProblem & problem, const Options::SearchOrder & search_order);
    static Solution Recursive(const PackagedProblem & problem, const Options & options);
};



class BranchAndBoundSolver{
    private:
    static int GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode, Solution currentSolution);

    static Solution DFSLateFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId);
    static Solution DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, const int & lower_bound);

    public:
    struct Options{
        enum class BoundingFunction { NONE, CONTINOUS, ACYCLIC, BASE_DYNAMIC };
        BoundingFunction bounding_function = BoundingFunction::NONE;
        bool late_fit = false;
    };

    static Solution BnBLateFitPath(const Problem & problem);
    static Solution BnBEarlyFitPath(const Problem & problem);

    static PackagedSolution Solve(const PackagedProblem & problem, const Options & options);
};



class GreedySolver{
    public:
    struct Options{
        Problem::SortMode sort_mode = Problem::SortMode::WEIGHT_VALUE_RATIO;
        int buffor = 1;
    };

    static PackagedSolution Solve(const PackagedProblem & problem, const Options & options);
    static Solution GreedyUniversal(const PackagedProblem & problem, const Options & options);
    static Solution GreedyIgnoreConnections(const Problem & problem, const Options & options);
    static Solution GreedyPath(const Problem & problem, const Options & options);
};



class FloydSolver{
    public:
    static PackagedSolution Solve(const PackagedProblem & problem);
    static Solution Connected(const Problem & problem);
};



} // namespace knapsack_solver



std::ostream& operator<<(std::ostream& os, const knapsack_solver::Solution& s);

std::ostream& operator<<(std::ostream& os, const knapsack_solver::PackagedSolution& ps);