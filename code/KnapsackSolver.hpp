#pragma once

#include <iostream> // std::ostream
#include <fstream> // std::ifstream
#include <vector>
#include <string>
#include <cfloat> // DBL_MAX
#include <functional> // std::function

// Batch Solve requires these three:
#include <filesystem> // std::filesystem::create_directory, std::filesystem::remove_all
#include "json.hpp"
#include "FileNameDefines.hpp"

#include "Problem.hpp"

namespace knapsack_solver {



class Solution {
    bool IsPathDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & length) const;
    bool IsCycleDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & start, const int & length) const;
    /// @brief Should not be called by anything else than IsCyclePossible()
    bool IsCyclePossibleDFS(const Problem & problem, std::vector<int> & visited, std::vector<int> _remaining_space, const int & current, const int & start) const;

    public:
    int max_value;
    std::vector<bool> selected;
    std::vector<int> remainingSpace;
    bool valid;

    Solution();
    Solution(const int & InstanceSize, const std::vector<int> & available_space);

    enum class FaultTreatment { INVALIDATE, THROW, IGNORE };

    void AddItem(const Problem & problem, const int & selected_item_id, const FaultTreatment & fit_fault = FaultTreatment::THROW/*, const FaultTreatment & structure_fault = FaultTreatment::IGNORE*/);
    void RemoveItem(const Problem & problem, const int & selected_item_id);
    /// @brief Checks if selected item would fit in the knapsack together with all selected items. Relays on remainingSpace vector to accelerate calculation, does not check if remainingSpace was calculated correctly.
    /// @param problem needed only to check weights of selected item
    /// @param selected_item_id index in the problem of the item to check
    bool Fits(const Problem & problem, const int & selected_item_id) const;
    
    /// @returns true if addition was succesfull, false if item did not fit or was already in solution
    bool AddItemIfFits(const Problem & problem, const int & selected_item_id);

    /// @brief Checks if solution fits in the knapsack only based on selected vector. Ignores remainingSpace vector.
    /// @param problem problem needed to define knapsack sizes
    bool IsFit(const Problem & problem) const;

    bool IsPath(const Problem & problem) const;
    bool IsCycle(const Problem & problem) const;
    bool IsTree(const Problem & problem) const;
    bool IsConnected(const Problem & problem) const;
    bool IsStructure(const PackagedProblem & problem) const;

    bool IsValid(const PackagedProblem & problem) const;

    /// @brief Check if cycle that would fit in the knapsack is possible, assuming that items can only be added and not removed.
    /// @param problem problem needed to define connections and knapsack sizes
    /// @throws invalid_argument - if problem.items is different size than this->selected
    bool IsCyclePossible(const Problem & problem) const;
};



class Validation{
    public:
    Validation() = delete;
    struct ValidationStatus{
        bool undergone = false;
        bool valid = true;
        bool value = true;
        bool remaining_space = true;
        bool fit = true;
        bool structure = true;
        bool self_valid = true;
        bool quality = true; // needs to be filled outside of validation method - by the solve method
     };
    static std::vector<int> CalculateRemainingSpaces(const std::vector<bool> & selected, const Problem & problem);
    static int CalculateMaxValue(const std::vector<bool> & selected, const Problem & problem);
    static ValidationStatus Validate(const Solution & solution, const PackagedProblem & problem);
};



class PackagedSolution {
    public:
        std::string algorithm;
        Solution solution;
        double quality;
        double solve_time;
        Validation::ValidationStatus validation_status;
    
    void ExportJSON(const std::string file_name) const;
};



class KnapsackSolver {
    public:
    static int GoalFunction(const Solution & solution, const PackagedProblem & problem);
};



class DynamicSolver{
    public:
    DynamicSolver() = delete;
    struct Options{

        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };
    static Solution Dynamic(const Problem & problem);

    static PackagedSolution Solve(PackagedProblem & problem, const Options & options);
    static std::string GetAlgorithmName(const Options & options);
    static bool expect_perfection;
};



class GreedySolver{
    public:
    GreedySolver() = delete;
    struct Options{
        Problem::SortMode sort_mode = Problem::SortMode::VALUE_WEIGHT_RATIO;
        bool multi_run = true;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    static Solution GreedyUniversal(const PackagedProblem & problem, const Options & options);
    static Solution GreedyIgnoreConnections(const Problem & problem, const Options & options);
    static Solution GreedyPath(const Problem & problem, const Options & options);

    static PackagedSolution Solve(PackagedProblem & problem, const Options & options);
    static std::string GetAlgorithmName(const Options & options);
    static bool expect_perfection;
};



class BruteForceSolver{
    public:
    BruteForceSolver() = delete;
    struct Options{
        enum class SearchOrder { ZERO_FIRST, ONE_FIRST, RANDOM, GRAY_CODE, UNCONSTRAINED };
        SearchOrder search_order = SearchOrder::UNCONSTRAINED;

        bool iterative = true;
        bool late_fit = true;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    private:
    static Solution Max(const Solution & a, const Solution & b, const PackagedProblem & problem);
    static Solution SolutionFromNumber(int num, const Problem & problem);
    static Solution DFS(const PackagedProblem & problem, Solution currentSolution, const Options::SearchOrder & search_order, const int & depth);
    static Solution DFS(const PackagedProblem & problem, Solution currentSolution, const Options::SearchOrder & search_order, const bool & add, const int & depth);

    public:
    static Solution Iterative(const PackagedProblem & problem, const Options::SearchOrder & search_order);
    static Solution Recursive(const PackagedProblem & problem, const Options & options);

    static PackagedSolution Solve(PackagedProblem & problem, const Options options);
    static std::string GetAlgorithmName(const Options & options);
    static bool expect_perfection;
};



class BranchAndBoundSolver{
    private:
    BranchAndBoundSolver() = delete;

    static Solution DFSLateFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId);
    static Solution DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId);

    public:
    //static int GreedyIgnoreConnections(const Problem & problem, Solution currentSolution);
    
    struct Options{
        enum class BoundingFunction { NONE, CONTINOUS, BASE_DYNAMIC };
        BoundingFunction bounding_function = BoundingFunction::NONE;
        //std::function<int(const Problem &,  Solution )> lower_bound = GreedyIgnoreConnections;
        //std::function<int(const Problem &, Solution )> upper_bound = GreedyIgnoreConnections;
        bool late_fit = true;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    // Under Construction
    //static Solution DFS(const PackagedProblem & problem, Solution current_solution, const bool & add, const int & depth, int lower_bound, const Options & options);

    static Solution BnBLateFitPath(const Problem & problem);
    static Solution BnBEarlyFitPath(const Problem & problem);

    static PackagedSolution Solve(PackagedProblem & problem, const Options & options);
    static std::string GetAlgorithmName(const Options & options);
    static bool expect_perfection;
};



class GreedyHeuristicSearchSolver{
    public:
    GreedyHeuristicSearchSolver() = delete;
    struct Options{
        Problem::SortMode sort_mode = Problem::SortMode::VALUE_WEIGHT_RATIO;
        double coverage = 0.25;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    /// @param current_solution for initial call pass `Solution(problem.items.size(), problem.knapsack_sizes)` - this will construct empty solution.
    /// @param sorted_item_ids for initial call pass `problem.GetSortedItemIds(options.sort_mode);`
    /// @param amount_to_visit for initial call pass `options.coverage * problem.items.size()`
    static Solution IgnoreConnectionsDFS(const PackagedProblem & problem, const Solution & current_solution, const std::vector<int> & sorted_item_ids, const int & amount_to_visit);
    /// @param current_solution for initial call pass `Solution(problem.items.size(), problem.knapsack_sizes)` - this will construct empty solution.
    /// @param previous_item_id for initial call pass `-1` - this will omit checking if previous item has connection to the next one in inial call.
    /// @param sorted_item_ids for initial call pass `problem.GetSortedItemIds(options.sort_mode);`
    /// @param amount_to_visit for initial call pass `options.coverage * problem.items.size()`
    static Solution PathDFS(const PackagedProblem & problem, const Solution & current_solution, const int & previous_item_id, const std::vector<int> & sorted_item_ids, const int & amount_to_visit);
    /// @param current_solution for initial call pass `Solution(problem.items.size(), problem.knapsack_sizes)` - this will construct empty solution.
    /// @param previous_item_id for initial call pass `-1` - this will omit checking if previous item has connection to the next one in inial call.
    /// @param start_item_id for initial call pass `-1`
    /// @param sorted_item_ids for initial call pass `problem.GetSortedItemIds(options.sort_mode);`
    /// @param amount_to_visit for initial call pass `options.coverage * problem.items.size()`
    static Solution CycleDFS(const PackagedProblem & problem, const Solution & current_solution, const int & previous_item_id, const int & start_item_id, const std::vector<int> & sorted_item_ids, const int & amount_to_visit);

    static Solution Universal(const PackagedProblem & problem, const Options & options);

    static PackagedSolution Solve(PackagedProblem & problem, const Options & options);
    static std::string GetAlgorithmName(const Options & options);
    static bool expect_perfection;
};



class FloydSolver{
    public:
    FloydSolver() = delete;
    static PackagedSolution Solve(PackagedProblem & problem);
    static Solution Connected(const Problem & problem);
};



template <typename T>
inline PackagedSolution Solve(PackagedProblem & problem, const typename T::Options & options){
    PackagedSolution ps = T::Solve(problem, options);

    // fill algorithm info / details
    ps.algorithm = T::GetAlgorithmName(options);
    ps.algorithm += "_" + ToString(problem.requirements.structureToFind);
    ps.algorithm += "_" + ToString(problem.requirements.weightTreatment);
    
    // validate
    ps.validation_status = Validation::Validate(ps.solution, problem);
    if (T::expect_perfection && ps.solution.max_value < problem.known_optimum) {
        ps.validation_status.quality = false;
        ps.validation_status.valid = false;
    }

    // optimum update
    if (ps.validation_status.undergone && ps.validation_status.valid && ps.solution.max_value > problem.known_optimum){
        problem.known_optimum = ps.solution.max_value;
        if (problem.associated_file != "") problem.ExportJSON(problem.associated_file);
    }
    if (ps.solution.max_value == 0) ps.quality = problem.known_optimum == 0.0 ? 1.0 : DBL_MAX;
    else ps.quality = static_cast<double>(problem.known_optimum) / static_cast<double>(ps.solution.max_value);
    return ps;
}


template <typename T>
inline void BatchSolve(const std::string & directory_path, const typename T::Options & options) {
    // get the amount of problems to be solved
    std::ifstream fin (directory_path + FND_BATCH_INFO_FILE);
    nlohmann::json data = nlohmann::json::parse(fin);
    int amount = data["amount"];
    fin.close();

    // create solutions
    for (int i = 0; i < amount; ++i){
        PackagedProblem problem(directory_path + FND_PROBLEMS_FOLDER + FND_PROBLEM_FILE + std::to_string(i) + ".json");
        PackagedSolution ps = Solve<T>(problem, options);
        if (i == 0) {
            std::filesystem::remove_all(directory_path + "/" + ps.algorithm);
            std::filesystem::create_directories(directory_path + "/" + ps.algorithm);
        }
        ps.ExportJSON(directory_path + "/" + ps.algorithm + FND_SOLUTION_FILE + std::to_string(i) + ".json");
    }
}



} // namespace knapsack_solver



std::ostream& operator<<(std::ostream& os, const knapsack_solver::Solution& s);
std::ostream& operator<<(std::ostream& os, const knapsack_solver::PackagedSolution& ps);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::BruteForceSolver::Options::SearchOrder & so);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction & bf);

std::string ToString(const knapsack_solver::BruteForceSolver::Options::SearchOrder & so);
std::string ToString(const knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction & bf);

knapsack_solver::BruteForceSolver::Options::SearchOrder ToSearchOrder(const std::string & str);
knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction ToBoundingFunction(const std::string & str);