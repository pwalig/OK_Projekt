#pragma once

#include <iostream> // std::ostream
#include <fstream> // std::ifstream
#include <vector>
#include <string>

// Batch Solve requires these three:
#include <filesystem> // std::filesystem::create_directory, std::filesystem::remove_all
#include "json.hpp"
#include "FileNameDefines.hpp"

#include "Problem.hpp"

namespace knapsack_solver {



class Solution {
    bool IsPathDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & length) const;
    bool IsCycleDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & start, const int & length) const;
    bool IsCyclePossibleDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & start) const;

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
    bool Fits(const Problem & problem, const int & selected_item_id) const;
    
    /// @returns true if addition was succesfull, false if item did not fit or was already in solution
    bool AddItemIfFits(const Problem & problem, const int & selected_item_id);

    bool IsFit(const Problem & problem) const;

    bool IsPath(const Problem & problem) const;
    bool IsCycle(const Problem & problem) const;
    bool IsTree(const Problem & problem) const;
    bool IsConnected(const Problem & problem) const;
    bool IsStructure(const PackagedProblem & problem) const;

    bool IsValid(const PackagedProblem & problem) const;

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
    static std::vector<int> CalculateRemainingSpaces(const Solution & solution, const Problem & problem);
    static int CalculateMaxValue(const Solution & solution, const Problem & problem);
    static ValidationStatus Validate(const Solution & solution, const PackagedProblem & problem);
    static int GoalFunction(const Solution & solution, const PackagedProblem & problem);
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
    static PackagedSolution Solve(PackagedProblem & problem, const Options options);
    static std::string GetAlgorithmName(const PackagedProblem & problem, const Options & options);
    
    static Solution Iterative(const PackagedProblem & problem, const Options::SearchOrder & search_order);
    static Solution Recursive(const PackagedProblem & problem, const Options & options);
};



class BranchAndBoundSolver{
    private:
    BranchAndBoundSolver() = delete;
    static int GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode, Solution currentSolution);

    static Solution DFSLateFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId);
    static Solution DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, int lower_bound);

    public:
    struct Options{
        enum class BoundingFunction { NONE, CONTINOUS, ACYCLIC, BASE_DYNAMIC };
        BoundingFunction bounding_function = BoundingFunction::NONE;
        bool late_fit = false;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    static Solution BnBLateFitPath(const Problem & problem);
    static Solution BnBEarlyFitPath(const Problem & problem);

    static PackagedSolution Solve(PackagedProblem & problem, const Options & options);
    static std::string GetAlgorithmName(const PackagedProblem & problem, const Options & options);
};



class GreedySolver{
    public:
    GreedySolver() = delete;
    struct Options{
        Problem::SortMode sort_mode = Problem::SortMode::VALUE_WEIGHT_RATIO;
        int buffor = 1;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    static PackagedSolution Solve(PackagedProblem & problem, const Options & options);
    static std::string GetAlgorithmName(const PackagedProblem & problem, const Options & options);

    static Solution GreedyUniversal(const PackagedProblem & problem, const Options & options);
    static Solution GreedyIgnoreConnections(const Problem & problem, const Options & options);
    static Solution GreedyPath(const Problem & problem, const Options & options);
};



class FloydSolver{
    public:
    FloydSolver() = delete;
    static PackagedSolution Solve(PackagedProblem & problem);
    static Solution Connected(const Problem & problem);
};



template <typename T>
inline PackagedSolution Solve(PackagedProblem & problem, const typename T::Options & options){
    return T::Solve(problem, options);
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
        PackagedSolution ps = T::Solve(problem, options);
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