#pragma once

#include <iostream> // std::ostream
#include <fstream> // std::ifstream
#include <vector>
#include <string>
#include <cfloat> // DBL_MAX
#include <functional> // std::function
#include <thread>

// Batch Solve requires these three:
#include <filesystem> // std::filesystem::create_directory, std::filesystem::remove_all
#include "json.hpp"
#include "FileNameDefines.hpp"

#include "Problem.hpp"
#include "Solution.hpp"

namespace knapsack_solver {




class KnapsackSolver {
    public:
    static int GoalFunction(const Solution & solution, const PackagedProblem & problem);
};

void RunMassTests(const std::vector<Problem::GenerationSettings> & gsv, const std::vector<Problem::Requirements> & rqv, const int & repeats);


class MassTestResult {
    struct ValidationStatus{
        double undergone = 0.0;
        double valid = 0.0;
        double value = 0.0;
        double remaining_space = 0.0;
        double fit = 0.0;
        double structure = 0.0;
        double self_valid = 0.0;
        double quality = 0.0;
    };

    double solve_time = 0.0;
    double quality = 0.0;
    int max_value = 0.0;
    std::vector<double> remaining;
    int opt_max_value_sum = 0;
    ValidationStatus validation_status;
    double overall_quality = 0.0;

    public:
    int amount = 0;
    void DivideByAmount();
    /// @brief division by ammount is handled by `MassTestResult::ExportJSON(const std::string & file_name)`
    void AddSolution(const PackagedSolution & ps, const int & optimum);
    /// @brief exports .json file and handles the division by amount
    void ExportJSON(const std::string & file_name);
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



class GreedySolver{
    public:
    GreedySolver() = delete;
    struct Options{
        Problem::SortMode sort_mode = Problem::SortMode::VALUE_WEIGHT_RATIO;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    static Solution NaiveUniversal(const PackagedProblem & problem, const Options & options);

    /// @deprecated replaced by `Solution` `GreedySolver::NaiveUniversal(const PackagedProblem & problem, const Options & options)`
    static Solution GreedyUniversal(const PackagedProblem & problem, const Options & options);
    /// @deprecated replaced by `Solution` `GreedySolver::NaiveUniversal(const PackagedProblem & problem, const Options & options)`
    static Solution GreedyIgnoreConnections(const Problem & problem, const Options & options);
    /// @deprecated replaced by `Solution` `GreedySolver::NaiveUniversal(const PackagedProblem & problem, const Options & options)`
    static Solution GreedyPath(const Problem & problem, const Options & options);

    static PackagedSolution Solve(PackagedProblem & problem, const Options & options);
    static std::string GetAlgorithmName(const Options & options);
    static bool expect_perfection;
};



class GRASPSolver{
    public:
    GRASPSolver() = delete;
    struct Options{
        Problem::SortMode sort_mode = Problem::SortMode::VALUE_WEIGHT_RATIO;
        int iterations = 1;
        double coverage = 0.0;
        double chose_from = 0.25;
        Options() = default;
        explicit Options(std::vector<std::string> & args);
    };

    static Solution Universal(const PackagedProblem & problem, const Options & options);
    static Solution MultiRun(const PackagedProblem & problem, const Options & options);

    static PackagedSolution Solve(PackagedProblem & problem, Options options);
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

    // create apropriate folder
    PackagedProblem problem0(directory_path + FND_PROBLEMS_FOLDER + FND_PROBLEM_FILE + "0.json");
    std::string algo_name = T::GetAlgorithmName(options);
    algo_name += "_" + ToString(problem0.requirements.structureToFind);
    algo_name += "_" + ToString(problem0.requirements.weightTreatment);
    std::filesystem::remove_all(directory_path + "/" + algo_name);
    std::filesystem::create_directories(directory_path + "/" + algo_name);

    // prepare threads
    struct Pair1 {
        PackagedSolution solution;
        int optimum;
    };
    std::vector<std::thread> threads;
    std::vector<Pair1> pairs(amount, Pair1());
    threads.reserve(amount);

    // run threads
    for (int t = 0; t < amount; ++t) {
        threads.push_back(std::thread([directory_path, options, &pairs, t](){

            PackagedProblem problem(directory_path + FND_PROBLEMS_FOLDER + FND_PROBLEM_FILE + std::to_string(t) + ".json");
            PackagedSolution ps = Solve<T>(problem, options);

            pairs[t].solution = ps;
            pairs[t].optimum = problem.known_optimum;
            
            ps.ExportJSON(directory_path + "/" + ps.algorithm + FND_SOLUTION_FILE + std::to_string(t) + ".json");
        }));
    }

    // calculate the averages
    MassTestResult mtr;
    mtr.amount = amount;
    for (int i = 0; i < amount; ++i) {
        threads[i].join();
        mtr.AddSolution(pairs[i].solution, pairs[i].optimum);
    }

    // save solution info to file
    mtr.ExportJSON(directory_path + "/" + algo_name + "/solve-info.json");
}



} // namespace knapsack_solver



std::ostream& operator<<(std::ostream & os, const knapsack_solver::BruteForceSolver::Options::SearchOrder & so);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction & bf);

std::string ToString(const knapsack_solver::BruteForceSolver::Options::SearchOrder & so);
std::string ToString(const knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction & bf);

knapsack_solver::BruteForceSolver::Options::SearchOrder ToSearchOrder(const std::string & str);
knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction ToBoundingFunction(const std::string & str);