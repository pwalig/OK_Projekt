#pragma once

#include <vector>
#include <iostream>
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
    void AddItemIfFits(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace);
};

struct PackagedSolution {
    std::string algorithm;
    Solution solution;
    double to_optimum_ratio;
    double solve_time;
    std::vector<int> remainingSpaces;
};

class KnapsackSolver{
    public:
    enum class BoundingFunction { VALUE_SUM, GREEDY_IGNORE_CONNECTIONS, GREEDY };

    private:
    static bool Fits(const std::vector<int> & weights, const std::vector<int> & remainingSpace);
    static std::vector<int> CalculateRemainingSpaces(const Solution & solution, const Problem & problem);

    static Solution DFSPath(const Problem & problem, Solution currentSolution, const int & currentItemId, std::vector<int> remainingSpace);

    static Solution BranchAndBoundDFSPath(const Problem & problem, Solution currentSolution, const int & currentItemId, std::vector<int> remainingSpace, const int & lower_bound);
    static int GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode, Solution currentSolution, std::vector<int> remainingSpace);
    static Solution GreedyPath(const Problem & problem, Solution currentSolution, const int & currentItemId, std::vector<int> remainingSpace);
    

    public:
    // has bugs
    static Solution FirstFit(const Problem & problem);
    // has bugs
    static Solution RandomFit(const Problem & problem);

    static Solution GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode);
    static Solution GreedyPath(const Problem & problem, const Problem::SortMode & sortMode);
    static Solution BruteForcePath(const Problem & problem);
    static Solution BranchAndBoundPath(const Problem & problem);

    static PackagedSolution Greedy(const Problem & problem, Requirements requirements);
    static PackagedSolution BruteForce(const Problem & problem, Requirements requirements);
    static PackagedSolution BranchAndBound(const Problem & problem, Requirements requirements);
};


} // namespace knapsack_solver

std::ostream& operator<<(std::ostream& os, const knapsack_solver::Solution& s);

std::ostream& operator<<(std::ostream& os, const knapsack_solver::PackagedSolution& ps);