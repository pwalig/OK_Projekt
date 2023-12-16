#pragma once

#include "Problem.hpp"


namespace knapsack_solver{


class Solution {
    bool IsPathDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & length) const;
    bool IsCycleDFS(const Problem & problem, std::vector<int> & visited, const int & current, const int & start, const int & length) const;
    /// @brief Should not be called by anything else than `IsCyclePossible()`
    bool IsCyclePossibleDFS(const Problem & problem, std::vector<int> & visited, const std::vector<int> & _remaining_space, const int & current, const int & start) const;
    /// @brief Should not be called by anything else than `IsPathPossible()`
    bool IsPathPossibleDFS(const Problem & problem, std::vector<int> & visited, const std::vector<int> & _remaining_space, const int & current) const;

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
    /// @brief Check if path that would fit in the knapsack is possible, assuming that items can only be added and not removed.
    /// @param problem problem needed to define connections and knapsack sizes
    /// @throws invalid_argument - if problem.items is different size than this->selected
    bool IsPathPossible(const Problem & problem) const;
    /// @brief Check if path that would fit in the knapsack is possible, assuming that items can only be added and not removed.
    /// @param problem problem needed to define connections and knapsack sizes and structure to find
    /// @throws invalid_argument - if problem.items is different size than this->selected
    bool IsStructurePossible(const PackagedProblem & problem) const;
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


} // namespace knapsack_solver


std::ostream& operator<<(std::ostream& os, const knapsack_solver::Solution& s);
std::ostream& operator<<(std::ostream& os, const knapsack_solver::PackagedSolution& ps);