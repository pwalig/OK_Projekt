#pragma once

#include <string>
#include <vector>
#include <map>
#include <functional> // std::function
#include "KnapsackSolver.hpp" // PackagedSolution

namespace knapsack_solver{



class CommandInterpreter{
    enum class Command { RUN_TASKS, PRINT, SOLVE, BATCH_SOLVE, GENERATE_PROBLEM };
    static std::map<std::string, Command> command_map;
    /// @brief Runs Solve<T>, based on given arguments and provided type
    /// @tparam T Solver Type (eg. GreedySolver)
    /// @param args vector of arguments. Requires arg[0] to be a file name
    template <typename T> static PackagedSolution SolveFromArgs(std::vector<std::string> & args);
    /// @brief Runs BatchSolve<T>, based on given arguments and provided type
    /// @tparam T Solver Type (eg. GreedySolver)
    /// @param args vector of arguments. Requires arg[0] to be a path to batch folder
    template <typename T> static void BatchSolveFromArgs(std::vector<std::string> & args);
    public:
    /// @brief Consume Argument: Finds value in args. Executes lambda. And erases value from args.
    /// @param args vector of arguments. WARNING - will be modified.
    /// @param value argument to find
    /// @param lambda action to perform if argument was found
    template <typename T> static void Consume(std::vector<std::string> & args, const std::string & value, const std::function<void(const std::string &, T &)> & lambda, T & to_modify);
    /// @brief Consume Argument: Finds value in args. Executes lambda. And erases value from args.
    /// @param args vector of arguments. WARNING - will be modified.
    /// @param value argument to find
    /// @param lambda action to perform if argument was found
    template <typename T> static void Consume(std::vector<std::string> & args, const std::string & value, const std::function<void(T &)> & lambda, T & to_modify);
    /// @brief Consume Argument: returns argument under arg_id and erases it from args.
    /// @param args vector of arguments. WARNING - will be modified.
    /// @param arg_id index of the argument
    /// @returns the argument
    static std::string Consume(std::vector<std::string> & args, const int & arg_id);
    static void InterpretCommand(const int & argc, char const * const * const argv);
    static void InterpretCommand(const std::string & command, std::vector<std::string> args);
    static void RunTasks(const std::string & file_name);
};



} // namespace knapsack_solver