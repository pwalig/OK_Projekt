#pragma once

#include <string>
#include <vector>
#include <map>
#include <functional> // std::function

namespace knapsack_solver{



class CommandInterpreter{
    enum class Command { RUN_TASKS, PRINT, SOLVE, BATCH_SOLVE, GENERATE_PROBLEM };
    static std::map<std::string, Command> command_map;
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
    static void InterpretCommand(const int & argc, char const * const * const argv);
    static void InterpretCommand(const std::string & command, std::vector<std::string> args);
    static void RunTasks(const std::string & file_name);
};



} // namespace knapsack_solver