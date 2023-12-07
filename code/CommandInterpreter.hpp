#pragma once

#include <string>
#include <vector>
#include <map>

namespace knapsack_solver{



class CommandInterpreter{
    enum class Command { RUN_TASKS, PRINT, SOLVE, GENERATE_PROBLEM };
    static std::map<std::string, Command> command_map;
    public:
    static void InterpretCommand(const int & argc, char const * const * const argv);
    static void InterpretCommand(const std::string & command, const std::vector<std::string> & args);
    static void RunTasks(const std::string & file_name);
};



} // namespace knapsack_solver