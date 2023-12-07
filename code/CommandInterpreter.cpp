#include "CommandInterpreter.hpp"

#include <fstream>
#include <iostream>
#include <sstream> // stringstream

#include "json.hpp"
#include "KnapsackSolver.hpp"

using std::endl;
using std::cout;
using std::string;
using std::vector;
using nlohmann::json;
using namespace knapsack_solver;

std::map<string, CommandInterpreter::Command> CommandInterpreter::command_map = 
{
    {"tasks", Command::RUN_TASKS}, 
    {"print", Command::PRINT}, 
    {"solve", Command::SOLVE},
    {"generate", Command::GENERATE_PROBLEM}
};

void CommandInterpreter::InterpretCommand(const int & argc, char const * const * const argv){
    vector<string> args;
    for (int i = 2; i < argc; ++i){
        args.push_back(argv[i]);
    }
    InterpretCommand(argv[1], args);
}

void CommandInterpreter::InterpretCommand(const string & command, const vector<string> & args) {
    if (command_map.find(command) == command_map.end())
        throw std::invalid_argument("" + command + " is not a recognised command");
        

    // declare iterator for browsing arguments
    auto it = args.begin();
    
    switch (command_map[command])
    {
    case Command::RUN_TASKS:{
        if (args.size() > 1) throw std::invalid_argument(args[1] + "is not a recognised argument for command " + command);
        RunTasks(args[0]);
        break;
    }

    case Command::SOLVE:{
        PackagedProblem pp(args[0]);
        PackagedSolution ps;

        if (args[1] == "greedy"){
            GreedySolver::Options op;
            
            for (it = args.begin() + 2; it != args.end(); ++it){
                if ((*it) == "-sort") {
                    ++it;
                    if (*it == "value/weight") op.sort_mode = Problem::SortMode::WEIGHT_VALUE_RATIO;
                    else if (*it == "value") op.sort_mode = Problem::SortMode::VALUE;
                    else if (*it == "weight") op.sort_mode = Problem::SortMode::WEIGHT;
                    else if (*it == "dont-sort") op.sort_mode = Problem::SortMode::DONT_SORT;
                    else if (*it == "random") op.sort_mode = Problem::SortMode::RANDOM;
                    else throw std::invalid_argument((*it) + " is not a valid sort method for greedy algorithm");
                }
                else if ((*it) == "-p") {}
                else if ((*it) == "-o") ++it;
                else throw std::invalid_argument((*it) + "is not a recognised argument for command " + command);
            }
            ps = GreedySolver::Solve(pp, op);
        }
        else if (args[1] == "brute-force"){
            BruteForceSolver::Options op;

            for (it = args.begin() + 2; it != args.end(); ++it){
                // recursive / iterative
                if ((*it) == "-recursive") op.iterative = false;
                // late fit / early fit
                else if ((*it) == "-early") op.late_fit = false;
                // search order
                else if ((*it) == "-order") {
                    ++it;
                    if (*it == "zero") op.search_order = BruteForceSolver::Options::SearchOrder::ZERO_FIRST;
                    else if (*it == "one") op.search_order = BruteForceSolver::Options::SearchOrder::ONE_FIRST;
                    else if (*it == "any") op.search_order = BruteForceSolver::Options::SearchOrder::UNCONSTRAINED;
                    else if (*it == "gray") op.search_order = BruteForceSolver::Options::SearchOrder::GRAY_CODE;
                    else if (*it == "random") op.search_order = BruteForceSolver::Options::SearchOrder::RANDOM;
                    else throw std::invalid_argument((*it) + " is not a valid visit order for brute force algorithm");
                }
                else if ((*it) == "-p") {}
                else if ((*it) == "-o") ++it;
                else throw std::invalid_argument((*it) + "is not a recognised argument for command " + command);
            }
            ps = BruteForceSolver::Solve(pp, op);
        }
        else if (args[1] == "branch-and-bound"){
            BranchAndBoundSolver::Options op;

            for (it = args.begin() + 2; it != args.end(); ++it){
                // late fit / early fit
                if ((*it) == "-early") op.late_fit = false;
                // search order
                else if ((*it) == "-bf") {
                    ++it;
                    if (*it == "dynamic") op.bounding_function = BranchAndBoundSolver::Options::BoundingFunction::BASE_DYNAMIC;
                    else if (*it == "continous") op.bounding_function = BranchAndBoundSolver::Options::BoundingFunction::CONTINOUS;
                    else if (*it == "none") op.bounding_function = BranchAndBoundSolver::Options::BoundingFunction::NONE;
                    else throw std::invalid_argument((*it) + " is not a valid bounding algorithm for branch and bound algorithm");
                }
                else if ((*it) == "-p") {}
                else if ((*it) == "-o") ++it;
                else throw std::invalid_argument((*it) + "is not a recognised argument for command " + command);
            }
            ps = BranchAndBoundSolver::Solve(pp, op);
        }
        else 
            throw std::invalid_argument(args[1] + " is not a recognised algorithm");

        // output file
        it = std::find(args.begin(), args.end(), "-o");
        if (it != args.end()) ps.ExportJSON(*(++it));
        
        // print on console
        if (std::find(args.begin(), args.end(), "-p") != args.end()) cout << ps;

        break;
    }

    case Command::PRINT:{
        if (args.size() > 1) throw std::invalid_argument(args[1] + "is not a recognised argument for command " + command);
        cout << args[0];
        break;
    }

    case Command::GENERATE_PROBLEM:{
        Problem::GenerationSettings gs;
        Problem::Requirements rq;
        int amount = 0;

        for (it = args.begin() + 1; it != args.end(); ++it){

            // instance size
            if ((*it) == "-isize") gs.instance_size = std::stoi(*(++it));
            // sub knackpacks
            else if ((*it) == "-subk") gs.sub_knapsacks = std::stoi(*(++it));
            // knapsack size limit
            else if ((*it) == "-maxks") gs.knapsack_size_limit_exclusive = std::stoi(*(++it));
            // value limit
            else if ((*it) == "-maxv") gs.value_limit_exclusive = std::stoi(*(++it));
            // weight limit
            else if ((*it) == "-maxw") gs.weight_limit_exclusive = std::stoi(*(++it));
            // connection density
            else if ((*it) == "-cd") gs.connection_density = std::stod(*(++it));

            // structure to find
            else if ((*it) == "-structure") {
                ++it;
                if (*it == "path") rq.structureToFind = Problem::Requirements::StructureToFind::PATH;
                else if (*it == "cycle") rq.structureToFind = Problem::Requirements::StructureToFind::CYCLE;
                else if (*it == "tree") rq.structureToFind = Problem::Requirements::StructureToFind::TREE;
                else if (*it == "connected") rq.structureToFind = Problem::Requirements::StructureToFind::CONNECTED_GRAPH;
                else if (*it == "ignore") rq.structureToFind = Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS;
                else throw std::invalid_argument((*it) + " is not a valid sort method for greedy algorithm");
            }

            // batch or not
            else if ((*it) == "-batch") {
                amount = std::stod(*(++it));
                if (amount <= 0) throw std::invalid_argument("in \"-batch\" amount should be more than 0");
            }
            
            // invalid argument
            else throw std::invalid_argument((*it) + "is not a recognised argument for command " + command);

        }
        
        // call generate
        if (amount > 0) PackagedProblem::BatchGeneratePackagedProblemsJSON(gs, rq, amount, args[0]);
        else PackagedProblem::GeneratePackagedProblemJSON(gs, rq, args[0]);
        
        break;
    }
    
    default:{
        throw std::invalid_argument("" + command + " is not implemented yet");
        break;
    }
    }
}

void CommandInterpreter::RunTasks(const std::string & file_name) {
    std::ifstream fin (file_name);
    json data = json::parse(fin);
    for (json command : data["commands"]){
        InterpretCommand(command["command"], command["args"]);
    }
    fin.close();
}