#include "CommandInterpreter.hpp"

#include <fstream>
#include <iostream>
#include <sstream> // stringstream

#include "json.hpp"
#include "KnapsackSolver.hpp"

/*
#define RUN_TASKS 0
#define GREEDY 10
#define BRUTE_FORCE 11
#define BRANCH_AND_BOUND 12
#define DYNAMIC 13
#define FLOYD 14
*/

using std::endl;
using std::cout;
using std::string;
using std::vector;
using nlohmann::json;
using namespace knapsack_solver;

std::map<string, CommandInterpreter::Command> CommandInterpreter::command_map = 
{
    {"tasks", Command::RUN_TASKS}, 
    {"test", Command::TEST_PRINT}, 
    {"solve", Command::SOLVE},
    {"generate", Command::GENERATE_PROBLEM}
};

void CommandInterpreter::InterpretCommand(const int & argc, char const * const * const argv){
    vector<string> args;
    for (int i = 2; i < argc; ++i){
        args.push_back(argv[i]);
    }
    //cout << argv[1] << " " << args[0] << endl;
    InterpretCommand(argv[1], args);
}

void CommandInterpreter::InterpretCommand(const string & command, const vector<string> & args) {
    if (command_map.find(command) == command_map.end())
        throw std::invalid_argument("" + command + " is not a recognised command");
        

    // declare iterator for browsing arguments
    auto it = args.begin();
    vector<bool> unrecognised(true, args.size());
    
    switch (command_map[command])
    {
    case Command::RUN_TASKS:{
        RunTasks(args[0]);
        break;
    }

    case Command::SOLVE:{
        cout << args[0] << endl;
        PackagedProblem pp(args[0]);
        PackagedSolution ps;

        if (args[1] == "greedy"){
            GreedySolver::Options op;
            it = std::find(args.begin(), args.end(), "-sort");
            if (it != args.end()) {
                ++it;
                if (*it == "value/weight") op.sort_mode = Problem::SortMode::WEIGHT_VALUE_RATIO;
                else if (*it == "value") op.sort_mode = Problem::SortMode::VALUE;
                else if (*it == "weight") op.sort_mode = Problem::SortMode::WEIGHT;
                else if (*it == "dont-sort") op.sort_mode = Problem::SortMode::DONT_SORT;
                else if (*it == "random") op.sort_mode = Problem::SortMode::RANDOM;
                else throw std::invalid_argument((*it) + " is not a valid sort method for greedy algorithm");
            }
            ps = GreedySolver::Solve(pp, op);
        }
        else if (args[1] == "brute-force"){
            BruteForceSolver::Options op;
            // recursive / iterative
            op.iterative = (std::find(args.begin(), args.end(), "-recursive") == args.end());
            // late fit / early fit
            op.late_fit = (std::find(args.begin(), args.end(), "-early") == args.end());
            // search order
            it = std::find(args.begin(), args.end(), "-order");
            if (it != args.end()) {
                ++it;
                if (*it == "zero") op.search_order = BruteForceSolver::Options::SearchOrder::ZERO_FIRST;
                else if (*it == "one") op.search_order = BruteForceSolver::Options::SearchOrder::ONE_FIRST;
                else if (*it == "any") op.search_order = BruteForceSolver::Options::SearchOrder::UNCONSTRAINED;
                else if (*it == "gray") op.search_order = BruteForceSolver::Options::SearchOrder::GRAY_CODE;
                else if (*it == "random") op.search_order = BruteForceSolver::Options::SearchOrder::RANDOM;
                else throw std::invalid_argument((*it) + " is not a valid sort method for greedy algorithm");
            }
            ps = BruteForceSolver::Solve(pp, op);
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

    case Command::TEST_PRINT:{
        cout << args[0];
        break;
    }

    case Command::GENERATE_PROBLEM:{
        // generation settings
        Problem::GenerationSettings gs;

        // instance size
        it = std::find(args.begin(), args.end(), "-isize");
        if (it != args.end()){
            gs.instance_size = std::stoi(*(++it));
        }
        // sub knackpacks
        it = std::find(args.begin(), args.end(), "-subk");
        if (it != args.end()){
            gs.sub_knapsacks = std::stoi(*(++it));
        }
        // knapsack size limit
        it = std::find(args.begin(), args.end(), "-maxks");
        if (it != args.end()){
             gs.knapsack_size_limit_exclusive = std::stoi(*(++it));
        }
        // value limit
        it = std::find(args.begin(), args.end(), "-maxv");
        if (it != args.end()){
            gs.value_limit_exclusive = std::stoi(*(++it));
        }
        // weight limit
        it = std::find(args.begin(), args.end(), "-maxw");
        if (it != args.end()){
            gs.weight_limit_exclusive = std::stoi(*(++it));
        }
        // connection density
        it = std::find(args.begin(), args.end(), "-cd");
        if (it != args.end()){
            gs.connection_density = std::stod(*(++it));
        }

        // requirements
        Problem::Requirements rq;

        // structure to find
        it = std::find(args.begin(), args.end(), "-structure");
        if (it != args.end()) {
            ++it;
            if (*it == "path") rq.structureToFind = Problem::Requirements::StructureToFind::PATH;
            else if (*it == "cycle") rq.structureToFind = Problem::Requirements::StructureToFind::CYCLE;
            else if (*it == "tree") rq.structureToFind = Problem::Requirements::StructureToFind::TREE;
            else if (*it == "connected") rq.structureToFind = Problem::Requirements::StructureToFind::CONNECTED_GRAPH;
            else if (*it == "ignore") rq.structureToFind = Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS;
            else throw std::invalid_argument((*it) + " is not a valid sort method for greedy algorithm");
        }
            
        // batch or not
        it = std::find(args.begin(), args.end(), "-batch");
        if (it != args.end()){
            int amount = std::stod(*(++it));
            PackagedProblem::BatchGeneratePackagedProblemsJSON(gs, rq, amount, args[0]);
        }
        else{
            PackagedProblem::GeneratePackagedProblemJSON(gs, rq, args[0]);
        }
        
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