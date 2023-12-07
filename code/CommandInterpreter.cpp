#include "CommandInterpreter.hpp"
//"CommandInterpreter.hpp" includes:
//#include <string>
//#include <vector>
//#include <map>
//#include <functional> // std::function

#include <fstream> // std::ifstream in RunTasks()
#include <iostream> // std::cout, std::endl
#include <stdexcept> // throw error types

#include "json.hpp" // required by RunTasks()
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
    {"batch-solve", Command::BATCH_SOLVE},
    {"generate", Command::GENERATE_PROBLEM}
};

template <typename T>
void CommandInterpreter::Consume(std::vector<std::string> & args, const std::string & value, const std::function<void(const std::string &, T &)> & lambda, T & to_modify){
    auto it = std::find(args.begin(), args.end(), value);
    if (it != args.end()) {
        it = args.erase(it);
        lambda((*it), to_modify);
        args.erase(it);
    }
}

template <typename T>
void CommandInterpreter::Consume(std::vector<std::string> & args, const std::string & value, const std::function<void(T &)> & lambda, T & to_modify){
    auto it = std::find(args.begin(), args.end(), value);
    if (it != args.end()) {
        it = args.erase(it);
        lambda(to_modify);
    }
}



// ---------- OPTIONS CONSTRUCTORS ----------

GreedySolver::Options::Options(std::vector<std::string> & args){
    // sort method
    CommandInterpreter::Consume<Problem::SortMode>(args, "-sort", [](const string & arg, Problem::SortMode & sm){
        if (arg == "value/weight") sm = Problem::SortMode::WEIGHT_VALUE_RATIO;
        else if (arg == "value") sm = Problem::SortMode::VALUE;
        else if (arg == "weight") sm = Problem::SortMode::WEIGHT;
        else if (arg == "dont-sort") sm = Problem::SortMode::DONT_SORT;
        else if (arg == "random") sm = Problem::SortMode::RANDOM;
        else throw std::invalid_argument(arg + " is not a valid visit order for brute force algorithm");
    }, this->sort_mode);
}

BruteForceSolver::Options::Options(std::vector<std::string> & args) : Options(){
    // recursive / iterative
    CommandInterpreter::Consume<bool>(args, "-recursive", [](bool & b){
        b = false;
    }, this->iterative);

    // late fit / early fit
    CommandInterpreter::Consume<bool>(args, "-early", [](bool & b){
        b = false;
    }, this->late_fit);

    // search order
    CommandInterpreter::Consume<Options::SearchOrder>(args, "-order", [](const string & arg, Options::SearchOrder & so){
        if (arg == "zero") so = BruteForceSolver::Options::SearchOrder::ZERO_FIRST;
        else if (arg == "one") so = BruteForceSolver::Options::SearchOrder::ONE_FIRST;
        else if (arg == "any") so = BruteForceSolver::Options::SearchOrder::UNCONSTRAINED;
        else if (arg == "gray") so = BruteForceSolver::Options::SearchOrder::GRAY_CODE;
        else if (arg == "random") so = BruteForceSolver::Options::SearchOrder::RANDOM;
        else throw std::invalid_argument(arg + " is not a valid visit order for brute force algorithm");
    }, this->search_order);
}

BranchAndBoundSolver::Options::Options(std::vector<std::string> & args) : Options(){
    // late fit / early fit
    CommandInterpreter::Consume<bool>(args, "-early", [](bool & b){
        b = false;
    }, this->late_fit);

    // bounding function
    CommandInterpreter::Consume<Options::BoundingFunction>(args, "-bf", [](const string & arg, Options::BoundingFunction & bf){
        if (arg == "dynamic") bf = BranchAndBoundSolver::Options::BoundingFunction::BASE_DYNAMIC;
        else if (arg == "continous") bf = BranchAndBoundSolver::Options::BoundingFunction::CONTINOUS;
        else if (arg == "none") bf = BranchAndBoundSolver::Options::BoundingFunction::NONE;
        else throw std::invalid_argument(arg + " is not a valid bounding algorithm for branch and bound algorithm");
    }, this->bounding_function);
}


// ---------- COMMAND INTERPRETING ----------

void CommandInterpreter::InterpretCommand(const int & argc, char const * const * const argv){
    vector<string> args;
    for (int i = 2; i < argc; ++i){
        args.push_back(argv[i]);
    }
    InterpretCommand(argv[1], args);
}

void CommandInterpreter::InterpretCommand(const string & command, vector<string> args) {
    if (command_map.find(command) == command_map.end())
        throw std::invalid_argument("" + command + " is not a recognised command");
        

    // declare iterator for browsing arguments
    auto it = args.begin();
    
    switch (command_map[command])
    {
    case Command::RUN_TASKS:{
        RunTasks(args[0]);
        args.erase(args.begin());
        break;
    }

    case Command::SOLVE:{
        PackagedProblem pp(args[0]);
        args.erase(args.begin());
        PackagedSolution ps;

        if (args[0] == "greedy"){
            GreedySolver::Options op(args);
            ps = GreedySolver::Solve(pp, op);
        }
        else if (args[0] == "brute-force"){
            BruteForceSolver::Options op(args);
            ps = BruteForceSolver::Solve(pp, op);
        }
        else if (args[0] == "branch-and-bound"){
            BranchAndBoundSolver::Options op(args);
            ps = BranchAndBoundSolver::Solve(pp, op);
        }
        else 
            throw std::invalid_argument(args[0] + " is not a recognised algorithm");
        args.erase(args.begin());

        // output file
        Consume<PackagedSolution>(args, "-o", [](const string & arg, PackagedSolution & el){
            el.ExportJSON(arg);
        }, ps);
        
        // print on console
        Consume<PackagedSolution>(args, "-p", [](PackagedSolution & el){
            cout << el;
        }, ps);

        break;
    }
    case Command::BATCH_SOLVE:{

        if (args[1] == "greedy"){
            GreedySolver::Options op(args);
            BatchSolve<GreedySolver>(args[0], op);
        }
        else if (args[1] == "brute-force"){
            BruteForceSolver::Options op(args);
            BatchSolve<BruteForceSolver>(args[0], op);
        }
        else if (args[1] == "branch-and-bound"){
            BranchAndBoundSolver::Options op(args);
            BatchSolve<BranchAndBoundSolver>(args[0], op);
        }
        else 
            throw std::invalid_argument(args[1] + " is not a recognised algorithm");
        args.erase(args.begin());
        args.erase(args.begin());

        break;
    }

    case Command::PRINT:{
        cout << args[0];
        args.erase(args.begin());
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
            else if ((*it) == "-maxks") { gs.knapsack_size_limit_exclusive = std::stoi(*(++it)); gs.randomize_knapsack_sizes = true; }
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

        args.clear(); // "generate" does not leave unhandled arguments - iterates through the whole thing and throws error on first unrecognised argument

        break;
    }
    
    default:{
        throw std::invalid_argument("" + command + " is not implemented yet");
        break;
    }
    }
    
    // unhandled arguments
    for (auto arg = args.begin(); arg != args.end(); ++arg)
        cout << "argument: " + (*arg) + " was not recognised and ignored.\n";
}

void CommandInterpreter::RunTasks(const std::string & file_name) {
    std::ifstream fin (file_name);
    json data = json::parse(fin);
    for (json command : data["commands"]){
        InterpretCommand(command["command"], command["args"]);
    }
    fin.close();
}