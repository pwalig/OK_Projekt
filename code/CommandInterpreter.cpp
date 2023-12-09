#include "CommandInterpreter.hpp"
//"CommandInterpreter.hpp" includes:
//#include <string>
//#include <vector>
//#include <map>
//#include <functional> // std::function
//#include "KnapsackSolver.hpp"

#include <fstream> // std::ifstream in RunTasks()
#include <iostream> // std::cout, std::endl
#include <stdexcept> // throw error types

#include "json.hpp" // required by RunTasks()

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

std::string CommandInterpreter::Consume(std::vector<std::string> & args, const int & arg_id){
    string arg = args[arg_id];
    args.erase(args.begin() + arg_id);
    return arg;
}



// ---------- OPTIONS CONSTRUCTORS ----------

GreedySolver::Options::Options(std::vector<std::string> & args){
    // multi-run / single-run
    CommandInterpreter::Consume<bool>(args, "-single-run", [](const string & arg, bool & sr){
        sr = false;
    }, this->multi_run);
    // sort method
    CommandInterpreter::Consume<Problem::SortMode>(args, "-sort", [](const string & arg, Problem::SortMode & sm){
        sm = ToSortMode(arg);
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
        so = ToSearchOrder(arg);
    }, this->search_order);
}

BranchAndBoundSolver::Options::Options(std::vector<std::string> & args) : Options(){
    // late fit / early fit
    CommandInterpreter::Consume<bool>(args, "-early", [](bool & b){
        b = false;
    }, this->late_fit);

    // bounding function
    CommandInterpreter::Consume<Options::BoundingFunction>(args, "-bf", [](const string & arg, Options::BoundingFunction & bf){
        bf = ToBoundingFunction(arg);
    }, this->bounding_function);
}


// ---------- COMMAND INTERPRETING ----------

template <typename T>
void CommandInterpreter::BatchSolveFromArgs(std::vector<std::string> & args){
    typename T::Options op(args);
    if (args.size() > 1) throw std::invalid_argument(args[1] + " is not a recognised argument for command " + "batch-solve");
    BatchSolve<T>(Consume(args, 0), op);
}

template <typename T>
PackagedSolution CommandInterpreter::SolveFromArgs(std::vector<std::string> & args){
    typename T::Options op(args);
    if (args.size() > 1) throw std::invalid_argument(args[1] + " is not a recognised argument for command " + "solve");
    PackagedProblem pp(Consume(args, 0));
    return Solve<T>(pp, op);
}

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
        if (args.size() > 1) throw std::invalid_argument("run tasks takes only one argument, provided: " + args.size());
        RunTasks(Consume(args, 0));
        break;
    }

    case Command::SOLVE:{
        // output file
        string export_file = "";
        Consume<string>(args, "-o", [](const string & arg, string & el){
            el = arg;
        }, export_file);
        
        // print on console
        bool print = false;
        Consume<bool>(args, "-p", [](bool & el){
            el = true;
        }, print);

        PackagedSolution ps;
        string algorithm = Consume(args, 1);
        if (algorithm == "greedy") ps = SolveFromArgs<GreedySolver>(args);
        else if (algorithm == "brute-force") ps = SolveFromArgs<BruteForceSolver>(args);
        else if (algorithm == "branch-and-bound") ps = SolveFromArgs<BranchAndBoundSolver>(args);
        else throw std::invalid_argument(algorithm + " is not a recognised algorithm");
            
        if (print) cout << ps;
        if (export_file != "") ps.ExportJSON(export_file);

        break;
    }
    case Command::BATCH_SOLVE:{
        string algorithm = Consume(args, 1);
        if (algorithm == "greedy") BatchSolveFromArgs<GreedySolver>(args);
        else if (algorithm == "brute-force") BatchSolveFromArgs<BruteForceSolver>(args);
        else if (algorithm == "branch-and-bound") BatchSolveFromArgs<BranchAndBoundSolver>(args);
        else throw std::invalid_argument(algorithm + " is not a recognised algorithm");
        break;
    }

    case Command::PRINT:{
        if (args.size() > 1) throw std::invalid_argument("print takes only one argument, provided: " + args.size());
        cout << Consume(args, 0);
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
            else if ((*it) == "-structure") rq.structureToFind = ToStructureToFind(*(++it));

            // batch or not
            else if ((*it) == "-batch") {
                amount = std::stod(*(++it));
                if (amount <= 0) throw std::invalid_argument("in \"-batch\" amount should be more than 0");
            }
            
            // invalid argument
            else throw std::invalid_argument((*it) + " is not a recognised argument for command " + command);

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