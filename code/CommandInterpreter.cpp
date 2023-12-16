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
    {"generate", Command::GENERATE_PROBLEM},
    {"gather", Command::GATHER}
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

void CommandInterpreter::Consume(std::vector<std::string> & args, const std::string & value, const std::function<void(const std::string &)> & lambda){
    auto it = std::find(args.begin(), args.end(), value);
    if (it != args.end()) {
        it = args.erase(it);
        lambda((*it));
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

void CommandInterpreter::Consume(std::vector<std::string> & args, const std::string & value, const std::function<void()> & lambda){
    auto it = std::find(args.begin(), args.end(), value);
    if (it != args.end()) {
        it = args.erase(it);
        lambda();
    }
}

std::string CommandInterpreter::Consume(std::vector<std::string> & args, const int & arg_id){
    string arg = args[arg_id];
    args.erase(args.begin() + arg_id);
    return arg;
}




// ---------- OPTIONS CONSTRUCTORS ----------

GreedySolver::Options::Options(std::vector<std::string> & args){
    // sort method
    CommandInterpreter::Consume<Problem::SortMode>(args, "-sort", [](const string & arg, Problem::SortMode & sm){
        sm = ToSortMode(arg);
    }, this->sort_mode);
}

GRASPSolver::Options::Options(std::vector<std::string> & args){
    // iterations
    CommandInterpreter::Consume<GRASPSolver::Options>(args, "-iterations", [](const string & arg, GRASPSolver::Options & op){
        op.iterations = std::stod(arg);
        if (op.iterations <= 0) throw std::invalid_argument("grasp iterations should be positive");
    }, *this);
    // to chose from
    CommandInterpreter::Consume<double>(args, "-chose-from", [](const string & arg, double & dbl){
        dbl = std::stod(arg);
    }, this->chose_from);
    // sort method
    CommandInterpreter::Consume<Problem::SortMode>(args, "-sort", [](const string & arg, Problem::SortMode & sm){
        sm = ToSortMode(arg);
    }, this->sort_mode);
}

GreedyHeuristicSearchSolver::Options::Options(std::vector<std::string> & args){
    // coverage
    CommandInterpreter::Consume<GreedyHeuristicSearchSolver::Options>(args, "-coverage", [](const string & arg, GreedyHeuristicSearchSolver::Options & op){
        op.to_visit = -1;
        op.coverage = std::stod(arg);
    }, *this);
    // to-visit
    CommandInterpreter::Consume<int>(args, "-to-visit", [](const string & arg, int & tvs){
        tvs = std::stod(arg);
    }, this->to_visit);
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
    CommandInterpreter::Consume<bool>(args, "-iterative", [](bool & b){
        if (b == false) throw std::invalid_argument("tried to use -recursive and -iterative simultaneously");
    }, this->iterative);

    // late fit / early fit
    CommandInterpreter::Consume<bool>(args, "-early", [](bool & b){
        b = false;
    }, this->late_fit);
    CommandInterpreter::Consume<bool>(args, "-late", [](bool & b){
        if (b == false) throw std::invalid_argument("tried to use -early and -late simultaneously");
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
    CommandInterpreter::Consume<bool>(args, "-late", [](bool & b){
        if (b == false) throw std::invalid_argument("tried to use -early and -late simultaneously");
    }, this->late_fit);

    // bounding function
    CommandInterpreter::Consume<Options::BoundingFunction>(args, "-bf", [](const string & arg, Options::BoundingFunction & bf){
        bf = ToBoundingFunction(arg);
    }, this->bounding_function);
}

DynamicSolver::Options::Options(std::vector<std::string> & args) : Options(){}




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
    PackagedProblem pp(Consume(args, 0));
    Consume<PackagedProblem>(args, "-structure", [](const string & arg, PackagedProblem & _p_p){
        if (_p_p.requirements.structureToFind == ToStructureToFind(arg)) return;
        _p_p.requirements.structureToFind = ToStructureToFind(arg);
        _p_p.associated_file = "";
        _p_p.known_optimum = -1;
    }, pp);
    Consume<PackagedProblem>(args, "-weight-treatment", [](const string & arg, PackagedProblem & _p_p){
        if (_p_p.requirements.weightTreatment == ToWeightTreatment(arg)) return;
        _p_p.requirements.weightTreatment = ToWeightTreatment(arg);
        _p_p.associated_file = "";
        _p_p.known_optimum = -1;
    }, pp);
    if (args.size() > 0) throw std::invalid_argument(args[0] + " is not a recognised argument for command " + "solve");
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
        Consume(args, "-o", [&export_file](const string & arg){
            export_file = arg;
        });
        
        // print on console
        bool print = false;
        Consume(args, "-p", [&print](){
            print = true;
        });

        PackagedSolution ps;
        string algorithm = Consume(args, 1);
        if (algorithm == "greedy") ps = SolveFromArgs<GreedySolver>(args);
        else if (algorithm == "greedy-heuristic-search") ps = SolveFromArgs<GreedyHeuristicSearchSolver>(args);
        else if (algorithm == "grasp") ps = SolveFromArgs<GRASPSolver>(args);
        else if (algorithm == "brute-force") ps = SolveFromArgs<BruteForceSolver>(args);
        else if (algorithm == "branch-and-bound") ps = SolveFromArgs<BranchAndBoundSolver>(args);
        else if (algorithm == "dynamic") ps = SolveFromArgs<DynamicSolver>(args);
        else throw std::invalid_argument(algorithm + " is not a recognised algorithm");
            
        if (print) cout << ps;
        if (export_file != "") ps.ExportJSON(export_file);

        break;
    }
    case Command::BATCH_SOLVE:{
        string algorithm = Consume(args, 1);
        
        /* would be nice
        bool fake;
        Consume(args, "-fake", [&fake](){
            fake = true;
        }); should create placeholder folder with algo name and solve-info.json
        */

        if (algorithm == "greedy") BatchSolveFromArgs<GreedySolver>(args);
        else if (algorithm == "greedy-heuristic-search") BatchSolveFromArgs<GreedyHeuristicSearchSolver>(args);
        else if (algorithm == "grasp") BatchSolveFromArgs<GRASPSolver>(args);
        else if (algorithm == "brute-force") BatchSolveFromArgs<BruteForceSolver>(args);
        else if (algorithm == "branch-and-bound") BatchSolveFromArgs<BranchAndBoundSolver>(args);
        else if (algorithm == "dynamic") BatchSolveFromArgs<DynamicSolver>(args);
        else throw std::invalid_argument(algorithm + " is not a recognised algorithm");
        break;
    }

    case Command::PRINT:{
        if (args.size() > 1) throw std::invalid_argument("print takes only one argument, provided: " + args.size());
        cout << Consume(args, 0);
        break;
    }

    case Command::GENERATE_PROBLEM:{
        string file_name;
        Problem::GenerationSettings gs;
        
        Consume(args, "-file", [&gs](const string & arg){
            gs = Problem::GenerationSettings(arg);
        });

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
            // weight treatment
            else if ((*it) == "-weight-treatment") rq.weightTreatment = ToWeightTreatment(*(++it));

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

    case Command::GATHER:{
        string field_name(Consume(args, 1));
        std::filesystem::path dir(Consume(args, 0));
        std::ofstream fout;
        if (std::filesystem::exists(dir.parent_path().string() + "/" + field_name + ".txt")) {
            fout.open(dir.parent_path().string() + "/" + field_name + ".txt", std::ios_base::app);
        }
        else {
            fout.open(dir.parent_path().string() + "/" + field_name + ".txt");
            fout << "^ N ^ ";
            for (auto const& dir_entry : std::filesystem::directory_iterator{dir}){
                if (!std::filesystem::is_directory(dir_entry)) continue;
                if (dir_entry.path().stem().string() == "problems") continue;
                fout << dir_entry.path().stem().string() << " ^ ";
            }
        }
        fout << endl << "| " << dir.stem().string() << " | ";
        for (auto const& dir_entry : std::filesystem::directory_iterator{dir}){
            if (!std::filesystem::is_directory(dir_entry)) continue;
            if (dir_entry.path().stem().string() == "problems") continue;
            for (auto const& file_entry : std::filesystem::directory_iterator{dir_entry.path()}){
                if (file_entry.path().filename().string() == "solve-info.json"){
                    std::ifstream fin (file_entry.path().string());
                    json data = json::parse(fin);
                    fout << data[field_name] << " | ";
                }
            }
        }
        fout.close();
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