#include "Problem.hpp"

#include <algorithm> // std::sort, std::find
#include <fstream> // std::ifstream, std::ofstream
#include <filesystem> // std::filesystem::create_directory, std::filesystem::remove_all

#include "UtilityFunctions.hpp"
#include "json.hpp"
#include "FileNameDefines.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using nlohmann::json;
using namespace knapsack_solver;


//--------- ITEM ----------

Item::Item(const int & _value, const vector<int> & _weights, const vector<int> & _connections) : value(_value), weights(_weights), connections(_connections) { }

bool Item::HasConnectionTo(const int & id) const{
    for (int i = 0; i < this->connections.size(); ++i){
        if (this->connections[i] == id) return true;
    }
    return false;
}

int Item::GetWeightSum() const{
    int sum = 0;
    for (int i = 0; i < this->weights.size(); ++i){
        sum += this->weights[i];
    }
    return sum;
}

ostream& operator<<(ostream & os, const Item & item){
    os << "item: " << &item << ", value: " << item.value /*<< ", weights: " << item.weights.size() << ", connections: " << item.connections.size()*/ << endl;
    os << "weights:";
    for (int i = 0; i < item.weights.size(); ++i) os << " " << item.weights[i];
    os << endl << "connections: ";
    for (int i = 0; i < item.connections.size(); ++i) os << " " << item.connections[i];
    return (os << endl);
}




//--------- PROBLEM ----------

Problem::Problem(const string & file_name){
    std::ifstream fin (file_name);
    json data = json::parse(fin);
    knapsack_sizes = data["knapsack_sizes"].get<vector<int>>();
    for (auto & j_item : data["items"]){
        Item i(j_item["value"], j_item["weights"].get<vector<int>>(), j_item["connections"].get<vector<int>>());
        this->items.push_back(i);
    }
    fin.close();
}

Problem::Problem(const GenerationSettings & gs) {
    if (gs.randomize_knapsack_sizes){
        for (int i = 0; i < gs.sub_knapsacks; ++i){
            this->knapsack_sizes.push_back(std::rand() % gs.knapsack_size_limit_exclusive);
        }
    }
    else
        this->knapsack_sizes = gs.fixed_knapsack_sizes;
    this->GenerateInnerItems(gs.instance_size, gs.sub_knapsacks, gs.value_limit_exclusive, gs.weight_limit_exclusive, gs.connection_density);
}

void Problem::ExportJSON(const string & file_name) const{
    json data;
    data["knapsack_sizes"] = knapsack_sizes;
    for (int i = 0; i < items.size(); ++i){
        data["items"][i] = { {"value", items[i].value}, {"weights", items[i].weights}, {"connections", items[i].connections} };
    }
    std::ofstream fout(file_name);
    fout << data.dump(4);
    fout.close();
}

void Problem::GenerateProblemJSON(const GenerationSettings & gs, const std::string file_name){
    Problem p(gs);
    p.ExportJSON(file_name);
}

void Problem::BatchGenerateProblemsJSON(const GenerationSettings & gs, const int & amount, const std::string & directory_path, const std::string & batch_name, const std::string & file_name){
    std::filesystem::remove_all(directory_path + "/" + batch_name);
    std::filesystem::create_directories(directory_path + "/" + batch_name + "/problems");
    json data;
    data["batch_name"] = batch_name;
    data["amount"] = amount;
    std::ofstream fout(directory_path + "/" + batch_name + "/" + batch_name + ".json");
    fout << data.dump(4);
    fout.close();
    for (int i = 0; i < amount; ++i){
        GenerateProblemJSON(gs, directory_path + "/" + batch_name + "/problems/" +  file_name + std::to_string(i) + ".json");
    }
}

vector<Item> Problem::GenerateItemsVector(const int & instance_size, const int & sub_knapsacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density){
    vector<Item> out;
    for(int i = 0; i < instance_size; ++i){
        vector<int> w, c;
        for (int j = 0; j < sub_knapsacks; ++j) w.push_back(rand() % weight_limit_exclusive);

        for (int j = 0; j < instance_size; ++j){
            if (j == i) continue;
            if (RandomT<double>(0.0, 1.0) <= connection_density) c.push_back(j);
        }
        Item item(rand() % value_limit_exclusive, w, c);
        out.push_back(item);
    }
    return out;
}

void Problem::GenerateInnerItems(const int & instance_size, const int & sub_knapsacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density) {
    items.clear();
    for(int i = 0; i < instance_size; ++i){
        vector<int> w, c;
        for (int j = 0; j < sub_knapsacks; ++j) w.push_back(rand() % weight_limit_exclusive);

        for (int j = 0; j < instance_size; ++j){
            if (j == i) continue;
            if (RandomT<double>(0.0, 1.0) <= connection_density) c.push_back(j);
        }
        Item item(rand() % value_limit_exclusive, w, c);
        items.push_back(item);
    }
}

vector<int> Problem::GetSortedItemIds(const SortMode & sortMode) const{
    struct el{
        int id;
        int value;
        int weight;
    };
    vector<el> toSort;
    for (int i = 0; i < this->items.size(); ++i){
        el e = {i, items[i].value, items[i].GetWeightSum()};
        toSort.push_back(e);
    }

    switch (sortMode)
    {
    case SortMode::VALUE_WEIGHT_RATIO:
        std::sort(toSort.begin(), toSort.end(), [](el a, el b){
            if (a.weight == 0 && b.weight == 0){
                return a.value > b.value;
            }
            else if (a.weight == 0) return true;
            else if (b.weight == 0) return false;
            return (double)a.value / (double)a.weight > (double)b.value / (double)b.weight;
        });
        break;
        
    case SortMode::WEIGHT:
        std::sort(toSort.begin(), toSort.end(), [](el a, el b){
            return a.weight < b.weight;
        });
        break;
        
    case SortMode::VALUE:
        std::sort(toSort.begin(), toSort.end(), [](el a, el b){
            return a.value > b.value;
        });
        break;
    
    case SortMode::RANDOM:
    {
        vector<int> sortedIds;
        while (sortedIds.size() < this->items.size())
        {
            int itemId = std::rand() % this->items.size();
            while (std::find(sortedIds.begin(), sortedIds.end(), itemId) != sortedIds.end()){
                itemId = std::rand() % this->items.size();
            }
            sortedIds.push_back(itemId);
        }
        return sortedIds;
        break;
    }

    default:
        break;
    }

    vector<int> sortedIds;
    for (int i = 0; i < toSort.size(); ++i){
        sortedIds.push_back(toSort[i].id);
    }
    return sortedIds;
}

int Problem::GetValueSum() const{
    int sum = 0;
    for (int i = 0; i < items.size(); ++i){
        sum += items[i].value;
    }
    return sum;
}

ostream& operator<<(ostream & os, const Problem & p) {
    int s = p.knapsack_sizes.size();
    os << "problem: " << &p << ", sub knapsacks: " << s << ", knapsack sizes:";
    for (int i = 0; i < s; ++i) os << " " << p.knapsack_sizes[i];
    s = p.items.size();
    os << endl << "items amount: " << s << ", items:\n";
    for (int i = 0; i < s; ++i) os << "[" << i << "] " << p.items[i];
    return os;
}




//--------- PACKAGED PROBLEM ----------

PackagedProblem::PackagedProblem(const string & file_name) : associated_file(file_name) {
    std::ifstream fin (file_name);
    json data = json::parse(fin);

    // problem
    problem.knapsack_sizes = data["knapsack_sizes"].get<vector<int>>();
    for (auto & j_item : data["items"]){
        Item i(j_item["value"], j_item["weights"].get<vector<int>>(), j_item["connections"].get<vector<int>>());
        problem.items.push_back(i);
    }

    // requirements
    requirements.structureToFind = data["structure_to_find"];
    requirements.weightTreatment = data["weight_treatment"];

    // known optimum
    known_optimum = data["optimum"];

    fin.close();
}

PackagedProblem::PackagedProblem(const Problem::GenerationSettings & gs, const Problem::Requirements & rq) : problem(gs), requirements(rq), known_optimum(-1), associated_file("") { }

void PackagedProblem::ExportJSON(const string & file_name) const{
    json data;
    data["knapsack_sizes"] = this->problem.knapsack_sizes;
    for (int i = 0; i < this->problem.items.size(); ++i){
        data["items"][i] = { {"value", this->problem.items[i].value}, {"weights", this->problem.items[i].weights}, {"connections", this->problem.items[i].connections} };
    }
    data["directed"] = this->problem.directed;
    data["cycle_guarantee"] = this->problem.cycleGuarantee;
    data["structure_to_find"] = this->requirements.structureToFind;
    data["weight_treatment"] = this->requirements.weightTreatment;
    data["optimum"] = this->known_optimum;
    std::ofstream fout(file_name);
    fout << data.dump(4);
    fout.close();
}

void PackagedProblem::GeneratePackagedProblemJSON(const Problem::GenerationSettings & gs, const Problem::Requirements & rq, const string file_name){
    PackagedProblem pp(gs, rq);
    pp.ExportJSON(file_name);
}

void PackagedProblem::BatchGeneratePackagedProblemsJSON(const Problem::GenerationSettings & gs, const Problem::Requirements & rq, const int & amount, const string & directory_path){
    std::filesystem::remove_all(directory_path);
    std::filesystem::create_directories(directory_path + FND_PROBLEMS_FOLDER);
    for (int i = 0; i < amount; ++i) {
        GeneratePackagedProblemJSON(gs, rq, directory_path + FND_PROBLEMS_FOLDER + FND_PROBLEM_FILE + std::to_string(i) + ".json");
    }
    
    // create info file in the directory
    json data;
    data["amount"] = amount;
    std::ofstream fout(directory_path + FND_BATCH_INFO_FILE);
    fout << data.dump(4);
    fout.close();
}

ostream& operator<<(ostream & os, const PackagedProblem & pp) {
    os << pp.problem;
    os << "structure to find: " << pp.requirements.structureToFind << endl;
    os << "weight treatment: " << pp.requirements.weightTreatment << endl;
    os << "optimum: ";
    pp.known_optimum >= 0 ? os << pp.known_optimum : os << "unknown";
    return (os << endl);
}




//--------- ENUMS ----------

namespace knapsack_solver{

NLOHMANN_JSON_SERIALIZE_ENUM( Problem::Requirements::StructureToFind, {
    {Problem::Requirements::StructureToFind::CYCLE, "cycle"},
    {Problem::Requirements::StructureToFind::PATH, "path"},
    {Problem::Requirements::StructureToFind::TREE, "tree"},
    {Problem::Requirements::StructureToFind::CONNECTED_GRAPH, "connected-graph"},
    {Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS, "ignore-connections"}
})

NLOHMANN_JSON_SERIALIZE_ENUM( Problem::CycleGuarantee, {
    {Problem::CycleGuarantee::CONTAINS_CYCLE, "contains-cycle"},
    {Problem::CycleGuarantee::NO_CYCLE, "no-cycle"},
    {Problem::CycleGuarantee::NO_GUARANTEES, "no-guarantees"}
})

NLOHMANN_JSON_SERIALIZE_ENUM( Problem::Requirements::WeightTreatment, {
    {Problem::Requirements::WeightTreatment::IGNORE_ALL, "ignore-all"},
    {Problem::Requirements::WeightTreatment::RESPECT_ALL, "respect-all"},
    {Problem::Requirements::WeightTreatment::RESPECT_FIRST_ONLY, "respect-first-only"},
    {Problem::Requirements::WeightTreatment::SET_ALL_TO_1, "set-all-to-1"}
})

}

std::string ToString(const Problem::CycleGuarantee & cg){
    switch(cg){
        case Problem::CycleGuarantee::CONTAINS_CYCLE:
        return "contains-cycle";
        break;
        case Problem::CycleGuarantee::NO_CYCLE:
        return "no-cycle";
        break;
        case Problem::CycleGuarantee::NO_GUARANTEES:
        return "no-guarantees";
        break;
        default:
        throw std::invalid_argument("unrecognised value");
        break;
    }

}
string ToString(const Problem::Requirements::StructureToFind & stf){
    switch(stf){
        case Problem::Requirements::StructureToFind::CYCLE:
        return "cycle";
        break;
        case Problem::Requirements::StructureToFind::PATH:
        return "path";
        break;
        case Problem::Requirements::StructureToFind::TREE:
        return "tree";
        break;
        case Problem::Requirements::StructureToFind::CONNECTED_GRAPH:
        return "connected-graph";
        break;
        case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
        return "ignore-connections";
        break;
        default:
        throw std::invalid_argument("unrecognised value");
        break;
    }
}
string ToString(const Problem::Requirements::WeightTreatment & wt){
    switch(wt){
        case Problem::Requirements::WeightTreatment::IGNORE_ALL:
        return "ignore-all";
        break;
        case Problem::Requirements::WeightTreatment::RESPECT_ALL:
        return "respect-all";
        break;
        case Problem::Requirements::WeightTreatment::RESPECT_FIRST_ONLY:
        return "respect-first_only";
        break;
        case Problem::Requirements::WeightTreatment::SET_ALL_TO_1:
        return "set-all-to-1";
        break;
        default:
        throw std::invalid_argument("unrecognised value");
        break;
    }
}
string ToString(const knapsack_solver::Problem::SortMode & sm){
    switch (sm)
    {
    case Problem::SortMode::VALUE_WEIGHT_RATIO:
        return "value/weight";
        break;
    case Problem::SortMode::VALUE:
        return "value";
        break;
    case Problem::SortMode::WEIGHT:
        return "weight";
        break;
    case Problem::SortMode::RANDOM:
        return "random";
        break;
    case Problem::SortMode::DONT_SORT:
        return "dont-sort";
        break;
        default:
        throw std::invalid_argument("unrecognised value");
        break;
    }
}

knapsack_solver::Problem::CycleGuarantee ToCycleGuarantee(const std::string & str){
    if (str == "contains-cycle") return Problem::CycleGuarantee::CONTAINS_CYCLE;
    else if (str == "no-cycle") return Problem::CycleGuarantee::NO_CYCLE;
    else if (str == "no-guarantees") return Problem::CycleGuarantee::NO_GUARANTEES;
    else throw std::invalid_argument("unrecognised value");

}
knapsack_solver::Problem::Requirements::StructureToFind ToStructureToFind(const std::string & str){
    if (str == "cycle") return Problem::Requirements::StructureToFind::CYCLE;
    else if (str == "path") return Problem::Requirements::StructureToFind::PATH;
    else if (str == "tree") return Problem::Requirements::StructureToFind::TREE;
    else if (str == "connected-graph") return Problem::Requirements::StructureToFind::CONNECTED_GRAPH;
    else if (str == "ignore-connections") return Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS;
    else throw std::invalid_argument("unrecognised value");

}
knapsack_solver::Problem::Requirements::WeightTreatment ToWeightTreatment(const std::string & str){
    if (str == "ignore-all") return Problem::Requirements::WeightTreatment::IGNORE_ALL;
    else if (str == "respect-all") return Problem::Requirements::WeightTreatment::RESPECT_ALL;
    else if (str == "respect-first_only") return Problem::Requirements::WeightTreatment::RESPECT_FIRST_ONLY;
    else if (str == "set-all-to-1") return Problem::Requirements::WeightTreatment::SET_ALL_TO_1;
    else throw std::invalid_argument("unrecognised value");

}
knapsack_solver::Problem::SortMode ToSortMode(const std::string & str){
    if (str == "value/weight") return Problem::SortMode::VALUE_WEIGHT_RATIO;
    else if (str == "value") return Problem::SortMode::VALUE;
    else if (str == "weight") return Problem::SortMode::WEIGHT;
    else if (str == "random") return Problem::SortMode::RANDOM;
    else if (str == "dont-sort") return Problem::SortMode::DONT_SORT;
    else throw std::invalid_argument("unrecognised value");
}

std::ostream& operator<<(std::ostream & os, const Problem::CycleGuarantee & cg){
    return os << ToString(cg);
}
std::ostream& operator<<(std::ostream & os, const Problem::Requirements::StructureToFind & stf){
    return os << ToString(stf);
}
std::ostream& operator<<(std::ostream & os, const Problem::Requirements::WeightTreatment & wt){
    return os << ToString(wt);
}
std::ostream& operator<<(std::ostream & os, const Problem::SortMode & sm){
    return os << ToString(sm);
}