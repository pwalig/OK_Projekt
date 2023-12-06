#include "Problem.hpp"
#include "UtilityFunctions.hpp"
#include "json.hpp"

#include <algorithm> // std::sort, std::find
#include <fstream> // std::ifstream, std::ofstream
#include <filesystem> // std::filesystem::create_directory

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
    if (gs.randomize_knapsack_sizes)
        for (int i = 0; i < gs.sub_knapsacks; ++i)
            this->knapsack_sizes.push_back(rand() % gs.knapsack_size_limit_exclusive);
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
    case SortMode::WEIGHT_VALUE_RATIO:
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

PackagedProblem::PackagedProblem(const string & file_name) : known_optimum(-1) {
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

    fin.close();
}

PackagedProblem::PackagedProblem(const Problem::GenerationSettings & gs, const Problem::Requirements & rq) : problem(gs), requirements(rq), known_optimum(-1) { }

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
    data["known_optimum"] = this->known_optimum;
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
    std::filesystem::create_directories(directory_path + "/problems");
    for (int i = 0; i < amount; ++i) {
        GeneratePackagedProblemJSON(gs, rq, directory_path + "/problems/packaged_problem_" + std::to_string(i) + ".json");
    }
    
    // create info file in the directory
    json data;
    data["amount"] = amount;
    std::ofstream fout(directory_path + "/batch-info.json");
    fout << data.dump(4);
    fout.close();
}

ostream& operator<<(ostream & os, const PackagedProblem & pp) {
    os << pp.problem;
    os << "structure to find: " << pp.requirements.structureToFind << endl;
    os << "weight treatment: " << pp.requirements.weightTreatment << endl;
    os << "known_optimum: " << pp.known_optimum << endl;
    return (os << endl);
}




//--------- ENUMS ----------

namespace knapsack_solver{

NLOHMANN_JSON_SERIALIZE_ENUM( Problem::Requirements::StructureToFind, {
    {Problem::Requirements::StructureToFind::CYCLE, "cycle"},
    {Problem::Requirements::StructureToFind::PATH, "path"},
    {Problem::Requirements::StructureToFind::TREE, "tree"},
    {Problem::Requirements::StructureToFind::CONNECTED_GRAPH, "connected_graph"},
    {Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS, "ignore_connections"}
})

NLOHMANN_JSON_SERIALIZE_ENUM( Problem::CycleGuarantee, {
    {Problem::CycleGuarantee::CONTAINS_CYCLE, "contains_cycle"},
    {Problem::CycleGuarantee::NO_CYCLE, "no_cycle"},
    {Problem::CycleGuarantee::NO_GUARANTEES, "no_guarantees"}
})

NLOHMANN_JSON_SERIALIZE_ENUM( Problem::Requirements::WeightTreatment, {
    {Problem::Requirements::WeightTreatment::IGNORE_ALL, "ignore_all"},
    {Problem::Requirements::WeightTreatment::RESPECT_ALL, "respect_all"},
    {Problem::Requirements::WeightTreatment::RESPECT_FIRST_ONLY, "respect_first_only"},
    {Problem::Requirements::WeightTreatment::SET_ALL_TO_1, "set_all_to_1"}
})

}

std::ostream& operator<<(std::ostream & os, const Problem::CycleGuarantee & cg){
    switch(cg){
        case Problem::CycleGuarantee::CONTAINS_CYCLE:
        os << "contains_cycle";
        break;
        case Problem::CycleGuarantee::NO_CYCLE:
        os << "no_cycle";
        break;
        case Problem::CycleGuarantee::NO_GUARANTEES:
        os << "no_guarantees";
        break;
    }
    return os;
}
std::ostream& operator<<(std::ostream & os, const Problem::Requirements::StructureToFind & stf){
    switch(stf){
        case Problem::Requirements::StructureToFind::CYCLE:
        os << "cycle";
        break;
        case Problem::Requirements::StructureToFind::PATH:
        os << "path";
        break;
        case Problem::Requirements::StructureToFind::TREE:
        os << "tree";
        break;
        case Problem::Requirements::StructureToFind::CONNECTED_GRAPH:
        os << "connected_graph";
        break;
        case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
        os << "ignore_connections";
        break;
    }
    return os;
}
std::ostream& operator<<(std::ostream & os, const Problem::Requirements::WeightTreatment & wt){
    switch(wt){
        case Problem::Requirements::WeightTreatment::IGNORE_ALL:
        os << "ignore_all";
        break;
        case Problem::Requirements::WeightTreatment::RESPECT_ALL:
        os << "respect_all";
        break;
        case Problem::Requirements::WeightTreatment::RESPECT_FIRST_ONLY:
        os << "respect_first_only";
        break;
        case Problem::Requirements::WeightTreatment::SET_ALL_TO_1:
        os << "set_all_to_1";
        break;
    }
    return os;
}