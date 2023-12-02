#include "Problem.hpp"
#include "UtilityFunctions.hpp"
#include "json.hpp"

#include <algorithm> // std::sort
#include <fstream> // std::ifstream, std::ofstream
#include <filesystem> // std::filesystem::create_directory

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using nlohmann::json;


//---------Item----------

Item::Item(const int & _value, const vector<int> & _weights, const vector<int> & _connections) : value(_value), weights(_weights), connections(_connections) { }

ostream& operator<<(ostream & os, const Item & item){
    os << "item: " << &item << ", value: " << item.value /*<< ", weights: " << item.weights.size() << ", connections: " << item.connections.size()*/ << endl;
    os << "weights:";
    for (int i = 0; i < item.weights.size(); ++i) os << " " << item.weights[i];
    os << endl << "connections: ";
    for (int i = 0; i < item.connections.size(); ++i) os << " " << item.connections[i];
    return (os << endl);
}

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


//---------Problem----------

//---------Constructors----------

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

Problem::Problem(const int & instance_size, const int & sub_knapsacks, const vector<int> fixed_knapsack_sizes, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density) : knapsack_sizes(fixed_knapsack_sizes) {
    this->GenerateInnerItems(instance_size, sub_knapsacks, value_limit_exclusive, weight_limit_exclusive, connection_density);
}

Problem::Problem(const int & instance_size, const int & sub_knapsacks, const int & knapsack_size_limit_exclusive, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density) {
    for (int i = 0; i < sub_knapsacks; ++i) {
        this->knapsack_sizes.push_back(rand() % knapsack_size_limit_exclusive);
    }
    this->GenerateInnerItems(instance_size, sub_knapsacks, value_limit_exclusive, weight_limit_exclusive, connection_density);
}


//---------Static Methods----------

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

//---------Methods----------

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

void Problem::SortItems(const SortMode & sortMode){
    switch (sortMode)
    {
    case SortMode::WEIGHT_VALUE_RATIO:
        std::sort(items.begin(), items.end(), [](Item a, Item b){
            double aw = a.GetWeightSum();
            double bw = b.GetWeightSum();
            if (aw == 0 && bw == 0){
                return a.value > b.value;
            }
            else if (aw == 0) return true;
            else if (bw == 0) return false;
            return (double)a.value / aw > (double)b.value / bw;
        });
        break;
        
    case SortMode::WEIGHT:
        std::sort(items.begin(), items.end(), [](Item a, Item b){
            return a.GetWeightSum() < b.GetWeightSum();
        });
        break;
        
    case SortMode::VALUE:
        std::sort(items.begin(), items.end(), [](Item a, Item b){
            return a.value > b.value;
        });
        break;
    
    default:
        break;
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


//---------Operators----------

ostream& operator<<(ostream & os, const Problem & p) {
    int s = p.knapsack_sizes.size();
    os << "problem: " << &p << ", sub knapsacks: " << s << ", knapsack sizes:";
    for (int i = 0; i < s; ++i) os << " " << p.knapsack_sizes[i];
    s = p.items.size();
    os << endl << "items amount: " << s << ", items:\n";
    for (int i = 0; i < s; ++i) os << "[" << i << "] " << p.items[i];
    return (os << endl);
}