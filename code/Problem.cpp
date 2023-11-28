#include "Problem.hpp"
#include "UtilityFunctions.hpp"

#include <algorithm> // sort

using std::vector;
using std::ostream;
using std::endl;


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

Problem::Problem(const int & instance_size, const int & sub_knackpacks, const vector<int> fixed_knapsack_sizes, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density) : knapsack_sizes(fixed_knapsack_sizes) {
    this->GenerateInnerItems(instance_size, sub_knackpacks, value_limit_exclusive, weight_limit_exclusive, connection_density);
}

Problem::Problem(const int & instance_size, const int & sub_knackpacks, const int & knapsack_size_limit_exclusive, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density) {
    for (int i = 0; i < sub_knackpacks; ++i) {
        this->knapsack_sizes.push_back(rand() % knapsack_size_limit_exclusive);
    }
    this->GenerateInnerItems(instance_size, sub_knackpacks, value_limit_exclusive, weight_limit_exclusive, connection_density);
}


//---------Methods----------

void Problem::GenerateInnerItems(const int & instance_size, const int & sub_knackpacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density) {
    items.clear();
    for(int i = 0; i < instance_size; ++i){
        vector<int> w, c;
        for (int j = 0; j < sub_knackpacks; ++j) w.push_back(rand() % weight_limit_exclusive);

        for (int j = 0; j < instance_size; ++j){
            if (j == i) continue;
            if (RandomT<double>(0.0, 1.0) <= connection_density) c.push_back(j);
        }
        Item item(rand() % value_limit_exclusive, w, c);
        items.push_back(item);
    }
}

vector<Item> Problem::GenerateItemsVector(const int & instance_size, const int & sub_knackpacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density){
    vector<Item> out;
    for(int i = 0; i < instance_size; ++i){
        vector<int> w, c;
        for (int j = 0; j < sub_knackpacks; ++j) w.push_back(rand() % weight_limit_exclusive);

        for (int j = 0; j < instance_size; ++j){
            if (j == i) continue;
            if (RandomT<double>(0.0, 1.0) <= connection_density) c.push_back(j);
        }
        Item item(rand() % value_limit_exclusive, w, c);
        out.push_back(item);
    }
    return out;
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
    os << "problem: " << &p << ", sub knackpacks: " << s << ", knapsack sizes:";
    for (int i = 0; i < s; ++i) os << " " << p.knapsack_sizes[i];
    s = p.items.size();
    os << endl << "items amount: " << s << ", items:\n";
    for (int i = 0; i < s; ++i) os << "[" << i << "] " << p.items[i];
    return (os << endl);
}