#pragma once

#include <vector>
#include <iostream>

class Item {
    public:
    int value;
    std::vector<int> weights;
    std::vector<int> connections;
    
    Item(const int & value, const std::vector<int> & weights, const std::vector<int> & connections);
    friend std::ostream& operator<<(std::ostream & os, const Item & item);

    bool HasConnectionTo(const int & id) const;
    int GetWeightSum() const;
};

struct Requirements {
    enum class StructureToFind { PATH, CYCLE, TREE, CONNECTED_GRAPH, IGNORE_CONNECTIONS };
    StructureToFind structureToFind = StructureToFind::PATH;
    enum class WeightTreatment { IGNORE_ALL, RESPECT_ALL, RESPECT_FIRST_ONLY, SET_ALL_TO_1 };
    WeightTreatment weightTreatment = WeightTreatment::RESPECT_ALL;
    bool directed = true;
    int known_optimum = -1;
};

class Problem {
    private:
    void GenerateInnerItems(const int & instance_size, const int & sub_knackpacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density);

    public:
    std::vector<int> knapsack_sizes;
    std::vector<Item> items;

    Problem(const int & instance_size, const int & sub_knackpacks, const std::vector<int> fixed_knapsack_sizes, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density);
    Problem(const int & instance_size, const int & sub_knackpacks, const int & knapsack_size_limit_exclusive, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density);

    static std::vector<Item> GenerateItemsVector(const int & instance_size, const int & sub_knackpacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density);

    enum class SortMode {WEIGHT_VALUE_RATIO, WEIGHT, VALUE};
    // DON'T USE
    // leaves connections unchanged and this messes up the problem
    void SortItems(const SortMode & sortMode);
    std::vector<int> GetSortedItemIds(const SortMode & sortMode) const;
    int GetValueSum() const;
};

std::ostream& operator<<(std::ostream & os, const Problem & p);