#pragma once

#include <iostream> // std::ostream
#include <vector>
#include <string>
#include <algorithm> // std::sort, std::find
#include "DebugDefines.hpp"


namespace knapsack_solver{

class Item {
    public:
    int value;
    std::vector<int> weights;
    std::vector<int> connections;
    
    Item(const int & value, const std::vector<int> & weights, const std::vector<int> & connections);

    bool HasConnectionTo(const int & id) const;
    int GetWeightSum() const;
};

class Problem {
    private:
    void GenerateInnerItems(const int & instance_size, const int & sub_knapsacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density);

    public:
    struct Requirements {
        enum class StructureToFind { CYCLE, PATH, TREE, CONNECTED_GRAPH, IGNORE_CONNECTIONS };
        StructureToFind structureToFind = StructureToFind::PATH;

        enum class WeightTreatment { IGNORE_ALL, RESPECT_ALL, RESPECT_FIRST_ONLY, SET_ALL_TO_1 };
        WeightTreatment weightTreatment = WeightTreatment::RESPECT_ALL;
    };

    std::vector<int> knapsack_sizes;
    std::vector<Item> items;

    bool directed = true;
    enum class CycleGuarantee { CONTAINS_CYCLE, NO_CYCLE, NO_GUARANTEES };
    CycleGuarantee cycleGuarantee = CycleGuarantee::NO_GUARANTEES;

    struct GenerationSettings{
        int instance_size;
        int sub_knapsacks;
        bool randomize_knapsack_sizes;
        int knapsack_size_limit_exclusive;
        int value_limit_exclusive;
        int weight_limit_exclusive;
        double connection_density;
        GenerationSettings() = default;
        GenerationSettings(const std::string & file_name);
    };

    Problem() = default;
    Problem(const std::string & file_name);
    Problem(const GenerationSettings & gs);

    void ExportJSON(const std::string & file_name) const;

    static std::vector<Item> GenerateItemsVector(const int & instance_size, const int & sub_knapsacks, const int & value_limit_exclusive, const int & weight_limit_exclusive, const double & connection_density);

    /// @brief creates one .json file named in directory and name specified in <file_name>
    /// @param gs GenerationSettings struct found in Problem::GenerationSettings
    /// @param file_name full path and file name with file extension
    /// @deprecated use PackagedProblem::GeneratePackagedProblemJSON()
    static void GenerateProblemJSON(const GenerationSettings & gs, const std::string file_name);

    /// @brief creates a folder named <batch_name> inside <directory_path> with a subfolder "problems"
    /// @param gs GenerationSettings struct found in Problem::GenerationSettings
    /// @param file_name general problem file name. x.json will be appended at the end, where x is the file number.
    /// @deprecated use PackagedProblem::BatchGeneratePackagedProblemsJSON()
    static void BatchGenerateProblemsJSON(const GenerationSettings & gs, const int & amount, const std::string & directory_path, const std::string & batch_name, const std::string & file_name);

    enum class SortMode {VALUE_WEIGHT_RATIO, WEIGHT, VALUE, RANDOM, DONT_SORT};
    
    template<typename T> T GetSortedItemIds(const SortMode & sortMode) const;
    int GetValueSum() const;
};

template<typename T>
inline T Problem::GetSortedItemIds(const SortMode & sortMode) const{
    struct el{
        int id;
        int value;
        int weight;
    };
    std::vector<el> toSort;
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
        T sortedIds;
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

    T sortedIds;
    for (int i = 0; i < toSort.size(); ++i){
        sortedIds.push_back(toSort[i].id);
    }
    return sortedIds;
}

class PackagedProblem{
    public:
    Problem problem;
    Problem::Requirements requirements;
    int known_optimum;
    std::string associated_file;
    
    PackagedProblem(const std::string & file_name);
    PackagedProblem(const Problem::GenerationSettings & gs, const Problem::Requirements & rq);

    void ExportJSON(const std::string & file_name) const;

    /// @brief creates one .json file named in directory and name specified in <file_name>
    /// @param file_name full path and file name with file extension
    static void GeneratePackagedProblemJSON(const Problem::GenerationSettings & gs, const Problem::Requirements & rq, const std::string file_name);
    static void BatchGeneratePackagedProblemsJSON(const Problem::GenerationSettings & gs, const Problem::Requirements & rq, const int & amount, const std::string & directory_path);
};

} // namesapace knapsack_solver

std::ostream& operator<<(std::ostream & os, const knapsack_solver::Item & item);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::Problem & p);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::PackagedProblem & pp);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::Problem::CycleGuarantee & cg);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::Problem::Requirements::StructureToFind & stf);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::Problem::Requirements::WeightTreatment & wt);
std::ostream& operator<<(std::ostream & os, const knapsack_solver::Problem::SortMode & sm);

std::string ToString(const knapsack_solver::Problem::CycleGuarantee & cg);
std::string ToString(const knapsack_solver::Problem::Requirements::StructureToFind & stf);
std::string ToString(const knapsack_solver::Problem::Requirements::WeightTreatment & wt);
std::string ToString(const knapsack_solver::Problem::SortMode & sm);

knapsack_solver::Problem::CycleGuarantee ToCycleGuarantee(const std::string & str);
knapsack_solver::Problem::Requirements::StructureToFind ToStructureToFind(const std::string & str);
knapsack_solver::Problem::Requirements::WeightTreatment ToWeightTreatment(const std::string & str);
knapsack_solver::Problem::SortMode ToSortMode(const std::string & str);