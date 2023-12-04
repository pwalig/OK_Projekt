#include "KnapsackSolver.hpp"
#include "json.hpp"

#include <chrono>
#include <fstream> // std::ifstream, std::ofstream
#include <cmath> // std::max
#include <stdexcept> // throw error types

using namespace knapsack_solver;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using nlohmann::json;



//----------SOLUTION----------

Solution::Solution() : max_value(0) {}

Solution::Solution(int InstanceSize) : max_value(0) {
    for (int i = 0; i < InstanceSize; ++i) this->selected.push_back(false);
}

ostream& operator<<(ostream& os, const Solution& s){
    os << "max value: " << s.max_value << "\nselected: ";
    int instance_size = s.selected.size();
    for (int i = 0; i < instance_size; ++i){
        os << (s.selected[i] ? 1 : 0) << " ";
    }
    return (os << endl);
}

void Solution::AddItem(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace){
    for (int j = 0; j < remainingSpace.size(); ++j){
        remainingSpace[j] -= problem.items[selected_item_id].weights[j];
        if (remainingSpace[j] < 0) {
            throw std::invalid_argument("item does not fit");
        }
    }
    this->max_value += problem.items[selected_item_id].value; // update max value
    if (this->selected[selected_item_id] == true) throw std::invalid_argument("item was already in the solution");
    this->selected[selected_item_id] = true; // add item to solution
}

void Solution::RemoveItem(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace){
    for (int j = 0; j < remainingSpace.size(); ++j){
        remainingSpace[j] += problem.items[selected_item_id].weights[j];
    }
    this->max_value -= problem.items[selected_item_id].value; // update max value
    if (this->selected[selected_item_id] == false) throw std::invalid_argument("item was not in the solution");
    this->selected[selected_item_id] = false; // remove item from solution
}

bool Solution::AddItemIfFits(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace){
    if (selected[selected_item_id]) return false; // don't add item if it is already in the solution

    for (int j = 0; j < remainingSpace.size(); ++j){
        if (remainingSpace[j] < problem.items[selected_item_id].weights[j]) return false;
    }

    for (int j = 0; j < remainingSpace.size(); ++j){
        remainingSpace[j] -= problem.items[selected_item_id].weights[j];
    }
    this->max_value += problem.items[selected_item_id].value; // update max value
    this->selected[selected_item_id] = true; // add item to solution
    return true;
}




//---------- PACKAGED SOLUTION ----------

void PackagedSolution::ExportJSON(const std::string file_name){
    json data;
    data["algorithm"] = algorithm;
    to_optimum_ratio >= 0 ? data["to_optimum_ratio"] = to_optimum_ratio : data["to_optimum_ratio"] = "optimum_unnown";
    data["solve_time"] = solve_time;
    data["remaining_spaces"] = remainingSpaces;
    data["solution"]["max_value"] = solution.max_value;
    for (int i = 0; i < solution.selected.size(); ++i){
        data["solution"]["selected"][i] = solution.selected[i];
    }
    std::ofstream fout(file_name);
    fout << data.dump(4);
    fout.close();
}

ostream& operator<<(ostream& os, const PackagedSolution& ps){
    os << "solved with: " << ps.algorithm << endl << ps.solution;
    os << "to_optimum_ratio: ";
    ps.to_optimum_ratio >= 0 ? os << ps.to_optimum_ratio : os << "optimum unknown";
    os << "\nsolve time: " << ps.solve_time << " s\nremaining spaces:";
    for (int i = 0; i < ps.remainingSpaces.size(); ++i){
        os << " " << ps.remainingSpaces[i];
    }
    return (os << endl);
}




//---------- KNAPSACK SOLVER ----------

bool KnapsackSolver::Fits(const std::vector<int> & weights, const std::vector<int> & remainingSpace){
    for (int j = 0; j < remainingSpace.size(); ++j){
        if (remainingSpace[j] < weights[j]){
            return false;
        }
    }
    return true;
}

vector<int> KnapsackSolver::CalculateRemainingSpaces(const Solution & solution, const Problem & problem){
    vector<int> remainigSpaces = problem.knapsack_sizes;
    for (int i = 0; i < solution.selected.size(); ++i){
        if (solution.selected[i]) {
            for (int j = 0; j < remainigSpaces.size(); ++j){
                remainigSpaces[j] -= problem.items[i].weights[j];
            }
        }
    }
    return remainigSpaces;
}




//---------- GREEDY SOLVER ----------

PackagedSolution GreedySolver::Solve(const PackagedProblem & problem, const Options & options){
    PackagedSolution ps;
    std::chrono::steady_clock::time_point start, end;
    Solution s;

    ps.algorithm = "";

    switch (options.sort_mode)
    {
    case Problem::SortMode::WEIGHT_VALUE_RATIO:
        ps.algorithm += "greedy sort by value/weight";
        break;
        
    case Problem::SortMode::VALUE:
        ps.algorithm += "greedy sort by value";
        break;
        
    case Problem::SortMode::WEIGHT:
        ps.algorithm += "greedy sort by weight";
        break;

    case Problem::SortMode::RANDOM:
        ps.algorithm += "random fit";
        break;
    
    case Problem::SortMode::DONT_SORT:
        ps.algorithm += "first fit";
        break;

    default:
        break;
    }

    switch (problem.requirements.structureToFind)
    {
    case Problem::Requirements::StructureToFind::PATH:
        ps.algorithm += " path";
        break;
    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
        ps.algorithm += " ignore connections";
        break;
    
    default:
        break;
    }
    
    start = std::chrono::steady_clock::now();
    s = GreedyUniversal(problem, options);
    end = std::chrono::steady_clock::now();

    // construct packaged solution
    ps.solution = s;
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = s.max_value / problem.known_optimum;
    ps.remainingSpaces = KnapsackSolver::CalculateRemainingSpaces(s, problem.problem);

    return ps;
}

Solution GreedySolver::GreedyUniversal(const PackagedProblem & problem, const Options & options) {
    int isiz = problem.problem.items.size();
    Solution globalSolution(isiz); // create empty solution
    vector<int> remainingSpace = problem.problem.knapsack_sizes;
    vector<int> sortedItemIds = problem.problem.GetSortedItemIds(options.sort_mode);

    switch (problem.requirements.structureToFind)
    {

    case Problem::Requirements::StructureToFind::PATH:
    {
        // Find best starting point
        int startItem = -1;
        for (int currentItemId : sortedItemIds) {
            if (KnapsackSolver::Fits(problem.problem.items[currentItemId].weights, remainingSpace)){
                globalSolution.AddItem(problem.problem, currentItemId, remainingSpace); // if current item fits add it to the current solution
                startItem = currentItemId;
                break;
            }
        }
        if (startItem == -1) return globalSolution; // nothing fits

        bool going = true;
        while (going){
            going = false;
            // find next item
            for (auto it = sortedItemIds.begin(); it != sortedItemIds.end(); ++it) {
                int currentItemId = *it;
                if (!globalSolution.selected[currentItemId] &&
                    problem.problem.items[startItem].HasConnectionTo(currentItemId) &&
                    KnapsackSolver::Fits(problem.problem.items[currentItemId].weights, remainingSpace))
                {
                    globalSolution.AddItem(problem.problem, currentItemId, remainingSpace);
                    startItem = currentItemId;
                    going = true;
                    break;
                }
            }
        }
        break;
    }

    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
    {
        for (int i = 0; i < isiz; ++i){
            int currentItemId = sortedItemIds[i];

            // if current item fits add it to the current solution
            globalSolution.AddItemIfFits(problem.problem, currentItemId, remainingSpace);
        }
        break;
    }
    
    default:
    {
        cout << "unimplemented yet\n";
        break;
    }
    }

    return globalSolution;
}

Solution GreedySolver::GreedyIgnoreConnections(const Problem & problem, const Options & options){
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution
    vector<int> remainingSpace = problem.knapsack_sizes;
    vector<int> sortedItemIds = problem.GetSortedItemIds(options.sort_mode);

    for (int i = 0; i < isiz; ++i){
        int currentItemId = sortedItemIds[i];

        // if current item fits add it to the current solution
        globalSolution.AddItemIfFits(problem, currentItemId, remainingSpace);
    }

    return globalSolution;
}

Solution GreedySolver::GreedyPath(const Problem & problem, const Options & options) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution
    vector<int> remainingSpace = problem.knapsack_sizes;
    vector<int> sortedItemIds = problem.GetSortedItemIds(options.sort_mode);

    for (int i = 0; i < isiz; ){
        int currentItemId = sortedItemIds[i];

        if (KnapsackSolver::Fits(problem.items[currentItemId].weights, remainingSpace))
            globalSolution.AddItem(problem, currentItemId, remainingSpace); // if current item fits add it to the current solution
        else return globalSolution; // if does not fit return what have we managed to fit so far

        i = -1;
        int nextItem;
        do{
            if (++i >= isiz) return globalSolution; // if can't find next possible item return what we had so far
            nextItem = sortedItemIds[i];
        } while (globalSolution.selected[nextItem] || !problem.items[currentItemId].HasConnectionTo(nextItem));
    }

    return globalSolution;
}




//---------- BRANCH AND BOUND SOLVER ----------

PackagedSolution BranchAndBoundSolver::Solve(const PackagedProblem & problem, const Options & options){
    PackagedSolution ps;
    Solution s;
    std::chrono::steady_clock::time_point start, end;
    ps.algorithm = "branch and bound";

    // solve with measured time
    if (options.late_fit){
        start = std::chrono::steady_clock::now();
        s = BnBLateFitPath(problem.problem);
        end = std::chrono::steady_clock::now();
        ps.algorithm += " late fit";
    }
    else{
        start = std::chrono::steady_clock::now();
        s = BnBEarlyFitPath(problem.problem);
        end = std::chrono::steady_clock::now();
        ps.algorithm += " early fit";
    }
    const std::chrono::duration<double> elapsed_seconds{end - start};

    ps.algorithm += " path";
    ps.solution = s;
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = s.max_value / problem.known_optimum;
    ps.remainingSpaces = KnapsackSolver::CalculateRemainingSpaces(s, problem.problem);

    return ps;
}

//----------Bounding Functions----------

int BranchAndBoundSolver::GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode, Solution currentSolution, std::vector<int> remainingSpace){
    int isiz = problem.items.size();
    vector<int> sortedItemIds = problem.GetSortedItemIds(sortMode);

    for (int i = 0; i < isiz; ++i){
        int currentItemId = sortedItemIds[i];

        if (currentSolution.selected[currentItemId]) continue; // don't add items that were already in the solution

        currentSolution.AddItemIfFits(problem, currentItemId, remainingSpace);
    }

    return currentSolution.max_value;
}

//---------- Early Fit ----------

Solution BranchAndBoundSolver::DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, vector<int> remainingSpace, const int & lower_bound) {

#ifdef DEBUG_BNB
    cout << "check for " << currentItemId << ", remaining weight: "<< remainingSpace[0] << endl;
#endif

    Solution globalSolution = currentSolution;

    // construct solution further with recursive DFSLateFitPath calls
    for (int nextItemId : problem.items[currentItemId].connections){
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        if (KnapsackSolver::Fits(problem.items[nextItemId].weights, remainingSpace)){
            Solution tempSolution = currentSolution;
            vector<int> tempRemainingSpace = remainingSpace;

            tempSolution.AddItem(problem, nextItemId, tempRemainingSpace);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::WEIGHT_VALUE_RATIO, tempSolution, remainingSpace) > lower_bound*/) {
                tempSolution = DFSEarlyFitPath(problem, tempSolution, nextItemId, tempRemainingSpace, lower_bound);
                if (tempSolution.max_value > globalSolution.max_value) globalSolution = tempSolution;
            }
        }
    }

#ifdef DEBUG_BNB
    cout << "exit with: " << globalSolution.max_value << endl;
#endif
    return globalSolution; // return the best solution from all paths starting from currentItemId
}

Solution BranchAndBoundSolver::BnBEarlyFitPath(const Problem & problem) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){
        vector<int> remainingSpace = problem.knapsack_sizes;

        if (KnapsackSolver::Fits(problem.items[i].weights, remainingSpace)){
            Solution currentSolution(isiz); // create empty solution
            currentSolution.AddItem(problem, i, remainingSpace);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::WEIGHT_VALUE_RATIO, currentSolution, remainingSpace) > globalSolution.max_value*/){
                currentSolution = DFSEarlyFitPath(problem, currentSolution, i, remainingSpace, globalSolution.max_value);
                if(currentSolution.max_value > globalSolution.max_value) globalSolution = currentSolution; // if newly found solution is better use it as you global solution
            }
        }
    }

    return globalSolution;
}

//---------- Late Fit ----------

Solution BranchAndBoundSolver::DFSLateFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, vector<int> remainingSpace) {

    if (KnapsackSolver::Fits(problem.items[currentItemId].weights, remainingSpace))
        currentSolution.AddItem(problem, currentItemId, remainingSpace); // if current item fits add it to the current solution
    else return currentSolution; // if does not fit return what have we managed to fit so far

#ifdef _DEBUG_MODE_
    cout << "check for " << currentItemId << ", remaining weight: "<< remainingSpace[0] << endl;
#endif

    
    Solution globalSolution(problem.items.size()); // create empty solution

    // construct solution further with recursive DFSLateFitPath calls
    for (int nextItemId : problem.items[currentItemId].connections){
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        Solution tempSolution = DFSLateFitPath(problem, currentSolution, nextItemId, remainingSpace);
        if (tempSolution.max_value > globalSolution.max_value) globalSolution = tempSolution;
    }

#ifdef _DEBUG_MODE_
    cout << "exit with: " << globalSolution.max_value << endl;
#endif
    return globalSolution; // return the best solution from all paths starting from currentItemId
}

Solution BranchAndBoundSolver::BnBLateFitPath(const Problem & problem) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){
        Solution currentSolution(isiz); // create empty solution 

        
        currentSolution = DFSLateFitPath(problem, currentSolution, i, problem.knapsack_sizes);
        if(currentSolution.max_value > globalSolution.max_value) globalSolution = currentSolution; // if newly found solution is better use it as you global solution
    }

    return globalSolution;
}




//---------- FLOYD SOLVER ----------

PackagedSolution FloydSolver::Solve(const PackagedProblem & problem){
    PackagedSolution ps;
    std::chrono::steady_clock::time_point start, end;
    Solution s;

    switch (problem.requirements.structureToFind)
    {
    case Problem::Requirements::StructureToFind::PATH:
        start = std::chrono::steady_clock::now();
        s = Connected(problem.problem);
        end = std::chrono::steady_clock::now();
        ps.algorithm = "floyd path";
        break;
    
    default:
        cout << "Not implemented yet\n";
        break;
    }

    // construct packaged solution
    ps.solution = s;
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = s.max_value / problem.known_optimum;
    ps.remainingSpaces = KnapsackSolver::CalculateRemainingSpaces(s, problem.problem);

    return ps;
}

//#define DEBUG_FLOYD_PATH

Solution FloydSolver::Connected(const Problem & problem){
    int isiz = problem.items.size(); // instance size

    struct SolutionProgress{
        Solution solution;
        vector<int> remainingSpace;
        SolutionProgress(const int & instance_size, const vector<int> & remaining) : solution(instance_size), remainingSpace(remaining){}
    };


    // Fill initial D vector
    vector<vector<SolutionProgress>> d;
    for (int i = 0; i < isiz; ++i) {
        vector<SolutionProgress> di;

        // if you can't fit element fill whole row with invalid answers
        if (!KnapsackSolver::Fits(problem.items[i].weights, problem.knapsack_sizes)) {
            SolutionProgress sp(isiz, problem.knapsack_sizes);
            sp.solution.max_value = -1;
            for (int j = 0; j < isiz; ++j) di.push_back(sp);
        }
        else {
            for (int j = 0; j < isiz; ++j){
                SolutionProgress sp(isiz, problem.knapsack_sizes);
                sp.solution.AddItem(problem, i, sp.remainingSpace); // i did fit, so add it

                if (!KnapsackSolver::Fits(problem.items[j].weights, problem.knapsack_sizes)) {
                    sp.solution.RemoveItem(problem, i, sp.remainingSpace);
                    sp.solution.max_value = -1; // item j does not fit - solution invalid
                }
                else if (problem.items[i].HasConnectionTo(j)) sp.solution.AddItem(problem, j, sp.remainingSpace);
                else {
                    sp.solution.RemoveItem(problem, i, sp.remainingSpace);
                    sp.solution.max_value = -2; // items lack connection - solution invalid
                }
                di.push_back(sp);
            }
        }
        d.push_back(di);
    }

#ifdef DEBUG_FLOYD_PATH
    for (int i = 0; i < isiz; ++i) {
        for (int j = 0; j < isiz; ++j){
            cout << d[i][j].solution.max_value << "\t";
        }
        cout << endl;
    }
#endif

    // Solve
    for (int k = 0; k < isiz; ++k){
        for (int i = 0; i < isiz; ++i){
            for (int j = 0; j < isiz; ++j){
                if (d[i][k].solution.max_value >= 0 && d[k][j].solution.max_value >= 0 && // if you can get from i to k and from k to j
                    !d[i][j].solution.selected[k] &&
                    d[i][j].solution.max_value != -1 && // items dont fit anyway
                    KnapsackSolver::Fits(problem.items[k].weights, d[i][j].remainingSpace) && // initial fit check
                    d[i][k].solution.max_value + d[k][j].solution.max_value > d[i][j].solution.max_value)
                {   
                    if (d[i][j].solution.max_value == -2){ // if items lacked connection but now have it through item k
                        d[i][j].solution.max_value = 0;
                        d[i][j].solution.AddItem(problem, i, d[i][j].remainingSpace);
                        if (i != j) d[i][j].solution.AddItem(problem, j, d[i][j].remainingSpace);
                        if (!KnapsackSolver::Fits(problem.items[k].weights, d[i][j].remainingSpace)) { // if k stopped to fit after adding i and j
                            d[i][j].solution.RemoveItem(problem, i, d[i][j].remainingSpace);
                            if (i != j) d[i][j].solution.RemoveItem(problem, j, d[i][j].remainingSpace);
                            d[i][j].solution.max_value = -2;
                            continue;
                        }
                        else d[i][j].solution.AddItem(problem, k, d[i][j].remainingSpace);
                    }
                    else d[i][j].solution.AddItem(problem, k, d[i][j].remainingSpace);
                }
            }
        }
    
#ifdef DEBUG_FLOYD_PATH
        cout << endl << "through: " << k << endl;
        for (int i = 0; i < isiz; ++i) {
            for (int j = 0; j < isiz; ++j){
                cout << d[i][j].solution.max_value << "\t";
            }
            cout << endl;
        }
#endif
    }

    // Find Best Solution
    int maxi, maxj, maxval = 0;
    for (int i = isiz - 1; i >= 0; --i) {
        for (int j = isiz - 1; j >= 0; --j){
            if (d[i][j].solution.max_value > maxval){
                maxval = d[i][j].solution.max_value;
                maxi = i;
                maxj = j;
            }
        }
    }

    return d[maxi][maxj].solution;
}