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

//----------Solution Definitions----------

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

void Solution::AddItemIfFits(const Problem & problem, const int & selected_item_id, std::vector<int> & remainingSpace){
    if (selected[selected_item_id]) return; // don't add item if it is already in the solution

    for (int j = 0; j < remainingSpace.size(); ++j){
        if (remainingSpace[j] < problem.items[selected_item_id].weights[j]) return;
    }

    for (int j = 0; j < remainingSpace.size(); ++j){
        remainingSpace[j] -= problem.items[selected_item_id].weights[j];
    }
    this->max_value += problem.items[selected_item_id].value; // update max value
    this->selected[selected_item_id] = true; // add item to solution
}


//----------Packaged Solution Definitions----------

void PackagedSolution::ExportJSON(const std::string file_name){
    json data;
    data["algorithm"] = algorithm;
    data["to_optimum_ratio"] = to_optimum_ratio;
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
    os << "to_optimum_ratio: " << ps.to_optimum_ratio << "\nsolve time: " << ps.solve_time << " s\nremaining spaces:";
    for (int i = 0; i < ps.remainingSpaces.size(); ++i){
        os << " " << ps.remainingSpaces[i];
    }
    return (os << endl);
}


//----------Solver Functions Definitions----------

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

PackagedSolution KnapsackSolver::Greedy(const Problem & problem, Requirements requirements){
    PackagedSolution ps;
    std::chrono::steady_clock::time_point start, end;
    Solution s;

    switch (requirements.structureToFind)
    {
    case Requirements::StructureToFind::PATH:
        start = std::chrono::steady_clock::now();
        s = GreedyPath(problem, Problem::SortMode::WEIGHT_VALUE_RATIO);
        end = std::chrono::steady_clock::now();
        ps.algorithm = "greedy path";
        break;
    case Requirements::StructureToFind::IGNORE_CONNECTIONS:
        start = std::chrono::steady_clock::now();
        s = GreedyIgnoreConnections(problem, Problem::SortMode::WEIGHT_VALUE_RATIO);
        end = std::chrono::steady_clock::now();
        ps.algorithm = "greedy ignore connections";
        break;
    
    default:
        cout << "Not implemented yet\n";
        break;
    }

    // construct packaged solution
    ps.solution = s;
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = s.max_value / requirements.known_optimum;
    ps.remainingSpaces = CalculateRemainingSpaces(s, problem);

    return ps;
}

PackagedSolution KnapsackSolver::BruteForce(const Problem & problem, Requirements requirements){

    // solve with measured time
    const auto start{std::chrono::steady_clock::now()};
    Solution s = BruteForcePath(problem);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};

    // construct packaged solution
    PackagedSolution ps;
    ps.algorithm = "brute force path";
    ps.solution = s;
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = s.max_value / requirements.known_optimum;
    ps.remainingSpaces = CalculateRemainingSpaces(s, problem);

    return ps;
}

PackagedSolution KnapsackSolver::BranchAndBound(const Problem & problem, Requirements requirements){

    // solve with measured time
    const auto start{std::chrono::steady_clock::now()};
    Solution s = BranchAndBoundPath(problem);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};

    // construct packaged solution
    PackagedSolution ps;
    ps.algorithm = "branch and bound path";
    ps.solution = s;
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = s.max_value / requirements.known_optimum;
    ps.remainingSpaces = CalculateRemainingSpaces(s, problem);

    return ps;
}


//----------GREEDY ALGORITHMS----------

Solution KnapsackSolver::GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode){
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution
    vector<int> remainingSpace = problem.knapsack_sizes;
    vector<int> sortedItemIds = problem.GetSortedItemIds(sortMode);

    for (int i = 0; i < isiz; ++i){
        int currentItemId = sortedItemIds[i];

        // if current item fits add it to the current solution
        globalSolution.AddItemIfFits(problem, currentItemId, remainingSpace);
    }

    return globalSolution;
}

Solution KnapsackSolver::GreedyPath(const Problem & problem, const Problem::SortMode & sortMode) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution
    vector<int> remainingSpace = problem.knapsack_sizes;
    vector<int> sortedItemIds = problem.GetSortedItemIds(sortMode);

    for (int i = 0; i < isiz; ){
        int currentItemId = sortedItemIds[i];

        if (Fits(problem.items[currentItemId].weights, remainingSpace))
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


//----------Brute Force Path----------

Solution KnapsackSolver::DFSPath(const Problem & problem, Solution currentSolution, const int & currentItemId, vector<int> remainingSpace) {

    if (Fits(problem.items[currentItemId].weights, remainingSpace))
        currentSolution.AddItem(problem, currentItemId, remainingSpace); // if current item fits add it to the current solution
    else return currentSolution; // if does not fit return what have we managed to fit so far

#ifdef _DEBUG_MODE_
    cout << "check for " << currentItemId << ", remaining weight: "<< remainingSpace[0] << endl;
#endif

    
    Solution globalSolution(problem.items.size()); // create empty solution

    // construct solution further with recursive DFSPath calls
    for (int i = 0; i < problem.items[currentItemId].connections.size(); ++i){
        int nextItemId = problem.items[currentItemId].connections[i];
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        Solution tempSolution = DFSPath(problem, currentSolution, nextItemId, remainingSpace);
        if (tempSolution.max_value > globalSolution.max_value) globalSolution = tempSolution;
    }

#ifdef _DEBUG_MODE_
    cout << "exit with: " << globalSolution.max_value << endl;
#endif
    return globalSolution; // return the best solution from all paths starting from currentItemId
}

Solution KnapsackSolver::BruteForcePath(const Problem & problem) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){
        Solution currentSolution(isiz); // create empty solution 

        
        currentSolution = DFSPath(problem, currentSolution, i, problem.knapsack_sizes);
        if(currentSolution.max_value > globalSolution.max_value) globalSolution = currentSolution; // if newly found solution is better use it as you global solution
    }

    return globalSolution;
}


//----------Branch and Bound----------

//----------Bounding Functions----------

int KnapsackSolver::GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode, Solution currentSolution, std::vector<int> remainingSpace){
    int isiz = problem.items.size();
    vector<int> sortedItemIds = problem.GetSortedItemIds(sortMode);

    for (int i = 0; i < isiz; ++i){
        int currentItemId = sortedItemIds[i];

        if (currentSolution.selected[currentItemId]) continue; // don't add items that were already in the solution

        currentSolution.AddItemIfFits(problem, currentItemId, remainingSpace);
    }

    return currentSolution.max_value;
}

//----------Branch and Bound----------

Solution KnapsackSolver::BranchAndBoundDFSPath(const Problem & problem, Solution currentSolution, const int & currentItemId, vector<int> remainingSpace, const int & lower_bound) {

#ifdef _DEBUG_MODE_
    cout << "check for " << currentItemId << ", remaining weight: "<< remainingSpace[0] << endl;
#endif

    Solution globalSolution = currentSolution;

    // construct solution further with recursive DFSPath calls
    for (int i = 0; i < problem.items[currentItemId].connections.size(); ++i){
        int nextItemId = problem.items[currentItemId].connections[i];
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        if (Fits(problem.items[nextItemId].weights, remainingSpace)){
            Solution tempSolution = currentSolution;
            vector<int> tempRemainingSpace = remainingSpace;

            tempSolution.AddItem(problem, nextItemId, tempRemainingSpace);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::WEIGHT_VALUE_RATIO, tempSolution, remainingSpace) > lower_bound*/) {
                tempSolution = BranchAndBoundDFSPath(problem, tempSolution, nextItemId, tempRemainingSpace, lower_bound);
                if (tempSolution.max_value > globalSolution.max_value) globalSolution = tempSolution;
            }
        }
    }

#ifdef _DEBUG_MODE_
    cout << "exit with: " << globalSolution.max_value << endl;
#endif
    return globalSolution; // return the best solution from all paths starting from currentItemId
}

Solution KnapsackSolver::BranchAndBoundPath(const Problem & problem) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){
        vector<int> remainingSpace = problem.knapsack_sizes;

        if (Fits(problem.items[i].weights, remainingSpace)){
            Solution currentSolution(isiz); // create empty solution
            currentSolution.AddItem(problem, i, remainingSpace);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::WEIGHT_VALUE_RATIO, currentSolution, remainingSpace) > globalSolution.max_value*/){
                currentSolution = BranchAndBoundDFSPath(problem, currentSolution, i, remainingSpace, globalSolution.max_value);
                if(currentSolution.max_value > globalSolution.max_value) globalSolution = currentSolution; // if newly found solution is better use it as you global solution
            }
        }
    }

    return globalSolution;
}


//---------- Floyd ----------

PackagedSolution Floyd::Solve(const Problem & problem, Requirements requirements){
    PackagedSolution ps;
    std::chrono::steady_clock::time_point start, end;
    Solution s;

    switch (requirements.structureToFind)
    {
    case Requirements::StructureToFind::PATH:
        start = std::chrono::steady_clock::now();
        s = Path(problem);
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
    ps.to_optimum_ratio = s.max_value / requirements.known_optimum;
    ps.remainingSpaces = KnapsackSolver::CalculateRemainingSpaces(s, problem);

    return ps;
}

//#define DEBUG_FLOYD_PATH

Solution Floyd::Path(const Problem & problem){
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