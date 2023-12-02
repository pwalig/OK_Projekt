#include "KnapsackSolver.hpp"
#include "json.hpp"

#include <chrono>
#include <fstream> // std::ifstream, std::ofstream

using namespace knapsack_solver;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using nlohmann::json;

//#define _DEBUG_MODE_

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
            cout << "pack full\n";
            throw;
        }
    }
    this->max_value += problem.items[selected_item_id].value; // update max value
    this->selected[selected_item_id] = true; // add item to solution
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

Solution KnapsackSolver::FirstFit(const Problem & problem){
    int subk = problem.knapsack_sizes.size();
    int isiz = problem.items.size();
    Solution solution(isiz);

    vector<int> fills = problem.knapsack_sizes;

    bool* checked = new bool[isiz];
    int checkedCount;
    for (checkedCount = 0; checkedCount < isiz; ++checkedCount) checked[checkedCount] = false;

    int testedItem = 0;
    checkedCount = 0;
    while(checkedCount < isiz) { // exit when all checked
        // item checked
        checked[testedItem] = true;
        checkedCount++;
        int nextItem = 0;

        bool fits = true;
        for (int j = 0; j < subk; ++j) {if(fills[j] - problem.items[testedItem].weights[j] < 0) {fits = false;}}
        if (fits){

            // add item to the answer
            solution.selected[testedItem] = true;
            solution.max_value += problem.items[testedItem].value;

            // fill knapsack - reduce remaining space
            for (int j = 0; j < subk; ++j) fills[j] -= problem.items[testedItem].weights[j];

            // choose next item
            int connectionsIterator = 0;
            do{
                nextItem = problem.items[testedItem].connections[connectionsIterator];
            }while(nextItem < isiz && checked[nextItem] && ++connectionsIterator < problem.items[testedItem].connections.size());
            if (nextItem == isiz) break;
        }
        else {
            // choose next item
            while(nextItem < isiz && checked[nextItem]){
                ++nextItem;
            }
            if (nextItem == isiz) break;
        }
        testedItem = nextItem;
    }
    delete[] checked;
    return solution;
}

Solution KnapsackSolver::RandomFit(const Problem & problem){
    int subk = problem.knapsack_sizes.size();
    int isiz = problem.items.size();
    Solution solution(isiz);

    vector<int> fills = problem.knapsack_sizes;

    bool* checked = new bool[isiz];
    int checkedCount = 0;
    for (int j = 0; j < isiz; ++j) checked[j] = false;

    int testedItem = rand() % isiz;
    while(checkedCount < isiz) { // exit when all checked
        ////////////////cout << "official check: " << testedItem << endl;
        // item checked
        checked[testedItem] = true;
        checkedCount++;

        int nextItem = 0;

        bool fits = true;
        for (int j = 0; j < subk; ++j) {if(fills[j] - problem.items[testedItem].weights[j] < 0) {fits = false;}}
        if (fits){
            /////////////cout << "fits!" << endl;

            // add item to the answer
            solution.selected[testedItem] = true;
            solution.max_value += problem.items[testedItem].value;

            // fill knapsack - reduce remaining space
            for (int j = 0; j < subk; ++j) fills[j] -= problem.items[testedItem].weights[j];


            // choose next item
            int connectionsAmount = problem.items[testedItem].connections.size();
            vector<bool> nextCheckedFit;
            int nextCheckedCount = 0;
            for (int j = 0; j < connectionsAmount; ++j) nextCheckedFit.push_back(false);

            do{
                nextItem = problem.items[testedItem].connections[rand() % connectionsAmount]; // try this one
                ///////////cout << "trying next: " << nextItem << endl;

                if (!nextCheckedFit[nextItem]){ // if this one was not tried before -> now it was -> note that
                    ++nextCheckedCount;
                    nextCheckedFit[nextItem] = true;
                }

            }while(nextCheckedCount < connectionsAmount && checked[nextItem]);
            if (nextCheckedCount >= connectionsAmount) break; // if all were tried -> stop
        }
        else {
            ///////////////cout << "no fit :(" << endl;
            // choose next item
            vector<bool> nextChecked;
            int nextCheckedCount = 0;
            for (int j = 0; j < isiz; ++j) nextChecked.push_back(false);
            do{
                nextItem = rand() % isiz; // try this one
                //cout << "trying next nf: " << nextItem << endl;

                if (!nextChecked[nextItem]){ // if this one was not tried before -> now it was -> note that
                    ++nextCheckedCount;
                    nextChecked[nextItem] = true;
                }
                ///////////////for (int j = 0; j < isiz; ++j) cout << " " << nextChecked[j];
                ////////////////cout << endl;

            }while(nextCheckedCount < isiz && checked[nextItem]);
            if (nextCheckedCount >= isiz) break; // if all were tried -> stop
        }
        testedItem = nextItem;
    }
    delete[] checked;
    return solution;
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