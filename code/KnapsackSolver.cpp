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

Solution::Solution() : max_value(0), valid(true) {}

Solution::Solution(const int & InstanceSize, const vector<int> & available_space) : max_value(0), valid(true), selected(InstanceSize, false), remainingSpace(available_space) { }

ostream& operator<<(ostream& os, const Solution& s){
    os << "max value: " << s.max_value << "\nselected: ";
    int instance_size = s.selected.size();
    for (int i = 0; i < instance_size; ++i){
        os << (s.selected[i] ? 1 : 0) << " ";
    }
    return (os << endl);
}

void Solution::AddItem(const Problem & problem, const int & selected_item_id){
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

void Solution::AddItemForce(const Problem & problem, const int & selected_item_id){
    for (int j = 0; j < remainingSpace.size(); ++j){
        remainingSpace[j] -= problem.items[selected_item_id].weights[j];
        if (remainingSpace[j] < 0) {
            this->valid = false;
        }
    }
    this->max_value += problem.items[selected_item_id].value; // update max value
    if (this->selected[selected_item_id] == true) throw std::invalid_argument("item was already in the solution");
    this->selected[selected_item_id] = true; // add item to solution
}

void Solution::RemoveItem(const Problem & problem, const int & selected_item_id){
    for (int j = 0; j < remainingSpace.size(); ++j){
        remainingSpace[j] += problem.items[selected_item_id].weights[j];
        if (remainingSpace[j] > problem.knapsack_sizes[j]) throw std::invalid_argument("remaining space exceeded problem capacity");
    }
    this->max_value -= problem.items[selected_item_id].value; // update max value
    if (this->selected[selected_item_id] == false) throw std::invalid_argument("item was not in the solution");
    this->selected[selected_item_id] = false; // remove item from solution
}

bool Solution::Fits(const Problem & problem, const int & selected_item_id){
    for (int j = 0; j < remainingSpace.size(); ++j){
        if (remainingSpace[j] < problem.items[selected_item_id].weights[j]){
            return false;
        }
    }
    return true;
}

bool Solution::AddItemIfFits(const Problem & problem, const int & selected_item_id){
    if (selected[selected_item_id]) throw std::invalid_argument("item was already in the solution"); // don't add item if it is already in the solution

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

int KnapsackSolver::GoalFunction(const Solution & solution, const PackagedProblem & problem){
    return 0;
}




//---------- BRUTE FORCE SOLVER ----------

PackagedSolution BruteForceSolver::Solve(const PackagedProblem & problem, const Options options){
    if (problem.requirements.structureToFind != Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS)
        throw std::logic_error("Structures other than ignore connections are uniplemented for brute force solver");
    
    PackagedSolution ps;
    std::chrono::steady_clock::time_point start, end;
    Solution s;

    ps.algorithm = "brute force; ignore connections; ";
    if (options.iterative) {
        ps.algorithm += "iterative; ";

        switch (options.search_order)
        {
        case BruteForceSolver::Options::SearchOrder::ZERO_FIRST:
            ps.algorithm += "zero first";
            break;

        case BruteForceSolver::Options::SearchOrder::ONE_FIRST:
            ps.algorithm += "one first";
            break;

        case BruteForceSolver::Options::SearchOrder::RANDOM:
            ps.algorithm += "random order";
            break;

        default:
            break;
        }
    }
    else {
        ps.algorithm += "recursive; ";
        ps.algorithm += options.late_fit ? "late_fit; " : "early_fit; ";
    }
    
    if (options.iterative) {
        start = std::chrono::steady_clock::now();
        s = Iterative(problem, options.search_order);
        end = std::chrono::steady_clock::now();
    }
    else {
        start = std::chrono::steady_clock::now();
        s = Recursive(problem, options);
        end = std::chrono::steady_clock::now();
    }

    // construct packaged solution
    ps.solution = s;
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = s.max_value / problem.known_optimum;
    ps.remainingSpaces = KnapsackSolver::CalculateRemainingSpaces(s, problem.problem);

    return ps;
}

Solution BruteForceSolver::Max(const Solution & a, const Solution & b){
    if (a.valid && !b.valid) return a;
    if (!a.valid && b.valid) return b;
    else if (a.max_value > b.max_value) return a;
    else return b;
}

Solution BruteForceSolver::SolutionFromNumber(int num, const Problem & problem) {
    // Vector to store the number
    Solution solution(problem.items.size(), problem.knapsack_sizes);
 
    int i = problem.items.size();
    while (num != 0) {
        --i;
        if (num % 2 == 1) solution.AddItemForce(problem, i);
 
        // Integer division
        // gives quotient
        num /= 2;
    }
    return solution;
}

Solution BruteForceSolver::Iterative(const PackagedProblem & problem, const Options::SearchOrder & search_order){
    int isiz = problem.problem.items.size();
    Solution global_solution(isiz, problem.problem.knapsack_sizes);

    switch (search_order)
    {
    case Options::SearchOrder::ZERO_FIRST:
        for (int i = 0; i < pow(2, isiz); ++i) {
            Solution temp_solution = SolutionFromNumber(i, problem.problem);
            if (temp_solution.valid && temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;
        }
        break;

    case Options::SearchOrder::ONE_FIRST:
        for (int i = pow(2, isiz) - 1; i >= 0; --i) {
            Solution temp_solution = SolutionFromNumber(i, problem.problem);
            if (temp_solution.valid && temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;
        }
        break;

    case Options::SearchOrder::RANDOM:
        throw std::logic_error("Random brute force order not implemented yet");
        break;
    }

    return global_solution;
}

Solution BruteForceSolver::DFS(const Problem & problem, Solution sol0, const Options::SearchOrder & search_order, const bool & add, const int & depth){
    if (add) sol0.AddItemForce(problem, depth);
    if (depth == problem.items.size() - 1) return sol0;

    return Max(DFS(problem, sol0, search_order, false, depth + 1), DFS(problem, sol0, search_order, true, depth + 1));
}

Solution BruteForceSolver::DFS(const Problem & problem, Solution sol0, const Options::SearchOrder & search_order, const int & depth){
    int isiz = problem.items.size();
    if (depth == isiz) return sol0;
    
    Solution sol1 = sol0;
    sol0 = DFS(problem, sol0, search_order, depth + 1);
    sol1.AddItemForce(problem, depth);
    sol1 = DFS(problem, sol1, search_order, depth + 1);

    return Max(sol0, sol1);
}

Solution BruteForceSolver::Recursive(const PackagedProblem & problem, const Options & options){
    int isiz = problem.problem.items.size();
    Solution sol0(isiz, problem.problem.knapsack_sizes);
    
    if (options.late_fit) {
        return Max(DFS(problem.problem, sol0, options.search_order, false, 0), DFS(problem.problem, sol0, options.search_order, true, 0));
    }
    else{
        Solution sol1 = sol0;
        sol0 = DFS(problem.problem, sol0, options.search_order, 1);
        sol1.AddItemForce(problem.problem, 0);
        sol1 = DFS(problem.problem, sol1, options.search_order, 1);
        return Max(sol0, sol1);
    }

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
    Solution globalSolution(isiz, problem.problem.knapsack_sizes); // create empty solution
    vector<int> sortedItemIds = problem.problem.GetSortedItemIds(options.sort_mode);

    switch (problem.requirements.structureToFind)
    {

    case Problem::Requirements::StructureToFind::PATH:
    {
        // Find best starting point
        int startItem = -1;
        for (int currentItemId : sortedItemIds) {
            if (globalSolution.Fits(problem.problem, currentItemId)){
                globalSolution.AddItem(problem.problem, currentItemId); // if current item fits add it to the current solution
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
                    globalSolution.Fits(problem.problem, currentItemId))
                {
                    globalSolution.AddItem(problem.problem, currentItemId);
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
            globalSolution.AddItemIfFits(problem.problem, currentItemId);
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
    Solution globalSolution(isiz, problem.knapsack_sizes); // create empty solution
    vector<int> sortedItemIds = problem.GetSortedItemIds(options.sort_mode);

    for (int i = 0; i < isiz; ++i){
        int currentItemId = sortedItemIds[i];

        // if current item fits add it to the current solution
        globalSolution.AddItemIfFits(problem, currentItemId);
    }

    return globalSolution;
}

Solution GreedySolver::GreedyPath(const Problem & problem, const Options & options) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz, problem.knapsack_sizes); // create empty solution
    vector<int> sortedItemIds = problem.GetSortedItemIds(options.sort_mode);

    for (int i = 0; i < isiz; ){
        int currentItemId = sortedItemIds[i];

        if (globalSolution.Fits(problem, currentItemId))
            globalSolution.AddItem(problem, currentItemId); // if current item fits add it to the current solution
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

int BranchAndBoundSolver::GreedyIgnoreConnections(const Problem & problem, const Problem::SortMode & sortMode, Solution currentSolution){
    int isiz = problem.items.size();
    vector<int> sortedItemIds = problem.GetSortedItemIds(sortMode);

    for (int i = 0; i < isiz; ++i){
        int currentItemId = sortedItemIds[i];

        if (currentSolution.selected[currentItemId]) continue; // don't add items that were already in the solution

        currentSolution.AddItemIfFits(problem, currentItemId);
    }

    return currentSolution.max_value;
}

//---------- Early Fit ----------

Solution BranchAndBoundSolver::DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, const int & lower_bound) {

#ifdef DEBUG_BNB
    cout << "check for " << currentItemId << ", remaining weight: "<< remainingSpace[0] << endl;
#endif

    Solution globalSolution = currentSolution;

    // construct solution further with recursive DFSLateFitPath calls
    for (int nextItemId : problem.items[currentItemId].connections){
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        if (currentSolution.Fits(problem, nextItemId)){
            Solution tempSolution = currentSolution;

            tempSolution.AddItem(problem, nextItemId);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::WEIGHT_VALUE_RATIO, tempSolution, remainingSpace) > lower_bound*/) {
                tempSolution = DFSEarlyFitPath(problem, tempSolution, nextItemId, lower_bound);
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
    Solution globalSolution(isiz, problem.knapsack_sizes); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){

        if (globalSolution.Fits(problem, i)){
            Solution currentSolution(isiz, problem.knapsack_sizes); // create empty solution
            currentSolution.AddItem(problem, i);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::WEIGHT_VALUE_RATIO, currentSolution, remainingSpace) > globalSolution.max_value*/){
                currentSolution = DFSEarlyFitPath(problem, currentSolution, i, globalSolution.max_value);
                if(currentSolution.max_value > globalSolution.max_value) globalSolution = currentSolution; // if newly found solution is better use it as you global solution
            }
        }
    }

    return globalSolution;
}

//---------- Late Fit ----------

Solution BranchAndBoundSolver::DFSLateFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId) {

    if (!currentSolution.AddItemIfFits(problem, currentItemId)) // if current item fits add it to the current solution
        return currentSolution; // if does not fit return what have we managed to fit so far

#ifdef _DEBUG_MODE_
    cout << "check for " << currentItemId << ", remaining weight: "<< currentSolution.remainingSpace[0] << endl;
#endif

    
    Solution globalSolution(problem.items.size(), problem.knapsack_sizes); // create empty solution

    // construct solution further with recursive DFSLateFitPath calls
    for (int nextItemId : problem.items[currentItemId].connections){
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        Solution tempSolution = DFSLateFitPath(problem, currentSolution, nextItemId);
        if (tempSolution.max_value > globalSolution.max_value) globalSolution = tempSolution;
    }

#ifdef _DEBUG_MODE_
    cout << "exit with: " << globalSolution.max_value << endl;
#endif
    return globalSolution; // return the best solution from all paths starting from currentItemId
}

Solution BranchAndBoundSolver::BnBLateFitPath(const Problem & problem) {
    int isiz = problem.items.size();
    Solution globalSolution(isiz, problem.knapsack_sizes); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){
        Solution currentSolution(isiz, problem.knapsack_sizes); // create empty solution 

        currentSolution = DFSLateFitPath(problem, currentSolution, i);
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

    // Fill initial D vector
    vector<vector<Solution>> d;
    for (int i = 0; i < isiz; ++i) {
        vector<Solution> di;
        Solution helper_solution (isiz, problem.knapsack_sizes);

        // if you can't fit element fill whole row with invalid answers
        if (!helper_solution.Fits(problem, i)) {
            Solution sp(isiz, problem.knapsack_sizes);
            sp.max_value = -1;
            for (int j = 0; j < isiz; ++j) di.push_back(sp);
        }
        else {
            for (int j = 0; j < isiz; ++j){
                Solution sp(isiz, problem.knapsack_sizes);
                sp.AddItem(problem, i); // i did fit, so add it

                if (!helper_solution.Fits(problem, j)) {
                    sp.RemoveItem(problem, i);
                    sp.max_value = -1; // item j does not fit - solution invalid
                }
                else if (problem.items[i].HasConnectionTo(j)) sp.AddItem(problem, j);
                else {
                    sp.RemoveItem(problem, i);
                    sp.max_value = -2; // items lack connection - solution invalid
                }
                di.push_back(sp);
            }
        }
        d.push_back(di);
    }

#ifdef DEBUG_FLOYD_PATH
    for (int i = 0; i < isiz; ++i) {
        for (int j = 0; j < isiz; ++j){
            cout << d[i][j].max_value << "\t";
        }
        cout << endl;
    }
#endif

    // Solve
    for (int k = 0; k < isiz; ++k){
        for (int i = 0; i < isiz; ++i){
            for (int j = 0; j < isiz; ++j){
                if (d[i][k].max_value >= 0 && d[k][j].max_value >= 0 && // if you can get from i to k and from k to j
                    !d[i][j].selected[k] &&
                    d[i][j].max_value != -1 && // items dont fit anyway
                    d[i][j].Fits(problem, k) && // initial fit check
                    d[i][k].max_value + d[k][j].max_value > d[i][j].max_value)
                {   
                    if (d[i][j].max_value == -2){ // if items lacked connection but now have it through item k
                        d[i][j].max_value = 0;
                        d[i][j].AddItem(problem, i);
                        if (i != j) d[i][j].AddItem(problem, j);
                        if (!d[i][j].Fits(problem, k)) { // if k stopped to fit after adding i and j
                            d[i][j].RemoveItem(problem, i);
                            if (i != j) d[i][j].RemoveItem(problem, j);
                            d[i][j].max_value = -2;
                            continue;
                        }
                        else d[i][j].AddItem(problem, k);
                    }
                    else d[i][j].AddItem(problem, k);
                }
            }
        }
    
#ifdef DEBUG_FLOYD_PATH
        cout << endl << "through: " << k << endl;
        for (int i = 0; i < isiz; ++i) {
            for (int j = 0; j < isiz; ++j){
                cout << d[i][j].max_value << "\t";
            }
            cout << endl;
        }
#endif
    }

    // Find Best Solution
    int maxi, maxj, maxval = 0;
    for (int i = isiz - 1; i >= 0; --i) {
        for (int j = isiz - 1; j >= 0; --j){
            if (d[i][j].max_value > maxval){
                maxval = d[i][j].max_value;
                maxi = i;
                maxj = j;
            }
        }
    }

    return d[maxi][maxj];
}