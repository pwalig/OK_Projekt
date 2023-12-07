#include "KnapsackSolver.hpp"
//"KnapsackSolver.hpp" already includes:
    //#include <iostream>  // std::ostream
    //#include <fstream> // std::ifstream
    //#include <vector>
    //#include <string>

    //Batch Solve requires these three:
    //#include <filesystem> // std::filesystem::create_directory
    //#include "json.hpp"
    //#include "FileNameDefines.hpp"

    //#include "Problem.hpp"

#include <chrono>
#include <fstream> // std::ofstream
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
    os << "max value: " << s.max_value << "\nselected:";
    int instance_size = s.selected.size();
    for (int i = 0; i < instance_size; ++i){
        os << " " << s.selected[i];
    }
    os << "\nremaining spaces:";
    for (int i = 0; i < s.remainingSpace.size(); ++i){
        os << " " << s.remainingSpace[i];
    }
    os << "\nvalid: " << s.valid;
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

bool Solution::IsPathDFS(const Problem & problem, vector<int> & visited, const int & current, const int & length) const{
    for (int next : problem.items[current].connections) {
        if (selected[next] && std::find(visited.begin(), visited.end(), next) == visited.end()){ // next item has to be selected and new
            if (visited.size() + 1 == length) return true; // path found
            if (visited.size() + 1 > length) return false; // path would have to be to long
            visited.push_back(next);
            if (IsPathDFS(problem, visited, next, length)) return true; // path found later
            visited.pop_back();
        }
    }
    return false; // path not found
}
bool Solution::IsPath(const Problem & problem) const{
    if (selected.size() != problem.items.size()) throw std::invalid_argument("amount of available items does not match");

    // calculate whats the length of the path
    int length = 0;
    for (int i = 0; i < selected.size(); ++i)
        if (selected[i]) ++length;
    
    if (length <= 1) return true;
    
    // check from each starting point
    vector<int> visited;
    for (int i = 0; i < selected.size(); ++i) {
        if (selected[i]){
            visited.push_back(i);
            if (IsPathDFS(problem, visited, i, length)) return true; // path found somewhere
            visited.pop_back();
        }
    }
    return false;
}

bool Solution::IsCycleDFS(const Problem & problem, vector<int> & visited, const int & current, const int & start, const int & length) const{
    for (int next : problem.items[current].connections) {
        if (selected[next] && std::find(visited.begin(), visited.end(), next) == visited.end()){ // next item has to be selected and new
            if (visited.size() + 1 == length && problem.items[next].HasConnectionTo(start)) return true; // cycle found
            if (visited.size() + 1 > length) return false; // cycle would have to be to long
            visited.push_back(next);
            if (IsCycleDFS(problem, visited, next, start, length)) return true; // cycle found later
            visited.pop_back();
        }
    }
    return false; // cycle not found
}
bool Solution::IsCycle(const Problem & problem) const{
    if (selected.size() != problem.items.size()) throw std::invalid_argument("amount of available items does not match");

    // calculate whats the length of the cycle
    int length = 0;
    for (int i = 0; i < selected.size(); ++i)
        if (selected[i]) ++length;
    
    if (length == 0) return true;
    
    // check from each starting point
    vector<int> visited;
    for (int i = 0; i < selected.size(); ++i) {
        if (selected[i]){
            if (length == 1) return problem.items[i].HasConnectionTo(i);
            visited.push_back(i);
            if (IsCycleDFS(problem, visited, i, i, length)) return true; // cycle found somewhere
            visited.pop_back();
        }
    }
    return false;
}

bool Solution::IsTree(const Problem & problem) const{
    throw std::logic_error("not implemented yet");
    return false;
}
bool Solution::IsConnected(const Problem & problem) const{
    throw std::logic_error("not implemented yet");
    return false;
}
bool Solution::IsCyclePossible(const Problem & problem) const{
    throw std::logic_error("not implemented yet");
    return false;
}




//---------- PACKAGED SOLUTION ----------

void PackagedSolution::ExportJSON(const std::string file_name){
    json data;

    // packaged soution
    data["algorithm"] = algorithm;
    to_optimum_ratio >= 0 ? data["to_optimum_ratio"] = to_optimum_ratio : data["to_optimum_ratio"] = "optimum_unnown";
    data["solve_time"] = solve_time;
    data["remaining_spaces"] = solution.remainingSpace;

    // solution
    data["solution"]["max_value"] = solution.max_value;
    for (int i = 0; i < solution.selected.size(); ++i){
        data["solution"]["selected"][i] = solution.selected[i];
    }

    // validation status
    if (validation_status.undergone) {
        data["validation_status"]["valid"] = validation_status.valid;
        data["validation_status"]["fit"] = validation_status.fit;
        data["validation_status"]["remaining_space"] = validation_status.remaining_space;
        data["validation_status"]["self_valid"] = validation_status.self_valid;
        data["validation_status"]["structure"] = validation_status.structure;
        data["validation_status"]["value"] = validation_status.value;
    }
    else data["validation_status"] = "did_not_undergo_validation";

    std::ofstream fout(file_name);
    fout << data.dump(4);
    fout.close();
}

ostream& operator<<(ostream& os, const PackagedSolution& ps){
    // packaged solution
    os << "solved with: " << ps.algorithm << endl << ps.solution;
    os << "to_optimum_ratio: ";
    ps.to_optimum_ratio >= 0 ? os << ps.to_optimum_ratio : os << "optimum unknown";
    os << "\nsolve time: " << ps.solve_time << endl;

    // validation status
    if (ps.validation_status.undergone) {
        os << "validaton status:\n";
        os << "\tvalid: " << ps.validation_status.valid << endl;
        os << "\tvalue: " << ps.validation_status.value << endl;
        os << "\tfit: " << ps.validation_status.fit << endl;
        os << "\tremaining_space: " << ps.validation_status.remaining_space << endl;
        os << "\tstructure: " << ps.validation_status.structure << endl;
        os << "\tself_valid: " << ps.validation_status.self_valid << endl;
    }
    else os << "did_not_undergo_validation";

    return (os << endl);
}




//---------- VALIDATION ----------

vector<int> Validation::CalculateRemainingSpaces(const Solution & solution, const Problem & problem){
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

int Validation::CalculateMaxValue(const Solution & solution, const Problem & problem){
    if (solution.selected.size() != problem.items.size()) throw std::invalid_argument("amount of available items does not mathch");
    int sum = 0;
    for (int i = 0; i < solution.selected.size(); ++i){
        if (solution.selected[i]) sum += problem.items[i].value;
    }
    return sum;
}

Validation::ValidationStatus Validation::Validate(const Solution & solution, const PackagedProblem & problem){
    ValidationStatus vs;
    if (solution.max_value != CalculateMaxValue(solution, problem.problem)) vs.value = false; // solution calculated its value wrong
    vector<int> calculatedRemainingSpace = CalculateRemainingSpaces(solution, problem.problem);
    if (solution.remainingSpace != calculatedRemainingSpace) vs.remaining_space = false; // solution calculated its remaining space wrong
    for (int weight : calculatedRemainingSpace) if (weight < 0) vs.fit = false; // items dont fit

    switch (problem.requirements.structureToFind)
    {
    case Problem::Requirements::StructureToFind::PATH:
        vs.structure = solution.IsPath(problem.problem);
        break;
        
    case Problem::Requirements::StructureToFind::CYCLE:
        vs.structure = solution.IsCycle(problem.problem);
        break;
    
    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
        vs.structure = true;
        break;

    default:
        break;
    }

    vs.valid = vs.value && vs.remaining_space && vs.fit && vs.structure; // set main valid
    if (solution.valid != vs.valid) {vs.self_valid = false; vs.valid = false;} // solution either: invalided itself for no reason, or: thought it was valid even though it was not

    vs.undergone = true;
    return vs;
}

int Validation::GoalFunction(const Solution & solution, const PackagedProblem & problem){

    switch (problem.requirements.structureToFind)
    {
    case Problem::Requirements::StructureToFind::PATH:
        if (!solution.IsPath(problem.problem)) return INT_MIN;
        break;
        
    case Problem::Requirements::StructureToFind::CYCLE:
        if (!solution.IsCycle(problem.problem)) return INT_MIN;
        break;

    default:
        break;
    }
    return CalculateMaxValue(solution, problem.problem);
}




//---------- BRUTE FORCE SOLVER ----------

PackagedSolution BruteForceSolver::Solve(const PackagedProblem & problem, const Options options){
    if (problem.requirements.structureToFind != Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS)
        throw std::logic_error("Structures other than ignore connections are uniplemented for brute force solver");
    
    PackagedSolution ps;
    std::chrono::steady_clock::time_point start, end;

    // fil algorithm info / details
    ps.algorithm = "brute-force_ignore-connections_";
    if (options.iterative) {
        ps.algorithm += "iterative_";

        switch (options.search_order)
        {
        case BruteForceSolver::Options::SearchOrder::ZERO_FIRST:
            ps.algorithm += "zero-first";
            break;
            
        case BruteForceSolver::Options::SearchOrder::UNCONSTRAINED:
            ps.algorithm += "unconstrained";
            break;

        case BruteForceSolver::Options::SearchOrder::ONE_FIRST:
            ps.algorithm += "one-first";
            break;

        default:
            throw std::logic_error("not implemented yet");
            break;
        }
    }
    else {
        ps.algorithm += "recursive_";
        ps.algorithm += options.late_fit ? "late-fit" : "early-fit";
    }
    
    // run with time measure
    if (options.iterative) {
        start = std::chrono::steady_clock::now();
        ps.solution = Iterative(problem, options.search_order);
        end = std::chrono::steady_clock::now();
    }
    else {
        start = std::chrono::steady_clock::now();
        ps.solution = Recursive(problem, options);
        end = std::chrono::steady_clock::now();
    }
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = 1.0;

    // validate
    ps.validation_status = Validation::Validate(ps.solution, problem);

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

    int n = pow(2, isiz);
    switch (search_order)
    {
    case Options::SearchOrder::UNCONSTRAINED:
    case Options::SearchOrder::ZERO_FIRST:
        for (int i = 0; i < n; ++i) {
            Solution temp_solution = SolutionFromNumber(i, problem.problem);
            if (temp_solution.valid && temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;
        }
        break;

    case Options::SearchOrder::ONE_FIRST:
        for (int i = n - 1; i >= 0; --i) {
            Solution temp_solution = SolutionFromNumber(i, problem.problem);
            if (temp_solution.valid && temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;
        }
        break;

    default:
        throw std::logic_error("not implemented yet");
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
    if (options.search_order != BruteForceSolver::Options::SearchOrder::UNCONSTRAINED)
        throw std::invalid_argument("recursive does not support searc order"); //it could be supported, but it slows algorithm down
    
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

    // fill algorithm info / details
    ps.algorithm = "";

    switch (options.sort_mode)
    {
    case Problem::SortMode::WEIGHT_VALUE_RATIO:
        ps.algorithm += "greedy-sort-by-value-weight";
        break;
        
    case Problem::SortMode::VALUE:
        ps.algorithm += "greedy-sort-by-value";
        break;
        
    case Problem::SortMode::WEIGHT:
        ps.algorithm += "greedy-sort-by-weight";
        break;

    case Problem::SortMode::RANDOM:
        ps.algorithm += "random-fit";
        break;
    
    case Problem::SortMode::DONT_SORT:
        ps.algorithm += "first-fit";
        break;

    default:
        break;
    }

    switch (problem.requirements.structureToFind)
    {
    case Problem::Requirements::StructureToFind::PATH:
        ps.algorithm += "_path";
        break;
    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
        ps.algorithm += "_ignore-connections";
        break;
    
    default:
        break;
    }
    
    // run with time measure
    std::chrono::steady_clock::time_point start, end;
    start = std::chrono::steady_clock::now();
    ps.solution = GreedyUniversal(problem, options);
    end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();

    ps.to_optimum_ratio = ps.solution.max_value / problem.known_optimum;

    // validate
    ps.validation_status = Validation::Validate(ps.solution, problem);

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
    ps.algorithm = "branch-and-bound";

    // solve with measured time
    if (options.late_fit){
        start = std::chrono::steady_clock::now();
        s = BnBLateFitPath(problem.problem);
        end = std::chrono::steady_clock::now();
        ps.algorithm += "_late-fit";
    }
    else{
        start = std::chrono::steady_clock::now();
        s = BnBEarlyFitPath(problem.problem);
        end = std::chrono::steady_clock::now();
        ps.algorithm += "_early-fit";
    }
    const std::chrono::duration<double> elapsed_seconds{end - start};

    ps.algorithm += "_path";
    ps.solution = s;
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();
    ps.to_optimum_ratio = 1.0;

    ps.validation_status = Validation::Validate(s, problem);

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

Solution BranchAndBoundSolver::DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & currentItemId, int lower_bound) {

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
                if (tempSolution.max_value > globalSolution.max_value) {
                    globalSolution = tempSolution;
                    lower_bound = tempSolution.max_value;
                }
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

#ifdef DEBUG_BNB_LF
    cout << "check for " << currentItemId << ", remaining weight: "<< currentSolution.remainingSpace[0] << endl;
#endif

    
    Solution globalSolution(problem.items.size(), problem.knapsack_sizes); // create empty solution

    // construct solution further with recursive DFSLateFitPath calls
    for (int nextItemId : problem.items[currentItemId].connections){
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        Solution tempSolution = DFSLateFitPath(problem, currentSolution, nextItemId);
        if (tempSolution.max_value > globalSolution.max_value) globalSolution = tempSolution;
    }

#ifdef DEBUG_BNB_LF
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