#include "KnapsackSolver.hpp"
//"KnapsackSolver.hpp" already includes:
    //#include <iostream>  // std::ostream
    //#include <fstream> // std::ifstream
    //#include <vector>
    //#include <string>
    //#include <cfloat> // DBL_MAX

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
using std::string;
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

void Solution::AddItem(const Problem & problem, const int & selected_item_id, const FaultTreatment & fit_fault){
    for (int j = 0; j < remainingSpace.size(); ++j){
        remainingSpace[j] -= problem.items[selected_item_id].weights[j];
        if (remainingSpace[j] < 0) {
            if (fit_fault == FaultTreatment::THROW) throw std::invalid_argument("item does not fit");
            else if (fit_fault == FaultTreatment::INVALIDATE) valid = false;
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

bool Solution::Fits(const Problem & problem, const int & selected_item_id) const{
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

bool Solution::IsFit(const Problem & problem) const{
    vector<int> remainigSpaces = problem.knapsack_sizes;
    for (int i = 0; i < this->selected.size(); ++i){
        if (this->selected[i]) {
            for (int j = 0; j < remainigSpaces.size(); ++j){
                remainigSpaces[j] -= problem.items[i].weights[j];
                if (remainigSpaces[j] < 0) return false;
            }
        }
    }
    return true;
}

bool Solution::IsPathDFS(const Problem & problem, vector<int> & visited, const int & current, const int & length) const{
    for (int next : problem.items[current].connections) {
        if (selected[next] && std::find(visited.begin(), visited.end(), next) == visited.end()){ // next item has to be selected and new
            visited.push_back(next);
            if (visited.size() == length) return true; // path found
            if (visited.size() > length) return false; // path would have to be to long
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
            visited.push_back(next);
            if (visited.size() == length && problem.items[next].HasConnectionTo(start)) return true; // cycle found
            if (visited.size() > length) return false; // cycle would have to be to long
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

bool Solution::IsStructure(const PackagedProblem & problem) const{
    switch (problem.requirements.structureToFind)
    {
    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
        return true;
        break;

    case Problem::Requirements::StructureToFind::CONNECTED_GRAPH:
        return IsConnected(problem.problem);
        break;

    case Problem::Requirements::StructureToFind::PATH:
        return IsPath(problem.problem);
        break;
        
    case Problem::Requirements::StructureToFind::CYCLE:
        return IsCycle(problem.problem);
        break;
        
    case Problem::Requirements::StructureToFind::TREE:
        return IsTree(problem.problem);
        break;
    
    default:
        throw std::invalid_argument("unnknown structure");
        break;
    }
}

bool Solution::IsValid(const PackagedProblem & problem) const{
    return IsFit(problem.problem) && IsStructure(problem);
}

bool Solution::IsCyclePossibleDFS(const Problem & problem, vector<int> & visited, vector<int> _remaining_space, const int & current, const int & start) const{
    for (int next : problem.items[current].connections) {
        if (std::find(visited.begin(), visited.end(), next) == visited.end()){ // next item has to be new (not visited yet)

            bool fit = true;
            for (int j = 0; j < _remaining_space.size(); ++j){
                if (_remaining_space[j] >= problem.items[next].weights[j]) _remaining_space[j] -= problem.items[next].weights[j];
                else { fit = false; break; }
            }
            if (!fit) continue;

            visited.push_back(next);
            if (problem.items[next].HasConnectionTo(start)){ // found some cycle lets check if it has all selected vertices
                bool _found = true;
                for (int i = 0; i < this->selected.size(); ++i){
                    if (this->selected[i] && std::find(visited.begin(), visited.end(), i) == visited.end()){
                        _found = false;
                        break;
                    }
                }
                if (_found) return true; // it has - cycle found
            }
            if (IsCyclePossibleDFS(problem, visited, _remaining_space, next, start)) return true; // cycle found later
            visited.pop_back();
        }
    }
    return false; // cycle not found
}
bool Solution::IsCyclePossible(const Problem & problem) const{
    if (selected.size() != problem.items.size()) throw std::invalid_argument("amount of available items does not match");
    vector<int> visited;
    vector<int> _remaining_space = problem.knapsack_sizes;
    for(int i = 0; i < selected.size(); ++i){
        bool fit = true;
        for (int j = 0; j < _remaining_space.size(); ++j){
            if (_remaining_space[j] >= problem.items[i].weights[j]) _remaining_space[j] -= problem.items[i].weights[j];
            else { fit = false; break; }
        } 
        if (!fit) continue;

        visited.push_back(i);
        if (IsCyclePossibleDFS(problem, visited, _remaining_space, i, i)) return true; // cycle found somewhere
        visited.pop_back();
    }
    return false;
}

bool Solution::IsPathPossibleDFS(const Problem & problem, vector<int> & visited, vector<int> _remaining_space, const int & current) const{
    for (int next : problem.items[current].connections) {
        if (std::find(visited.begin(), visited.end(), next) == visited.end()){ // next item has to be new (not visited yet)

            bool fit = true;
            for (int j = 0; j < _remaining_space.size(); ++j){
                if (_remaining_space[j] >= problem.items[next].weights[j]) _remaining_space[j] -= problem.items[next].weights[j];
                else { fit = false; break; }
            }
            if (!fit) continue;

            visited.push_back(next);
            // found some path lets check if it has all selected vertices
            bool _found = true;
            for (int i = 0; i < this->selected.size(); ++i){
                if (this->selected[i] && std::find(visited.begin(), visited.end(), i) == visited.end()){
                    _found = false;
                    break;
                }
            }
            if (_found) return true; // it has - path found
            if (IsPathPossibleDFS(problem, visited, _remaining_space, next)) return true; // path found later
            visited.pop_back();
        }
    }
    return false; // cycle not found
}
bool Solution::IsPathPossible(const Problem & problem) const{
    if (selected.size() != problem.items.size()) throw std::invalid_argument("amount of available items does not match");
    vector<int> visited;
    vector<int> _remaining_space = problem.knapsack_sizes;
    for(int i = 0; i < selected.size(); ++i){
        bool fit = true;
        for (int j = 0; j < _remaining_space.size(); ++j){
            if (_remaining_space[j] >= problem.items[i].weights[j]) _remaining_space[j] -= problem.items[i].weights[j];
            else { fit = false; break; }
        } 
        if (!fit) continue;

        visited.push_back(i);
        if (IsPathPossibleDFS(problem, visited, _remaining_space, i)) return true; // cycle found somewhere
        visited.pop_back();
    }
    return false;
}

bool Solution::IsStructurePossible(const PackagedProblem & problem) const{
    switch (problem.requirements.structureToFind)
    {
    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS:
        return true;
        break;

    case Problem::Requirements::StructureToFind::CONNECTED_GRAPH:
    throw std::logic_error("not implemented yet");
        break;

    case Problem::Requirements::StructureToFind::PATH:
        return IsPathPossible(problem.problem);
        break;
        
    case Problem::Requirements::StructureToFind::CYCLE:
        return IsCyclePossible(problem.problem);
        break;
        
    case Problem::Requirements::StructureToFind::TREE:
    throw std::logic_error("not implemented yet");
        break;
    
    default:
        throw std::invalid_argument("unnknown structure");
        break;
    }
}





//---------- PACKAGED SOLUTION ----------

void PackagedSolution::ExportJSON(const std::string file_name) const{
    json data;

    // packaged soution
    data["algorithm"] = algorithm;
    quality >= 0 ? data["to_optimum_ratio"] = quality : data["to_optimum_ratio"] = "optimum_unnown";
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
        data["validation_status"]["quality"] = validation_status.quality;
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
    ps.quality >= 0 ? os << ps.quality : os << "optimum unknown";
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
        os << "\tquality: " << ps.validation_status.quality << endl;
    }
    else os << "did_not_undergo_validation";

    return (os << endl);
}




//---------- VALIDATION ----------

vector<int> Validation::CalculateRemainingSpaces(const std::vector<bool> & selected, const Problem & problem){
    vector<int> remainigSpaces = problem.knapsack_sizes;
    for (int i = 0; i < selected.size(); ++i){
        if (selected[i]) {
            for (int j = 0; j < remainigSpaces.size(); ++j){
                remainigSpaces[j] -= problem.items[i].weights[j];
            }
        }
    }
    return remainigSpaces;
}

int Validation::CalculateMaxValue(const std::vector<bool> & selected, const Problem & problem){
    if (selected.size() != problem.items.size()) throw std::invalid_argument("amount of available items does not mathch");
    int sum = 0;
    for (int i = 0; i < selected.size(); ++i){
        if (selected[i]) sum += problem.items[i].value;
    }
    return sum;
}

Validation::ValidationStatus Validation::Validate(const Solution & solution, const PackagedProblem & problem){
    ValidationStatus vs;
    if (solution.max_value != CalculateMaxValue(solution.selected, problem.problem)) vs.value = false; // solution calculated its value wrong
    vector<int> calculatedRemainingSpace = CalculateRemainingSpaces(solution.selected, problem.problem);
    if (solution.remainingSpace != calculatedRemainingSpace) vs.remaining_space = false; // solution calculated its remaining space wrong
    for (int weight : calculatedRemainingSpace) if (weight < 0) vs.fit = false; // items dont fit

    vs.structure = solution.IsStructure(problem);

    vs.valid = vs.value && vs.remaining_space && vs.fit && vs.structure; // set main valid
    if (solution.valid != vs.valid) {vs.self_valid = false; vs.valid = false;} // solution either: invalided itself for no reason, or: thought it was valid even though it was not

    vs.undergone = true;
    return vs;
}




//---------- KNAPSACK SOLVER ----------

int KnapsackSolver::GoalFunction(const Solution & solution, const PackagedProblem & problem){
    if (!solution.IsValid(problem)) return INT_MIN;
    return Validation::CalculateMaxValue(solution.selected, problem.problem);
}




//---------- DYNAMIC SOLVER ----------

bool DynamicSolver::expect_perfection = true;

string DynamicSolver::GetAlgorithmName(const Options & options){
    string name = "dynamic-programming";
    return name;
}

PackagedSolution DynamicSolver::Solve(PackagedProblem & problem, const Options & options){
    if (problem.requirements.weightTreatment != Problem::Requirements::WeightTreatment::RESPECT_FIRST_ONLY && problem.problem.knapsack_sizes.size() > 1)
        throw std::invalid_argument("dynamic programming can't solve for multiple weights");
    if (problem.requirements.structureToFind != Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS)
        throw std::invalid_argument("dynamic programming can't solve for required structure");
    
    PackagedSolution ps;
    
    // run with time measure
    std::chrono::steady_clock::time_point start, end;
    start = std::chrono::steady_clock::now();
    ps.solution = Dynamic(problem.problem);
    end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();

    return ps;
}

Solution DynamicSolver::Dynamic(const Problem & problem){

    // Solution
    int isiz = problem.items.size();
    vector<vector<int>> dp(isiz+1, vector<int>(problem.knapsack_sizes[0]+1, 0)); // initialise dynamic programming 2D NxP array with zeros
    for (int i = 1; i <= isiz; ++i){
        for (int j = 1; j <= problem.knapsack_sizes[0]; ++j){
            if (problem.items[i-1].weights[0] <= j)
                dp[i][j] = std::max(problem.items[i-1].value + dp[i-1][j - problem.items[i-1].weights[0]], dp[i-1][j]);
            else
                dp[i][j] = dp[i-1][j];
        }
    }

    // Back tracking
    Solution s(isiz, problem.knapsack_sizes);
    int res = dp[isiz][problem.knapsack_sizes[0]];

    int w = problem.knapsack_sizes[0];
    for (int i = isiz; i > 0 && res > 0; i--) {
         
        // either the result comes from the top
        // (K[i-1][w]) or from (val[i-1] + K[i-1]
        // [w-wt[i-1]]) as in Knapsack table. If
        // it comes from the latter one/ it means
        // the item is included.
        if (res == dp[i - 1][w])
            continue;    
        else {
            // This item is included.
            s.AddItem(problem, i-1);
             
            // Since this weight is included its
            // value is deducted
            res -= problem.items[i - 1].value;
            w -= problem.items[i - 1].weights[0];
        }
    }
#ifdef KSS_DEVELOPMENT
    if (s.max_value != dp[isiz][problem.knapsack_sizes[0]])
        throw std::logic_error("back tracking went wrong! Got: " + std::to_string(s.max_value) + ", expected: " + std::to_string(dp[isiz][problem.knapsack_sizes[0]]));
#endif

    return s;
}




//---------- BRUTE FORCE SOLVER ----------

bool BruteForceSolver::expect_perfection = true;

string BruteForceSolver::GetAlgorithmName(const Options & options){
    string name = "brute-force_";
    if (options.iterative) {
        name += "iterative_";
        name += ToString(options.search_order);
    }
    else {
        name += "recursive_";
        name += options.late_fit ? "late-fit" : "early-fit";
    }
    return name;
}

PackagedSolution BruteForceSolver::Solve(PackagedProblem & problem, const Options options){    
    PackagedSolution ps;
    
    // run with time measure
    std::chrono::steady_clock::time_point start, end;
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

    return ps;
}

Solution BruteForceSolver::SolutionFromNumber(int num, const Problem & problem) {
    // Vector to store the number
    Solution solution(problem.items.size(), problem.knapsack_sizes);
 
    int i = problem.items.size();
    while (num != 0) {
        --i;
        if (num % 2 == 1) solution.AddItem(problem, i, Solution::FaultTreatment::INVALIDATE);
 
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
            if (temp_solution.valid && temp_solution.IsStructure(problem) && temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;
        }
        break;

    case Options::SearchOrder::ONE_FIRST:
        for (int i = n - 1; i >= 0; --i) {
            Solution temp_solution = SolutionFromNumber(i, problem.problem);
            if (temp_solution.valid && temp_solution.IsStructure(problem) && temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;
        }
        break;

    default:
        throw std::logic_error("not implemented yet");
        break;
    }

    return global_solution;
}

Solution BruteForceSolver::Max(const Solution & a, const Solution & b, const PackagedProblem & problem){
    bool va = a.valid && a.IsStructure(problem);
    bool vb = b.valid && b.IsStructure(problem);
    if (va && !vb) return a;
    if (!va && vb) return b;
    else if (a.max_value > b.max_value) return a;
    else return b;
}

Solution BruteForceSolver::DFS(const PackagedProblem & problem, Solution sol0, const Options::SearchOrder & search_order, const bool & add, const int & depth){
    if (add) sol0.AddItem(problem.problem, depth, Solution::FaultTreatment::INVALIDATE);
    if (depth == problem.problem.items.size() - 1) return sol0;

    return Max(DFS(problem, sol0, search_order, false, depth + 1), DFS(problem, sol0, search_order, true, depth + 1), problem);
}

Solution BruteForceSolver::DFS(const PackagedProblem & problem, Solution sol0, const Options::SearchOrder & search_order, const int & depth){
    if (depth == problem.problem.items.size()) return sol0;
    
    Solution sol1 = sol0;
    sol0 = DFS(problem, sol0, search_order, depth + 1);
    sol1.AddItem(problem.problem, depth, Solution::FaultTreatment::INVALIDATE);
    sol1 = DFS(problem, sol1, search_order, depth + 1);

    return Max(sol0, sol1, problem);
}

Solution BruteForceSolver::Recursive(const PackagedProblem & problem, const Options & options){
    if (options.search_order != BruteForceSolver::Options::SearchOrder::UNCONSTRAINED)
        throw std::invalid_argument("recursive does not support search order"); //it could be supported, but it slows algorithm down
    
    int isiz = problem.problem.items.size();
    Solution sol0(isiz, problem.problem.knapsack_sizes);
    
    if (options.late_fit) {
        return Max(DFS(problem, sol0, options.search_order, false, 0), DFS(problem, sol0, options.search_order, true, 0), problem);
    }
    else{
        Solution sol1 = sol0;
        sol0 = DFS(problem, sol0, options.search_order, 1);
        sol1.AddItem(problem.problem, 0, Solution::FaultTreatment::INVALIDATE);
        sol1 = DFS(problem, sol1, options.search_order, 1);
        return Max(sol0, sol1, problem);
    }

}



//---------- GREEDY SOLVER ----------

bool GreedySolver::expect_perfection = false;

string GreedySolver::GetAlgorithmName(const Options & options){
    string name = "";
    switch (options.sort_mode)
    {
    case Problem::SortMode::RANDOM:
        name += "random-fit";
        break;
    
    case Problem::SortMode::DONT_SORT:
        name += "first-fit";
        break;

    default:
        name += "greedy-sort-by-" + ToString(options.sort_mode);
        break;
    }
    name += "_";
    name += options.multi_run ? "multi-run" : "single-run";
    return name;
}

PackagedSolution GreedySolver::Solve(PackagedProblem & problem, const Options & options){
    PackagedSolution ps;
    
    // run with time measure
    std::chrono::steady_clock::time_point start, end;
    start = std::chrono::steady_clock::now();
    ps.solution = NaiveUniversal(problem, options);
    end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();

    return ps;
}

Solution GreedySolver::NaiveUniversal(const PackagedProblem & problem, const Options & options) {
    int isiz = problem.problem.items.size();
    Solution global_solution(isiz, problem.problem.knapsack_sizes); // create empty solution
    vector<int> sortedItemIds = problem.problem.GetSortedItemIds(options.sort_mode);
    
    for (int i = 0; i < isiz; ++i){
        int current_item_id = sortedItemIds[i];
        
        if (global_solution.Fits(problem.problem, current_item_id)){
            global_solution.AddItem(problem.problem, current_item_id); // if current item fits add it to the current solution
            if (!global_solution.IsStructurePossible(problem)) global_solution.RemoveItem(problem.problem, current_item_id); // if adding it would make it impossible to construct a valid answer then remove it
        }
    }
    
    return global_solution;
}

Solution GreedySolver::GreedyUniversal(const PackagedProblem & problem, const Options & options) {
    int isiz = problem.problem.items.size();
    Solution global_solution(isiz, problem.problem.knapsack_sizes); // create empty solution
    vector<int> sortedItemIds = problem.problem.GetSortedItemIds(options.sort_mode);

    switch (problem.requirements.structureToFind)
    {

    case Problem::Requirements::StructureToFind::PATH:
    {
        if (options.multi_run){ // ---------- PATH MULTI RUN ----------
            // Find best starting point
            int startItem = -1;
            for (int current_item_id : sortedItemIds) {
                if (global_solution.Fits(problem.problem, current_item_id)){
                    global_solution.AddItem(problem.problem, current_item_id); // if current item fits add it to the current solution
                    startItem = current_item_id;
                    break;
                }
            }
            if (startItem == -1) return global_solution; // nothing fits

            bool going = true;
            while (going){
                going = false;
                // find next item
                for (auto it = sortedItemIds.begin(); it != sortedItemIds.end(); ++it) {
                    int current_item_id = *it;
                    if (!global_solution.selected[current_item_id] &&
                        problem.problem.items[startItem].HasConnectionTo(current_item_id) &&
                        global_solution.Fits(problem.problem, current_item_id))
                    {
                        global_solution.AddItem(problem.problem, current_item_id);
                        startItem = current_item_id;
                        going = true;
                        break;
                    }
                }
            }
        }
        else { // ---------- PATH SINGLE RUN ----------
            for (int i = 0; i < isiz; ){
                int current_item_id = sortedItemIds[i];

                if (global_solution.Fits(problem.problem, current_item_id))
                    global_solution.AddItem(problem.problem, current_item_id); // if current item fits add it to the current solution
                else return global_solution; // if does not fit return what have we managed to fit so far

                i = -1;
                int nextItem;
                do{
                    if (++i >= isiz) return global_solution; // if can't find next possible item return what we had so far
                    nextItem = sortedItemIds[i];
                } while (global_solution.selected[nextItem] || !problem.problem.items[current_item_id].HasConnectionTo(nextItem));
            }
        }
        break;
    }

    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS: {
        for (int i = 0; i < isiz; ++i){
            int current_item_id = sortedItemIds[i];

            // if current item fits add it to the current solution
            global_solution.AddItemIfFits(problem.problem, current_item_id);
        }
        break;
    }
    
    case Problem::Requirements::StructureToFind::CYCLE: {
        for (int i = 0; i < isiz; ++i){
            int current_item_id = sortedItemIds[i];
            
            if (global_solution.Fits(problem.problem, current_item_id)){
                global_solution.AddItem(problem.problem, current_item_id); // if current item fits add it to the current solution
                if (!global_solution.IsCyclePossible(problem.problem)) global_solution.RemoveItem(problem.problem, current_item_id); // if adding it would make it impossible to construct a valid answer then remove it
            }
        }
        break;
    }
    
    default: {
        throw std::logic_error("not implemented yet");
        break;
    }
    }

    return global_solution;
}

Solution GreedySolver::GreedyIgnoreConnections(const Problem & problem, const Options & options){
    int isiz = problem.items.size();
    Solution global_solution(isiz, problem.knapsack_sizes); // create empty solution
    vector<int> sortedItemIds = problem.GetSortedItemIds(options.sort_mode);

    for (int i = 0; i < isiz; ++i){
        int current_item_id = sortedItemIds[i];

        // if current item fits add it to the current solution
        global_solution.AddItemIfFits(problem, current_item_id);
    }

    return global_solution;
}

Solution GreedySolver::GreedyPath(const Problem & problem, const Options & options) {
    int isiz = problem.items.size();
    Solution global_solution(isiz, problem.knapsack_sizes); // create empty solution
    vector<int> sortedItemIds = problem.GetSortedItemIds(options.sort_mode);

    for (int i = 0; i < isiz; ){
        int current_item_id = sortedItemIds[i];

        if (global_solution.Fits(problem, current_item_id))
            global_solution.AddItem(problem, current_item_id); // if current item fits add it to the current solution
        else return global_solution; // if does not fit return what have we managed to fit so far

        i = -1;
        int nextItem;
        do{
            if (++i >= isiz) return global_solution; // if can't find next possible item return what we had so far
            nextItem = sortedItemIds[i];
        } while (global_solution.selected[nextItem] || !problem.items[current_item_id].HasConnectionTo(nextItem));
    }

    return global_solution;
}




//---------- BRANCH AND BOUND SOLVER ----------

bool BranchAndBoundSolver::expect_perfection = true;

string BranchAndBoundSolver::GetAlgorithmName(const Options & options){
    string name = "branch-and-bound_";
    if (options.late_fit) name += "late-fit";
    else name += "early-fit";
    return name;
}

PackagedSolution BranchAndBoundSolver::Solve(PackagedProblem & problem, const Options & options){
    if (problem.requirements.structureToFind != Problem::Requirements::StructureToFind::PATH)
        throw std::logic_error("branch and bound does not support structures other than path");

    PackagedSolution ps;

    // solve with measured time
    std::chrono::steady_clock::time_point start, end;
    if (options.late_fit){
        start = std::chrono::steady_clock::now();
        ps.solution = BnBLateFitPath(problem.problem);
        end = std::chrono::steady_clock::now();
    }
    else{
        start = std::chrono::steady_clock::now();
        ps.solution = BnBEarlyFitPath(problem.problem);
        end = std::chrono::steady_clock::now();
    }
    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();

    return ps;
}

/*
Solution BranchAndBoundSolver::DFS(const PackagedProblem & problem, Solution sol0, const bool & add, const int & depth, int lower_bound, const Options & options){
    if (add){
        if (!sol0.Fits(problem.problem, depth)) return sol0;
        sol0.AddItem(problem.problem, depth);
    }
    if (depth == problem.problem.items.size() - 1) return sol0; // just added (or not) last item return what we have there will be no further DFS calls

    int up = options.upper_bound(problem.problem, sol0); // we just added an item (or not) to the solution what is possible theoretically to find in this barnch
    if (up <= lower_bound) return sol0; // if we cant get more than what we found already (or know we can find elsewhere) dont chceck further

    Solution sol1 = DFS(problem, sol0, true, depth + 1, options.lower_bound(problem.problem, sol0), options);
    lower_bound = std::max(lower_bound, options.lower_bound(problem.problem, sol1)); // what can we find goind sol1 route

    if (up <= lower_bound) return sol1; // if we cant match that 
    Solution sol2 = DFS(problem, sol0, false, depth + 1, options.lower_bound(problem.problem, sol0), options);

    if (sol1.max_value > sol2.max_value) return sol1;
    return sol2;
}*/

//---------- Early Fit ----------

Solution BranchAndBoundSolver::DFSEarlyFitPath(const Problem & problem, Solution currentSolution, const int & current_item_id) {

#ifdef DEBUG_BNB
    cout << "check for " << current_item_id << ", remaining weight: "<< remainingSpace[0] << endl;
#endif

    Solution global_solution = currentSolution; // in worst case we will return what we have so far (it is a valid solution since we adding items in path fashion)

    // construct solution further with recursive DFSLateFitPath calls
    for (int nextItemId : problem.items[current_item_id].connections){
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        if (currentSolution.Fits(problem, nextItemId)){
            Solution tempSolution = currentSolution;

            tempSolution.AddItem(problem, nextItemId);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::VALUE_WEIGHT_RATIO, tempSolution, remainingSpace) > lower_bound*/) {
                tempSolution = DFSEarlyFitPath(problem, tempSolution, nextItemId);
                if (tempSolution.max_value > global_solution.max_value)  global_solution = tempSolution;
            }
        }
    }

#ifdef DEBUG_BNB
    cout << "exit with: " << global_solution.max_value << endl;
#endif
    return global_solution; // return the best solution from all paths starting from current_item_id
}

Solution BranchAndBoundSolver::BnBEarlyFitPath(const Problem & problem) {
    int isiz = problem.items.size();
    Solution global_solution(isiz, problem.knapsack_sizes); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){

        if (global_solution.Fits(problem, i)){
            Solution currentSolution(isiz, problem.knapsack_sizes); // create empty solution
            currentSolution.AddItem(problem, i);
            if (true/*GreedyIgnoreConnections(problem, Problem::SortMode::VALUE_WEIGHT_RATIO, currentSolution, remainingSpace) > global_solution.max_value*/){
                currentSolution = DFSEarlyFitPath(problem, currentSolution, i);
                if(currentSolution.max_value > global_solution.max_value) global_solution = currentSolution; // if newly found solution is better use it as you global solution
            }
        }
    }

    return global_solution;
}

//---------- Late Fit ----------

Solution BranchAndBoundSolver::DFSLateFitPath(const Problem & problem, Solution currentSolution, const int & current_item_id) {

    if (!currentSolution.AddItemIfFits(problem, current_item_id)) // if current item fits add it to the current solution
        return currentSolution; // if does not fit return what have we managed to fit so far

#ifdef DEBUG_BNB_LF
    cout << "check for " << current_item_id << ", remaining weight: "<< currentSolution.remainingSpace[0] << endl;
#endif

    
    Solution global_solution = currentSolution; // in worst case we will return what we have so far (it is a valid solution since we adding items in path fashion)

    // construct solution further with recursive DFSLateFitPath calls
    for (int nextItemId : problem.items[current_item_id].connections){
        if (currentSolution.selected[nextItemId]) continue; //don't try adding if next item is already in the solution

        Solution tempSolution = DFSLateFitPath(problem, currentSolution, nextItemId);
        if (tempSolution.max_value > global_solution.max_value) global_solution = tempSolution;
    }

#ifdef DEBUG_BNB_LF
    cout << "exit with: " << global_solution.max_value << endl;
#endif
    return global_solution; // return the best solution from all paths starting from current_item_id
}

Solution BranchAndBoundSolver::BnBLateFitPath(const Problem & problem) {
    int isiz = problem.items.size();
    Solution global_solution(isiz, problem.knapsack_sizes); // create empty solution

    // check possible solutions from every starting point
    for (int i = 0; i < isiz; ++i){
        Solution currentSolution(isiz, problem.knapsack_sizes); // create empty solution 

        currentSolution = DFSLateFitPath(problem, currentSolution, i);
        if(currentSolution.max_value > global_solution.max_value) global_solution = currentSolution; // if newly found solution is better use it as you global solution
    }

    return global_solution;
}




//---------- GRASP SOLVER ----------

bool GRASPSolver::expect_perfection = false;

string GRASPSolver::GetAlgorithmName(const Options & options){
    string name = "grasp-" + ToString(options.sort_mode);
    name += "_chose-from-" + std::to_string(options.chose_from);

    if (options.iterations > 0) {
        if (options.coverage > 0.0) throw std::invalid_argument("both iteratations and coverage are positive");
        else name += "_iterations-" + std::to_string(options.iterations);
    }
    else {
        if (options.coverage > 0.0) name += "_coverage-" + std::to_string(options.coverage);
        else throw std::invalid_argument("both iteratations and coverage are non positive");
    }

    return name;
}

PackagedSolution GRASPSolver::Solve(PackagedProblem & problem, Options options){
    PackagedSolution ps;

    if (options.iterations > 0) {
        if (options.coverage > 0.0) throw std::invalid_argument("both iteratations and coverage are positive");
    }
    else {
        if (options.coverage > 0.0) options.iterations = static_cast<double>(pow(2, problem.problem.items.size())) * options.coverage;
        else throw std::invalid_argument("both iteratations and coverage are non positive");
    }
    
    // run with time measure
    std::chrono::steady_clock::time_point start, end;

    if (options.iterations == 1){
        start = std::chrono::steady_clock::now();
        ps.solution = Universal(problem, options);
        end = std::chrono::steady_clock::now();
    }
    else {
        start = std::chrono::steady_clock::now();
        ps.solution = MultiRun(problem, options);
        end = std::chrono::steady_clock::now();
    }

    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();

    return ps;
}

Solution GRASPSolver::Universal(const PackagedProblem & problem, const Options & options) {
    Solution global_solution(problem.problem.items.size(), problem.problem.knapsack_sizes);
    vector<int> sorted_item_ids = problem.problem.GetSortedItemIds(options.sort_mode);
    int amount_to_chose_from = problem.problem.items.size() * options.chose_from;
    while (sorted_item_ids.size() > 0) {
        int pick = std::rand() % std::min(amount_to_chose_from, static_cast<int>(sorted_item_ids.size()));
        int current_item_id = sorted_item_ids[pick];

        if (global_solution.Fits(problem.problem, current_item_id)){
            global_solution.AddItem(problem.problem, current_item_id); // if current item fits add it to the current solution
            if (!global_solution.IsStructurePossible(problem)) global_solution.RemoveItem(problem.problem, current_item_id); // if adding it would make it impossible to construct a valid answer then remove it
        }

        sorted_item_ids.erase(sorted_item_ids.begin() + pick);
    }
    return global_solution;
}
Solution GRASPSolver::MultiRun(const PackagedProblem & problem, const Options & options) {
    Solution global_solution(problem.problem.items.size(), problem.problem.knapsack_sizes);

    for (int i = 0; i < options.iterations; ++i){
        Solution temp_solution = Universal(problem, options);
        if (temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;
    }
    return global_solution;
}




//---------- GREEDY HEURISTIC SEARCH SOLVER ----------

bool GreedyHeuristicSearchSolver::expect_perfection = false;

string GreedyHeuristicSearchSolver::GetAlgorithmName(const Options & options){
    string name = "greedy-heuristic-search-sort-by-" + ToString(options.sort_mode);
    name += "_coverage-" + std::to_string(options.coverage);
    return name;
}

PackagedSolution GreedyHeuristicSearchSolver::Solve(PackagedProblem & problem, const Options & options){
    PackagedSolution ps;
    
    // run with time measure
    std::chrono::steady_clock::time_point start, end;
    switch (problem.requirements.structureToFind) {

    case Problem::Requirements::StructureToFind::IGNORE_CONNECTIONS: {
        start = std::chrono::steady_clock::now();
        ps.solution = IgnoreConnectionsDFS(problem, Solution(problem.problem.items.size(), problem.problem.knapsack_sizes), problem.problem.GetSortedItemIds(options.sort_mode), problem.problem.items.size() * options.coverage);
        end = std::chrono::steady_clock::now();
        break;
    }
    case Problem::Requirements::StructureToFind::PATH: {
        start = std::chrono::steady_clock::now();
        ps.solution = PathDFS(problem, Solution(problem.problem.items.size(), problem.problem.knapsack_sizes), -1, problem.problem.GetSortedItemIds(options.sort_mode), problem.problem.items.size() * options.coverage);
        end = std::chrono::steady_clock::now();
        break;
    }
    case Problem::Requirements::StructureToFind::CYCLE: {
        start = std::chrono::steady_clock::now();
        ps.solution = CycleDFS(problem, Solution(problem.problem.items.size(), problem.problem.knapsack_sizes), -1, -1, problem.problem.GetSortedItemIds(options.sort_mode), problem.problem.items.size() * options.coverage);
        end = std::chrono::steady_clock::now();
        break;
    }
    
    default: {
        throw std::logic_error("not implemented yet");
        break;
    }
    }

    const std::chrono::duration<double> elapsed_seconds{end - start};
    ps.solve_time = std::chrono::duration<double>(elapsed_seconds).count();

    return ps;
}

Solution GreedyHeuristicSearchSolver::PathDFS(const PackagedProblem & problem, const Solution & current_solution, const int & previous_item_id, const vector<int> & sorted_item_ids, const int & amount_to_visit) {
    Solution global_solution(current_solution);
    int counted_visits = 0;
    int total_visits = 0;
    while (counted_visits < amount_to_visit && total_visits < sorted_item_ids.size()) {
        ++total_visits;

        int current_item_id = sorted_item_ids[total_visits-1];
        if (current_solution.selected[current_item_id]) continue; // if item is already selected dont count as visit otherwise algorithm would add at most amount_to_visit items to every solution it finds
        if (previous_item_id >= 0 && !problem.problem.items[previous_item_id].HasConnectionTo(current_item_id)) continue;
        if (!current_solution.Fits(problem.problem, current_item_id)) continue; // if item does not fit dont count as visit otherwise algorithm would give empty answer if first amout_to_visit items dont fit

        Solution temp_solution(current_solution);
        temp_solution.AddItem(problem.problem, current_item_id);
        vector<int> new_sorted_item_ids(sorted_item_ids);
        new_sorted_item_ids.erase(new_sorted_item_ids.begin() + total_visits - 1);
        temp_solution = PathDFS(problem, temp_solution, current_item_id, new_sorted_item_ids, amount_to_visit);
        if (temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;

        ++counted_visits;
    }
    return global_solution;
}

Solution GreedyHeuristicSearchSolver::CycleDFS(const PackagedProblem & problem, const Solution & current_solution, const int & previous_item_id, const int & start_item_id, const vector<int> & sorted_item_ids, const int & amount_to_visit) {
    Solution global_solution(current_solution);
    int counted_visits = 0;
    int total_visits = 0;
    while (counted_visits < amount_to_visit && total_visits < sorted_item_ids.size()) {
        ++total_visits;

        int current_item_id = sorted_item_ids[total_visits-1];
        if (current_solution.selected[current_item_id]) continue; // if item is already selected dont count as visit otherwise algorithm would add at most amount_to_visit items to every solution it finds
        if (!current_solution.Fits(problem.problem, current_item_id)) continue; // if item does not fit dont count as visit otherwise algorithm would give empty answer if first amout_to_visit items dont fit

        Solution temp_solution(current_solution);
        temp_solution.AddItem(problem.problem, current_item_id);
        if (!temp_solution.IsCyclePossible(problem.problem)) continue;
        vector<int> new_sorted_item_ids(sorted_item_ids);
        new_sorted_item_ids.erase(new_sorted_item_ids.begin() + total_visits - 1);
        temp_solution = CycleDFS(problem, temp_solution, current_item_id, start_item_id >= 0 ? start_item_id : current_item_id, new_sorted_item_ids, amount_to_visit);
        if (temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;

        ++counted_visits;
    }
    return global_solution;
}

Solution GreedyHeuristicSearchSolver::IgnoreConnectionsDFS(const PackagedProblem & problem, const Solution & current_solution, const vector<int> & sorted_item_ids, const int & amount_to_visit) {
    Solution global_solution(current_solution);
    int counted_visits = 0;
    int total_visits = 0;
    while (counted_visits < amount_to_visit && total_visits < sorted_item_ids.size()) {
        ++total_visits;

        int current_item_id = sorted_item_ids[total_visits-1];
        if (current_solution.selected[current_item_id]) continue; // if item is already selected dont count as visit otherwise algorithm would add at most amount_to_visit items to every solution it finds
        if (!current_solution.Fits(problem.problem, current_item_id)) continue; // if item does not fit dont count as visit otherwise algorithm would give empty answer if first amout_to_visit items dont fit

        Solution temp_solution(current_solution);
        temp_solution.AddItem(problem.problem, current_item_id);
        vector<int> new_sorted_item_ids(sorted_item_ids);
        new_sorted_item_ids.erase(new_sorted_item_ids.begin() + total_visits - 1);
        temp_solution = IgnoreConnectionsDFS(problem, temp_solution, new_sorted_item_ids, amount_to_visit);
        if (temp_solution.max_value > global_solution.max_value) global_solution = temp_solution;

        ++counted_visits;
    }
    return global_solution;
}




//---------- FLOYD SOLVER ----------

PackagedSolution FloydSolver::Solve(PackagedProblem & problem){
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

    // validate
    ps.validation_status = Validation::Validate(ps.solution, problem);

    // optimum update
    if (ps.validation_status.undergone && ps.validation_status.valid && ps.solution.max_value > problem.known_optimum){
        problem.known_optimum = ps.solution.max_value;
        if (problem.associated_file != "") problem.ExportJSON(problem.associated_file);
    }
    if (ps.solution.max_value == 0) ps.quality = DBL_MAX;
    else ps.quality = static_cast<double>(problem.known_optimum) / static_cast<double>(ps.solution.max_value);

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



//---------- ENUMS ----------

std::string ToString(const knapsack_solver::BruteForceSolver::Options::SearchOrder & so){
    switch (so)
    {
    case BruteForceSolver::Options::SearchOrder::ZERO_FIRST:
        return "zero-first";
        break;
    case BruteForceSolver::Options::SearchOrder::ONE_FIRST:
        return "one-first";
        break;
    case BruteForceSolver::Options::SearchOrder::UNCONSTRAINED:
        return "unconstrained";
        break;
    case BruteForceSolver::Options::SearchOrder::GRAY_CODE:
        return "gray-code";
        break;   
    case BruteForceSolver::Options::SearchOrder::RANDOM:
        return "random";
        break;   
    default:
        throw std::invalid_argument("unrecognised value");
        break;
    }

}
std::string ToString(const knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction & bf){
    switch (bf)
    {
    case BranchAndBoundSolver::Options::BoundingFunction::NONE:
        return "none";
        break;
    case BranchAndBoundSolver::Options::BoundingFunction::BASE_DYNAMIC:
        return "dynamic";
        break;
    case BranchAndBoundSolver::Options::BoundingFunction::CONTINOUS:
        return "continous";
        break;    
    default:
        throw std::invalid_argument("unrecognised value");
        break;
    }
}

knapsack_solver::BruteForceSolver::Options::SearchOrder ToSearchOrder(const std::string & str){
    if (str == "zero-first") return BruteForceSolver::Options::SearchOrder::ZERO_FIRST;
    else if (str == "one-first") return BruteForceSolver::Options::SearchOrder::ONE_FIRST;
    else if (str == "unconstrained") return BruteForceSolver::Options::SearchOrder::UNCONSTRAINED;
    else if (str == "gray-code") return BruteForceSolver::Options::SearchOrder::GRAY_CODE;
    else if (str == "random") return BruteForceSolver::Options::SearchOrder::RANDOM;
    else throw std::invalid_argument("unrecognised value");
}
knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction ToBoundingFunction(const std::string & str){
    if (str == "none") return BranchAndBoundSolver::Options::BoundingFunction::NONE;
    else if (str == "dynamic") return BranchAndBoundSolver::Options::BoundingFunction::BASE_DYNAMIC;
    else if (str == "continous") return BranchAndBoundSolver::Options::BoundingFunction::CONTINOUS;
    else throw std::invalid_argument("unrecognised value");
}

std::ostream& operator<<(std::ostream & os, const knapsack_solver::BruteForceSolver::Options::SearchOrder & so){
    return os << ToString(so);
}
std::ostream& operator<<(std::ostream & os, const knapsack_solver::BranchAndBoundSolver::Options::BoundingFunction & bf){
    return os << ToString(bf);
}