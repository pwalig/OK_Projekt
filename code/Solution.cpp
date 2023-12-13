#include "Solution.hpp"
#include "json.hpp"

#include <fstream>
#include <iostream>

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