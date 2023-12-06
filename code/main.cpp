#include <iostream> // cout, endl
#include <vector>

#include <cstdlib> // std::srand()
#include <ctime> // std::time()

#include "Problem.hpp"
#include "KnapsackSolver.hpp"
#include "CommandInterpreter.hpp"

#define KSS_DEVELOPMENT

using namespace knapsack_solver;
using std::vector;
using std::cout;
using std::endl;

int main(int argc, char * argv[]) {
    
#ifdef KSS_DEVELOPMENT

    /*std::srand(std::time(0) * 1000);
    cout << "Random value on [0, " << RAND_MAX << "]: " << std::rand() << '\n';*/
    if (argc == 1) CommandInterpreter::RunTasks("../tests/commands/tasks1.json");
    else CommandInterpreter::InterpretCommand(argc, argv);

#else

    CommandInterpreter::InterpretCommand(argc, argv);

#endif

    return 0;
}