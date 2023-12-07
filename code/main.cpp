#include <iostream> // cout, endl

#include <cstdlib> // std::srand()
#include <ctime> // std::time()

#include "CommandInterpreter.hpp"

#define KSS_DEVELOPMENT

using namespace knapsack_solver;
using std::cout;
using std::endl;

int main(int argc, char * argv[]) {
    std::srand(std::time(0) * 1000);
    
#ifdef KSS_DEVELOPMENT
    /*cout << "Random value on [0, " << RAND_MAX << "]: " << std::rand() << '\n';*/
    if (argc == 1) CommandInterpreter::RunTasks("../tests/commands/tasks1.json");
    else CommandInterpreter::InterpretCommand(argc, argv);

#else

    CommandInterpreter::InterpretCommand(argc, argv);

#endif

    return 0;
}