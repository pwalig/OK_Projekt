# OK_Projekt
Solver for multiple Knapsack problem variants.

> [!WARNING]
> Project discontinued!
> 
> See the successor: [graphsack](https://github.com/pwalig/graphsack)

## Status

Project discontinued in favour of [graphsack](https://github.com/pwalig/graphsack).

### Done
- Brute Force Methods
- Primitive Branch and Bound Methods
- Greedy algorithms
- Dynamic Programming
- measure calculation time
- json format support
- Heuristics

### Missing
- Multi threading
- GPU support

## Problems solved by OK_Project
- Multi dimentional knapsack problem, but items are given as a graph, where each vertex in the graph is an item to pack. The problem is to find either path, cycle or tree in this graph, that fits in the knapsack and has maximum value.
- [Knapsack Problem](https://en.wikipedia.org/wiki/Knapsack_problem)
- [Multi Dimentional Knapsack Problem](https://en.wikipedia.org/wiki/Knapsack_problem#Multi-dimensional_knapsack_problem)

## Using

OK_Project uses:
- C++ language
- [nlohmann/json](https://github.com/nlohmann/json) JSON parser for Modern C++

## Usage

For full usage and building instructions checkout [wiki](https://github.com/pwalig/OK_Projekt/wiki).

## Building

```
git clone https://github.com/pwalig/OK_Project.git
cd OK_Project
g++ -std=c++17 -g code/*cpp -o knapsack-solver.exe
```
