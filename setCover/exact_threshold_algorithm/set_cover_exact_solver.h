#ifndef SCES_H
#define SCES_H

#include <string>
#include <vector>
#include <stack>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream
#include <math.h>
#include <algorithm>
#include "util_functions.h"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

class SetCoverExactSolver
{
// note that X = {0, 1, ..., n_sets - 1} 
private:
    int n_sets; // number of sets
    IloEnv  env;        // store simplex environment
    IloNumVarArray x1;  // indicator for set cover 
    
    IloConstraintArray cover_constraints;
    
    IloModel model;
    IloCplex cplex;
    IloObjective obj;  // objective function.
    
public:
    // constructor 
    SetCoverExactSolver (int n_sets);
    // destructor
    ~SetCoverExactSolver();
    
    // create a element x set matrix where each entry (=1) correspond to if element is belong to a set 
    std::vector<std::vector<int>> create_set_covers_mat(std::vector<std::vector<int>> & cover_set, int n_elem);

    void solveSetCover(std::vector<std::vector<int>> & cover_set, int n_elem);

    double getObjValue();
    
    std::vector<int> getIndicatorSolution();
    
    void printSolution();
    
};


#endif 
