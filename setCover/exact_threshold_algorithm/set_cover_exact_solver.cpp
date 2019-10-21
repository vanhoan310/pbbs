#include "set_cover_exact_solver.h"
using namespace std;

SetCoverExactSolver::SetCoverExactSolver(int n_sets): model(env), cplex(model), x1(env, n_sets, 0, 1, ILOINT), cover_constraints(env), n_sets(n_sets) 
{
    //obj 
    IloExpr cost_sum(env);
    for (int i=0; i < n_sets; i++)
    {
        cost_sum += x1[i];
    }
    IloObjective obj = IloMinimize(env, cost_sum);
    model.add(obj);
    cost_sum.end();
}

std::vector<std::vector<int>> SetCoverExactSolver::create_set_covers_mat(std::vector<std::vector<int>> & cover_set, int n_elem)
{
   vector<vector<int>> results(n_elem, vector<int>(n_sets));
   int idx = 0;
   for (auto & v: cover_set){
       for (auto & j: v)
       {
           results[j][idx] = 1;
       }
       idx++;
   }
   return results;
}
void SetCoverExactSolver::solveSetCover(std::vector<std::vector<int>> & cover_set, int n_elem)
{
    // remove existing cover constraints (if exists)
    for (IloInt i = 0; i < cover_constraints.getSize(); i++)
    {
        model.remove(cover_constraints[i]);
    }
    
    // Add cover constraints 
    // get elements x sets matrix where (i,j)=1 iff element i belongs set C_j, 0 otherwise
    vector<vector<int>> set_cover_mat = create_set_covers_mat(cover_set, n_elem);
    // for (auto v: set_cover_mat)
    //     print_vec(v);

    for (auto &v: set_cover_mat){
        IloExpr cv_cons(env);
        for (int i=0; i < n_sets; ++i){
            cv_cons += v[i]*x1[i];
        }
        IloConstraint cons = cv_cons >= 1;
        model.add(cons); 
        cover_constraints.add(cons);
    }

    // Solve model 
    cplex.setOut(env.getNullStream());
    try 
    {
        if ( !cplex.solve() ) {
            env.error() << "Failed to optimize ILP" << endl;
            throw(-1);
        }
    }   catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
}

double SetCoverExactSolver::getObjValue()
{
    return cplex.getObjValue();
}

std::vector<int> SetCoverExactSolver::getIndicatorSolution()
{
    std::vector<int> v;
    for (int i=0; i < n_sets; i++)
    {
        if (cplex.getValue(x1[i]) > 0)
        {
            v.push_back(1);
        }
        else 
        {
            v.push_back(0);
        }
    }
    return v;
}

void SetCoverExactSolver::printSolution()
{
    cout << "Min cover size = " << cplex.getObjValue() << endl;
    for (int i=0; i < n_sets; i++)
    {
        cout << cplex.getValue(x1[i]) << " ";
    }
    cout << endl;
}

SetCoverExactSolver::~SetCoverExactSolver()
{
	env.end();
}
