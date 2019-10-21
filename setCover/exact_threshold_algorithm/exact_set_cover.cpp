#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using namespace std;
#include <string>
#include <vector>
#include <stack>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream
#include <math.h>
#include <algorithm>

int main() {
    
    
    // cover instance
    int n_sets = 10;
    vector<vector<int>> covers;
    vector<int> v1 = {0,1,2,4};
    vector<int> v2 = {2,3,5,7};
    covers.push_back(v1);
    covers.push_back(v2);
    
    
    // set cover model
    IloEnv env;
    IloModel model(env);
    // variables
    IloNumVarArray x1(env, n_sets, 0, 1, ILOINT);
    // constraints
    for (auto &v: covers){
        IloExpr cv_cons(env);
        for (auto &i: v){
            cv_cons += x1[i];
        }
        model.add(cv_cons >= 1);
    }
    
    //obj 
    IloExpr cost_sum(env);
    for (int i=0; i < n_sets; i++)
    {
        cost_sum += x1[i];
    }
    IloObjective obj = IloMinimize(env, cost_sum);
    model.add(obj);
    cost_sum.end();
    
    //solve model 
    IloCplex cplex(model);
    cplex.solve();
    cout << "Min = " << cplex.getObjValue() << endl;
    for (int i=0; i < n_sets; i++)
    {
        cout << cplex.getValue(x1[i]) << " ";
    }
    cout << endl;
    
    env.end();
    
    return 0;
}
