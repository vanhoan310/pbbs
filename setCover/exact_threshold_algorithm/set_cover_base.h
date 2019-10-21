#ifndef SCB_H
#define SCB_H

#include <string>
#include <vector>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream
#include <math.h>
#include <algorithm>
using namespace std;

class SetCoverBase
{
private:
    int n_sets; // number of sets
    int n_clusters;
    vector<vector<int>> hit_sets; // store set of {A_j | x_i \in A_j}
    
public:
    int max_elem; // number of elements in the universe
    // constructor 
    SetCoverBase (std::vector<std::vector<int>> & clusters, int _clusters);
    // destructor
    ~SetCoverBase();

    vector<vector<int>> thresholdGraph(vector<vector<double>> & distance_mat, double L, bool only_kcenter=false); // return set cover 
    
};


#endif 
