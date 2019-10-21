#include "set_cover_base.h"
#include "util_functions.h"
#include<algorithm>
using namespace std;

SetCoverBase::SetCoverBase (std::vector<std::vector<int>> & clusters, int _clusters)
{
    n_clusters = _clusters;
    n_sets = (int) clusters[0].size();
    // convert cluster indicator to cluster of objects 
    vector<vector<int>> item_clusters = clusters; // note that each row of 'clusters' corresponds to a clustering, e.g clusters = [[1,1,2,2], [1,1,1,2]]
    int cluster_id = n_sets;
    for (auto &v: item_clusters)
    {
        for (int k = 1; k <= n_clusters; k++)
        {
            if ( std::find(v.begin(), v.end(), k) != v.end())
            {
                std::replace(v.begin(), v.end(), k, cluster_id);
                cluster_id++;
            }
        }
    }

    max_elem = cluster_id;
    hit_sets = transpose(item_clusters);
}

vector<vector<int>> SetCoverBase::thresholdGraph(vector<vector<double>> & distance_mat, double L, bool only_kcenter)
{
    if (only_kcenter)
    {
        cout << "Clear hit_sets (pure k center without using hitting constraints)\n";
        max_elem = n_sets;
        for (int i=0; i < n_sets; ++i)
        {
            hit_sets[i].clear();
            assert(hit_sets[i].size() == 0);
        }
    }
    vector<vector<int>> cover_set = hit_sets; //take init of {A_j| x_i \in A_j}
    for (int i=0; i < n_sets; ++i)
    {
        for (int j=0; j < n_sets; ++j)
        {
            if (distance_mat[i][j] <= L)
            {
                cover_set[i].push_back(j); // [x_i] take x_i and its neighbors 
            }
        }
    }
    return cover_set;
}

SetCoverBase::~SetCoverBase(){}
