#include<iostream>
#include<vector>
#include "util_functions.h"
#include "set_cover_base.h"
#include "set_cover_exact_solver.h"
#include <chrono>

using namespace std;
int main(int argc, char *argv[])
{
    std::string dataname = argv[1]; // zeisel grun
    vector<int> samples_size = {1, 2, 4, 6, 8, 10};
    // vector<int> samples_size = {50};
    std::string clusterfilename = "../../data/"+dataname+"_clusters.csv";
    vector<vector<int>> clusters = parse2DCsvFile(clusterfilename);

    std::string dfilename = "../../data/"+dataname+"_dmat.csv";
    vector<vector<double>> distance_mat = parse2DCsvFile2Double(dfilename);
    int n_sets = (int) distance_mat[0].size();
    vector<double> L_vec;
    cout << "scale by 100 \n";
    for (int i=0; i < n_sets; ++i)
    {
        for (int j=i+1; j<n_sets; ++j)
        {
            L_vec.push_back(std::ceil(100*distance_mat[i][j])/100);
        }
    }
    std::sort(L_vec.begin(), L_vec.end());
    auto ip = std::unique(L_vec.begin(), L_vec.end());
    L_vec.resize(std::distance(L_vec.begin(), ip));
    // print_vec(L_vec);
    
    auto start = chrono::steady_clock::now();
    // find number of clusters 
    int n_clusters = 0;
    for (auto & v: clusters) {
        for (auto & i: v) {
            if (i > n_clusters)
                n_clusters = i;
        }
    }

    SetCoverExactSolver solver(n_sets);
    bool only_kcenter = true; //TODO: use hitting condition or not is here
    int argv2 = stoi(argv[2]);
    if(argv2 == 0)
    {
        only_kcenter = false; //TODO: use hitting condition or not is here
    }
    vector<vector<int>> indicator_cover_vec;
    for (auto & kidx: samples_size)
    {
        int k_target = kidx*n_sets/100;
        vector<vector<int>> feasible_clusters = clusters;
        // find feasible clusters (because there might be no solution) 
        int bestK = 1000000;
        if (only_kcenter == false)
        {
            while (bestK > k_target && feasible_clusters.size() > 1)
            {
                SetCoverBase feasible_scb(feasible_clusters, n_clusters);
                vector<vector<int>> cover_set = feasible_scb.thresholdGraph(distance_mat, L_vec.back());
                int n_elem = feasible_scb.max_elem;
                solver.solveSetCover(cover_set, n_elem);
                bestK = (int) solver.getObjValue();
                cout << "bestK: " << bestK << ", k_target = " << k_target << endl;
                if (bestK > k_target){
                    int L_reduce = (int) feasible_clusters.size()/2;
                    feasible_clusters = vector<vector<int>>(feasible_clusters.begin(), feasible_clusters.begin()+L_reduce);
                }
            }
        }

        SetCoverBase scb(feasible_clusters, n_clusters);

        int min_idx = 0;
        int max_idx = (int) L_vec.size();
        int k = -1;
        int space = 1;
        cout << "TODO: stop when k==k_target conditon is meet \n";
        while (max_idx > min_idx + space ) //  && k != k_target
        {
            int mid_idx = ceil((max_idx+min_idx)/2);
            double L = L_vec[mid_idx];
            vector<vector<int>> cover_set = scb.thresholdGraph(distance_mat, L, only_kcenter);
            int n_elem = scb.max_elem;
            // print 
            // for (auto &v: cover_set)
            //     print_vec(v);
            // cout << endl;
            
            // solve set cover 
            solver.solveSetCover(cover_set, n_elem);
            k = (int) solver.getObjValue();
            if (k > k_target){
                min_idx = mid_idx;
            }
            else {
                max_idx = mid_idx;
            }
            cout << "mid_idx = " << mid_idx << ", L = " << L << ", obj = " << k << ", min: " << min_idx << ", max: " << max_idx<< endl;
            if (max_idx <= min_idx + space)
            {
                double L = L_vec[max_idx];
                cover_set = scb.thresholdGraph(distance_mat, L);
                solver.solveSetCover(cover_set, n_elem);
            }
        }
        k = (int) solver.getObjValue();
        assert(k==k_target);
        cout << "n_clusters = " << n_clusters << ", k_target = " << k_target << endl;
        solver.printSolution();
        indicator_cover_vec.push_back(solver.getIndicatorSolution());
    }
    // TODO: change file name accordingly  
    if (only_kcenter == false){
        writeVec2File("/data/hoan/spectral_clustering_matlab/main/cpp/output/" +dataname + "_indicator_sets.csv", indicator_cover_vec);
    }
    else { 
        writeVec2File("/data/hoan/spectral_clustering_matlab/main/cpp/output/" +dataname + "_indicator_sets_kcenter.csv", indicator_cover_vec);
    }
    // writeVec2File("/data/hoan/spectral_clustering_matlab/main/cpp/output/" + dataname + "_indicator_sets_hitk.csv", indicator_cover_vec);

    // for (auto &v : clusters)
    //     print_vec(v);
    // cout << "distance_mat" << endl;
    // for (auto &v: distance_mat)
    //     print_vec(v);
    //
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    return 0;
}
