#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <iomanip>      // std::setprecision
#include <omp.h>
#include "util_functions.h"
#include <chrono>
using namespace std;

template <class T>
vector<T> HDscore(const vector<vector<T>> & distance_mat, const vector<int> & setcover, int n_sets, int n_dim, int k)
{
    vector<T> scores(n_sets);
    int i,j;
    int cover_size = (int) setcover.size();
    #pragma omp parallel for private(j)
    for (i=0; i < n_sets; ++i){
        T dmin = correlationDistance(distance_mat.at(i), distance_mat.at(setcover[0]), n_dim);
        for (j=1; j< cover_size; ++j){
            T dd = correlationDistance(distance_mat.at(i), distance_mat.at(setcover[j]), n_dim);
            if (dd < dmin){
                dmin = dd;
            }
        }
        scores[i] = dmin;
    }

    sort(scores.begin(), scores.end());
    vector<T> newscore(scores.end()-k, scores.end());
    return newscore;

}

int main(int argc, char* argv[])
{
    std::string dataname = argv[1]; // zeisel grun
    std::string dfilename = "../../../spectral_clustering_matlab/data/"+dataname+"_pca.csv";
    auto start = chrono::steady_clock::now();
    vector<vector<double>> distance_mat = parse2DCsvFile2Double(dfilename);
    auto end = chrono::steady_clock::now();
    cout << "read file in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    int n_sets = (int) distance_mat.size();
    int n_dim =  (int) distance_mat[0].size();
    cout << "# of items: " << n_sets << ", n_dim: " << n_dim << endl;

    string solfilename = dataname+"_setcover_solutions.csv";
    vector<vector<int>> solutions = parse2DCsvFile(solfilename);
    int q = max(1, (int) (1e-4* (float) n_sets));
    for (auto sol: solutions){
        vector<double> hdscore = HDscore(distance_mat, sol, n_sets, n_dim, q);
        for (size_t i=0; i < hdscore.size(); ++i){
            cout << std::fixed << std::setprecision(2) << hdscore.at(i) << "  ";
        }
        cout << "\n";
    }

    return 0;
}

