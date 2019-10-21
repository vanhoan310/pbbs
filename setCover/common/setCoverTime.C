// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "parallel.h"
#include "IO.h"
#include "graphIO.h"
#include "parseCommandLine.h"
#include "setCover.h"
#include <omp.h>
#include "util_functions.h"
#include "grid_cover.h"
#include <chrono>
using namespace std;
using namespace benchIO;

template <class T>
void create_threshold_data(vector<vector<T>> & Dist, T LL, int set_size, vector<int> & instance){
    int n;
    n=(int) set_size;
    int total_elem = 0;
    int first_elem_set_idx = 0;
    instance.clear(); // clear the vector 
    instance.push_back(n);
    instance.push_back(total_elem); // will change later 
    for (int i=0; i < n; ++i){
        instance.push_back(i);
    }

    for(int i = 0; i<n; i++){
        instance[i+2] = first_elem_set_idx;
        for(int j = 0; j < n; j++){
            if(Dist[i][j] <= LL)
            {
                instance.push_back(j);
                first_elem_set_idx++;
            }
        }
    }
    instance[1] = first_elem_set_idx;
}

void timesetCover(graph<intT> G, int rounds, char* outFile) {
  _seq<intT> Sets = _seq<intT>(NULL,0);
  for (int i=0; i < rounds; i++) {
    if (Sets.A != NULL) free(Sets.A);
    graph<intT> GG = G.copy();
    startTime();
    Sets = setCover(GG);
    nextTimeN();
    GG.del();
  }
  cout << endl;

  if (outFile != NULL)
    writeIntArrayToFile(Sets.A, Sets.n, outFile);

  Sets.del();
  G.del();
}

template <class T>
vector<int> run_setcover_binarysearch(vector<vector<T>> & Dist, T d_min, T d_max, int n_points, int n_iters, string dataname, int n_sets)
{
    // create the instance 
    for (int iter=0; iter < n_iters; ++iter)
    {
        T L = (d_min + d_max)/2.0;
        vector<int> instance;
        create_threshold_data(Dist, L, n_sets, instance);
        // run setcover
        int rounds = 1;
        graph<intT> G = readGraphFromVector<intT>(instance);
        string outputFile = dataname+"_setcover_sol__.csv";
        char * oFile = &outputFile[0];
        timesetCover(G, rounds, oFile);
        vector<int> set_cover =  parse1DcsvFile2Int(outputFile, 1);
        std::remove(oFile);
        int cover_size = (int) set_cover.size();
        cout << "iter: " << iter << ", n_points: " << n_points << ", set_cover size: " << cover_size << "\n";
        if (cover_size < n_points){
            d_max = L;
        }
        else{
            d_min = L;
        }
        if(iter == n_iters - 1){ return set_cover; }
    }
}

vector<int> fill_in(vector<int> & solution, int n_sets, int n_points){
    int sol_size = (int) solution.size();
    if (sol_size != n_points)
        cout << "fill in solution \n";
    if (n_points > sol_size){
       for (int i=0; i < n_points-sol_size; ++i){
          solution.push_back(rand() % (n_sets-1));
       }
      return solution; 
    }
    else
    {
        vector<int> newsolution(solution.begin(), solution.begin()+n_points);
        return newsolution;
    }
}

int parallel_main(int argc, char* argv[]) {
  // commandLine P(argc, argv, "[-o <outFile>] [-r <rounds>] <inFile>");
  // char* iFile = P.getArgument(0);
  // char* oFile = P.getOptionValue("-o");
  // int rounds = P.getOptionIntValue("-r",1);

    std::string dataname = argv[1]; // zeisel grun
    std::string dfilename = "../../../spectral_clustering_matlab/data/"+dataname+"_pca.csv";
    // TODO: need to remove the labels 
    // std::string dfilename = "../../../data/"+dataname+"-prepare-log_count_pca.csv";
    // if (!is_file_exist(dfilename)){
    //     dfilename = "../../../data/"+dataname+"_pca.csv";
    // }
    // if (!is_file_exist(dfilename)){
    //     dfilename = "../../../data/"+dataname+"-prepare-log_count_pca2000.csv";
    // }
    auto start = chrono::steady_clock::now();
    vector<vector<double>> distance_mat = parse2DCsvFile2Double(dfilename);
    auto end = chrono::steady_clock::now();
    cout << "read file in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    int n_sets = (int) distance_mat.size();
    int n_dim =  (int) distance_mat[0].size();
    cout << "# of items: " << n_sets << ", n_dim: " << n_dim << endl;

    int n_boxes = 80000;
    vector<vector<int>> sc_solutions;
    if ( n_sets < n_boxes)
    {
        start = chrono::steady_clock::now();
        vector<vector<float>> Dist(n_sets, vector<float>(n_sets,0.0));
        int i,j; 
        float d_max = 0.0;
        #pragma omp parallel for private(i,j)
        for (i=0; i < n_sets; ++i){
            for (j=0; j<n_sets; ++j){
                Dist[i][j] = correlationDistance(distance_mat.at(i), distance_mat.at(j), n_dim);
                if (Dist[i][j] > d_max){
                    d_max = Dist[i][j];
                }
            }
        }
        end = chrono::steady_clock::now();
        cout << "compute distance matrix in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;

        start = chrono::steady_clock::now();
        d_max = d_max/2.0;
        float d_min = d_max/100.0;
        int n_points = min(1000, n_sets);
        int n_iters = 14;
        vector<int> setcover = run_setcover_binarysearch(Dist, d_min, d_max, n_points, n_iters, dataname, n_sets);
        vector<int> fill_setcover = fill_in(setcover, n_sets, n_points);
        sc_solutions.push_back(fill_setcover);
        end = chrono::steady_clock::now();
        cout << "setcover in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    }
    else 
    {
        // run grid box algorithm 
        auto start2 = chrono::steady_clock::now();
        int dim_x = 10; double N = 17.0; int L_box = n_boxes;
        int L_min = min(n_sets/3, L_box-2000); int L_max = min(n_sets/2, L_box+2000);
        vector<int> box_idx = run_box_binarysearch(distance_mat, dim_x, N, L_min, L_max);
        auto end2 = chrono::steady_clock::now();
        cout << "run box algorithm in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;

        // compute distance matrix
        start = chrono::steady_clock::now();
        n_boxes = (int) box_idx.size();
        vector<vector<float>> Dist(n_boxes, vector<float>(n_boxes,0.0));
        int i,j; 
        float d_max = 0.0;
        #pragma omp parallel for private(i,j)
        for (i=0; i < n_boxes; ++i){
            for (j=0; j< n_boxes; ++j){
                Dist[i][j] = correlationDistance(distance_mat.at(box_idx[i]), distance_mat.at(box_idx[j]), n_dim);
                if (Dist[i][j] > d_max){
                    d_max = Dist[i][j];
                }
            }
        }
        end = chrono::steady_clock::now();
        cout << "compute distance matrix in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;

        // compute set cover 
        start = chrono::steady_clock::now();
        d_max = d_max/2.0;
        float d_min = d_max/100.0;
        int n_points = 10000;
        int n_iters = 14;
        vector<int> setcover = run_setcover_binarysearch(Dist, d_min, d_max, n_points, n_iters, dataname, n_boxes);
        vector<int> convert_boxid2realid;
        for(size_t k=0; k < setcover.size(); ++k){
            convert_boxid2realid.push_back(box_idx[setcover[k]]);
        }
        vector<int> fill_setcover = fill_in(convert_boxid2realid, n_sets, n_points);
        sc_solutions.push_back(fill_setcover);
        end = chrono::steady_clock::now();
        cout << "setcover in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    }
    //write solution to file
    string fileName = dataname + "_setcover_solutions.csv";
    writeVec2File(fileName, sc_solutions);

}
