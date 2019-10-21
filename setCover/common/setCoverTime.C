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
    
    auto start2 = chrono::steady_clock::now();
    vector<vector<float>> Dist(n_sets, vector<float>(n_sets,0.0));
    // #pragma omp parallel for collapse(2)
    int i,j;
    #pragma omp parallel for private(i,j)
    for (i=0; i < n_sets; ++i){
        for (j=0; j<n_sets; ++j){
            Dist[i][j] = correlationDistance(distance_mat.at(i), distance_mat.at(j), n_dim);
        }
    }
    auto end2 = chrono::steady_clock::now();
    cout << "compute distance matrix in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;

    // create the instance 
    start2 = chrono::steady_clock::now();
    vector<int> instance;
    float L = 0.5;
    create_threshold_data(Dist, L, n_sets, instance);
    end2 = chrono::steady_clock::now();
    cout << "create instance in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;

    // solve set cover 
    start2 = chrono::steady_clock::now();
    int rounds = 1;
    graph<intT> G = readGraphFromVector<intT>(instance);
    string outputFile = dataname+"_setcover_sol__.csv";
    char * oFile = &outputFile[0];
    timesetCover(G, rounds, oFile);
    end2 = chrono::steady_clock::now();
    cout << "solve setcover in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;


  // graph<intT> G = readGraphFromFile<intT>(iFile);
  // timesetCover(G, rounds, oFile);
}
