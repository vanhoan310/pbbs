#include<iostream>
#include<vector>
#include <Eigen/Dense>
#include "util_functions.h"
#include <chrono>

using namespace std;

float correlationDistance(vector<double> & X, vector<double> &Y, int n) 
{ 

    float sum_X = 0.0, sum_Y = 0.0, sum_XY = 0.0; 
    float squareSum_X = 0.0, squareSum_Y = 0.0; 

    for (int i = 0; i < n; i++) 
    { 
        // sum of elements of array X. 
        sum_X = sum_X + X[i]; 

        // sum of elements of array Y. 
        sum_Y = sum_Y + Y[i]; 

        // sum of X[i] * Y[i]. 
        sum_XY = sum_XY + X[i] * Y[i]; 

        // sum of square of array elements. 
        squareSum_X = squareSum_X + X[i] * X[i]; 
        squareSum_Y = squareSum_Y + Y[i] * Y[i]; 
    } 

    // use formula for calculating correlation coefficient. 
    float corr = (n * sum_XY - sum_X * sum_Y) 
                / sqrt((n * squareSum_X - sum_X * sum_X) 
                    * (n * squareSum_Y - sum_Y * sum_Y)); 

    return 1.0 - corr; 
} 

int main(int argc, char *argv[])
{
    std::string dataname = argv[1]; // zeisel grun
    float scale_factor = 100.0;
    cout << "Distance is scale by " << scale_factor << endl;
    std::string dfilename = "../../data/"+dataname+"_pca.csv";
    auto start = chrono::steady_clock::now();
    vector<vector<double>> distance_mat = parse2DCsvFile2Double(dfilename);
    auto end = chrono::steady_clock::now();
    cout << "Elapsed time in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    int n_sets = (int) distance_mat.size();
    int n_dim =  (int) distance_mat[0].size();
    
    auto start2 = chrono::steady_clock::now();
    Eigen::MatrixXi Dist = Eigen::MatrixXi::Zero(n_sets,n_sets);
    #pragma omp parallel for collapse(2)
    for (int i=0; i < n_sets; ++i){
        for (int j=1; j<n_sets; ++j){
            Dist(i,j) = (int) (scale_factor * correlationDistance(distance_mat.at(i), distance_mat.at(j), n_dim));
        }
    }

    // cout << Dist;

    int max_value = Dist.maxCoeff();

    auto end2 = chrono::steady_clock::now();
    cout << "Elapsed time in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;

    vector<int> test;
    #pragma omp parallel for 
    for (int i=0; i < 10; ++i){
    #pragma omp critical
        test.push_back(i);
    }
    cout << test.size() << endl;
    // print_vec(test);

    return 0;
}

