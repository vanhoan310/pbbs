#include <iostream>
#include <vector>
#include <random>
#include <set>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <omp.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::MatrixXf;
using Eigen::MatrixXi;
int main(int argc, char* argv[])
{
    // cout << m.colwise().sum() << endl;
    // MatrixXf::Index   maxIndex;
    // float maxNorm = m.colwise().sum().maxCoeff(&maxIndex);

    int n_elem = atoi(argv[1]);
    // Generate random set cover matrix
    MatrixXi M = MatrixXi::Zero(n_elem, n_elem);
    std::mt19937 e(std::random_device{}());
    std::bernoulli_distribution d;
    for (int i=0; i < n_elem; i=i+1){
        for (int j=0; j < n_elem; j=j+n_elem/8){
            M(i,j) = d(e);

        }
    }
    for (int i=0; i < n_elem; ++i){
        M(i,i) = 1;
    }
    MatrixXi N = M;
    
    // auto start = std::chrono::steady_clock::now();
    // set<int> cover_items; //store cover elements 
    // vector<int> solution;
    //
    // set<int> CC; // store uncovered elements
    // for (int i=0; i < n_elem; ++i)
    //     CC.insert(i);
    // while(CC.size() > 0) {
    //     MatrixXi::Index maxIndex;
    //     auto maxNorm = N.colwise().sum().maxCoeff(&maxIndex);
    //     solution.push_back(maxIndex);
    //     set<int> select_col;
    //     vector<int> select_col_vec;
    //     for (auto idx: CC){
    //         if (N(idx, maxIndex) == 1){
    //             select_col.insert(idx);
    //             select_col_vec.push_back(idx);
    //         }
    //     }
    //     std::set<int>::iterator it;
    //     // #pragma omp parallel for private(it, N)
    //     for (it = select_col.begin(); it != select_col.end(); ++it){
    //         N.row(*it).setZero();
    //     }
    //     set<int> result;
    //     std::set_difference(CC.begin(), CC.end(), select_col.begin(), select_col.end(),
    //                 std::inserter(result, result.end()));
    //     CC = result;
    // }
    // auto finish = std::chrono::steady_clock::now();
    // std::cout << "\nTime difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()/1000.0 << "[s]" << std::endl;
    // for(auto v: solution)
    //     cout << v +1 << " ";
    // cout << endl;
   
    cout << "Start my algorithm: n = " << n_elem << "\n";
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    set<int> C; // store uncovered elements
    for (int i=0; i < n_elem; ++i)
        C.insert(i);

    vector<int> cover_set; // store cover set 
    int count = 0;
    while (C.size() > 0)
    {
        VectorXi::Index maxS;
        if(count == 0){
            auto maxNorm = N.colwise().sum().maxCoeff(&maxS);
        } else {
            Eigen::VectorXi s = Eigen::VectorXi::Zero(n_elem);
            std::set<int>::iterator it;
            for (it = C.begin(); it != C.end(); ++it) {
                s += M.row(*it);
            }
            int maxValue = s.maxCoeff(&maxS);
        }
        count++;
        cover_set.push_back(maxS);
        set<int> select_col;
        for (auto idx: C){
            if (M(idx, maxS) == 1){
                select_col.insert(idx);
            }
        }
        set<int> result;
        std::set_difference(C.begin(), C.end(), select_col.begin(), select_col.end(),
                    std::inserter(result, result.end()));
        C = result;
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "\nTime difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()/1000.0 << "[s]" << std::endl;
    cout << "\ncover set:";
    for(auto v: cover_set)
        cout << v + 1<< " ";
    cout << endl;
    cout << "cover size:" << cover_set.size() << endl;
    cout << "pct size:" << (double) 100.0*cover_set.size()/n_elem << endl;


}
