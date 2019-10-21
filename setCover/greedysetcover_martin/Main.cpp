#include <iostream>
#include "LinearMultiArray.h"
#include "SetCover.h"
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <istream>
#include <cstring>
#include <stdlib.h> 
#include <Eigen/Dense>
#include <chrono>
#include <omp.h>
using namespace Eigen;
using namespace std;

struct set_data {
    unsigned int ** sets;
    unsigned short ** weights;
    unsigned int * set_sizes;
    unsigned int * element_size_lookup;
    unsigned int set_count;
    unsigned int uniqu_element_count;
    unsigned int all_element_count;
    unsigned int max_weight;
};

set_data create_threshold_data(Eigen::MatrixXi & Dist, int LL, int set_size, int elment_count ){
    int n,m;
    m=set_size;
    n=elment_count;
    set_data ret_struct;
    ret_struct.set_count = m;
    ret_struct.uniqu_element_count = n;
    unsigned int * element_buffer=(unsigned int *)malloc(n * sizeof(unsigned int));
    ret_struct.element_size_lookup = (unsigned int *)malloc((n+2) * sizeof(unsigned int));
    unsigned int * element_size=ret_struct.element_size_lookup;
    memset(&element_size[0], 0, sizeof(unsigned int)*n+2);
    
    ret_struct.set_sizes = (unsigned int *)malloc((m+1) * sizeof(unsigned int));
    ret_struct.set_sizes[m]=0;
    unsigned int * set_sizes=ret_struct.set_sizes;
    
    ret_struct.sets = (unsigned int **)malloc(m * sizeof(unsigned int*));
    unsigned int ** sets   = ret_struct.sets;

    ret_struct.weights = (unsigned short **)malloc(m * sizeof(unsigned short*));
    unsigned short ** weights = ret_struct.weights;
    
    ret_struct.max_weight = 0;
    ret_struct.all_element_count=0;
    for(int i = 0; i<m; i++){
        unsigned int element_counter=0;
        for(int j = 0; j < m-1; j++){ //TODO: why fail here if I use j = m-1
            int curr_element=j + 1;
            if(Dist(i,j) <= LL)
            {
                element_buffer[element_counter++]=curr_element;
                element_size[curr_element]++;
                ret_struct.all_element_count++;   
            }
        }
        if(ret_struct.max_weight < element_counter){
            ret_struct.max_weight = element_counter;
        }
        
        unsigned int * elements = new unsigned int[element_counter];
        memcpy(elements,element_buffer,sizeof(int)*element_counter);
        unsigned short * weight = new unsigned short[element_counter];
        std::fill_n(weight, element_counter, 1);
        
        weights[i] = weight;
        sets[i] = elements;
        set_sizes[i]=element_counter;
        
    }
    free(element_buffer);
    return ret_struct;
}


set_data create_random_set_data(int set_size, int elment_count ){
    int n,m;
    m=set_size;
    n=elment_count;
    set_data ret_struct;
    ret_struct.set_count = m;
    ret_struct.uniqu_element_count = n;
    
    
    unsigned int * element_buffer=(unsigned int *)malloc(n * sizeof(unsigned int));
    ret_struct.element_size_lookup = (unsigned int *)malloc((n+2) * sizeof(unsigned int));
    unsigned int * element_size=ret_struct.element_size_lookup;
    memset(&element_size[0], 0, sizeof(unsigned int)*n+2);
    
    ret_struct.set_sizes = (unsigned int *)malloc((m+1) * sizeof(unsigned int));
    ret_struct.set_sizes[m]=0;
    unsigned int * set_sizes=ret_struct.set_sizes;
    
    ret_struct.sets = (unsigned int **)malloc(m * sizeof(unsigned int*));
    unsigned int ** sets   = ret_struct.sets;

    
    ret_struct.weights = (unsigned short **)malloc(m * sizeof(unsigned short*));
    unsigned short ** weights = ret_struct.weights;
    
    ret_struct.max_weight = 0;
    ret_struct.all_element_count=0;
    for(int i = 0; i<m; i++){
        unsigned int element_counter=0;
        for(int j = 0; j < 10; j++)        {
            int curr_element=rand()%(elment_count-2 + 1) + 1;
            if(curr_element == 0)
                continue;
            element_buffer[element_counter++]=curr_element;
            element_size[curr_element]++;
            ret_struct.all_element_count++;
        }
        if(ret_struct.max_weight < element_counter){
            ret_struct.max_weight = element_counter;
        }
        
        unsigned int * elements = new unsigned int[element_counter];
        memcpy(elements,element_buffer,sizeof(int)*element_counter);
        unsigned short * weight = new unsigned short[element_counter];
        std::fill_n(weight, element_counter, 1);
        
        weights[i] = weight;
        sets[i] = elements;
        set_sizes[i]=element_counter;
        
    }
    free(element_buffer);
    
    return ret_struct;
}

set_data read_in_set_data(char * dbname){
    // Open the stream
    std::ifstream file(dbname);
    if (!file )
    {
        std::cout << "File dosen't exist" << std::endl;
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();
    set_data ret_struct;
    std::string line;
    std::getline(buffer, line);
    std::istringstream line_stream(line);
    int n,m;
    line_stream >> n >> m;
    ret_struct.set_count = m;
    ret_struct.uniqu_element_count = n;   
    unsigned int * element_buffer=(unsigned int *)malloc(n * sizeof(unsigned int));
    unsigned int * element_size=new unsigned int[n+2];
    memset(&element_size[0], 0, sizeof(unsigned int)*(n+2));
    ret_struct.element_size_lookup = element_size;
    unsigned int * set_size=new unsigned int[m];
    
    ret_struct.set_sizes = set_size;
    unsigned int ** sets   = new unsigned int*[m];
    ret_struct.sets = sets;
    unsigned short ** weights = new unsigned short*[m];
    ret_struct.weights = weights;
    ret_struct.max_weight = 0;
    ret_struct.all_element_count=0;

    for(int i = 0; i<m; i++){
        std::getline(buffer, line);
        std::string token;
        unsigned int element_counter=0;
        line_stream.clear();
        line_stream.str(line);
        while (std::getline(line_stream, token, ' '))
        {
            int curr_element=atoi(token.c_str());
            if(curr_element == 0)
                continue;
            element_buffer[element_counter++]=curr_element;
            element_size[curr_element]++;
            ret_struct.all_element_count++;
        }
        if(ret_struct.max_weight < element_counter){
            ret_struct.max_weight = element_counter;
        }

        unsigned int * elements = new unsigned int[element_counter];
        memcpy(elements,element_buffer,sizeof(int)*element_counter);
        unsigned short * weight = new unsigned short[element_counter];
        std::fill_n(weight, element_counter, 1);

        weights[i] = weight;
        sets[i] = elements;
        set_size[i]=element_counter;

    }
    return ret_struct;
}

bool validate_result(std::list<set *> * ret,unsigned int uniqu_element_count){
    std::list<set *>::const_iterator iterator;
    unsigned int result_element_count=0;
    for (iterator = ret->begin(); iterator != ret->end(); ++iterator) {
        set::element * element =(*iterator)->elements;
        do{
            result_element_count++;
        }while((element=element->next)!=NULL);
    }
    return (uniqu_element_count==result_element_count)? true : false;
}

vector<vector<double>> parse2DCsvFile2Double(string inputFileName) {
    vector<vector<double> > data;
    ifstream inputFile(inputFileName);
    int l = 0; // 0
    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<double> record;
            int removefist = 0;
 
            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    if (removefist >= 0)    
                    record.push_back(std::stod(line));
                    else
                        removefist++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }
            data.push_back(record);
        }
    }
 
    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }
 
    return data;
}

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

float euclideanDistance(vector<double> & X, vector<double> &Y, int n) 
{ 

    float sum_XY = 0.0;
    for (int i = 0; i < n; i++) 
    { 
        sum_XY += (X[i]-Y[i])*(X[i]-Y[i]);
    } 
    return sqrt(sum_XY);
} 

bool is_file_exist(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char* argv[])
{

    std::string dataname = argv[1]; // zeisel grun
    float scale_factor = 100.0;
    cout << "ELEMENT INDEX START FROM 0 -----------------------------\n";
    cout << "Distance is scale by " << scale_factor << endl;
    std::string dfilename = "../../../data/"+dataname+"-prepare-log_count_pca.csv";
    if (!is_file_exist(dfilename)){
        dfilename = "../../../data/"+dataname+"_pca.csv";
    }
    if (!is_file_exist(dfilename)){
        dfilename = "../../../data/"+dataname+"-prepare-log_count_pca2000.csv";
    }
    auto start = chrono::steady_clock::now();
    vector<vector<double>> distance_mat = parse2DCsvFile2Double(dfilename);
    auto end = chrono::steady_clock::now();
    cout << "read file in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;
    int n_sets = (int) distance_mat.size();
    int n_dim =  (int) distance_mat[0].size();
    cout << "# of items: " << n_sets << ", n_dim: " << n_dim << endl;
    
    auto start2 = chrono::steady_clock::now();
    Eigen::MatrixXi Dist = Eigen::MatrixXi::Zero(n_sets,n_sets);
    // #pragma omp parallel for collapse(2)
    int i,j;
    #pragma omp parallel for private(i,j)
    for (i=0; i < n_sets; ++i){
        for (j=0; j<n_sets; ++j){
            Dist(i,j) = (int) (scale_factor * correlationDistance(distance_mat.at(i), distance_mat.at(j), n_dim));
        }
    }
    for (int i=n_sets-10; i < n_sets; ++i){
        for (int j=n_sets- 10; j<n_sets; ++j){
            cout << Dist(i,j) << "\t";
        }
        cout << endl;
    }
    auto end2 = chrono::steady_clock::now();
    cout << "compute distance matrix in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;
    
    int L_max = Dist.maxCoeff()/2;
    int L_min = L_max/8;
    int pct = atoi(argv[2]);
    int k_target = round(pct*n_sets/100);
    cout << "\nL_max: "<< L_max << ", pct: " << pct << ", k_target: " << k_target << endl;
    int n_iterations = 6;
    start = chrono::steady_clock::now();
    for (int iter=0; iter < n_iterations; iter++){
        cout << "---------------------------------------------------------------\n";
        start2 = chrono::steady_clock::now();
        std::cout << "start readin";
        std::cout.flush();
        // set_data set_data = read_in_set_data("testdatafile");
        // set_data set_data = create_random_set_data(n_sets, n_sets);
        int LL = (L_min + L_max)/2;
        set_data set_data = create_threshold_data(Dist, LL, n_sets, n_sets);
        std::cout << " --- Done" << std::endl;
        std::cout << "init setcover";
        std::cout.flush();
        set_cover setcover(set_data.set_count, set_data.uniqu_element_count, set_data.max_weight,
                           set_data.all_element_count, set_data.element_size_lookup);
        for(unsigned int i = 0; i < set_data.set_count; i++){
                setcover.add_set(i+1, set_data.set_sizes[i] ,(const unsigned int*)set_data.sets[i],
                          (const unsigned short*)set_data.weights[i], set_data.set_sizes[i]); free(set_data.sets[i]); free(set_data.weights[i]);
        }
        std::cout << " --- Done" << std::endl;
        end2 = chrono::steady_clock::now();
        cout << "create the instance in seconds : " << chrono::duration_cast<chrono::seconds>(end2 - start2).count() << " sec" << endl;

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        std::cout << "execute setcover" ; std::cout.flush();
        std::list<set *> ret= setcover.execute_set_cover();
        std::cout << " --- Done" << std::endl;
        // std::cout << "validate result ";
        // if(validate_result(&ret,set_data.uniqu_element_count))
        //     std::cout << " VALID ";
        // else
        //     std::cout << " NOT VALID ";
        // std::cout << " --- Done" << std::endl;
        
        std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
        std::cout << "L_min: " << L_min << ", L: "<<LL<<", L_max: " << L_max << endl;
        std::cout << "execute setcover in " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin).count()/1000.0 << "[s]" << std::endl;
        int cover_size = (int) ret.size();
        if (cover_size > k_target){
            L_min = LL;
        }else{
            L_max = LL;
        }
        std::cout << "size of the set cover: " << ret.size() << std::endl;
        // std::cout << "pct: " << (double) 100*ret.size()/n_sets << std::endl;
    }
    end = chrono::steady_clock::now();
    cout << "---------------------------------------------------------------\n";
    cout << "---------------------------------------------------------------\n";
    cout << "total runtime in seconds : " << chrono::duration_cast<chrono::seconds>(end - start).count() << " sec" << endl;

//     std::list<set *>::const_iterator iterator;
//     for (iterator = ret.begin(); iterator != ret.end(); ++iterator) {
//         std::cout << "set id " << (*iterator)->set_id << " Element: ";
//         set::element * element =(*iterator)->elements;
//         do{
//             std::cout << element->element_id << ", ";
//         }while((element=element->next)!=NULL);
//         std::cout << std::endl;
//     }

}

