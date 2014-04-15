#include "kongrad.hh"
#include <vector>
#include <iostream>

///Constructor
KonGrad::KonGrad(vector< vector<double> > Matrix, vector<double> vec){
    KonGrad::A=Matrix;
    KonGrad::b=vec;
}


void KonGrad::testmv(const vector<double> vecin){
    vector<double> vecout;
    KonGrad::matrixVector(this->A,vecin, vecout);
}


/**
* @brief Mulitply a matrix and a vector
* 
*/
void KonGrad::matrixVector(vector< vector<double> > Matrix, vector<double> vecin, vector<double> vecout){
    //check if dimensions are correct
    const int vecinDim = vecin.size();
    
    std::cout << "blub" << std::endl;
}
