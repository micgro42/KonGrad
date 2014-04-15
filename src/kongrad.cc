#include "kongrad.hh"
#include <vector>
#include <iostream>
#include <cassert>

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
void KonGrad::matrixVector(const vector< vector<double> > &matrix, const vector<double> &vecin, vector<double> &vecout){
    //check if dimensions are correct
    const int vecinDim = vecin.size();
    const int matrixlineDim = matrix.size();
    assert(vecinDim==matrixlineDim);
    vecout.assign(vecinDim,0);
    for (int i=0;i<vecinDim;++i){
        for (int j=0;j<vecinDim;++j){
            vecout.at(j)+=matrix.at(i).at(j)*vecin.at(j);
        }
        
    }
    this->printVector(vecout);
}


void KonGrad::printVector (const vector<double> &vec){
    cout << "vector ";
    for( vector<double>::const_iterator i = vec.begin(); i != vec.end(); ++i){
        cout << *i << ' ';
    }
    cout << endl;
}