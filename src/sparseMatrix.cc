//sparseMatrix.cc
#include "sparseMatrix.hh"

#include <cassert>
#include <boost/log/trivial.hpp>


sparseMatrix::sparseMatrix(vector<vector<double> > matrix){
    _A=matrix;
}

sparseMatrix::sparseMatrix(){
    vector<double> line;
    for (int i=0;i<10;++i){
        line.assign(10,0);
        line.at(i)=1;
        _A.push_back(line);
    }
}



void sparseMatrix::matrixVector(const vector<double> &vecin, vector<double> &vecout){
    //check if dimensions are correct
    const int vecinDim = vecin.size();
    const int matrixlineDim = _A.size();
    assert(vecinDim==matrixlineDim);
    vecout.assign(vecinDim,0);
//     BOOST_LOG_TRIVIAL(trace) << "matrixVector: vecin " << printVector(vecin);
    for (int i=0;i<vecinDim;++i){
        for (int j=0;j<vecinDim;++j){
            vecout.at(j)+=_A.at(i).at(j)*vecin.at(j);
        }
        
    }
//     BOOST_LOG_TRIVIAL(trace) << "matrixVector: vecout " << printVector(vecout);
}