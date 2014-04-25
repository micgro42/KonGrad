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
    BOOST_LOG_TRIVIAL(trace) << "matrixVector: vecin " << printVector(vecin);
    for (int i=0;i<vecinDim;++i){
        for (int j=0;j<vecinDim;++j){
            vecout.at(j)+=_A.at(i).at(j)*vecin.at(j);
        }
        
    }
    BOOST_LOG_TRIVIAL(trace) << "matrixVector: vecout " << printVector(vecout);
}


int sparseMatrix::printVector (const vector<double> &vec){
    // this prints out the vector if the size is smaller than 20.
    bool isSmallerThan20 = ( vec.size() < 20 );
    if ( isSmallerThan20 ){
        for( vector<double>::const_iterator i = vec.begin(); i != vec.end(); ++i){
            cout << *i << ' ';
        }
    }
    cout << endl;
    return !isSmallerThan20;
}

int sparseMatrix::printMatrix (){
    // this prints out the matrix if the size is smaller than 20.
    bool isSmallerThan20 = ( _A.size() < 20 );
    if ( isSmallerThan20 ){
        for ( unsigned int i = 0; i<_A.size(); ++i){
            for( vector<double>::const_iterator j = _A.at(i).begin(); j != _A.at(i).end(); ++j){
                cout << *j << ' ';
            }
            cout << endl;
        }
    }
    return !isSmallerThan20;
}











