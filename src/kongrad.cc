#include "kongrad.hh"
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>


KonGrad::KonGrad(vector< vector<double> > Matrix, vector<double> vec){
    _A=Matrix;
    _b=vec;
}

KonGrad::KonGrad(){
    vector<double> line;
    for (int i=0;i<3;++i){
        line.assign(3,0);
        line.at(i)=1;
        _A.push_back(line);
    }
    
    _b.assign(3,0);
}




void KonGrad::testmv(const vector<double> vecin){
    vector<double> vecout;
    KonGrad::matrixVector(this->_A,vecin, vecout);
}



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

void KonGrad::diffVector(const vector<double> &vecin1, const vector<double> &vecin2, vector<double> &vecout){
    const int vecin1Dim = vecin1.size();
    const int vecin2Dim = vecin2.size();
    assert(vecin1Dim == vecin2Dim);
    
    vecout.clear();
    
    for (int i=0;i<vecin1Dim;++i){
        vecout.push_back(vecin1.at(i)-vecin2.at(i));
    }
}


void KonGrad::sumVector(const vector<double> &vecin1, const vector<double> &vecin2, vector<double> &vecout){
    const int vecin1Dim = vecin1.size();
    const int vecin2Dim = vecin2.size();
    assert(vecin1Dim == vecin2Dim);
    
    vecout.clear();
    
    for (int i=0;i<vecin1Dim;++i){
        vecout.push_back(vecin1.at(i)+vecin2.at(i));
    }
}


void KonGrad::skalarVector(const double alpha, const vector<double> &vecin, vector<double> &vecout){
    const int vecinDim = vecin.size();
    vecout.clear();
    
    for (int i=0;i<vecinDim;++i){
        vecout.push_back(vecin.at(i)*alpha);
    }
}

void KonGrad::printVector (const vector<double> &vec){
    cout << "vector ";
    for( vector<double>::const_iterator i = vec.begin(); i != vec.end(); ++i){
        cout << *i << ' ';
    }
    cout << endl;
}

void KonGrad::solve (){
    ///@todo start mit Nullvektor einbauen
}

double KonGrad::skalarProd(const vector<double> &vecin1, const vector<double> &vecin2){
    const int vecin1Dim = vecin1.size();
    const int vecin2Dim = vecin2.size();
    assert(vecin1Dim == vecin2Dim);
    double skalarProd=0;
    for (int i=0;i<vecin1Dim;++i){
        skalarProd+=vecin1.at(i)*vecin2.at(i);
    }
    return skalarProd;
}


void KonGrad::solve (const vector<double> &startvec){
    const double tol=pow(10,-8);
    const double bnorm=sqrt(skalarProd(_b,_b));
    double alpha;
    double beta;
    vector<double> r;
    vector<double> rnew;
    vector<double> p;
    vector<double> pnew;
    vector<double> s;
    vector<double> x=startvec;
    vector<double> xnew;
    
    //temporäre Vektoren
    vector<double> tmpvec;
    vector<double> tmpvarMatvec;
    vector<double> tmpvarvecdiff;
    
    matrixVector(_A,startvec,tmpvarMatvec);
    diffVector(_b,tmpvarMatvec, r);
    if(sqrt(skalarProd(r,r))/bnorm < tol){
        cout << "done" << endl;
        return; /// @todo write better exit at start
    }
    p=r;
    
    bool notYetSmallEnough=true;
    int iternum=0;
    while(notYetSmallEnough){
        ++iternum;
        matrixVector(_A,p,s);
        alpha=skalarProd(p,r)/skalarProd(p,s);
        skalarVector(alpha,p,tmpvec);
        sumVector(x,tmpvec,xnew);
        skalarVector(alpha,s,tmpvec);
        diffVector(r,tmpvec,rnew);
        
        if(sqrt(skalarProd(rnew,rnew))/bnorm < tol){
            cout << "done, iterations: " << iternum << endl;
            notYetSmallEnough=false; 
        }
        
        beta=sqrt(skalarProd(rnew,rnew))/sqrt(skalarProd(r,r));
        skalarVector(beta,p,tmpvec);
        sumVector(rnew,tmpvec,pnew);
        p=pnew;
        r=rnew;
    }
    
}