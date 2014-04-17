#include "kongrad.hh"
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <boost/log/trivial.hpp>


KonGrad::KonGrad(vector< vector<double> > matrix, vector<double> vec) : _A(matrix), _b(vec) {
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
    KonGrad::matrixVector(_A,vecin, vecout);
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


double KonGrad::getRandomUni(double seed){
    /// @todo zufallsgenerator so umschreiben, dass der generator nur einmal gestartet wird und danach nur weitere Zahlen geholt werden.
    // generators for non-zero elements of matrix. for now only diagonal matrices are created and so the position is fixed.
    double trueSeed=time(NULL)+seed;
    default_random_engine generator(trueSeed);
    uniform_real_distribution<double> distributiond(0.0,1.0); // to generate the values
    // uniform_int_distribution<int> distributioni(1,fMatDim); // to generate the position of the nonzero elements
    
    return distributiond(generator);
}

void KonGrad::createRandomSparseSymmetricMatrix(const int dim, const int seed, vector< vector<double> > &matrixout){
    for (int i=0;i<dim;++i){
        vector<double> line(dim,0);
        int numNonZero=3*getRandomUni(seed);
        for (int j=0;j<numNonZero;++j){
            int posInLine;
            do{
                posInLine=dim*getRandomUni(seed);
            }while(line.at(posInLine)!=0);
            line.at(posInLine)=getRandomUni(seed);
        }
        matrixout.push_back(line);
    }
}


void KonGrad::solve (const vector<double> &startvec){
    const double tol=pow(10,-8);
    const double bnorm=sqrt(skalarProd(_b,_b));
    vector<double> r;
    vector<double> rnew;
    vector<double> p;
    vector<double> pnew;
    vector<double> s;
    vector<double> x=startvec;
    vector<double> xnew;
    
    //temporäre Vektoren
    vector<double> tmpvec;
    
    matrixVector(_A,startvec,tmpvec);
    diffVector(_b,tmpvec, r);
    if(sqrt(skalarProd(r,r))/bnorm < tol){
        cout << "done" << endl;
        return; /// @todo write better exit at start
    }
    p=r;
    
    bool converged=false;
    int iternum=0;
    while(!converged){
        double alpha;
        double beta;
        ++iternum;
        matrixVector(_A,p,s);
        alpha=skalarProd(p,r)/skalarProd(p,s);
        skalarVector(alpha,p,tmpvec);
        sumVector(x,tmpvec,xnew);
        skalarVector(alpha,s,tmpvec);
        diffVector(r,tmpvec,rnew);
        
        if(sqrt(skalarProd(rnew,rnew))/bnorm < tol){
            BOOST_LOG_TRIVIAL(info) << "The algorithm converged. Iterations: " << iternum;
            converged=true;
        }
        
        beta=sqrt(skalarProd(rnew,rnew))/sqrt(skalarProd(r,r));
        skalarVector(beta,p,tmpvec);
        sumVector(rnew,tmpvec,pnew);
        p=pnew;
        r=rnew;
        x=xnew;
    }
    
    printVector(xnew);
    
}