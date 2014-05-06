#include "kongrad.hh"
#include "global.h"
//#include "geom_pbc.c"
#include <vector>
#include <random>
#include <cassert>
#include <cmath>
#include <boost/log/trivial.hpp>

/**
 * @file kongrad.cc
 * 
 * @brief sourcecode for the KonGrad class
 * 
 * 
 */


KonGrad::KonGrad(vector< vector<double> > matrix, vector<double> vec) : _A(matrix), _b(vec) {
    _mass=0.1;
}

KonGrad::KonGrad(){
    vector<double> line;
    for (int i=0;i<3;++i){
        line.assign(3,0);
        line.at(i)=1;
        _A.push_back(line);
    }
    _mass=0.1;
    _b.assign(3,0);
}




void KonGrad::testmv(const vector<double> vecin){
    vector<double> vecout;
    KonGrad::matrixVector(vecin, vecout);
}



void KonGrad::matrixVector(const vector<double> &vecin, vector<double> &vecout){
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


void KonGrad::matrixVectorLaplace(const vector<double> &vecin, vector<double> &vecout){
	const int vecinDim = vecin.size();
	vecout.assign(vecinDim,0);
	for (int i=0; i<vecinDim;++i){
        vecout.at(i)=(2*ndim+_mass*_mass)*vecin.at(i);
        for (int k=1;k<=ndim;++k){
        	vecout.at(i)-=(vecin.at(nn[k][i])+vecin.at(nn[k+ndim][i]));
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

int KonGrad::printVector (const vector<double> &vec){
    // this prints out the vector if the size is smaller than 20.
    bool isSmallerThan20 = ( vec.size() < 20 );
    if ( isSmallerThan20 ){
        for( vector<double>::const_iterator i = vec.begin(); i != vec.end(); ++i){
            cout << *i << ' ';
        }
    }
    cout << endl;
    return 0;
}

int KonGrad::printMatrix (const  vector< vector<double> > &matrix){
    // this prints out the matrix if the size is smaller than 20.
    bool isSmallerThan20 = ( matrix.size() < 20 );
    if ( isSmallerThan20 ){
        for ( unsigned int i = 0; i<matrix.size(); ++i){
            for( vector<double>::const_iterator j = matrix.at(i).begin(); j != matrix.at(i).end(); ++j){
                cout << *j << ' ';
            }
            cout << endl;
        }
    }
    return !isSmallerThan20;
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


double KonGrad::getRandomUni(){

    uniform_real_distribution<double> distributiond(0.0,1.0); // to generate the values
    return distributiond(_randGenerator);
}

void KonGrad::startRandomGenerator (double seed){
    _randGenerator.seed(seed);
}



void KonGrad::createRandomSparseSymmetricMatrix(const int dim, vector< vector<double> > &matrixout){
    for (int i=0;i<dim;++i){
        vector<double> line(dim,0);
        int numNonZero=3*getRandomUni();
        for (int j=0;j<numNonZero;++j){
            int posInLine;
            do{
                posInLine=dim*getRandomUni();
            }while(line.at(posInLine)!=0);
            line.at(posInLine)=getRandomUni();
        }
        matrixout.push_back(line);
    }
}

void KonGrad::solve (const string method, const vector< vector<double> > &matrixin, const vector<double> &knownRightSide, const vector<double> &startvec, vector<double> &vecout){
    _A=matrixin;
    _b=knownRightSide;
    solve(method, startvec, vecout);
}


void KonGrad::solve (const string method, const vector<double> &startvec, vector<double> &vecout){
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
    
    BOOST_LOG_TRIVIAL(debug) << "solve: start vector " << printVector(startvec);
    BOOST_LOG_TRIVIAL(debug) << "solve: known right side " << printVector(_b);
    
    if (method=="sparseMatrix"){
        BOOST_LOG_TRIVIAL(debug) << "solve: matrix " << printMatrix(_A);
    }

    applyA(method,startvec,tmpvec);
    diffVector(_b,tmpvec, r);
    BOOST_LOG_TRIVIAL(trace) << "solve: r " << printVector(r);
    if(sqrt(skalarProd(r,r))/bnorm < tol){
        cout << "done" << endl;
        vecout=startvec;
        return; /// @todo write better exit at start
    }
    p=r;
    BOOST_LOG_TRIVIAL(trace) << "solve: p " << printVector(p);
    /// @todo weitere rechnungen in dezidierten Funktionen zusammenfassen
    bool converged=false;
    unsigned int iternum=0;
    while(!converged){
        double alpha;
        double beta;
        double relrest;
        ++iternum;
        applyA(method,p,s);
        BOOST_LOG_TRIVIAL(trace) << "solve: s " << printVector(s);
        
        alpha=skalarProd(p,r)/skalarProd(p,s);
        BOOST_LOG_TRIVIAL(trace) << "solve: alpha " << alpha;
        
        skalarVector(alpha,p,tmpvec);
        sumVector(x,tmpvec,xnew);
        BOOST_LOG_TRIVIAL(trace) << "solve: xnew " << printVector(xnew);
        
        skalarVector(alpha,s,tmpvec);
        diffVector(r,tmpvec,rnew);
        BOOST_LOG_TRIVIAL(trace) << "solve: rnew " << printVector(rnew);
        
        relrest=sqrt(skalarProd(rnew,rnew))/bnorm;
        BOOST_LOG_TRIVIAL(debug) << "relrest: " << relrest;
        if( relrest < tol){
            BOOST_LOG_TRIVIAL(info) << "The algorithm converged. Iterations: " << iternum;
            converged=true;
            BOOST_LOG_TRIVIAL(info) << "resultvector: " << printVector(xnew);
        }
        
        if ( iternum > 2*_A.size() ){
            BOOST_LOG_TRIVIAL(error) << "The algorithm did not converge. Aborted. Iterations: " << iternum;
            break;
        }
        
        BOOST_LOG_TRIVIAL(debug) << "resultvector x at iteration " << iternum << ": " << printVector(xnew);
        
        beta=skalarProd(rnew,rnew)/skalarProd(r,r);
        BOOST_LOG_TRIVIAL(trace) << "solve: beta " << beta;
        
        
        skalarVector(beta,p,tmpvec);
        sumVector(rnew,tmpvec,pnew);
        p=pnew;
        r=rnew;
        x=xnew;
    }
    if (converged){
        vecout=xnew;
    }
    
}



void KonGrad::applyA(const string method, const vector<double> &vecin, vector<double> &vecout){
	if (method=="sparseMatrix"){
		matrixVector(vecin, vecout);
	}
	if (method=="Laplace"){
	    matrixVectorLaplace(vecin, vecout);
	}
}




