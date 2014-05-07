#include "kongrad.hh"
#include "global.h"
#include "timedif.h"
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


KonGrad::~KonGrad(){

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

/// Performance: (1 + ndim*2) Flops
void KonGrad::matrixVectorLaplace(const vector<double> &vecin, vector<double> &vecout){
	const int vecinDim = vecin.size();
	vecout.assign(vecinDim,0);
	const double phivar=2*ndim+_mass*_mass;
	int k;
#pragma omp parallel for shared(vecout,vecin,nn,ndim) private(k)
	for (int i=0; i<vecinDim;++i){
        vecout.at(i)=phivar*vecin.at(i); //1 Flop
        for (k=1;k<=ndim;++k){ //ndim times
        	vecout.at(i)-=(vecin.at(nn[k][i])+vecin.at(nn[k+ndim][i])); //2 Flops
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

int KonGrad::calculateKonRate(){
    double gamma=_mass/sqrt(ndim);
    double tol=pow(10,-8);
    double steps = -log(tol)/gamma;
    cout << "gamma " << gamma << endl;
    cout << "steps " << steps << endl;
    return ceil(steps);
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


void KonGrad::addVector(const double alpha, const vector<double> &vecin1, const double beta, const vector<double> &vecin2, vector<double> &vecout){
	const int vecin1Dim = vecin1.size();
	const int vecin2Dim = vecin2.size();
	assert(vecin1Dim == vecin2Dim);
	vecout.assign(vecin1Dim,0);
#pragma omp parallel for shared(vecout, vecin1, vecin2)
	for (int i=0;i<vecin1Dim;++i){
        vecout.at(i)=alpha*vecin1.at(i)+beta*vecin2.at(i); //3 Flops per lattice-dot
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


///Performance: 2 Flops
double KonGrad::skalarProd(const vector<double> &vecin1, const vector<double> &vecin2){
    const int vecin1Dim = vecin1.size();
    const int vecin2Dim = vecin2.size();
    assert(vecin1Dim == vecin2Dim);
    double skalarProd=0;
#pragma omp parallel for shared(vecin1, vecin2) reduction(+: skalarProd)
    for (int i=0;i<vecin1Dim;++i){
        skalarProd+=vecin1.at(i)*vecin2.at(i);//2 Flops
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

void KonGrad::createRandomVector(const int dim, vector<double> &vecout){
	vecout.clear();
    for (int i=0;i<dim;++i){
    	vecout.push_back(getRandomUni());
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
    const int bsize=_b.size();
    vector<double> r;
    vector<double> rnew;
    vector<double> p;
    vector<double> pnew;
    vector<double> s;
    vector<double> x=startvec;
    vector<double> xnew;
    vecout.clear();
    
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
    float totalcputime=0;
    float itercputime=0;
    float totalclocktime=0;
    float iterclocktime=0;
    double alpha;
    double beta;
    double rnewnorm,rnorm;
    rnorm=skalarProd(r,r);
    itercputime=cpudif();
    iterclocktime=clkdif();
    while(!converged){
        ++iternum;
        applyA(method,p,s); //1+ndim*2
        alpha=rnorm/skalarProd(p,s); //2+1=3
        
        addVector(1,x,alpha,p,xnew); //3
        
        addVector(1,r,-alpha,s,rnew); //3
        
        rnewnorm=skalarProd(rnew,rnew); //2
        if( sqrt(rnewnorm)/bnorm < tol){
            BOOST_LOG_TRIVIAL(info) << "The algorithm converged. Iterations: " << iternum;
            converged=true;
            BOOST_LOG_TRIVIAL(info) << "resultvector: " << printVector(xnew);
        }
        
        if ( iternum > 2*bsize ){
            BOOST_LOG_TRIVIAL(error) << "The algorithm did not converge. Aborted. Iterations: " << iternum;
            break;
        }

        beta=0.5*rnewnorm/rnorm; //3
        addVector(1,rnew,beta,p,pnew); //3
        p=pnew;
        r=rnew;
        x=xnew;
        rnorm=rnewnorm;
    }
    itercputime=cpudif();
    iterclocktime=clkdif();
    totalcputime+=itercputime;
    totalclocktime+=iterclocktime;
    int NumberOfFlops = 1+ndim*2 +3+3+3+2+3+3;
    BOOST_LOG_TRIVIAL(info) << "total cpu time: " << totalcputime << " s";
    BOOST_LOG_TRIVIAL(info) << "cpu time per iteration: " << totalcputime/iternum << " s";
    BOOST_LOG_TRIVIAL(info) << "cpu time per iteration and lattice-dot: " << totalcputime/iternum/nvol*pow(10,9) << " ns";
    BOOST_LOG_TRIVIAL(info) << "total clock time: " << totalclocktime << " s";
    BOOST_LOG_TRIVIAL(info) << "clock time per iteration: " << totalclocktime/iternum << " s";
    BOOST_LOG_TRIVIAL(info) << "clock time per iteration and lattice-dot: " << totalclocktime/iternum/nvol*pow(10,9) << " ns";
    BOOST_LOG_TRIVIAL(info) << "cputime/clktime: " << (double)totalcputime/totalclocktime;
    BOOST_LOG_TRIVIAL(info) << "Flops per iteration: " << NumberOfFlops;
    BOOST_LOG_TRIVIAL(info) << "Total: " << NumberOfFlops*iternum*nvol/pow(10,9) << " GFlops";
    BOOST_LOG_TRIVIAL(info) << "Performance: " << NumberOfFlops*iternum*nvol/pow(10,9)/totalcputime << " GFlops/s";
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




