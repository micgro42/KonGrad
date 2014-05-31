#include "linsysequ.hh"
#include "global.h"
#include "timedif.h"
#include "eivtris.h"
//#include "geom_pbc.c"
#include <vector>
#include <random>
#include <cassert>
#include <cmath>
#include <boost/log/trivial.hpp>


#ifndef DEBUG
#define at(x) operator[](x)
#endif

/**
 * @file linsysequ.cc
 * 
 * @brief sourcecode for the LinSysEqu class
 * 
 * 
 */


LinSysEqu::LinSysEqu(vector< vector<double> > matrix, vector<double> vec) : _A(matrix), _b(vec) {
    _mass=0.1;
}

LinSysEqu::LinSysEqu(){
    vector<double> line;
    for (int i=0;i<3;++i){
        line.assign(3,0);
        line.at(i)=1;
        _A.push_back(line);
    }
    _mass=0.1;
    _b.assign(3,0);
}


LinSysEqu::~LinSysEqu(){

}


void LinSysEqu::testmv(const vector<double> vecin){
    vector<double> vecout;
    LinSysEqu::matrixVector(vecin, vecout);
}



void LinSysEqu::matrixVector(const vector<double> &vecin, vector<double> &vecout){
    //check if dimensions are correct
    const int vecinDim = vecin.size();
    const int matrixlineDim = _A.size();
    assert(vecinDim==matrixlineDim);
    vecout.assign(vecinDim,0);
    BOOST_LOG_TRIVIAL(trace) << "matrixVector: vecin " << printVector(vecin);
    ///@todo: make parallel
    for (int i=0;i<vecinDim;++i){
        for (int j=0;j<vecinDim;++j){
            vecout.at(j)+=_A.at(i).at(j)*vecin.at(j);
        }
        
    }
    BOOST_LOG_TRIVIAL(trace) << "matrixVector: vecout " << printVector(vecout);
}

/// Performance: (1 + ndim*2) Flops
void LinSysEqu::matrixVectorLaplace(const vector<double> &vecin, vector<double> &vecout){
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

///@deprecated use LinSysEqu::addVector with a negative scalar instead
void LinSysEqu::diffVector(const vector<double> &vecin1, const vector<double> &vecin2, vector<double> &vecout){
    const int vecin1Dim = vecin1.size();
    const int vecin2Dim = vecin2.size();
    assert(vecin1Dim == vecin2Dim);
    
    vecout.clear();
    
    for (int i=0;i<vecin1Dim;++i){
        vecout.push_back(vecin1.at(i)-vecin2.at(i));
    }
}

int LinSysEqu::calculateKonRate(){
    double gamma=_mass/sqrt(ndim);
    double tol=pow(10,-8);
    double steps = -log(tol)/gamma;
    cout << "gamma " << gamma << endl;
    cout << "steps " << steps << endl;
    return ceil(steps);
}

///@deprecated use LinSysEqu::addVector with 1 as scalar instead
void LinSysEqu::sumVector(const vector<double> &vecin1, const vector<double> &vecin2, vector<double> &vecout){
    const int vecin1Dim = vecin1.size();
    const int vecin2Dim = vecin2.size();
    assert(vecin1Dim == vecin2Dim);
    
    vecout.clear();
    
    for (int i=0;i<vecin1Dim;++i){
        vecout.push_back(vecin1.at(i)+vecin2.at(i));
    }
}


void LinSysEqu::addVector(const double alpha, const vector<double> &vecin1, const double beta, const vector<double> &vecin2, vector<double> &vecout){
	const int vecin1Dim = vecin1.size();
	const int vecin2Dim = vecin2.size();
	const int vecoutDim = vecout.size();
	assert(vecin1Dim == vecin2Dim);
	if (vecoutDim!=vecin1Dim){
		vecout.assign(vecin1Dim,0);
	}
#pragma omp parallel for shared(vecout, vecin1, vecin2)
	for (int i=0;i<vecin1Dim;++i){
		vecout.at(i)=alpha*vecin1.at(i)+beta*vecin2.at(i); //3 Flops per lattice-dot
    }
}

void LinSysEqu::skalarVector(const double alpha, const vector<double> &vecin, vector<double> &vecout){
    const unsigned int vecinDim = vecin.size();
    if (vecout.size()!=vecinDim){
    	vecout.assign(vecinDim,0);
    }
    ///@todo make parallel
    for (unsigned int i=0;i<vecinDim;++i){
        vecout.at(i)=vecin.at(i)*alpha;
    }
}

int LinSysEqu::printVector (const vector<double> &vec){
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

int LinSysEqu::printMatrix (const  vector< vector<double> > &matrix){
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


///Performance: 2 Flops per element of vecin
double LinSysEqu::skalarProd(const vector<double> &vecin1, const vector<double> &vecin2){
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


double LinSysEqu::getRandomUni(){

    uniform_real_distribution<double> distributiond(0.0,1.0); // to generate the values
    return distributiond(_randGenerator);
}

void LinSysEqu::startRandomGenerator (double seed){
    _randGenerator.seed(seed);
}



void LinSysEqu::createRandomSparseSymmetricMatrix(const int dim, vector< vector<double> > &matrixout){
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

void LinSysEqu::createRandomVector(const int dim, vector<double> &vecout){
	vecout.clear();
    for (int i=0;i<dim;++i){
    	vecout.push_back(getRandomUni());
    }
}

void LinSysEqu::solveLSE (const string method, const vector< vector<double> > &matrixin, const vector<double> &knownRightSide, const vector<double> &startvec, vector<double> &vecout){
    _A=matrixin;
    _b=knownRightSide;
    solveLSE(method, startvec, vecout);
}

///@return 0 if everything went well, 80 if there are too many interations, and 81 if the startvector already solves the system
int LinSysEqu::solveLSE (const string method, const vector<double> &startvec, vector<double> &vecout){
	int exitcode=1;
    const double tol=pow(10,-8);
    const double bnorm=sqrt(skalarProd(_b,_b));
    const unsigned int bsize=_b.size();
    vector<double> r,p,s,xnew;
    vector<double> x=startvec;
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
    	BOOST_LOG_TRIVIAL(warning) << "The input vector already solves the system. Exiting.";
        vecout=startvec;
        exitcode = 81;
        return exitcode; /// @todo write better exit at start
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
        alpha=rnorm/skalarProd(p,s); //2
        
        addVector(1,x,alpha,p,xnew); //2
        
        addVector(1,r,-alpha,s,r); //2
        
        rnewnorm=skalarProd(r,r); //2
        BOOST_LOG_TRIVIAL(debug) << "iteration: " << iternum << " norm of the rest: " << rnewnorm;
        if( sqrt(rnewnorm)/bnorm < tol){
            BOOST_LOG_TRIVIAL(info) << "The algorithm converged. Iterations: " << iternum;
            exitcode=0;
            converged=true;
            BOOST_LOG_TRIVIAL(info) << "resultvector: " << printVector(xnew);
        }
        
        if ( iternum > 2*bsize ){
            BOOST_LOG_TRIVIAL(error) << "The algorithm did not converge. Aborted. Iterations: " << iternum;
            exitcode=80;
            break;
        }

        beta=0.5*rnewnorm/rnorm;
        addVector(1,r,beta,p,p); //2
        x=xnew;///@todo: zuweisung parallelisieren bzw. prüfen ob x überschrieben werden können
        rnorm=rnewnorm;
    }
    if (converged){
    	itercputime=cpudif();
    	iterclocktime=clkdif();
    	totalcputime+=itercputime;
    	totalclocktime+=iterclocktime;
    	int NumberOfFlops = 1+ndim*2 +2+2+2+2+2;
    	BOOST_LOG_TRIVIAL(info) << "total cpu time: " << totalcputime << " s";
    	BOOST_LOG_TRIVIAL(info) << "cpu time per iteration: " << totalcputime/iternum << " s";
    	BOOST_LOG_TRIVIAL(info) << "cpu time per iteration and lattice-dot: " << totalcputime/iternum/nvol*pow(10,9) << " ns";
    	BOOST_LOG_TRIVIAL(info) << "total clock time: " << totalclocktime << " s";
    	BOOST_LOG_TRIVIAL(info) << "clock time per iteration: " << totalclocktime/iternum << " s";
    	BOOST_LOG_TRIVIAL(info) << "clock time per iteration and lattice-dot: " << totalclocktime/iternum/nvol*pow(10,9) << " ns";
    	BOOST_LOG_TRIVIAL(info) << "cputime/clktime: " << (double)totalcputime/totalclocktime;
    	BOOST_LOG_TRIVIAL(info) << "Flops per iteration: " << NumberOfFlops;
    	BOOST_LOG_TRIVIAL(info) << "Total: " << (double)NumberOfFlops*iternum*nvol/pow(10,9) << " GFlops";
    	BOOST_LOG_TRIVIAL(info) << "Performance: " << (double)NumberOfFlops*iternum*nvol/pow(10,9)/totalcputime << " GFlops/s";
        vecout=xnew;
    }
    return exitcode;
    
}

void LinSysEqu::eigenvLanczos(const string method, const vector<double> &startvec, vector<double> &vecout){
	double s = skalarProd(startvec,startvec);
	const double tol=pow(10,-8);
	if (s==0){
		vecout=startvec;
		return;
	}
	vector<double> v,q,alpha,beta,vnp1;
	skalarVector(1/s,startvec,v);
	applyA(method,v,q);
	bool converged=false;
	unsigned int iternum=0;
	while (!converged){
		alpha.push_back(skalarProd(v,q));
		addVector(1,q,-alpha.at(iternum),v,q);
		beta.push_back(skalarProd(q,q));
		if(beta.at(iternum)<tol){
			BOOST_LOG_TRIVIAL(info) << "The algorithm converged. Iterations: " << iternum;
			converged=true;
		}
		if(iternum>startvec.size()){
			BOOST_LOG_TRIVIAL(error) << "The algorithm did not converge. Aborted. Iterations: " << iternum;
			break;
		}
		skalarVector(1/beta.at(iternum),q,vnp1);
		applyA(method,vnp1,q);
		addVector(1,q,-beta.at(iternum),v,q);
		++iternum;
		v=vnp1;
	}

}



void LinSysEqu::applyA(const string method, const vector<double> &vecin, vector<double> &vecout){
	if (method=="sparseMatrix"){
		matrixVector(vecin, vecout);
	}
	if (method=="Laplace"){
	    matrixVectorLaplace(vecin, vecout);
	}
}




