#define BOOST_TEST_MODULE kongrad_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include "kongrad.hh"
#define DEFINE_GLOBAL
//#include "global.h"
#include "geom_pbc.c"
#include <vector>
/**
 * @file test.cc
 * 
 * @brief this file contains the unit test suite for the KonGrad class
 * 
 * 
 * 
 * 
 */
 

/**
 * 
 * @class F
 * 
 * @brief Fixture to setup and teardown the class instances for the unit tests where neccesary
 * 
 * 
 */

namespace logging = boost::log;

struct F {

//     F() : i( 0 ) { std::cout << "setup" << std::endl; }
    F(){
        KonGrad testSLE; // SLE = system of linear equations
        logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::info);
    }
    ~F()          {  }
    
    KonGrad testSLE;

};

struct scalarfield {
	scalarfield(){
		KonGrad testSLE; // SLE = system of linear equations
		logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::trace);
		ndim=2;

		lsize = (int *) malloc((ndim+1) * sizeof(int));
        lsize[1]=3;
	    lsize[2]=3;
	    geom_pbc();

	}
	~scalarfield() {   }

	KonGrad testSLE;
};


BOOST_AUTO_TEST_SUITE (kongrad_test) 




/**
 * @brief test if that the the Laplaceoperator is correctly applied to a small 2-D field
 *
 * @details The Input field looks like:<br>
 * 1 2 1 <br>
 * 2 3 2<br>
 * 1 2 1
 *
 * The output field should look like:<br>
 * -2 1 -2<br>
 * 1 4 1<br>
 * -2 1 -2
 *
 * This unittest tests the function KonGrad::matrixVectorLaplace
 */
BOOST_FIXTURE_TEST_CASE( matrixVectorLaplace_m0, scalarfield ){
	vector <double> vecin;
	vector <double> vecout;
	testSLE.setmass(0);

	vecin.push_back(1);
	vecin.push_back(2);
	vecin.push_back(1);
	vecin.push_back(2);
	vecin.push_back(3);
	vecin.push_back(2);
	vecin.push_back(1);
	vecin.push_back(2);
	vecin.push_back(1);
	testSLE.matrixVectorLaplace(vecin, vecout);
    BOOST_CHECK_EQUAL(vecout.at(0),-2);
    BOOST_CHECK_EQUAL(vecout.at(1),1);
    BOOST_CHECK_EQUAL(vecout.at(2),-2);
    BOOST_CHECK_EQUAL(vecout.at(3),1);
    BOOST_CHECK_EQUAL(vecout.at(4),4);
    BOOST_CHECK_EQUAL(vecout.at(5),1);
    BOOST_CHECK_EQUAL(vecout.at(6),-2);
    BOOST_CHECK_EQUAL(vecout.at(7),1);
    BOOST_CHECK_EQUAL(vecout.at(8),-2);

}




/**
 * test if that the skalar product of (2;2) and (2;2) is indeed 8.
 *
 * This unittest tests the function KonGrad::skalarProd
 */
BOOST_FIXTURE_TEST_CASE( skalarProd, F ){
    
    vector<double> b(2,2);
    
    BOOST_CHECK_EQUAL( testSLE.skalarProd(b,b) , 8 );
}


/**
 * test if some small integers are subtracted from each other correctly.
 *
 * This unittest tests the function KonGrad::diffVector
 */
BOOST_FIXTURE_TEST_CASE( diffVector, F ){
    
    vector<double> a,b,diff;
    
    for (int i=1;i<=5;++i){
        a.push_back(i);
        b.push_back(i+1);
    }
    
    BOOST_REQUIRE(a.size()==b.size());
    
    testSLE.diffVector(a, b, diff);
    
    for (int i=0;i<5;++i){
        BOOST_CHECK_EQUAL(diff.at(i),-1);
    }
}


/**
 * test if some small integers are summed correctly.
 *
 * This unittest tests the function KonGrad::sumVector
 */
BOOST_FIXTURE_TEST_CASE( sumVector, F ){
    
    vector<double> a,b,sum;
    
    for (int i=1;i<=5;++i){
        a.push_back(i);
        b.push_back(i+1);
    }
    
    BOOST_REQUIRE(a.size()==b.size());
    
    testSLE.sumVector(a, b, sum);
    
    for (int i=0;i<5;++i){
        BOOST_CHECK_EQUAL(sum.at(i),(i+i+3));
    }
}


/**
 * test if a vector of precise doubles is correctly multplied by another double
 *
 * This unittest tests the function KonGrad::skalarVector .
 */
BOOST_FIXTURE_TEST_CASE( skalarVector, F ){
    vector<double> vecin,prod;
    double scalar=2.5;
    
    for (int i=0;i<5;++i){
        vecin.push_back(i+0.5);
    }
    
    testSLE.skalarVector(scalar, vecin, prod);
    
    double truth;
    for (int i=0;i<5;++i){
        truth=(i+0.5)*2.5;
        BOOST_CHECK_EQUAL(prod.at(i),truth);
    }
}


/**
 * @brief test that a diagonal matrix and a vector are multiplied correctly
 *
 * @details This unittest tests the function KonGrad::matrixVector
 */
BOOST_FIXTURE_TEST_CASE( matrixVector, F ){
    
    vector<double> vecin(3,1);
    vector<double> vecout;
    
    vector< vector<double> > matrix;
    vector<double> line;
    for (int i=0;i<3;++i){
        line.assign(3,0);
        line.at(i)=i;
        matrix.push_back(line);
    }
    testSLE.setMatrix(matrix);
    testSLE.matrixVector(vecin, vecout);
    
    for (int i=0;i<3;++i){
        BOOST_CHECK_EQUAL(vecout.at(i),i);
    }
    

}





/**
 * @brief test that a diagonal matrix and a vector are multiplied correctly
 *
 * @details This unittest tests the function KonGrad::solve(const vector<double> &startvec, vector<double> &vecout)
 */
BOOST_FIXTURE_TEST_CASE(solve_sparseMatrix, F){
    
    
    //create 10x10 unit matrix
    vector< vector<double> > matrix;
    vector<double> line,b,startvector;
    for (int i=0;i<10;++i){
        line.assign(10,0);
        line.at(i)=i+1;
        matrix.push_back(line);
        b.push_back(i+1);
        startvector.push_back(2);
    }
    
    testSLE.setMatrix(matrix);
    testSLE.setb(b);
    vector<double> resultvector;
    testSLE.solve("sparseMatrix",startvector, resultvector);
    
    for( vector<double>::const_iterator i = resultvector.begin(); i != resultvector.end(); ++i){
        BOOST_CHECK_CLOSE(*i,1,0.000001); //tolerance 10^-8
    }
    
}


/**
 *
 *
 * @details this test goes the opposite direction of the test BOOST_FIXTURE_TEST_CASE( matrixVectorLaplace_m0, scalarfield )
 * the startvector is the "known right side"
 */
BOOST_FIXTURE_TEST_CASE(solve_LaplaceOp_periodic_1, scalarfield){
	vector <double> b,startvector,resultvector;
	testSLE.setmass(0);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(4);
	b.push_back(1);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(-2);
	testSLE.setb(b);
	testSLE.solve("Laplace",b, resultvector);
	BOOST_CHECK_EQUAL(resultvector.at(0),1);
	BOOST_CHECK_EQUAL(resultvector.at(1),2);
	BOOST_CHECK_EQUAL(resultvector.at(2),1);
	BOOST_CHECK_EQUAL(resultvector.at(3),2);
	BOOST_CHECK_EQUAL(resultvector.at(4),3);
	BOOST_CHECK_EQUAL(resultvector.at(5),2);
	BOOST_CHECK_EQUAL(resultvector.at(6),1);
	BOOST_CHECK_EQUAL(resultvector.at(7),2);
	BOOST_CHECK_EQUAL(resultvector.at(8),1);
}


/// the startvector is always the same number
BOOST_FIXTURE_TEST_CASE(solve_LaplaceOp_periodic_2, scalarfield){
	vector <double> b,startvector,resultvector;
	testSLE.setmass(0);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(4);
	b.push_back(1);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(-2);
	testSLE.setb(b);
	startvector.assign(9,5);
	testSLE.solve("Laplace",startvector, resultvector);
	BOOST_CHECK_EQUAL(resultvector.at(0),1);
	BOOST_CHECK_EQUAL(resultvector.at(1),2);
	BOOST_CHECK_EQUAL(resultvector.at(2),1);
	BOOST_CHECK_EQUAL(resultvector.at(3),2);
	BOOST_CHECK_EQUAL(resultvector.at(4),3);
	BOOST_CHECK_EQUAL(resultvector.at(5),2);
	BOOST_CHECK_EQUAL(resultvector.at(6),1);
	BOOST_CHECK_EQUAL(resultvector.at(7),2);
	BOOST_CHECK_EQUAL(resultvector.at(8),1);
}


///this unittest tests if solve correctly returns the start vector if it already solves the problem.
BOOST_FIXTURE_TEST_CASE(solve_LaplaceOp_periodic_3, scalarfield){
	vector <double> b,startvector,resultvector;
	testSLE.setmass(0);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(4);
	b.push_back(1);
	b.push_back(-2);
	b.push_back(1);
	b.push_back(-2);
	testSLE.setb(b);
	startvector.push_back(1);
	startvector.push_back(2);
	startvector.push_back(1);
	startvector.push_back(2);
	startvector.push_back(3);
	startvector.push_back(2);
	startvector.push_back(1);
	startvector.push_back(2);
	startvector.push_back(1);
	testSLE.solve("Laplace",startvector, resultvector);
	BOOST_CHECK_EQUAL(resultvector.at(0),1);
	BOOST_CHECK_EQUAL(resultvector.at(1),2);
	BOOST_CHECK_EQUAL(resultvector.at(2),1);
	BOOST_CHECK_EQUAL(resultvector.at(3),2);
	BOOST_CHECK_EQUAL(resultvector.at(4),3);
	BOOST_CHECK_EQUAL(resultvector.at(5),2);
	BOOST_CHECK_EQUAL(resultvector.at(6),1);
	BOOST_CHECK_EQUAL(resultvector.at(7),2);
	BOOST_CHECK_EQUAL(resultvector.at(8),1);
}



/**
 * @brief test that a two consecutive random numbers are not the same
 *
 * @details This unittest tests the function KonGrad::getRandomUni
 */
BOOST_FIXTURE_TEST_CASE(getRandomUni, F){
    testSLE.startRandomGenerator(10);
    for (int i=0; i<100000;++i){
        BOOST_CHECK(testSLE.getRandomUni()!=testSLE.getRandomUni());
    }

}

/**
 * @brief test that a the created matrix has the dimensions it is supposed to have.
 *
 * @details This unittest tests the function KonGrad::createRandomSparseSymmetricMatrix
 */
BOOST_FIXTURE_TEST_CASE(createRandomSparseSymmetricMatrix_dim, F){
    
    
    vector< vector<double> > matrix;
    int dim = 1000;
    
    testSLE.createRandomSparseSymmetricMatrix(dim, matrix);
    
    BOOST_CHECK_EQUAL(matrix.size(),dim);
    
    for (int i=0;i<dim;++i){
        BOOST_CHECK_EQUAL(matrix.at(i).size(),dim);
    }
}

/**
 * @brief test that a the created matrix has at most 2 non-zero entries per row
 *
 * @details This unittest tests the function KonGrad::createRandomSparseSymmetricMatrix
 */
BOOST_FIXTURE_TEST_CASE(createRandomSparseSymmetricMatrix_numNonZero, F){
    
    vector< vector<double> > matrix;
    int dim = 1000;
    
    testSLE.createRandomSparseSymmetricMatrix(dim, matrix);
    int nonZeroTotal=0;
    for (int i=0;i<dim;++i){
        int nonZero=0;
        for (int j=0;j<dim;++j){
            if (matrix.at(i).at(j)!=0){
                ++nonZero;
            }
        }
        BOOST_CHECK(nonZero<3);
        nonZeroTotal+=nonZero;
    }
    BOOST_LOG_TRIVIAL(info) << "The dim of testmatrix: " << dim << " nonZeroTotal: " << nonZeroTotal;
}


BOOST_AUTO_TEST_SUITE_END( )













