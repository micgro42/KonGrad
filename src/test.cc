#define BOOST_TEST_MODULE kongrad_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include "kongrad.hh"
#include <vector>
/**
 * @file test.cc
 * 
 * @brief this file contains the unit test suite for the conjugate gradiant project
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

BOOST_AUTO_TEST_SUITE (kongrad_test) 


BOOST_FIXTURE_TEST_CASE( skalarProd, F ){
    
    vector<double> b(2,2);
    
    BOOST_CHECK_EQUAL( testSLE.skalarProd(b,b) , 8 );
}

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



BOOST_FIXTURE_TEST_CASE(solve1, F){
    
    
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
    testSLE.solve(startvector, resultvector);
    
    for( vector<double>::const_iterator i = resultvector.begin(); i != resultvector.end(); ++i){
        BOOST_CHECK_CLOSE(*i,1,0.000001); //tolerance 10^-8
    }
    
}


BOOST_FIXTURE_TEST_CASE(getRandomUni, F){
    testSLE.startRandomGenerator(10);
    for (int i=0; i<100000;++i){
        BOOST_CHECK(testSLE.getRandomUni()!=testSLE.getRandomUni());
    }

}


BOOST_FIXTURE_TEST_CASE(createRandomSparseSymmetricMatrix_dim, F){
    
    
    vector< vector<double> > matrix;
    int dim = 1000;
    
    testSLE.createRandomSparseSymmetricMatrix(dim, matrix);
    
    BOOST_CHECK_EQUAL(matrix.size(),dim);
    
    for (int i=0;i<dim;++i){
        BOOST_CHECK_EQUAL(matrix.at(i).size(),dim);
    }
}

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













