#define BOOST_TEST_MODULE kongrad_test
#include <boost/test/unit_test.hpp>
#include <boost/log/trivial.hpp>
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
struct F {
//     F() : i( 0 ) { std::cout << "setup" << std::endl; }
    F(){
        KonGrad testSLE; // SLE = system of linear equations
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
    
    testSLE.matrixVector(matrix, vecin, vecout);
    
    for (int i=0;i<3;++i){
        BOOST_CHECK_EQUAL(vecout.at(i),i);
    }
    

}


/// @todo solve1 testcase fertig schreiben
BOOST_FIXTURE_TEST_CASE(solve1, F){
    
    
    //create 10x10 unit matrix
    vector< vector<double> > matrix;
    vector<double> line,b,startvector;
    for (int i=0;i<10;++i){
        line.assign(10,0);
        line.at(i)=1;
        matrix.push_back(line);
        b.push_back(i);
        startvector.push_back(0);
    }
    
    testSLE.setMatrix(matrix);
    testSLE.setb(b);
    
    testSLE.solve(startvector);
    
}

BOOST_FIXTURE_TEST_CASE(getRandomUni, F){
    cout << testSLE.getRandomUni(0)<<endl;
    cout << testSLE.getRandomUni(0)<<endl;
    cout << testSLE.getRandomUni(0)<<endl;
    cout << testSLE.getRandomUni(0)<<endl;
    cout << testSLE.getRandomUni(0)<<endl;
    double seed=time(NULL);
    cout << "seed " << seed << endl;
    cout << testSLE.getRandomUni(seed)<<endl;
    cout << testSLE.getRandomUni(seed)<<endl;
    cout << testSLE.getRandomUni(seed)<<endl;
    cout << testSLE.getRandomUni(seed)<<endl;
}


BOOST_FIXTURE_TEST_CASE(createRandomSparseSymmetricMatrix_dim, F){
    
    
    vector< vector<double> > matrix;
    int dim = 1000;
    int seed =0;
    
    testSLE.createRandomSparseSymmetricMatrix(dim, seed, matrix);
    
    BOOST_CHECK_EQUAL(matrix.size(),dim);
    
    for (int i=0;i<dim;++i){
        BOOST_CHECK_EQUAL(matrix.at(i).size(),dim);
    }
}

BOOST_FIXTURE_TEST_CASE(createRandomSparseSymmetricMatrix_numNonZero, F){
    
    vector< vector<double> > matrix;
    int dim = 1000;
    int seed =0;
    
    testSLE.createRandomSparseSymmetricMatrix(dim, seed, matrix);
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













