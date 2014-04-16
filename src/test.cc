#define BOOST_TEST_MODULE kongrad_test
#include <boost/test/included/unit_test.hpp>
#include "kongrad.hh"
#include <vector>
#include <iostream>
 
struct F {
//     F() : i( 0 ) { std::cout << "setup" << std::endl; }
    F(){
        std::cout << "setup" << std::endl; 
        KonGrad testSLE; // SLE = system of linear equations
    }
    ~F()          { std::cout << "teardown" << std::endl; }
    
    KonGrad testSLE;

};

BOOST_AUTO_TEST_SUITE (kongrad_test) 


BOOST_FIXTURE_TEST_CASE( skalarProd, F ){
    
    vector<double> b;
    b.push_back(2);
    b.push_back(2);
    
//     KonGrad LGS01;
    
    BOOST_CHECK_EQUAL( testSLE.skalarProd(b,b) , 8 );
}

BOOST_AUTO_TEST_SUITE_END( )