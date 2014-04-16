#define BOOST_TEST_MODULE "kongrad test"
#include <boost/test/included/unit_test.hpp>
#include "kongrad.hh"

using namespace std;

BOOST_AUTO_TEST_CASE( skalarProd ){
    
    //erschaffe 2*Einheitsmatrix 3x3
    vector< vector<double> > A;
    vector<double> line;
    for (int i=0;i<2;++i){
        line.assign(2,0);
        line.at(i)=2;
        A.push_back(line);
    }
    
    vector<double> b;
    b.push_back(2);
    b.push_back(2);
    
//     KonGrad LGS01(A, b);
    
//     BOOST_CHECK_EQUAL( LGS01.skalarProd(b,b) , 8 );
    BOOST_CHECK_EQUAL(0,0);
    BOOST_CHECK_EQUAL(0,1);
}