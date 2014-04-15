#include "kongrad.hh"
#include <iostream>

using namespace std;


int main(){
    
    //erschaffe 2*Einheitsmatrix 3x3
    vector< vector<double> > A;
    vector<double> line;
    for (int i=0;i<3;++i){
        line.assign(3,0);
        line.at(i)=2;
        A.push_back(line);
    }
    
    vector<double> b;
    b.push_back(1);
    b.push_back(2);
    b.push_back(3);
    
    KonGrad LGS01(A, b);
    LGS01.testmv(b);
    
    vector<double> c;
    c.push_back(1);
    c.push_back(1);
    c.push_back(1);
    LGS01.solve(c);
    
    return 0;
}


