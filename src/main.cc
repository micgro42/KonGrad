#include "kongrad.hh"
#include <iostream>

using namespace std;


int main(){
    cout << "hello world" << endl;
    
    //erschaffe Einheitsmatrix 3x3
    vector< vector<double> > A;
    vector<double> line;
    for (int i=0;i<3;++i){
        line.assign(3,0);
        line.at(i)=1;
        A.push_back(line);
    }
    
    vector<double> b;
    b.push_back(1);
    b.push_back(2);
    b.push_back(3);
    
    KonGrad LGS01(A, b);
    LGS01.testmv(b);
    
    return 0;
}


