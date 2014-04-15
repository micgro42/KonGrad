#ifndef KONGRAD_HH
#define KONGRAD_HH
#include <vector>

using namespace std;

class KonGrad{
    public:
        KonGrad(vector<vector<double> > Matrix, vector<double> b);
        void testmv (vector<double> vecin);
        void printVector (const vector<double> &vec);
    private:
        void matrixVector(const vector< vector<double> > &Matrix, const vector<double> &vecin, vector<double> &vecout);

        vector<vector<double> > A;
        vector<double> b;
};

#endif
