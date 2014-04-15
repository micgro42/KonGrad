#ifndef KONGRAD_HH
#define KONGRAD_HH
#include <vector>

using namespace std;

class KonGrad{
    public:
        ConGrad(vector<vector<double> > Matrix, vector<double> b);
        void testmv (vector<double> vecin);
    private:
        void matrixVector(vector< vector<double> > Matrix, vector<double> vecin, vector<double> vecout);

        vector<vector<double> > Matrix;
        vector<double> b;
};

#endif
