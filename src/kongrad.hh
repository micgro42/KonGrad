#ifndef KONGRAD_HH
#define KONGRAD_HH
#include <vector>

using namespace std;

///@author Michael Große

/**
 * @class KonGrad
 * 
 * @brief Class to store a linear system of equations
 * 
 */
class KonGrad{
    public:
        ///The Constructor
        KonGrad(vector<vector<double> > Matrix, vector<double> b);
        
        /**
         * @brief test the matrix multiplication routine
         * 
         */
        void testmv (vector<double> vecin);
    private:
        
        /**
         * @brief print the vector to stdout
         */
        void printVector (const vector<double> &vec);
        
        /**
         * @brief Mulitply a matrix and a vector
         * 
         * @param[in] matrix ...
         * 
         * @param[out] vecout ...
         * 
         */
        void matrixVector(const vector< vector<double> > &Matrix, const vector<double> &vecin, vector<double> &vecout);
        
        /// (sparse) matrix
        vector<vector<double> > A;
        
        /// "known right side"
        vector<double> b;
};

#endif
