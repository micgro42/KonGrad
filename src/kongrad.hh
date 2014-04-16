#ifndef KONGRAD_HH
#define KONGRAD_HH
#include <vector>

using namespace std;

///@author Michael Große

/**
 * @class KonGrad
 * 
 * @brief Class to store and solve linear system of equations
 * 
 */
class KonGrad{
    public:
        ///The Constructor
        KonGrad(vector<vector<double> > Matrix, vector<double> vec);
        
        /**
         * @brief test the matrix multiplication routine
         * 
         */
        void testmv (vector<double> vecin);
        
        double skalarProd(const vector<double> &vecin1, const vector<double> &vecin2);
        
        
        /**
         * @brief subtracts one vector from another and writes the result to a third
         */
        void diffVector(const vector<double> &vecin1, const vector<double> &vecin2, vector<double> &vecout);
        
        /// sum of two vectors
        void sumVector(const vector<double> &vecin1, const vector<double> &vecin2, vector<double> &vecout);
        
        ///skalar times vector
        void skalarVector(const double alpha, const vector<double> &vecin, vector<double> &vecout);
        
        /**
         * @brief Mulitply a matrix and a vector
         * 
         * @param[in] matrix ...
         * 
         * @param[in] vecin ...
         * 
         * @param[out] vecout ...
         * 
         */
        void matrixVector(const vector< vector<double> > &matrix, const vector<double> &vecin, vector<double> &vecout);
        
        void solve (const vector<double> &startvec);
        void solve ();
        
    private:
        
        /**
         * @brief print the vector to stdout
         */
        void printVector (const vector<double> &vec);
        
        
        
        
        /// (sparse) matrix
        vector<vector<double> > _A;
        
        /// "known right side"
        vector<double> _b;
};

#endif
