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
        
        ///The default constructor, creates 3x3 unitmatirx and 3-dim 0-vector
        KonGrad();
        
        void setMatrix(const vector< vector<double> > &matrix) {_A=matrix;};
        void getMatrix(vector< vector<double> > &matrix) {matrix=_A;};
        
        void setb(const vector<double> &bvec) {_b=bvec;};
        void getb(vector<double> &bvec) {bvec=_b;};
        
        
        
        /**
         * @brief test the matrix multiplication routine
         * 
         */
        void testmv (vector<double> vecin);
        
        /**
         * @brief create the skalar product of two vectors
         * 
         */
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
        
        /// check if a matrix is symmetric
        bool matrixIsSymmetric(const vector< vector<double> > &matrix);
        
        /// make a matrix is symmetric
        void makeMatrixSymmetric(const vector< vector<double> > &matrixin, vector< vector<double> > &matrixout);
        
        
        void createRandomSparseSymmetricMatrix(int seed, vector< vector<double> > &matrixout);
        
        
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
