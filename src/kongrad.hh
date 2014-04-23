#ifndef KONGRAD_HH
#define KONGRAD_HH
#include <vector>
#include <random>

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
        
        /// set the Matrix of the system of linear equations
        void setMatrix(const vector< vector<double> > &matrix) {_A=matrix;};
        
        /// set the Matrix of the object
        void getMatrix(vector< vector<double> > &matrix) {matrix=_A;};
        
        /// set the "known right side"
        void setb(const vector<double> &bvec) {_b=bvec;};
        
        /// set the "known right side"
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
         * @param[in] matrix this vector remains unchanged
         * 
         * @param[in] vecin this vector remains unchanged
         * 
         * @param[out] vecout the content of this vector is removed and the new vector is written to it
         * 
         */
        void matrixVector(const vector< vector<double> > &matrix, const vector<double> &vecin, vector<double> &vecout);
        
        /// check if a matrix is symmetric
        bool matrixIsSymmetric(const vector< vector<double> > &matrix);
        
        /// make a matrix is symmetric
        void makeMatrixSymmetric(const vector< vector<double> > &matrixin, vector< vector<double> > &matrixout);
        
        
        ///starts random generator
        void startRandomGenerator(double seed);
        
        /// returns a random double number between 0 and 1 in an uniform distribution
        double getRandomUni();
        
        /// creates a sparse matrix
        void createRandomSparseSymmetricMatrix(const int dim, vector< vector<double> > &matrixout);
        
        /**
         * 
         * @brief function to solve the system of linear equations
         * 
         * @param[in] startvec vector from which to start
         * 
         * @param[out] vecout vector to which the result is written. potential content will be deleted.
         * 
         * 
         */
        void solve (const vector<double> &startvec, vector<double> &vecout);
        
        /**
         * @brief solve with option to handover the matrix and known right side to the solve function
         * 
         * 
         * 
         * 
         */
        void solve (const vector< vector<double> > &matrixin, const vector<double> &knownRightSide, const vector<double> &startvec, vector<double> &vecout);
        
        
        void solve ();
        
    private:
        
        /**
         * @brief print the vector to stdout
         */
        int printVector (const vector<double> &vec);
        
        /**
         * @brief print the matrix to stdout if the dim is < 20
         */
        int printMatrix (const  vector< vector<double> > &matrix);
        
        
        
        
        /**
         *
         *@brief a (sparse) matrix
         *
         * @todo make a new class for A
         * 
         */
        vector<vector<double> > _A;
        
        /// "known right side"
        vector<double> _b;
        
        ///random generator
        default_random_engine _randGenerator;
};

#endif
