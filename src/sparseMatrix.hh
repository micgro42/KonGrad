#ifndef SPARSEMATRIX_HH
#define SPARSEMATRIX_HH
#include <vector>
/**
 * @file sparseMatrix.hh
 * 
 * @brief header for the sparseMatrix class
 * 
 * 
 */

using namespace std;

/**
 * @class sparseMatrix
 * 
 * @brief Class to efficently store a (sparse) Matrix and perform some operations of it/with it.
 * 
 * @details right now it is still only a vector of vectors of double, but I will add a more efficent format later. 
 * 
 * @todo implement more efficent storage format
 * 
 */
class sparseMatrix{
    public:
    
        sparseMatrix(vector<vector<double> > matrix);
    
        sparseMatrix();
        
        unsigned int size();
        
        void pushLine(vector<double> &line);
        
        /// set the Matrix of the system of linear equations
        void setMatrix(const vector< vector<double> > &matrix) {_A=matrix;};
        
        /// set the Matrix of the object
        void getMatrix(vector< vector<double> > &matrix) {matrix=_A;};
        
        /**
         * @brief Mulitply the matrix and a vector
         * 
         * @param[in] vecin this vector remains unchanged and should have the same dimension as the matrix in this instance.
         * 
         * @param[out] vecout the content of this vector is removed and the new vector is written to it
         * 
         */
        void matrixVector(const vector<double> &vecin, vector<double> &vecout);
        
        /// check if a matrix is symmetric
        bool matrixIsSymmetric();
        
        /// make a matrix is symmetric
        void makeMatrixSymmetric();
        
        /**
         * @brief print the vector to stdout if the dim is < 20
         */
        int printVector (const vector<double> &vec);
        
        /**
         * @brief print the matrix to stdout if the dim is < 20
         */
        int printMatrix ();
    
    
    private:
        
        
        vector<vector<double> > _A;
};


#endif
