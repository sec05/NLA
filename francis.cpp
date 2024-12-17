#include "matrixAlgorithms.h"

using namespace NLA;

const double tolerance = pow(10,-15);

/*
An implementation of Francis's algorithm with and without shifts.
Returns a matrix Q of vectors and A which is the eigenvalues.
*/

Matrix** francis(Matrix* A, int iterations){
    Matrix* vectors  = new Matrix(A->rows,A->rows,"identity");
    Matrix* values = A->copyMatrix();

    for(int i = 0; i < iterations; i++){
        Matrix** QR = modifiedGramSchmidt(values);
        *vectors = *vectors * *QR[0];
        values = &(*QR[1] * *QR[0]);
        delete QR;
    }
    Matrix** arr = (Matrix**) malloc(sizeof(Matrix*)*2);
    arr[0] = values;
    arr[1] = vectors;
    delete vectors;
    delete values;
    return arr;
}

/*
An implementation of Francis's algorithm with and without shifts.
Returns a matrix Q of vectors and A which is the eigenvalues.
*/

Matrix** francisWithShifts(Matrix* A, int iterations){
    Matrix* T = householderUpperHessenberg(A);
    Matrix* vectors  = new Matrix(T->rows,T->rows,"identity");
    Matrix* values = T;
    for(int i = 0; i < iterations; i++){
        double mu = (*values)(values->rows,values->rows);
        Matrix* I = new Matrix(T->rows,T->rows,"identity");
        Matrix** Recursed = (Matrix**) malloc(sizeof(Matrix*)*2);
        *I = *I * mu;
        *T = *T - *I;
        Matrix** QR = modifiedGramSchmidt(T);
        *T = *QR[1] * *QR[0];
        *T = *T + *I;
        for(int j = 0; j < T->rows-1; j++){
            if((*T)(j,j+1) <= tolerance || (*T)(j+1,j) <= tolerance){
                (*T)(j,j+1) = 0;
                (*T)(j+1,j) = 0;
               // Recursed = francisRecurse(T,iterations);
                break;
            }
        }

    }
}

Matrix** francisRecurse(Matrix* A, int iterations){


}