#ifndef MATRIX_H
#define MATRIX_H
#include "vector.h"
#include <string>

namespace NLA{
    class Matrix{
        public:
        
        enum class MAT{
            IDENTITY,
            RANDOM,
            RANDOMSYM,
            HILBERT,
        };

        Matrix(int,int);
        Matrix(int,int,MAT);
        ~Matrix();

        int rows, columns;
        double** data;

        // arithemtic operations
        Matrix& operator+(const Matrix&);
        Matrix& operator-(const Matrix&);
        Matrix& operator*(const double);
        Matrix& operator*(const Matrix&);
        Vector& operator*(const Vector&);
        double& operator()(const int,const int);

        // munipulation algorithms
        void transpose();
        bool equals(Matrix*);
        NLA::Vector* getRow(int);
        NLA::Vector* getColumn(int);
        

        //utility algorithms
        bool outputToFile(std::string);
        Matrix* copyMatrix();

        // math algorithms
        double frobeniusNorm();

        private:
        void scale(double);
        void add(const Matrix&);
        void subtract(const Matrix&);
        Matrix& multiply(const Matrix&);
        Vector& multiply(const Vector&);
    };
}
#endif