#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include "matrix.h"

size_t luDecomposition(const Matrix& A, Matrix& P, Matrix& L, Matrix& U);

Vector Forward_Substitution(const Matrix& L, const Vector& Pb);

Vector Backward_Substitution(const Matrix& U, const Vector& y);

Matrix Inverse(const Matrix& A, const Matrix& L, const Matrix& U, const Matrix& P);

Vector E_Column(size_t n, size_t ind);

Matrix Single_Matrix(size_t n);

double Determinant(const Matrix& U, size_t count_perm);

void Print_Result(const Matrix& A, const Vector& b);

#endif // lu_decomposition