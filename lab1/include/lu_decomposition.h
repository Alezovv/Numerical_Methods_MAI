#ifndef LU_DECOMPOSITION_H
#define LU_DECOMPOSITION_H

#include <iomanip>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include "matrix.h"

size_t luDecomposition(const Matrix &A, Matrix &P, Matrix &L, Matrix &U);

Vector Forward_Substitution(const Matrix &L, const Vector &Pb);

Vector Backward_Substitution(const Matrix &U, const Vector &y);

Matrix Inverse(const Matrix &A, const Matrix &L, const Matrix &U, const Matrix &P);

Vector E_Column(size_t n, size_t ind);

Matrix Single_Matrix(size_t n);

double Determinant(const Matrix &U, size_t count_perm);

void Solve_Task_1_1(const std::string &filename);

#endif // lu_decomposition