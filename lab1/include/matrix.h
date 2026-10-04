#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix Multiply_Matrix(const Matrix &A, const Matrix &B);

void Print_Matrix(const Matrix &A);

void Print_Vector(const Vector &V);

bool Equal_Matrix(const Matrix &A, const Matrix &B);

Vector Multiply_Vector(const Matrix &A, const std::vector<double> &V);

std::string To_Mixed_Fraction(double val, double eps, int max_den);

#endif // MATRIX_H