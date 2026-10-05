#include "matrix.h"

std::string To_Mixed_Fraction(double val, double eps, int max_den)
{
    if (std::abs(val) < eps)
        return "0";

    bool negative = val < 0;
    val = std::abs(val);

    long long h0 = 0, h1 = 1, h2 = 0;
    long long k0 = 1, k1 = 0, k2 = 0;
    double x = val;

    for (int iter = 0; iter < 15; ++iter)
    {
        long long a = static_cast<long long>(std::floor(x));
        h2 = a * h1 + h0;
        k2 = a * k1 + k0;

        if (k2 > max_den)
            break;

        h0 = h1;
        h1 = h2;
        k0 = k1;
        k1 = k2;

        double frac = static_cast<double>(h1) / k1;
        if (std::abs(frac - val) < eps)
            break;

        double rem = x - a;
        if (std::abs(rem) < 1e-10)
            break;
        x = 1.0 / rem;
    }

    std::ostringstream oss;
    if (negative)
        oss << "-";

    long long whole = h1 / k1;
    long long num = h1 % k1;

    if (k1 == 1)
    {
        oss << whole;
    }
    else if (whole == 0)
    {
        oss << num << "/" << k1;
    }
    else
    {
        oss << whole << " " << num << "/" << k1;
    }

    return oss.str();
}

Matrix Multiply_Matrix(const Matrix &A, const Matrix &B)
{
    if (A.empty() || B.empty())
        return {};

    size_t rows_A = A.size();
    size_t cols_A = A[0].size();
    size_t cols_B = B[0].size();

    Matrix res(rows_A, std::vector<double>(cols_B, 0.0));

    for (size_t i = 0; i < rows_A; i++)
    {
        for (size_t j = 0; j < cols_B; j++)
        {
            for (size_t k = 0; k < cols_A; k++)
            {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

Vector Multiply_Vector(const Matrix &A, const std::vector<double> &V)
{
    if (A.empty())
        return {};

    size_t rows = A.size();
    size_t cols = A[0].size();

    std::vector<double> res(rows, 0.0);

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t k = 0; k < cols; k++)
        {
            res[i] += A[i][k] * V[k];
        }
    }

    return res;
}

void Print_Matrix(const Matrix &A)
{
    for (const auto &row : A)
    {
        for (double value : row)
        {
            std::cout << std::setw(14) << To_Mixed_Fraction(value, 1e-5, 1000);
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

void Print_Vector(const Vector &V)
{
    for (double value : V)
    {
        std::cout << std::setw(14) << To_Mixed_Fraction(value, 1e-5, 1000);
    }
    std::cout << '\n';
}

bool Equal_Matrix(const Matrix &A, const Matrix &B)
{
    if (A.size() != B.size())
        return false;
    if (A.empty())
        return true;
    if (A[0].size() != B[0].size())
        return false;

    size_t rows = A.size();
    size_t cols = A[0].size();
    const double EPS = 1e-9;

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            if (std::abs(A[i][j] - B[i][j]) > EPS)
            {
                return false;
            }
        }
    }

    return true;
}