#include "lu_decomposition.h"

size_t luDecomposition(const Matrix& A, Matrix& P, Matrix& L, Matrix& U) {
    size_t n = A.size();
    size_t count_perm = 0;
    U = A;

    L = Matrix(n, std::vector<double>(n, 0));
    P = Matrix(n, std::vector<double>(n, 0));

    for (size_t i = 0; i < n; i++) {
        L[i][i] = 1;
        P[i][i] = 1;
    }

    for (size_t i = 0; i < n - 1; i++) {
        double mmax = std::abs(U[i][i]);
        size_t n_piv = i;

        for (size_t j = i + 1; j < n; j++) {
            if (mmax < std::abs(U[j][i])) {
                mmax = std::abs(U[j][i]);
                n_piv = j;
            }
        }

        if (n_piv != i) {
            std::swap(U[i], U[n_piv]);
            std::swap(P[i], P[n_piv]);
            count_perm++;

            for (size_t j = 0; j < i; j++) {
                std::swap(L[i][j], L[n_piv][j]);
            }
        }

        for (size_t j = i + 1; j < n; j++) {
            L[j][i] = U[j][i] / U[i][i];

            for (size_t k = i; k < n; k++) {
                U[j][k] -= L[j][i] * U[i][k];
            }
        }
    }
    return count_perm;
}

Vector Forward_Substitution(const Matrix& L, const Vector& Pb) {
    size_t n = L.size();
    Vector y = Vector(n);
    for (int i = 0; i < n; i++) {
        double summary = 0;
        for (int j = 0; j < i; j++) {
            summary += L[i][j] * y[j];
        }
        y[i] = Pb[i] - summary;
    }
    return y;
}

Vector Backward_Substitution(const Matrix& U, const Vector& y) {
    size_t n = U.size();
    Vector x = Vector(n);
    for (int i = n - 1; i >= 0; i--) {
        double summary = 0;
        for (int j = i + 1; j < n; j++) {
            summary += U[i][j] * x[j];
        }
        x[i] = (y[i] - summary) / U[i][i];
    }
    return x;
}

Vector E_Column(size_t n, size_t ind) {
    Vector e = Vector(n, 0);
    e[ind] = 1;
    return e;
}

Matrix Inverse(const Matrix& A, const Matrix& L, const Matrix& U, const Matrix& P) {
    size_t n = A.size();
    Matrix A_opp(n, Vector(n, 0));
    for (size_t i = 0; i < n; i++) {
        Vector e = E_Column(n, i);
        Vector Pe = Multiply_Vector(P, e);

        Vector y = Forward_Substitution(L, Pe);
        Vector x = Backward_Substitution(U, y);

        for (size_t j = 0; j < n; j++) {
            A_opp[j][i] = x[j];
        }
    }
    return A_opp;
}

Matrix Single_Matrix(size_t n) {
    Matrix E = Matrix(n, Vector(n, 0));
    for (int i = 0; i < n; i++) {
        E[i][i] = 1;
    }

    return E;
}

double Determinant(const Matrix& U, size_t count_perm) {
    size_t n = U.size();
    double det = 1;

    for (size_t i = 0; i < n; i++) {
        det *= U[i][i];
    }

    if (count_perm % 2 != 0) {
        det = -det;
    }

    return det;
}

void Print_Result(const Matrix& A, const Vector& b) {
    // 1.1
    //{{2, 5, 1}, {8, 3, 4}, {6, 7, 9}};

    Matrix P, L, U;

    size_t count_perm = luDecomposition(A, P, L, U);

    std::cout << "Матрица A:\n";
    Print_Matrix(A);

    std::cout << "Матрица перестановок P:\n";
    Print_Matrix(P);

    std::cout << "Нижняя треугольная матрица L:\n";
    Print_Matrix(L);

    std::cout << "Верхняя треугольная матрица U:\n";
    Print_Matrix(U);

    Matrix LU = Multiply_Matrix(L, U);
    Matrix PA = Multiply_Matrix(P, A);

    std::cout << "L * U:\n";
    Print_Matrix(LU);

    std::cout << "P * A:\n";
    Print_Matrix(PA);

    std::cout << "Проверка L * U = P * A: ";

    if (Equal_Matrix(LU, PA)) {
        std::cout << "PASSED\n";
    } else {
        std::cout << "FAILED\n";
    }

    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "Ax = b\n";
    std::cout << "\n\n";

    Vector Pb = Multiply_Vector(P, b);
    std::cout << "P * b:\n";
    Print_Vector(Pb);
    std::cout << std::endl;

    std::cout << "Прямой ход Ly = b:\n";
    std::cout << "Матрица y:\n";
    Vector y = Forward_Substitution(L, Pb);
    Print_Vector(y);
    std::cout << std::endl;

    std::cout << "Обратный ход Ux = y:\n";
    std::cout << "Матрица x:\n";
    Vector x = Backward_Substitution(U, y);
    Print_Vector(x);
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "A^(-1) = \n";
    // AA^(-1) = E
    // PAx(i) = Px(i)
    // LUx(i) = Pe(i)
    // Ux(i) = y(i)
    // Ly(i) = Pe(i)
    // A ^ (-1) = x

    Matrix A_opp = Inverse(A, L, U, P);
    Print_Matrix(A_opp);

    std::cout << "A * A ^ (-1):\n";
    Matrix mult = Multiply_Matrix(A, A_opp);
    Print_Matrix(mult);

    Matrix E = Single_Matrix(A.size());
    std::cout << "Проверка A * A ^ (-1) = E: ";
    if (Equal_Matrix(mult, E)) {
        std::cout << "PASSED\n";
    } else {
        std::cout << "FAILED\n";
    }

    std::cout << "det(A) = " << Determinant(U, count_perm) << std::endl;
}
