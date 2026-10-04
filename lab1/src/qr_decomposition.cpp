#include "qr_decomposition.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Функция для выполнения одного шага QR-разложения (Метод Грама-Шмидта)
void Perform_QR(const Matrix &A, Matrix &Q, Matrix &R)
{
    size_t n = A.size();
    Q = Matrix(n, Vector(n, 0.0));
    R = Matrix(n, Vector(n, 0.0));

    for (size_t j = 0; j < n; ++j)
    {
        Vector v(n);
        for (size_t i = 0; i < n; ++i)
            v[i] = A[i][j];

        for (size_t i = 0; i < j; ++i)
        {
            double sum = 0.0;
            for (size_t k = 0; k < n; ++k)
                sum += Q[k][i] * A[k][j];
            R[i][j] = sum;

            for (size_t k = 0; k < n; ++k)
                v[k] -= R[i][j] * Q[k][i];
        }

        double norm = 0.0;
        for (size_t k = 0; k < n; ++k)
            norm += v[k] * v[k];
        R[j][j] = std::sqrt(norm);

        if (R[j][j] > 1e-9)
        {
            for (size_t k = 0; k < n; ++k)
                Q[k][j] = v[k] / R[j][j];
        }
    }
}

void Solve_Task_1_5(const std::string &filename)
{
    std::cout << "--- Запуск Лабораторной 1.5 (QR-разложение) ---\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    size_t n;
    file >> n;
    Matrix A(n, Vector(n));
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            file >> A[i][j];
        }
    }
    file.close();

    // 1. Исходная матрица A
    std::cout << "1. Исходная матрица A:\n";
    for (const auto &row : A)
    {
        for (double val : row)
            std::cout << std::setw(12) << std::fixed << std::setprecision(5) << val;
        std::cout << '\n';
    }
    std::cout << '\n';

    // 2. Матрицы Q и R (для первой итерации)
    Matrix Q, R;
    Perform_QR(A, Q, R);

    std::cout << "2. Матрица Q (ортогональная):\n";
    for (const auto &row : Q)
    {
        for (double val : row)
            std::cout << std::setw(12) << std::fixed << std::setprecision(5) << val;
        std::cout << '\n';
    }
    std::cout << "\nМатрица R (верхнетреугольная):\n";
    for (const auto &row : R)
    {
        for (double val : row)
            std::cout << std::setw(12) << std::fixed << std::setprecision(5) << val;
        std::cout << '\n';
    }
    std::cout << '\n';

    // 3. Результат умножения Q * R
    Matrix QR_mult = Multiply_Matrix(Q, R);
    std::cout << "3. Проверка умножения Q * R (должно равняться A):\n";
    for (const auto &row : QR_mult)
    {
        for (double val : row)
            std::cout << std::setw(12) << std::fixed << std::setprecision(5) << val;
        std::cout << '\n';
    }
    std::cout << '\n';

    // 4. Поиск собственных значений (Итерационный процесс)
    Matrix Ak = A;
    double eps = 1e-6;
    bool converged = false;
    int max_iters = 1000;

    for (int iter = 0; iter < max_iters; ++iter)
    {
        Matrix Qk, Rk;
        Perform_QR(Ak, Qk, Rk);
        Ak = Multiply_Matrix(Rk, Qk); // A_{k+1} = R_k * Q_k

        // Проверка сходимости (поддиагональные элементы стремятся к нулю)
        double max_subdiag = 0.0;
        for (size_t i = 1; i < n; ++i)
        {
            if (std::abs(Ak[i][i - 1]) > max_subdiag)
            {
                max_subdiag = std::abs(Ak[i][i - 1]);
            }
        }

        if (max_subdiag < eps)
        {
            converged = true;
            break;
        }
    }

    std::cout << "4. Найденные собственные значения (вещественные):\n";
    if (converged)
    {
        for (size_t i = 0; i < n; ++i)
        {
            std::cout << "Lambda_" << i + 1 << " = " << std::fixed << std::setprecision(5) << Ak[i][i] << "\n";
        }
    }
    else
    {
        std::cout << "ВНИМАНИЕ: Метод не сошелся. Возможно, матрица имеет комплексные собственные значения (см. Python скрипт).\n";
    }
}