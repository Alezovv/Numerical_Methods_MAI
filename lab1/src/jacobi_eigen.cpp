#define _USE_MATH_DEFINES
#include "jacobi_eigen.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Вспомогательная функция для поиска максимального внедиагонального элемента
void FindMaxOffDiagonal(const Matrix &A, size_t &max_i, size_t &max_j, double &max_val)
{
    size_t n = A.size();
    max_val = 0.0;
    for (size_t i = 0; i < n - 1; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            if (std::abs(A[i][j]) > max_val)
            {
                max_val = std::abs(A[i][j]);
                max_i = i;
                max_j = j;
            }
        }
    }
}

void Solve_Task_1_4(const std::string &filename)
{
    std::cout << "--- Запуск Лабораторной 1.4 (Метод вращений Якоби) ---\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    size_t n;
    double eps;
    file >> n >> eps;

    Matrix A(n, Vector(n));
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            file >> A[i][j];
        }
    }
    file.close();

    Matrix A_orig = A; // Сохраняем исходную матрицу для проверки в конце

    // Матрица собственных векторов (изначально единичная)
    Matrix V(n, Vector(n, 0.0));
    for (size_t i = 0; i < n; ++i)
        V[i][i] = 1.0;

    size_t i_max = 0, j_max = 0;
    double max_val = 0.0;

    // Основной цикл метода Якоби
    while (true)
    {
        FindMaxOffDiagonal(A, i_max, j_max, max_val);

        // Условие выхода: максимальный элемент меньше погрешности
        if (max_val < eps)
            break;

        // Вычисление угла поворота phi
        double phi;
        if (std::abs(A[i_max][i_max] - A[j_max][j_max]) < 1e-9)
        {
            phi = M_PI / 4.0;
        }
        else
        {
            phi = 0.5 * std::atan(2.0 * A[i_max][j_max] / (A[i_max][i_max] - A[j_max][j_max]));
        }

        double c = std::cos(phi);
        double s = std::sin(phi);

        // Обновляем матрицу A (пересчитываем только изменившиеся элементы)
        for (size_t k = 0; k < n; ++k)
        {
            if (k != i_max && k != j_max)
            {
                double aki = A[k][i_max];
                double akj = A[k][j_max];
                A[k][i_max] = A[i_max][k] = c * aki + s * akj;
                A[k][j_max] = A[j_max][k] = -s * aki + c * akj;
            }
        }

        double aii = A[i_max][i_max];
        double ajj = A[j_max][j_max];
        double aij = A[i_max][j_max];

        A[i_max][i_max] = c * c * aii + 2.0 * s * c * aij + s * s * ajj;
        A[j_max][j_max] = s * s * aii - 2.0 * s * c * aij + c * c * ajj;
        A[i_max][j_max] = A[j_max][i_max] = 0.0; // Принудительно зануляем

        // Обновляем матрицу собственных векторов V
        for (size_t k = 0; k < n; ++k)
        {
            double vki = V[k][i_max];
            double vkj = V[k][j_max];
            V[k][i_max] = c * vki + s * vkj;
            V[k][j_max] = -s * vki + c * vkj;
        }
    }

    // ВЫВОД РЕЗУЛЬТАТОВ (По ТЗ)
    std::cout << "1. Заданная точность вычислений (eps): " << eps << "\n\n";

    std::cout << "2. Найденные собственные значения (Lambda):\n";
    Vector eigenvalues(n);
    for (size_t i = 0; i < n; ++i)
        eigenvalues[i] = A[i][i];

    for (double val : eigenvalues)
    {
        std::cout << std::fixed << std::setprecision(5) << val << "  ";
    }
    std::cout << "\n\n";

    std::cout << "3. Матрица собственных векторов V (столбцы - векторы):\n";
    for (const auto &row : V)
    {
        for (double val : row)
        {
            std::cout << std::setw(12) << std::fixed << std::setprecision(5) << val;
        }
        std::cout << '\n';
    }
    std::cout << '\n';

    // 4. Проверка A * V = V * Lambda
    Matrix Lambda(n, Vector(n, 0.0));
    for (size_t i = 0; i < n; ++i)
        Lambda[i][i] = eigenvalues[i];

    Matrix AV = Multiply_Matrix(A_orig, V);
    Matrix VLambda = Multiply_Matrix(V, Lambda);

    std::cout << "4. Проверка A * V = V * Lambda:\n";
    std::cout << "A * V:\n";
    for (const auto &row : AV)
    {
        for (double val : row)
        {
            std::cout << std::setw(12) << std::fixed << std::setprecision(5) << val;
        }
        std::cout << '\n';
    }
    std::cout << "\nV * Lambda:\n";
    for (const auto &row : VLambda)
    {
        for (double val : row)
        {
            std::cout << std::setw(12) << std::fixed << std::setprecision(5) << val;
        }
        std::cout << '\n';
    }

    // Собственная проверка на равенство с небольшим допуском
    bool passed = true;
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            // Используем eps из ввода, немного расширив коридор погрешности
            if (std::abs(AV[i][j] - VLambda[i][j]) > eps * 10.0)
                passed = false;
        }
    }

    if (passed)
    {
        std::cout << "\nРезультат проверки: PASSED\n";
    }
    else
    {
        std::cout << "\nРезультат проверки: FAILED (Слишком большая погрешность)\n";
    }
}