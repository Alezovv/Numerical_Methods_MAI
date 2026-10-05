#include "iterative_methods.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

// Вспомогательная функция для оценки нормы вектора (максимальный по модулю элемент)
double Vector_Norm(const Vector &v1, const Vector &v2)
{
    double max_diff = 0.0;
    for (size_t i = 0; i < v1.size(); ++i)
    {
        max_diff = std::max(max_diff, std::abs(v1[i] - v2[i]));
    }
    return max_diff;
}

// Эталонное решение методом Гаусса
Vector Gauss_Exact_Solution(Matrix A, Vector b)
{
    size_t n = A.size();
    for (size_t i = 0; i < n; ++i)
    {
        // Поиск главного элемента
        size_t max_row = i;
        for (size_t k = i + 1; k < n; ++k)
        {
            if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
                max_row = k;
        }
        std::swap(A[i], A[max_row]);
        std::swap(b[i], b[max_row]);

        for (size_t k = i + 1; k < n; ++k)
        {
            double factor = A[k][i] / A[i][i];
            for (size_t j = i; j < n; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }
    Vector x(n, 0.0);
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = 0.0;
        for (size_t j = i + 1; j < n; ++j)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }
    return x;
}

void Solve_Task_1_3(const std::string &filename)
{
    std::cout << "--- Запуск Лабораторной 1.3 (Метод простых итераций и Зейделя) ---\n";

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
    Vector b(n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            file >> A[i][j];
        }
    }
    for (size_t i = 0; i < n; ++i)
    {
        file >> b[i];
    }
    file.close();
    
    std::cout << "1. Заданная точность вычислений (eps): " << eps << "\n\n";

    // Приведение к эквивалентному виду x = alpha * x + beta
    Matrix alpha(n, Vector(n, 0.0));
    Vector beta(n, 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        beta[i] = b[i] / A[i][i];
        for (size_t j = 0; j < n; ++j)
        {
            if (i != j)
            {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
    }

    // --- Метод простых итераций ---
    Vector x_simple = beta; // Начальное приближение
    Vector x_simple_new(n, 0.0);
    int iters_simple = 0;
    double current_eps = eps + 1; // Чтобы зайти в цикл

    while (current_eps > eps && iters_simple < 10000)
    {
        x_simple_new = Multiply_Vector(alpha, x_simple);
        for (size_t i = 0; i < n; i++)
            x_simple_new[i] += beta[i];

        current_eps = Vector_Norm(x_simple_new, x_simple);
        x_simple = x_simple_new;
        iters_simple++;
    }

    // --- Метод Зейделя ---
    Vector x_seidel = beta; // Начальное приближение
    Vector x_seidel_new = x_seidel;
    int iters_seidel = 0;
    current_eps = eps + 1;

    while (current_eps > eps && iters_seidel < 10000)
    {
        for (size_t i = 0; i < n; ++i)
        {
            double sum = beta[i];
            // Используем уже обновленные значения (x_seidel_new)
            for (size_t j = 0; j < i; ++j)
                sum += alpha[i][j] * x_seidel_new[j];
            // Используем старые значения (x_seidel)
            for (size_t j = i; j < n; ++j)
                sum += alpha[i][j] * x_seidel[j];

            x_seidel_new[i] = sum;
        }
        current_eps = Vector_Norm(x_seidel_new, x_seidel);
        x_seidel = x_seidel_new;
        iters_seidel++;
    }

    // Эталонное решение
    Vector exact_x = Gauss_Exact_Solution(A, b);
    
    std::cout << "--- Решение методом простых итераций ---\n";
    std::cout << "Вектор решения x:\n";
    Print_Vector(x_simple);
    std::cout << "Количество итераций: " << iters_simple << "\n\n";

    std::cout << "--- Решение методом Зейделя ---\n";
    std::cout << "Вектор решения x:\n";
    Print_Vector(x_seidel);
    std::cout << "Количество итераций: " << iters_seidel << "\n\n";

    std::cout << "--- Эталонное решение (Гаусс) для проверки ---\n";
    std::cout << "Вектор точного решения:\n";
    Print_Vector(exact_x);
    std::cout << "\n";

    // 4. Сравнение методов по скорости сходимости
    std::cout << "--- Сравнение скорости сходимости ---\n";
    if (iters_seidel < iters_simple)
    {
        std::cout << "Метод Зейделя сошелся быстрее на " << (iters_simple - iters_seidel) << " итераций.\n";
        std::cout << "Это типично, так как он использует уточненные значения корней в рамках одной итерации.\n";
    }
    else if (iters_seidel > iters_simple)
    {
        std::cout << "Метод простых итераций сошелся быстрее (что бывает редко).\n";
    }
    else
    {
        std::cout << "Оба метода сошлись за одинаковое количество итераций.\n";
    }
}