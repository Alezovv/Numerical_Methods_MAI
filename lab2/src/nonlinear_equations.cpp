#include "nonlinear_equations.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

// ФУНКЦИИ ДЛЯ 2.1
double f_1(double x)
{
    return std::pow(2.0, x) + x * x - 2.0;
}

double df_1(double x)
{
    return std::pow(2.0, x) * std::log(2.0) + 2.0 * x;
}

double phi_1(double x)
{
    return std::sqrt(2.0 - std::pow(2.0, x));
}

// ФУНКЦИИ ДЛЯ 2.2
double sys_f1(double x1, double x2)
{
    return x1 * x1 - 2.0 * std::log10(x2) - 1.0;
}

double sys_f2(double x1, double x2)
{
    return x1 * x1 - 2.0 * x1 * x2 + 2.0;
}

// Простые итерации для системы
double sys_phi1(double x2)
{
    return std::sqrt(2.0 * std::log10(x2) + 1.0);
}

double sys_phi2(double x1)
{
    return (x1 * x1 + 2.0) / (2.0 * x1);
}

void get_jacobian(double x1, double x2, double J[2][2])
{
    J[0][0] = 2.0 * x1;
    J[0][1] = -2.0 / (x2 * std::log(10.0));
    J[1][0] = 2.0 * x1 - 2.0 * x2;
    J[1][1] = -2.0 * x1;
}

void Solve_Task_2_1(const std::string &filename)
{
    std::cout << "\n--- Задача 2.1: Нелинейное уравнение ---\n";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    double eps, x0;
    file >> eps >> x0;
    file.close();

    std::cout << "Точность: " << eps << ", Начальное приближение: " << x0 << "\n\n";

    std::cout << ">>> МЕТОД ПРОСТОЙ ИТЕРАЦИИ <<<\n";
    std::cout << std::setw(5) << "k" << " | " << std::setw(12) << "x_k" << " | " << std::setw(12) << "f(x_k)" << " | " << std::setw(12) << "|x_k - x_{k-1}|\n";
    std::cout << std::string(50, '-') << "\n";

    double x_curr = x0, x_prev = x0;
    int iter = 0;
    double diff = eps + 1.0;

    while (diff > eps && iter < 1000)
    {
        x_prev = x_curr;
        x_curr = phi_1(x_prev);
        diff = std::abs(x_curr - x_prev);
        iter++;
        std::cout << std::setw(5) << iter << " | " << std::setw(12) << std::fixed << std::setprecision(6) << x_curr << " | " << std::setw(12) << f_1(x_curr) << " | " << std::setw(12) << diff << "\n";
    }
    std::cout << "Корень: " << x_curr << " (Итераций: " << iter << ")\n\n";

    std::cout << ">>> МЕТОД НЬЮТОНА <<<\n";
    std::cout << std::setw(5) << "k" << " | " << std::setw(12) << "x_k" << " | " << std::setw(12) << "f(x_k)" << " | " << std::setw(12) << "|x_k - x_{k-1}|\n";
    std::cout << std::string(50, '-') << "\n";

    x_curr = x0;
    iter = 0;
    diff = eps + 1.0;

    while (diff > eps && iter < 1000)
    {
        x_prev = x_curr;
        x_curr = x_prev - f_1(x_prev) / df_1(x_prev);
        diff = std::abs(x_curr - x_prev);
        iter++;
        std::cout << std::setw(5) << iter << " | " << std::setw(12) << std::fixed << std::setprecision(6) << x_curr << " | " << std::setw(12) << f_1(x_curr) << " | " << std::setw(12) << diff << "\n";
    }
    std::cout << "Корень: " << x_curr << " (Итераций: " << iter << ")\n\n";
}

void Solve_Task_2_2(const std::string &filename)
{
    std::cout << "\n--- Задача 2.2: Система нелинейных уравнений ---\n";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    double eps, x1_0, x2_0;
    file >> eps >> x1_0 >> x2_0;
    file.close();

    std::cout << "Точность: " << eps << "\nНачальное приближение: (" << x1_0 << ", " << x2_0 << ")\n\n";

    std::cout << ">>> МЕТОД ПРОСТОЙ ИТЕРАЦИИ <<<\n";
    std::cout << std::setw(5) << "k" << " | " << std::setw(12) << "x1_k" << " | " << std::setw(12) << "x2_k" << " | " << std::setw(12) << "max|diff|\n";
    std::cout << std::string(50, '-') << "\n";

    double x1_c = x1_0, x2_c = x2_0;
    double x1_p = x1_0, x2_p = x2_0;
    int iter = 0;
    double diff = eps + 1.0;

    while (diff > eps && iter < 1000)
    {
        x1_p = x1_c;
        x2_p = x2_c;

        x1_c = sys_phi1(x2_p);
        x2_c = sys_phi2(x1_p);

        diff = std::max(std::abs(x1_c - x1_p), std::abs(x2_c - x2_p));
        iter++;
        std::cout << std::setw(5) << iter << " | " << std::setw(12) << std::fixed << std::setprecision(6) << x1_c << " | " << std::setw(12) << x2_c << " | " << std::setw(12) << diff << "\n";
    }
    std::cout << "Решение: x1 = " << x1_c << ", x2 = " << x2_c << " (Итераций: " << iter << ")\n\n";

    std::cout << ">>> МЕТОД НЬЮТОНА <<<\n";
    std::cout << std::setw(5) << "k" << " | " << std::setw(12) << "x1_k" << " | " << std::setw(12) << "x2_k" << " | " << std::setw(12) << "max|diff|\n";
    std::cout << std::string(50, '-') << "\n";

    x1_c = x1_0;
    x2_c = x2_0;
    iter = 0;
    diff = eps + 1.0;

    while (diff > eps && iter < 1000)
    {
        x1_p = x1_c;
        x2_p = x2_c;

        double J[2][2];
        get_jacobian(x1_p, x2_p, J);

        double f1 = sys_f1(x1_p, x2_p);
        double f2 = sys_f2(x1_p, x2_p);

        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        // Решение системы J * dx = -F через обратную матрицу
        double dx1 = (-f1 * J[1][1] - (-f2) * J[0][1]) / det;
        double dx2 = (J[0][0] * (-f2) - J[1][0] * (-f1)) / det;

        x1_c = x1_p + dx1;
        x2_c = x2_p + dx2;

        diff = std::max(std::abs(dx1), std::abs(dx2));
        iter++;

        std::cout << std::setw(5) << iter << " | " << std::setw(12) << std::fixed << std::setprecision(6) << x1_c << " | " << std::setw(12) << x2_c << " | " << std::setw(12) << diff << "\n";
    }
    std::cout << "Решение: x1 = " << x1_c << ", x2 = " << x2_c << " (Итераций: " << iter << ")\n\n";
}