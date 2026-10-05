#include "tridiagonal_sweep.h"
#include <iostream>
#include <fstream>
#include <vector>

void Solve_Task_1_2(const std::string &filename)
{
    std::cout << "--- Запуск Лабораторной 1.2 (Метод прогонки) ---\n";

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << "\n";
        return;
    }

    size_t n;
    file >> n;

    // Векторы для трех диагоналей и вектора свободных членов
    // a - нижняя диагональ (индексы 1..n-1, a[0] = 0)
    // b - главная диагональ (индексы 0..n-1)
    // c - верхняя диагональ (индексы 0..n-2, c[n-1] = 0)
    Vector a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);

    // Считываем нижнюю диагональ (n - 1 элементов)
    for (size_t i = 1; i < n; ++i)
        file >> a[i];
    // Считываем главную диагональ (n элементов)
    for (size_t i = 0; i < n; ++i)
        file >> b[i];
    // Считываем верхнюю диагональ (n - 1 элементов)
    for (size_t i = 0; i < n - 1; ++i)
        file >> c[i];
    // Считываем вектор свободных членов (n элементов)
    for (size_t i = 0; i < n; ++i)
        file >> d[i];

    file.close();

    // Векторы для прогоночных коэффициентов
    Vector P(n, 0.0);
    Vector Q(n, 0.0);

    // Прямой ход
    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];

    for (size_t i = 1; i < n; ++i)
    {
        double denominator = b[i] + a[i] * P[i - 1];

        if (i < n - 1)
        {
            P[i] = -c[i] / denominator;
        }
        else
        {
            P[i] = 0.0; // Для последней строки P_n не вычисляется (или равен 0)
        }

        Q[i] = (d[i] - a[i] * Q[i - 1]) / denominator;
    }

    // Пункт 1: Все прогоночные коэффициенты
    std::cout << "Прогоночные коэффициенты P:\n";
    Print_Vector(P);

    std::cout << "Прогоночные коэффициенты Q:\n";
    Print_Vector(Q);

    // Обратный ход
    Vector x(n, 0.0);
    x[n - 1] = Q[n - 1];

    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = P[i] * x[i + 1] + Q[i];
    }

    // Пункт 2: Итоговое решение
    std::cout << "Итоговое решение системы уравнений x:\n";
    Print_Vector(x);
}