#include "interpolation.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// Локальная функция для решения СЛАУ методом Гаусса
std::vector<double> solve_gauss(std::vector<std::vector<double>> A, std::vector<double> B)
{
    int n = A.size();
    for (int i = 0; i < n; ++i)
    {
        int max_row = i;
        for (int k = i + 1; k < n; ++k)
            if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
                max_row = k;
        std::swap(A[i], A[max_row]);
        std::swap(B[i], B[max_row]);

        for (int k = i + 1; k < n; ++k)
        {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j)
                A[k][j] -= factor * A[i][j];
            B[k] -= factor * B[i];
        }
    }
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = B[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

void Solve_Task_3_3(const std::string &filename)
{
    std::cout << "--- Задача 3.3: Метод наименьших квадратов ---\n";
    std::ifstream file(filename);
    if (!file.is_open())
        return;
    int n;
    file >> n;
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
        file >> x[i] >> y[i];
    file.close();

    std::vector<std::vector<double>> coeffs(3);
    std::vector<double> mse(3, 0.0);

    for (int m = 1; m <= 3; ++m)
    { // степени многочлена 1, 2, 3
        int dim = m + 1;
        std::vector<std::vector<double>> A(dim, std::vector<double>(dim, 0.0));
        std::vector<double> B(dim, 0.0);

        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                for (int k = 0; k < n; ++k)
                    A[i][j] += std::pow(x[k], i + j);
            }
            for (int k = 0; k < n; ++k)
                B[i] += y[k] * std::pow(x[k], i);
        }

        coeffs[m - 1] = solve_gauss(A, B);

        // Расчет MSE
        double sum_err = 0.0;
        for (int k = 0; k < n; ++k)
        {
            double P = 0.0;
            for (int i = 0; i <= m; ++i)
                P += coeffs[m - 1][i] * std::pow(x[k], i);
            sum_err += std::pow(P - y[k], 2);
        }
        mse[m - 1] = sum_err / n;

        // Вывод результатов
        std::cout << "\nМногочлен " << m << "-й степени:\n";
        std::cout << "P" << m << "(x) = ";
        for (int i = 0; i <= m; ++i)
        {
            if (i > 0 && coeffs[m - 1][i] >= 0)
                std::cout << " + ";
            std::cout << std::fixed << std::setprecision(4) << coeffs[m - 1][i];
            if (i > 0)
                std::cout << "x^" << i;
        }
        std::cout << "\nMSE" << m << " = " << std::scientific << mse[m - 1] << std::fixed << "\n";
    }

    int best_deg = 1;
    if (mse[1] < mse[0])
        best_deg = 2;
    if (mse[2] < mse[best_deg - 1])
        best_deg = 3;
    std::cout << "\nМногочлен, лучше всего приближающий данные: P" << best_deg << "(x)\n";

    // Экспорт для Python
    std::ofstream csv("data/lab3/plot_3_3.csv");
    csv << "x,P1,P2,P3\n";
    double step = (x.back() - x.front()) / 50.0;
    for (double xi = x.front(); xi <= x.back() + 1e-9; xi += step)
    {
        csv << xi;
        for (int m = 1; m <= 3; ++m)
        {
            double P = 0.0;
            for (int i = 0; i <= m; ++i)
                P += coeffs[m - 1][i] * std::pow(xi, i);
            csv << "," << P;
        }
        csv << "\n";
    }
    csv.close();

    std::ofstream meta("data/lab3/meta_3_3.txt");
    for (int i = 0; i < n; ++i)
        meta << x[i] << " " << y[i] << "\n";
    meta.close();
    std::cout << "\n[Данные для графика выгружены. Запустите plotter_3_3.py]\n";
}