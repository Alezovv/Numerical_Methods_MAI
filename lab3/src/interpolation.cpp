#include "interpolation.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>

double exact_function(double x)
{
    if (x < 0)
        return 0.0; // Защита от корня из отрицательного числа при отрисовке левого края
    return std::sqrt(x);
}

// Форматирование чисел для красивого вывода многочленов
std::string fmt(double val)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(4) << val;
    return out.str();
}

// Вычисление значения полинома Лагранжа
double calc_lagrange(const std::vector<double> &x, const std::vector<double> &y, double x_star)
{
    double result = 0.0;
    int n = x.size();
    for (int i = 0; i < n; i++)
    {
        double term = y[i];
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                term *= (x_star - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }
    return result;
}

// Построение строки полинома Лагранжа
std::string build_lagrange_string(const std::vector<double> &x, const std::vector<double> &y)
{
    std::string poly = "L(x) = ";
    int n = x.size();
    for (int i = 0; i < n; i++)
    {
        if (i > 0 && y[i] >= 0)
            poly += " + ";
        poly += fmt(y[i]);

        double denominator = 1.0;
        std::string numerator = "";
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                numerator += "(x - " + fmt(x[j]) + ")";
                denominator *= (x[i] - x[j]);
            }
        }
        poly += " * [" + numerator + " / " + fmt(denominator) + "]";
    }
    return poly;
}

// Вычисление таблицы разделенных разностей для Ньютона
std::vector<std::vector<double>> calc_divided_differences(const std::vector<double> &x, const std::vector<double> &y)
{
    int n = x.size();
    std::vector<std::vector<double>> table(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
        table[i][0] = y[i];

    for (int j = 1; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
        }
    }
    return table;
}

// Вычисление значения полинома Ньютона
double calc_newton(const std::vector<double> &x, const std::vector<std::vector<double>> &table, double x_star)
{
    int n = x.size();
    double result = table[0][0];
    double term = 1.0;
    for (int i = 1; i < n; i++)
    {
        term *= (x_star - x[i - 1]);
        result += table[0][i] * term;
    }
    return result;
}

// Построение строки полинома Ньютона
std::string build_newton_string(const std::vector<double> &x, const std::vector<std::vector<double>> &table)
{
    std::string poly = "P(x) = " + fmt(table[0][0]);
    int n = x.size();
    std::string term = "";
    for (int i = 1; i < n; i++)
    {
        term += "(x - " + fmt(x[i - 1]) + ")";
        double coef = table[0][i];
        if (coef >= 0)
            poly += " + " + fmt(coef) + "*" + term;
        else
            poly += " - " + fmt(std::abs(coef)) + "*" + term;
    }
    return poly;
}

void Solve_Task_3_1(const std::string &filename)
{
    std::cout << "\n--- Задача 3.1: Интерполяция (Лагранж и Ньютон) ---\n";
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть " << filename << "\n";
        return;
    }

    int n;
    file >> n;
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
    {
        file >> x[i] >> y[i];
    }
    double x_star;
    file >> x_star;
    file.close();

    // 1. Исходные узлы
    std::cout << "1. Исходные узлы:\n";
    for (int i = 0; i < n; i++)
    {
        std::cout << "x[" << i << "] = " << std::setw(6) << x[i] << ", y[" << i << "] = " << y[i] << "\n";
    }

    // 2. Многочлен Лагранжа
    std::cout << "\n2. Многочлен Лагранжа:\n"
              << build_lagrange_string(x, y) << "\n";

    // 3. Таблица разделенных разностей
    auto diff_table = calc_divided_differences(x, y);
    std::cout << "\n3. Таблица разделенных разностей:\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n - i; j++)
        {
            std::cout << std::setw(10) << std::fixed << std::setprecision(4) << diff_table[i][j] << " ";
        }
        std::cout << "\n";
    }

    // 4. Многочлен Ньютона
    std::cout << "\n4. Многочлен Ньютона:\n"
              << build_newton_string(x, diff_table) << "\n";

    // 5-9. Расчеты в точке x*
    double lagrange_val = calc_lagrange(x, y, x_star);
    double newton_val = calc_newton(x, diff_table, x_star);
    double exact_val = exact_function(x_star);

    std::cout << "\n5. L(x*) = " << lagrange_val;
    std::cout << "\n6. P(x*) = " << newton_val;
    std::cout << "\n7. f(x*) = " << exact_val;
    std::cout << "\n8. Абсолютная погрешность:\n";
    std::cout << "   |L(x*) - f(x*)| = " << std::abs(lagrange_val - exact_val) << "\n";
    std::cout << "   |P(x*) - f(x*)| = " << std::abs(newton_val - exact_val) << "\n";
    std::cout << "9. Разность |L(x*) - P(x*)| = " << std::abs(lagrange_val - newton_val) << "\n";

    // Экспорт данных для Python-графика
    std::ofstream csv("data/lab3/plot_3_1.csv");
    csv << "x,f(x),lagrange,newton\n";
    double step = (x.back() - x.front()) / 100.0;
    for (double xi = x.front() - step * 10; xi <= x.back() + step * 10; xi += step)
    {
        csv << xi << "," << exact_function(xi) << ","
            << calc_lagrange(x, y, xi) << "," << calc_newton(x, diff_table, xi) << "\n";
    }
    csv.close();

    std::ofstream meta("data/lab3/meta_3_1.txt");
    for (int i = 0; i < n; ++i)
        meta << x[i] << " " << y[i] << "\n";
    meta << "STAR\n"
         << x_star << " " << lagrange_val << "\n";
    meta.close();

    std::cout << "\n[Данные для графика выгружены. График будет доступен после запуска Python-скрипта.]\n";
}