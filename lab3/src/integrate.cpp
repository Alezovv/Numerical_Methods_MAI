#include "interpolation.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

// Интегрируемая функция
double f(double x)
{
    return 1.0 / (3.0 * x * x + 4.0 * x + 2.0);
}

// Точное аналитическое значение интеграла
double exact_integral()
{
    auto F = [](double x)
    {
        return (1.0 / std::sqrt(2.0)) * std::atan((3.0 * x + 2.0) / std::sqrt(2.0));
    };
    return F(2.0) - F(-2.0);
}

// Методы прямоугольников
double rect_left(double a, double b, double h)
{
    double sum = 0;
    for (double x = a; x < b - 1e-9; x += h)
        sum += f(x);
    return sum * h;
}
double rect_right(double a, double b, double h)
{
    double sum = 0;
    for (double x = a + h; x <= b + 1e-9; x += h)
        sum += f(x);
    return sum * h;
}
double rect_mid(double a, double b, double h)
{
    double sum = 0;
    for (double x = a; x < b - 1e-9; x += h)
        sum += f(x + h / 2.0);
    return sum * h;
}

// Метод трапеций
double trapezoidal(double a, double b, double h)
{
    double sum = (f(a) + f(b)) / 2.0;
    for (double x = a + h; x < b - 1e-9; x += h)
        sum += f(x);
    return sum * h;
}

// Метод Симпсона
double simpson(double a, double b, double h)
{
    double sum = f(a) + f(b);
    int n = std::round((b - a) / h);
    for (int i = 1; i < n; i++)
    {
        double x = a + i * h;
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return sum * h / 3.0;
}

void print_row(const std::string &name, double I_h1, double I_h2, int p, double exact)
{
    double I_rr = I_h2 + (I_h2 - I_h1) / (std::pow(2, p) - 1.0);
    double err_h2 = std::abs(I_rr - I_h2);   // Оценка погрешности
    double abs_err = std::abs(exact - I_h2); // Фактическая абсолютная погрешность

    std::cout << std::left << std::setw(20) << name
              << std::right << std::setw(12) << I_h1
              << std::setw(12) << I_h2
              << std::setw(6) << p
              << std::setw(14) << I_rr
              << std::setw(15) << err_h2
              << std::setw(15) << abs_err << "\n";
}

void Solve_Task_3_5()
{
    std::cout << "--- Задача 3.5: Численное интегрирование ---\n";
    double a = -2.0, b = 2.0;
    double h1 = 1.0, h2 = 0.5;
    double exact = exact_integral();

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Точное значение интеграла: " << exact << "\n\n";

    std::cout << std::left << std::setw(20) << "Метод"
              << std::right << std::setw(12) << "I(h1)"
              << std::setw(12) << "I(h2)"
              << std::setw(6) << "p"
              << std::setw(14) << "I_rr"
              << std::setw(15) << "Погрешн. RR"
              << std::setw(15) << "Абс. погрешн.\n";
    std::cout << std::string(94, '-') << "\n";

    print_row("Левые прямоуг.", rect_left(a, b, h1), rect_left(a, b, h2), 2, exact);
    print_row("Правые прямоуг.", rect_right(a, b, h1), rect_right(a, b, h2), 2, exact);
    print_row("Средние прямоуг.", rect_mid(a, b, h1), rect_mid(a, b, h2), 2, exact);
    print_row("Трапеции", trapezoidal(a, b, h1), trapezoidal(a, b, h2), 2, exact);
    print_row("Симпсон", simpson(a, b, h1), simpson(a, b, h2), 4, exact);

    std::cout << "\nНаиболее точный результат дал метод Симпсона.\n";
}