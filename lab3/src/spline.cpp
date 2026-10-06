#include "interpolation.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

struct SplineTuple
{
    double a, b, c, d, x;
};

void Solve_Task_3_2(const std::string &filename)
{
    std::cout << "\n--- Задача 3.2: Кубический сплайн ---\n";
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

    int intervals = n - 1;
    std::vector<double> h(n);
    for (int i = 1; i < n; i++)
        h[i] = x[i] - x[i - 1];

    // Инициализация СЛАУ для коэффициентов c_i
    std::vector<double> A(n - 1, 0.0), B(n - 1, 0.0), C(n - 1, 0.0), D(n - 1, 0.0);
    for (int i = 1; i < n - 1; i++)
    {
        A[i] = h[i];
        B[i] = 2.0 * (h[i] + h[i + 1]);
        C[i] = h[i + 1];
        D[i] = 3.0 * ((y[i + 1] - y[i]) / h[i + 1] - (y[i] - y[i - 1]) / h[i]);
    }

    // Метод прогонки
    std::vector<double> P(n - 1, 0.0), Q(n - 1, 0.0), c(n, 0.0);
    P[1] = -C[1] / B[1];
    Q[1] = D[1] / B[1];
    for (int i = 2; i < n - 1; i++)
    {
        double denom = B[i] + A[i] * P[i - 1];
        P[i] = -C[i] / denom;
        Q[i] = (D[i] - A[i] * Q[i - 1]) / denom;
    }
    for (int i = n - 2; i > 0; i--)
    {
        c[i] = P[i] * c[i + 1] + Q[i];
    }
    // Естественные граничные условия: c[0] = 0, c[n-1] = 0

    // Расчет a, b, d
    std::vector<SplineTuple> splines(n);
    for (int i = 1; i < n; i++)
    {
        splines[i].a = y[i];
        splines[i].c = c[i];
        splines[i].d = (c[i] - c[i - 1]) / (3.0 * h[i]);
        splines[i].b = (y[i] - y[i - 1]) / h[i] + (2.0 * c[i] + c[i - 1]) * h[i] / 3.0;
        splines[i].x = x[i];
    }

    // 1-3. Вывод коэффициентов и участков
    std::cout << "Коэффициенты сплайна на интервалах:\n";
    std::cout << std::setw(15) << "Интервал" << " | "
              << std::setw(10) << "a" << " | "
              << std::setw(10) << "b" << " | "
              << std::setw(10) << "c" << " | "
              << std::setw(10) << "d" << "\n";
    std::cout << std::string(65, '-') << "\n";

    for (int i = 1; i < n; i++)
    {
        std::cout << "[" << std::fixed << std::setprecision(2) << x[i - 1] << ", " << x[i] << "] | "
                  << std::setw(10) << std::setprecision(4) << splines[i].a << " | "
                  << std::setw(10) << splines[i].b << " | "
                  << std::setw(10) << splines[i].c << " | "
                  << std::setw(10) << splines[i].d << "\n";
    }

    // 4-5. Поиск интервала и значения в x*
    int target_interval = 1;
    for (int i = 1; i < n; i++)
    {
        if (x_star >= x[i - 1] && x_star <= x[i])
        {
            target_interval = i;
            break;
        }
    }
    double dx = x_star - splines[target_interval].x;
    double s_star = splines[target_interval].a + splines[target_interval].b * dx +
                    splines[target_interval].c * dx * dx + splines[target_interval].d * dx * dx * dx;

    std::cout << "\nТочка x* = " << x_star << " попадает в интервал " << target_interval << "\n";
    std::cout << "S(x*) = " << s_star << "\n";

    // 6. Проверка прохождения через узлы
    std::cout << "\nПроверка узлов (S(x_i) == y_i):\n";
    for (int i = 1; i < n; i++)
    {
        double dx_node = x[i] - splines[i].x; // = 0
        double s_node = splines[i].a + splines[i].b * dx_node + splines[i].c * dx_node * dx_node + splines[i].d * dx_node * dx_node * dx_node;
        std::cout << "Узел " << i << ": S(" << x[i] << ") = " << s_node << " (Факт: " << y[i] << ")\n";
    }

    // 7. Проверка второй производной на концах S''(x) = 2c + 6d(x - x_i)
    double S_prime2_start = 2.0 * splines[1].c + 6.0 * splines[1].d * (x[0] - splines[1].x);
    double S_prime2_end = 2.0 * splines[n - 1].c + 6.0 * splines[n - 1].d * (x[n - 1] - splines[n - 1].x); // dx = 0, остается 2c

    std::cout << "\nГраничные условия (вторые производные на концах):\n";
    std::cout << "S''(x_0) = " << S_prime2_start << " (Ожидается 0)\n";
    std::cout << "S''(x_n) = " << S_prime2_end << " (Ожидается 0)\n";

    // Экспорт для графика
    std::ofstream csv("data/lab3/plot_3_2.csv");
    csv << "x,S\n";
    for (int i = 1; i < n; i++)
    {
        double step = (x[i] - x[i - 1]) / 20.0;
        for (double xi = x[i - 1]; xi <= x[i]; xi += step)
        {
            double dxi = xi - splines[i].x;
            double S_val = splines[i].a + splines[i].b * dxi + splines[i].c * dxi * dxi + splines[i].d * dxi * dxi * dxi;
            csv << xi << "," << S_val << "\n";
        }
    }
    csv.close();

    std::ofstream meta("data/lab3/meta_3_2.txt");
    for (int i = 0; i < n; ++i)
        meta << x[i] << " " << y[i] << "\n";
    meta << "STAR\n"
         << x_star << " " << s_star << "\n";
    meta.close();
    std::cout << "\n[Данные для графика выгружены. Запустите plotter_3_2.py]\n";
}