#include "interpolation.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

void Solve_Task_3_4(const std::string &filename)
{
    std::cout << "--- Задача 3.4: Численное дифференцирование ---\n";
    std::ifstream file(filename);
    if (!file.is_open())
        return;

    int n;
    file >> n;
    std::vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
        file >> x[i] >> y[i];
    double x_star;
    file >> x_star;
    file.close();

    // 1-2. Исходные данные
    std::cout << "Исходная таблица:\n";
    for (int i = 0; i < n; ++i)
        std::cout << "x[" << i << "]=" << std::setw(5) << x[i] << ", y[" << i << "]=" << y[i] << "\n";
    std::cout << "\nТочка для вычисления X* = " << x_star << "\n";

    // Поиск индекса точки X*
    int idx = -1;
    for (int i = 0; i < n; i++)
    {
        if (std::abs(x[i] - x_star) < 1e-6)
        {
            idx = i;
            break;
        }
    }

    if (idx <= 0 || idx >= n - 1)
    {
        std::cout << "Ошибка: точка X* должна находиться внутри таблицы (не на краях) для центральных разностей.\n";
        return;
    }

    double h = x[idx] - x[idx - 1]; // предполагается равномерная сетка

    // 3. Используемые узлы
    std::cout << "Для вычислений использованы узлы: x[" << idx - 1 << "], x[" << idx << "], x[" << idx + 1 << "]\n";

    // 4-5. Производные (формулы центральных разностей)
    double dy = (y[idx + 1] - y[idx - 1]) / (2.0 * h);
    double d2y = (y[idx + 1] - 2.0 * y[idx] + y[idx - 1]) / (h * h);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nПервая производная y'(X*) = " << dy << "\n";
    std::cout << "Вторая производная y''(X*) = " << d2y << "\n";
}