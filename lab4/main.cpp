#include <iostream>
#include <string>
#include "ode_solvers.h"

#ifdef _WIN32
#include <windows.h>
#endif

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    const std::string PATH_PREFIX = "data/lab4/";
    int choice = -1;

    while (choice != 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "      Лабораторная работа №4\n";
        std::cout << "  Численные методы решения ОДУ (Вариант 7)\n";
        std::cout << "========================================\n";
        std::cout << "1. Решение задачи Коши (4.1)\n";
        std::cout << "2. Решение краевой задачи (4.2)\n";
        std::cout << "0. Выход\n";
        std::cout << "========================================\n";
        std::cout << "Выберите задание (0-2): ";

        if (!(std::cin >> choice))
        {
            break;
        }
        std::cout << "\n";

        switch (choice)
        {
        case 1:
            Solve_Task_4_1(PATH_PREFIX + "input_4_1.txt");
            break;
        case 2:
            Solve_Task_4_2(PATH_PREFIX + "input_4_2.txt");
            break;
        case 0:
            std::cout << "Завершение работы программы.\n";
            break;
        default:
            std::cout << "Неверный ввод. Пожалуйста, выберите 0, 1 или 2.\n";
            break;
        }
    }
    return 0;
}
