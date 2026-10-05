#include <iostream>
#include <string>
#include "nonlinear_equations.h"

int main()
{
    const std::string PATH_PREFIX = "lab2/data/";
    int choice = -1;

    while (choice != 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "      Лабораторная работа №2\n";
        std::cout << "========================================\n";
        std::cout << "1. Решение нелинейных уравнений (2.1)\n";
        std::cout << "2. Решение систем нелинейных уравнений (2.2)\n";
        std::cout << "0. Выход\n";
        std::cout << "========================================\n";
        std::cout << "Выберите задание (0-2): ";

        std::cin >> choice;
        std::cout << "\n";

        switch (choice)
        {
        case 1:
            Solve_Task_2_1(PATH_PREFIX + "input_2_1.txt");
            break;
        case 2:
            Solve_Task_2_2(PATH_PREFIX + "input_2_2.txt");
            break;
        case 0:
            std::cout << "Завершение работы программы.\n";
            break;
        default:
            std::cout << "Неверный ввод.\n";
            break;
        }
    }
    return 0;
}