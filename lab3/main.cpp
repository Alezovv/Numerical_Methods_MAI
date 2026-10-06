#include <iostream>
#include <string>
#include <windows.h>
#include "interpolation.h"

int main()
{
    SetConsoleOutputCP(CP_UTF8);

    const std::string PATH_PREFIX = "data/lab3/";
    int choice = -1;

    while (choice != 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "      Лабораторная работа №3\n";
        std::cout << "========================================\n";
        std::cout << "1. Интерполяция Лагранжа и Ньютона (3.1)\n";
        std::cout << "2. Кубический сплайн (3.2)\n";
        std::cout << "3. Метод наименьших квадратов (3.3)\n";
        std::cout << "4. Численное дифференцирование (3.4)\n";
        std::cout << "5. Численное интегрирование (3.5)\n";
        std::cout << "0. Выход\n";
        std::cout << "========================================\n";
        std::cout << "Выберите задание (0-5): ";

        std::cin >> choice;
        std::cout << "\n";

        switch (choice)
        {
        case 1:
            Solve_Task_3_1(PATH_PREFIX + "input_3_1.txt");
            break;
        case 2:
            Solve_Task_3_2(PATH_PREFIX + "input_3_2.txt");
            break;
        case 3:
            Solve_Task_3_3(PATH_PREFIX + "input_3_3.txt");
            break;
        case 4:
            Solve_Task_3_4(PATH_PREFIX + "input_3_4.txt");
            break;
        case 5:
            Solve_Task_3_5();
            break;
        case 0:
            std::cout << "Завершение работы.\n";
            break;
        default:
            std::cout << "Неверный ввод.\n";
            break;
        }
    }
    return 0;
}