#include <iostream>
#include <string>
#include <cstdlib> // Для вызова Python скрипта через system()

// Подключаем все наши модули
#include "lu_decomposition.h"
#include "tridiagonal_sweep.h"
#include "iterative_methods.h"
#include "jacobi_eigen.h"
#include "qr_decomposition.h"

int main()
{
    const std::string PATH_PREFIX = "data/lab1/";

    int choice = -1;

    while (choice != 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "      Лабораторная работа №1\n";
        std::cout << "========================================\n";
        std::cout << "1. LU-разложение (1.1)\n";
        std::cout << "2. Метод прогонки (1.2)\n";
        std::cout << "3. Метод простых итераций и Зейделя (1.3)\n";
        std::cout << "4. Метод вращений Якоби (1.4)\n";
        std::cout << "5. QR-разложение (1.5)\n";
        std::cout << "6. Проверка комплексных корней (Python скрипт для 1.5)\n";
        std::cout << "0. Выход\n";
        std::cout << "========================================\n";
        std::cout << "Выберите задание (0-6): ";

        std::cin >> choice;

        std::cout << "\n";

        switch (choice)
        {
        case 1:
            Solve_Task_1_1(PATH_PREFIX + "input_1_1.txt");
            break;
        case 2:
            Solve_Task_1_2(PATH_PREFIX + "input_1_2.txt");
            break;
        case 3:
            Solve_Task_1_3(PATH_PREFIX + "input_1_3.txt");
            break;
        case 4:
            Solve_Task_1_4(PATH_PREFIX + "input_1_4.txt");
            break;
        case 5:
            Solve_Task_1_5(PATH_PREFIX + "input_1_5.txt");
            break;
        case 6:
            std::cout << "Запуск скрипта check_complex.py...\n\n";
            // В зависимости от системы, команда может быть "python", "python3" или "py"
            std::system("python lab1/check_complex.py");
            break;
        case 0:
            std::cout << "Завершение работы программы.\n";
            break;
        default:
            std::cout << "Неверный ввод. Пожалуйста, выберите число от 0 до 6.\n";
            break;
        }
    }

    return 0;
}