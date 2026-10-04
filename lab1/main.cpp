#include <iostream>
#include "file_io.h"
#include "lu_decomposition.h"
#include "matrix.h"

int main() {
    Matrix A;
    Vector b;

    if (Read_System_From_File("input_1_1.txt", A, b)) {
        Print_Result(A, b);
    } else {
        std::cerr << "Ошибка: не удалось открыть файл input_1_1.txt\n";
    }

    return 0;
}