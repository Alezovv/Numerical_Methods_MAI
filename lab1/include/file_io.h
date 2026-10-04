#ifndef FILE_IO_H
#define FILE_IO_H

#include <string>   // Обязательно для std::string
#include "matrix.h" // Чтобы компилятор знал типы Matrix и Vector

bool Read_System_From_File(const std::string &filename, Matrix &A, Vector &b);

#endif // FILE_IO_H