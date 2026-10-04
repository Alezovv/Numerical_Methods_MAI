#include <fstream>
#include "matrix.h"

bool Read_System_From_File(const std::string& filename, Matrix& A, Vector& b) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    size_t n;
    file >> n;

    A.assign(n, Vector(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }

    b.assign(n, 0);
    for (size_t i = 0; i < n; ++i) {
        file >> b[i];
    }
    return true;
}