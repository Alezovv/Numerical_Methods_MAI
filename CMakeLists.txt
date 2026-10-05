cmake_minimum_required(VERSION 3.10)

project(Numerical_Methods_MAI VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Включаем строгие проверки кода (поможет избежать глупых опечаток)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra)
endif()

# Складываем все .exe файлы в общую папку build/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(lab1)
add_subdirectory(lab2)
# add_subdirectory(lab3)
# add_subdirectory(lab4)