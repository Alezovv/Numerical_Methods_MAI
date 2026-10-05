add_executable(lab2
    main.cpp
    src/nonlinear_equations.cpp
)

# Подключаем папку с заголовочными файлами
target_include_directories(lab2 PRIVATE include)