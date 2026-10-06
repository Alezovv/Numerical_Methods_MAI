add_executable(lab3
    main.cpp
    src/interpolation.cpp
    src/spline.cpp
    src/ols.cpp
    src/diff.cpp
    src/integrate.cpp
)

target_include_directories(lab3 PRIVATE include)