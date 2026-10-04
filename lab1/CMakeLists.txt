add_executable(lab1
    main.cpp
    src/functions.cpp
    src/matrix.cpp
)

target_include_directories(lab1
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set_target_properties(lab1 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)