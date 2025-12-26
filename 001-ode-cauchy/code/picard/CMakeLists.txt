cmake_minimum_required(VERSION 3.16)

project(picard LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

# --- OpenMP (macOS, AppleClang) ---
find_library(OMP_LIBRARY omp
    HINTS
        /usr/local/opt/libomp/lib
        /opt/homebrew/opt/libomp/lib
)

find_path(OMP_INCLUDE_DIR omp.h
    HINTS
        /usr/local/opt/libomp/include
        /opt/homebrew/opt/libomp/include
)

if (NOT OMP_LIBRARY OR NOT OMP_INCLUDE_DIR)
    message(FATAL_ERROR "libomp not found. Install with: brew install libomp")
endif()

add_executable(picard
    main.c
    picard.c
    generate.c
    matrix_tools.c
    study.c
)

target_include_directories(picard PRIVATE ${OMP_INCLUDE_DIR})
target_compile_options(picard PRIVATE -Xpreprocessor -fopenmp)
target_link_libraries(picard PRIVATE ${OMP_LIBRARY})
