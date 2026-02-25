cmake_minimum_required(VERSION 3.16)

project(src LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

find_package(MPI REQUIRED)

add_executable(seidel
    main.c
    seidel.c
    tools.c
    functions.c
    study.c
)

target_link_libraries(seidel MPI::MPI_C)