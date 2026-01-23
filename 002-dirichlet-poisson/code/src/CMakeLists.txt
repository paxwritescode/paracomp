cmake_minimum_required(VERSION 3.16)

project(src LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

add_executable(seidel
    main.c
    seidel.c
    tools.c
)