project(picard LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug) 

add_executable(picard main.c picard.c picard.h generate.c generate.h matrices.c matrices.h)