cmake_minimum_required(VERSION 3.22)
project(numerical_task_9)

set(CMAKE_CXX_COMPILER g++-12)
add_compile_options(-Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -pedantic)
set(CMAKE_CXX_STANDARD 20)

add_executable(numerical_task_9 main.cpp Context.cpp Context.h Integrator.cpp Integrator.h)
