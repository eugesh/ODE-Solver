cmake_minimum_required(VERSION 3.16)
project(numerical_task_9)

set(CMAKE_CXX_COMPILER g++-10)
add_compile_options(-Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -pedantic -pthread)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
find_package (Threads REQUIRED)

add_executable(numerical_task_9 src/main.cpp
        src/Context.cpp src/Context.h
        src/ODESolver.cpp src/ODESolver.h
        src/Algebra.cpp src/Algebra.h
        src/NewtonSolver.cpp src/NewtonSolver.h
        src/Utils.cpp src/Utils.h)


    target_link_libraries(numerical_task_9
        ${CMAKE_THREAD_LIBS_INIT}
    )
