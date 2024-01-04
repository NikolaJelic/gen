cmake_minimum_required(VERSION 3.24)

project(genetic)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(${PROJECT_NAME}  src/main.cpp
                                src/gene.hpp
                                src/gene.cpp
                                src/population.cpp
                                src/population.hpp
                                src/simulation.cpp
                                src/simulation.hpp
                                )

