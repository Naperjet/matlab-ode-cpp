cmake_minimum_required(VERSION 3.10)
project(matlab-ode-cpp CXX)

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the linear algebra library and the ODE suite as subdirectories
add_subdirectory(lin-alg)
add_subdirectory(ode-solvers)
