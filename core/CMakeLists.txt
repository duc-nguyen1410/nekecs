cmake_minimum_required(VERSION 3.12)
project(nekecs_core CXX)

# Header-only, backend-agnostic Newton-Krylov-hookstep library.
add_library(nekecs_core INTERFACE)
target_include_directories(nekecs_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(nekecs_core INTERFACE cxx_std_17)

option(NEKECS_CORE_TESTS "Build nekecs core tests" ON)
if(NEKECS_CORE_TESTS)
    enable_testing()
    foreach(t test_linalg test_lorenz test_arnoldi test_continuation)
        add_executable(${t} tests/${t}.cpp)
        target_link_libraries(${t} PRIVATE nekecs_core)
        target_compile_options(${t} PRIVATE -Wall -Wextra -Wpedantic)
        add_test(NAME ${t} COMMAND ${t})
    endforeach()
endif()
