cmake_minimum_required(VERSION 3.12)
project(nekecs CXX)

# Backend-agnostic algorithms (header-only) and their tests (ctest).
enable_testing()
add_subdirectory(core)

# --- Nektar++ backend -------------------------------------------------------
# Configure with -DNektar++_DIR=<prefix>/lib64/nektar++/cmake, or with
# -DNEKECS_WITH_NEKTAR=OFF to build and test only the core library.
option(NEKECS_WITH_NEKTAR "Build the Nektar++ backend and applications" ON)
if(NOT NEKECS_WITH_NEKTAR)
    return()
endif()

find_package(Nektar++ REQUIRED)

add_library(nekecs_nektar INTERFACE)
target_include_directories(nekecs_nektar INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/nektar/include
    ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
target_compile_definitions(nekecs_nektar INTERFACE ${NEKTAR++_DEFINITIONS})
# Solver classes register with Nektar++'s factory from static initialisers;
# nothing references their symbols, so --as-needed would drop the library.
target_link_options(nekecs_nektar INTERFACE "LINKER:--no-as-needed")
target_link_libraries(nekecs_nektar INTERFACE nekecs_core
    ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})

# Applications: the solver, the DNS helper, and diagnostics.
foreach(app nekecs shed diag_startup diag_param)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE nekecs_nektar)
endforeach()
install(TARGETS nekecs shed RUNTIME DESTINATION bin)

# Adapter test: needs a run directory with a mesh (see README).
add_executable(test_flowmap nektar/tests/test_flowmap.cpp)
target_link_libraries(test_flowmap PRIVATE nekecs_nektar)
