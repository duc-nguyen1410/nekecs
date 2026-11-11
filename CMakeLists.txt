cmake_minimum_required(VERSION 3.12)
project(nekecs CXX)

# Backend-agnostic algorithms (header-only) and their tests.
enable_testing()
add_subdirectory(core)

# --- Nektar++ backend -------------------------------------------------------
# Configure with -DNektar++_DIR=<prefix>/lib64/nektar++/cmake
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

add_executable(nekecs apps/nekecs.cpp)
target_link_libraries(nekecs PRIVATE nekecs_nektar)

add_executable(test_flowmap nektar/tests/test_flowmap.cpp)
target_link_libraries(test_flowmap PRIVATE nekecs_nektar)

add_executable(flowmap_probe apps/probe/flowmap_probe.cpp)
target_include_directories(flowmap_probe PRIVATE
    ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
target_compile_definitions(flowmap_probe PRIVATE ${NEKTAR++_DEFINITIONS})
target_link_options(flowmap_probe PRIVATE "LINKER:--no-as-needed")
target_link_libraries(flowmap_probe PRIVATE ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})

install(TARGETS nekecs RUNTIME DESTINATION bin)

add_executable(diag_startup apps/diag_startup.cpp)
target_link_libraries(diag_startup PRIVATE nekecs_nektar)

add_executable(shed apps/shed.cpp)
target_link_libraries(shed PRIVATE nekecs_nektar)

add_executable(diag_param apps/diag_param.cpp)
target_link_libraries(diag_param PRIVATE nekecs_nektar)
