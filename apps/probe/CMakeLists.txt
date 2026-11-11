cmake_minimum_required(VERSION 3.12)
project(nekecs_probe CXX)

# Point CMAKE_PREFIX_PATH (or Nektar++_DIR) at the Nektar++ install, e.g.
#   cmake -DNektar++_DIR=$HOME/opt/nektar++/lib64/nektar++/cmake ..
find_package(Nektar++ REQUIRED)

add_executable(flowmap_probe flowmap_probe.cpp)
set_target_properties(flowmap_probe PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_include_directories(flowmap_probe PRIVATE ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
target_compile_definitions(flowmap_probe PRIVATE ${NEKTAR++_DEFINITIONS})
# Solver classes (UnsteadyNavierStokes, ...) register themselves with the
# EquationSystem factory from static initialisers. Nothing references their
# symbols directly, so Ubuntu's default --as-needed would drop the library.
target_link_options(flowmap_probe PRIVATE "LINKER:--no-as-needed")
target_link_libraries(flowmap_probe PRIVATE ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
