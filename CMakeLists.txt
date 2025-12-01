cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


include_directories(src nanoblas/src)

add_subdirectory (src)
add_subdirectory (nanoblas)


add_executable (test_ode demos/test_ode.cpp)
target_link_libraries (test_ode PUBLIC nanoblas)

add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries (demo_autodiff PUBLIC nanoblas)

add_executable(demo_rc_circuit demos/demo_rc_circuit.cpp)
target_link_libraries(demo_rc_circuit PUBLIC nanoblas)

add_executable(demo_pendulum demos/demo_pendulum.cpp)
target_link_libraries(demo_pendulum PUBLIC nanoblas)

add_executable(demo_legendre demos/demo_legendre.cpp)
target_link_libraries(demo_legendre PUBLIC nanoblas)
