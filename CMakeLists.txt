cmake_minimum_required(VERSION 3.0.0)
project(lek12 VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()
include_directories("/home/emil/cpp/Numerical/NR_C301/code")
add_executable(lek12 main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
