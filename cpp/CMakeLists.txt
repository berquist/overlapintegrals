cmake_minimum_required(VERSION 3.11)
project(overlapintegrals CXX)

find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})

enable_testing()
# find_package(GTest REQUIRED)

add_executable(main.x main.cpp)
# target_link_libraries(main.x ${GTEST_BOTH_LIBRARIES})
# add_test(main main.x)
