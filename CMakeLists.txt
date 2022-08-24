cmake_minimum_required(VERSION 3.14)
project(hungarian_algorithm)


set(CMAKE_CXX_STANDARD 17)


#===============
#  Google GTest
#===============

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()


#==============
#  Build Setup
#==============

set( INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hungarian_algorithm/include
)
include_directories(
    ${INCLUDE_DIRS}
)


set( SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hungarian_algorithm/src/hungarian_assigner.cc
)

add_library( hungarian_algorithm
    ${SOURCE}
)


#==============
#  Test Setup
#==============

add_executable( test_hungarian_algorithm
    test/test.cc
)

target_link_libraries( test_hungarian_algorithm
    hungarian_algorithm
    gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_hungarian_algorithm)
