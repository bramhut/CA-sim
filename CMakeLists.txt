# CMakeList.txt : CMake project for sim-aos, include source and define
# project specific logic here.
#

# Minimum version required for standard 23
cmake_minimum_required (VERSION 3.8)

project ("sim-aos")

set (CMAKE_CXX_STANDARD 20)

if (MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
else()
    add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors)
endif()

# Add source to this project's executable.
add_executable (sim-aos "sim-aos.cpp" "sim-aos.h" "object.h")

# TODO: Add tests and install targets if needed.
