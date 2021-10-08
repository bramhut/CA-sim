# CMakeList.txt : CMake project for sim-aos, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("sim-aos")

if (MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors)
endif()

# Add source to this project's executable.
add_executable (sim-aos "sim-aos.cpp" "sim-aos.h" "object.h")

# TODO: Add tests and install targets if needed.
