cmake_minimum_required(VERSION "4.1.2")
project("lto-test" LANGUAGES Fortran)

cmake_policy(SET CMP0069 NEW)

add_library(foo foo.f)
add_executable(boo main.f)
target_link_libraries(boo PUBLIC foo)
