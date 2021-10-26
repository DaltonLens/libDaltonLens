cmake_minimum_required(VERSION 3.16)

project(libDaltonLens)

include_directories(${CMAKE_SOURCE_DIR})

if (UNIX)
    add_compile_options(-Wall -pedantic)
endif()

add_definitions(-DTEST_IMAGES_DIR="${CMAKE_SOURCE_DIR}/tests/images/")

# Can't believe I had to go through that, but test aren't built automatically by cmake,
# so we need a special dependency hack. See https://stackoverflow.com/a/10824578
function (add_dl_test)
    add_test (NAME "${ARGV0}_BUILD" COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target "${ARGV0}" "--config" "$<CONFIG>" )
    add_test ("${ARGV0}" "${ARGV0}")
    add_executable("${ARGV0}" 
        "tests/${ARGV0}.c" 
        libDaltonLens.c
        libDaltonLens.h
    )
    set_tests_properties("${ARGV0}" PROPERTIES DEPENDS "${ARGV0}_BUILD")
    # Common lib for all tests.
    if (UNIX)
        target_link_libraries("${ARGV0}" m)
    endif()
endfunction()

add_dl_test (test_simulation)
enable_testing()
