cmake_minimum_required(VERSION 3.15)
project(dep-libcint)

include(ExternalProject)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

set(CINT_SRC "" CACHE STRING "Directory of source code (libcint.a)")
set(CINT_VER "v6.1.2" CACHE STRING "Version of library (libcint.a)")
if(NOT CINT_SRC AND (NOT $ENV{CINT_SRC} EQUAL ""))
    set(CINT_SRC $ENV{CINT_SRC})
endif()

if(CINT_SRC)
    message(STATUS "User specified directory of source code (libcint.a): ${CINT_SRC}")
else()
    set(CINT_SRC https://github.com/sunqm/libcint.git)
    message(STATUS "Download source code (libcint.a): ${CINT_SRC}")
endif()


# build with git source
ExternalProject_Add(libcint
    GIT_REPOSITORY ${CINT_SRC}
    GIT_TAG ${CINT_VER}
    PREFIX ${PROJECT_BINARY_DIR}/deps
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=RELEASE
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DWITH_RANGE_COULOMB:STRING=1
        -DWITH_FORTRAN:STRING=0
        -DWITH_CINT2_INTERFACE:STRING=1
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
        -DENABLE_STATIC=1
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/deps/lib DESTINATION ${CMAKE_INSTALL_PREFIX})

