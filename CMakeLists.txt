#
# CMake file for creating build system for diffusion problem
#

cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(diffusion_fdm_2d VERSION 1.0 LANGUAGES CXX)


#
# Testing and configuration part
#

# Testing for operating system
if( CMAKE_SYSTEM_NAME STREQUAL "Linux" )
    message( STATUS "Compiling for Linux ..." )
else()
    message( STATUS "Configuring on/for ${CMAKE_SYSTEM_NAME}" )
endif()

# Testing for compiler ID, compiler, etc.
message( STATUS "Is the C++ compiler loaded? ${CMAKE_CXX_COMPILER_LOADED}" )
if( CMAKE_CXX_COMPILER_LOADED )
  message( STATUS "The C++ compiler ID is: ${CMAKE_CXX_COMPILER_ID}" )
  message( STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}" )
endif()

# Which kind of host processor (e.g. x86_64)
message( STATUS "The host processor is: ${CMAKE_HOST_SYSTEM_PROCESSOR}" ) 

add_subdirectory(src)


