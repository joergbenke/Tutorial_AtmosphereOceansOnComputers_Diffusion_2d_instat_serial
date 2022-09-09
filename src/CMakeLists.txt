#
# CMake file for creating build system for diffusion problem
#


# Define required minimum CMake version
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)


# Define project name
project(diffusion_fdm_1d LANGUAGES CXX)


#
# Testing and configuration part
#

# Testing for operating system
if( CMAKE_SYSTEM_NAME STREQUAL "Linux" )
    message( STATUS "Compiling for Linux ...." )
else()
    message( STATUS "Configuring on/for ${CMAKE_SYSTEM_NAME}" )
endif()

# Testing for compiler ID, compiler, etc.
message( STATUS "Is the C++ compiler loaded? ${CMAKE_CXX_COMPILER_LOADED}" )
if( CMAKE_CXX_COMPILER_LOADED )
  message( STATUS "The C++ compiler ID is: ${CMAKE_CXX_COMPILER_ID}" )
#  message( STATUS "Is the C++ from GNU? ${CMAKE_COMPILER_IS_GNUCXX}" )
  message( STATUS "The C++ compiler version is: ${CMAKE_CXX_COMPILER_VERSION}" )
endif()

# Which kind of host processor (e.g. x86_64)
message( STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}" ) 

#
# Compile part
#

#
# Add library and set properties (C++ standard) for library files

# Add library files and create static library
add_library( diffusion_1d STATIC Diffusion_FDM_1d.h Diffusion_FDM_1d.cc )

# Set properties of C++ (Standrad, Extensions, ...)
set_target_properties( diffusion_1d 
		      PROPERTIES 
		      CXX_STANDARD 14 
		      CXX_EXTENSIONS OFF 
		      CXX_STANDARD_REQUIRED ON ) 

# Test for Compiler (GCC or Intel) and set appropriate flags for library file compilation
if( CMAKE_CXX_COMPILER_ID MATCHES GNU )
    list( APPEND flags_diffusion_1d "-O2" "-ftree-vectorize" "-Wall" )
endif()

if( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    list( APPEND flags_diffusion_1d "-O2" "-xhost" "-vec-report" "-Wall" )
endif()

# Set the compile options for the library files and compile with those flags
target_compile_options( diffusion_1d PUBLIC ${flags_diffusion_1d} ) 



#
# Add executable and set properties (C++ standard) for executable
# (analogous to library files) 

add_executable( diffusion diffusion.cc )

set_target_properties( diffusion 
		      PROPERTIES 
		      CXX_STANDARD 14 
		      CXX_EXTENSIONS OFF 
		      CXX_STANDARD_REQUIRED ON ) 

# or e.g. set( CMAKE_CXX_STANDARD 14)

if( CMAKE_CXX_COMPILER_ID MATCHES GNU )
    list( APPEND flags_diffusion "-O2" "-ftree-vectorize" "-Wall" )
endif()

if( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    list( APPEND flags_diffusion "-O2" "-xhost" "-vec-report" "-Wall" )
endif()

target_compile_options( diffusion PUBLIC ${flags_diffusion} ) 

# Output of Compile flags for library and Executables
message( STATUS "C++ compiler flags Library: ${flags_diffusion_1d}")
message( STATUS "C++ compiler flags Executable: ${flags_diffusion}")


#
# Linkage part
#

# Linkage
target_link_libraries( diffusion diffusion_1d )



