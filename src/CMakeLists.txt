# Add library and set properties (C++ standard) for library files
# Set the compile options for the library files and compile with those flags

# Test for Compiler (GCC or Intel) and set appropriate flags for library file compilation
if( CMAKE_CXX_COMPILER_ID MATCHES GNU )
    list( APPEND flags_diffusion_2d "-O2" "-ftree-vectorize" "-Wall" )
endif()

if( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    list( APPEND flags_diffusion_2d "-O2" "-xhost" "-vec-report" "-Wall" )
endif()


# Test for Compiler (GCC or Intel) and set appropriate flags for library file compilation
if( CMAKE_CXX_COMPILER_ID MATCHES GNU )
    list( APPEND flags_diffusion "-O2" "-ftree-vectorize" "-Wall" "-Wpedantic" )
endif()

if( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    list( APPEND flags_diffusion "-O2" "-xhost" "-vec-report" "-Wall" )
endif()


#
# Add library files and create static library

add_library( diffusion_2d STATIC Diffusion_FDM_2d.h Diffusion_FDM_2d.cc )
target_compile_options( diffusion_2d PUBLIC ${flags_diffusion} ) 

#
# Add executable and set properties (C++ standard) for executable
# (analogous to library files) 

add_executable( diffusion diffusion.cc )
target_compile_options( diffusion PUBLIC ${flags_diffusion} ) 

#
# Output of Compile flags for library and Executables

message( STATUS "C++ compiler flags Library: ${flags_diffusion_2d}")
message( STATUS "C++ compiler flags Executable: ${flags_diffusion}")


#
# Set properties of C++ (Standrad, Extensions, ...)

set_target_properties( diffusion_2d 
		      PROPERTIES 
		      CXX_STANDARD 14 
		      CXX_EXTENSIONS OFF 
		      CXX_STANDARD_REQUIRED ON ) 

set_target_properties( diffusion 
		      PROPERTIES 
		      CXX_STANDARD 14 
		      CXX_EXTENSIONS OFF 
		      CXX_STANDARD_REQUIRED ON ) 

#
# Linkage part

target_link_libraries( diffusion diffusion_2d )

