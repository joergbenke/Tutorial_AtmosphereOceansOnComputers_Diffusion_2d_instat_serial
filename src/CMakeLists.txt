cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(diffusion_1d_ftcs LANGUAGES CXX)

add_executable(Diffusion diffusion.cc)
set_target_properties( Diffusion 
		      PROPERTIES 
		      CXX_STANDARD 14 
		      CXX_EXTENSIONS OFF 
		      CXX_STANDARD_REQUIRED ON ) 

# or e.g. set( CMAKE_CXX_STANDARD 14)

#
if( CMAKE_SYSTEM_NAME STREQUAL "Linux" )
    message( "Compiling for Linux ...." )
else()
    message(STATUS "Configuring on/for ${CMAKE_SYSTEM_NAME}")
endif()

#
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
message( "GNU" )
			 list(APPEND CMAKE_CXX_FLAGS "-fno-rtti -fno-exceptions")
  list(APPEND CMAKE_CXX_FLAGS_DEBUG "-Wsuggest-final-types" "-Wsuggest-final-methods" "-Wsuggest-override")
  list(APPEND CMAKE_CXX_FLAGS_RELEASE "-O3" "-Wno-unused")
endif()

#list( APPEND flags "-O1" "Wall" )
#if( NOT WIN32 )
#    list( APPEND flags "-Wextra" "-Wpedantic" )
#endif()

target_compile_options( Diffusion PUBLIC "-O1" "-Wall" "-fPIC" )
message( "C++ compiler flags: ${CMAKE_CXX_FLAGS}" )

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
#endif()


#
if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "i386")
  message(STATUS "i386 architecture detected")
elseif(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "i686")
  message(STATUS "i686 architecture detected")
elseif(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "x86_64")
  message(STATUS "x86_64 architecture detected")
else()
  message(STATUS "host processor architecture is unknown")
endif()

#
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

message(STATUS "C++ flags, Debug configuration: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "C++ flags, Release configuration: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "C++ flags, Release configuration with Debug info: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "C++ flags, minimal Release configuration: ${CMAKE_CXX_FLAGS_MINSIZEREL}")

