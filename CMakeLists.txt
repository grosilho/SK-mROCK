cmake_minimum_required(VERSION 3.17...3.24)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
if(APPLE)
	message(STATUS "Platform: MacOSX")
	set(Eigen3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/eigen-3.4.0-MacOSX")
	set(autodiff_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/autodiff-0.5.0-MacOSX")
	set(Gmsh_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/gmsh-4.10.5-MacOSX")
	set(FFTW3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/fftw-3.3.10-MacOSX")
else()
	message(STATUS "Platform: Linux")
	set(Eigen3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/eigen-3.4.0")
	set(autodiff_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/autodiff-0.5.0")
	set(Gmsh_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/gmsh-4.10.5")
	set(FFTW3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/fftw-3.3.10")
endif()

project(MultirateIntegrators
        VERSION 1.0 
        DESCRIPTION "Implmementatiion of explicit stabilized multirate methods,  and other methods as well for comparison." 
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_OPENMP "Enables OpenMP" ON)
option(BUILD_SHARED_LIBS "Build shared (ON) or static (OFF) libraries" OFF)
option(GPT_PROFILE "Enables gperf-tools profiling by linking the appropriate libraries" OFF)
option(GNU_PROFILE "Enables GNU profiling by linking the appropriate libraries" OFF)
option(SCOREP_PROFILE "Enables ScoreP profiling by prepending scorep instrumentation" OFF)

if(GNU_PROFILE AND APPLE)
    message(FATAL_ERROR "Cannot profile with GNU Profile on MacOSX.")
endif()
if(GPT_PROFILE AND APPLE)
    message(FATAL_ERROR "Cannot profile with gperftools on MacOSX.")
endif()

find_package(Eigen3 REQUIRED CONFIG) # Eigen3::Eigen
find_package(autodiff REQUIRED CONFIG) # autodiff::autodiff
find_package(Gmsh REQUIRED) # ${GMSH_LIB}, ${GMSH_INC}
find_package(FFTW3 REQUIRED CONFIG) #${FFTW3_INCLUDE_DIRS}, ${FFTW3_LIBRARY_DIRS}, ${FFTW3_LIBRARIES}

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        link_libraries(OpenMP::OpenMP_CXX)
    endif()
endif()

if(GPT_PROFILE)
    link_libraries(-Wl,--no-as-needed,-lprofiler,--as-needed)
endif()
if(GNU_PROFILE)
    link_libraries(-pg)
endif()
if(SCOREP_PROFILE)
    if(ENABLE_OPENMP)
        message(AUTHOR_WARNING "Remember to set this environment variable before building:\n"
                 "export SCOREP_WRAPPER_INSTRUMENTER_FLAGS=--thread=omp")
    endif()
    message(AUTHOR_WARNING "It is useful to set some environment variables to limit profiling data, as\n"
            "export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=10")
    if(APPLE)
    	set(CMAKE_CXX_COMPILER "/opt/scorep/bin/scorep-clang++")
    else()
	set(CMAKE_CXX_COMPILER "scorep-g++")
    endif()
endif()

include_directories(BEFORE include)

add_subdirectory(src)

add_executable(main src/main.cpp)
target_link_libraries(main PUBLIC MainHeader parameters monte_carlo)

