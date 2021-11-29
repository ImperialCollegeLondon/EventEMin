cmake_minimum_required(VERSION 3.16)

set(LIB_NAME EventEMin)
project(${LIB_NAME} LANGUAGES CXX)

# C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(${LIB_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(${LIB_NAME}_SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(${LIB_NAME}_TEST_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/test)

# CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -march=native -mtune=native")

# Default to Release
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Dependencies

# OpenMP
option(USE_OPENMP "If available, use OpenMP for parallelization." ON)
if(USE_OPENMP)
  find_package(OpenMP)
endif()

if(OPENMP_FOUND)
  message("OpenMP found")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

# Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# OpenCV
find_package(OpenCV REQUIRED)

# Build library
add_library(${LIB_NAME} INTERFACE)

# Build batch mode
option(${LIB_NAME}_BATCH_MODE "Use batch mode." OFF)
if(${LIB_NAME}_BATCH_MODE)
  # gsl
  find_package(GSL REQUIRED)

  # interface batch lib
  add_library(${LIB_NAME}_BATCH_LIB INTERFACE)
  target_include_directories(${LIB_NAME}_BATCH_LIB INTERFACE ${${LIB_NAME}_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GSL_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
  target_link_libraries(${LIB_NAME}_BATCH_LIB INTERFACE ${GSL_LIBRARIES} ${OpenCV_LIBS})
  
  target_compile_definitions(${LIB_NAME} INTERFACE ${LIB_NAME}_BATCH_MODE)
  
  target_link_libraries(${LIB_NAME} INTERFACE ${LIB_NAME}_BATCH_LIB)
endif()

# Build incremental mode
option(${LIB_NAME}_INCREMENTAL_MODE "Use incremental mode." OFF)
if(${LIB_NAME}_INCREMENTAL_MODE)
  # interface incremental lib
  add_library(${LIB_NAME}_INCREMENTAL_LIB INTERFACE)
  target_include_directories(${LIB_NAME}_INCREMENTAL_LIB INTERFACE ${${LIB_NAME}_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS})
  target_link_libraries(${LIB_NAME}_INCREMENTAL_LIB INTERFACE ${OpenCV_LIBS})
  
  target_compile_definitions(${LIB_NAME} INTERFACE ${LIB_NAME}_INCREMENTAL_MODE)
  
  target_link_libraries(${LIB_NAME} INTERFACE ${LIB_NAME}_INCREMENTAL_LIB)
endif()

# Option for fast exp computation
option(${LIB_NAME}_FAST_EXP "Enable fast exp" ON)
if(${LIB_NAME}_FAST_EXP)
  target_compile_definitions(${LIB_NAME} INTERFACE ${LIB_NAME}_FAST_EXP)
endif()

add_subdirectory(${${LIB_NAME}_SOURCE_DIRS})
add_subdirectory(${${LIB_NAME}_TEST_DIRS})
