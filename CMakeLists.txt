cmake_minimum_required(VERSION 3.0)

set(LIB_NAME EventEMin)
project(${LIB_NAME} LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(${LIB_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(${LIB_NAME}_SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(${LIB_NAME}_TEST_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/test)

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

find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -march=native -std=c++11")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

# interface lib
add_library(${LIB_NAME}_lib INTERFACE)
target_include_directories(${LIB_NAME}_lib INTERFACE ${${LIB_NAME}_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${GSL_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
target_link_libraries(${LIB_NAME}_lib INTERFACE ${GSL_LIBRARIES} ${OpenCV_LIBS})
# export as ${LIB_NAME}_lib
set_property(TARGET ${LIB_NAME}_lib PROPERTY EXPORT_NAME ${LIB_NAME}_lib)

add_subdirectory(${${LIB_NAME}_SOURCE_DIRS})

# build example binaries
option(EVENT_MODEL_BUILD_EXAMPLES "Build event model examples" ON)
if(EVENT_MODEL_BUILD_EXAMPLES)
  add_subdirectory(${${LIB_NAME}_TEST_DIRS})
endif()
