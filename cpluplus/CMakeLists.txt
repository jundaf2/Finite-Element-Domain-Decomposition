cmake_minimum_required (VERSION 3.10)
project(cs484project VERSION 1.0.0 LANGUAGES C CXX)

# Location of additional cmake modules

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake
    )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_DATA_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/data)

include(FindOpenMP)
if(OPENMP_FOUND)
		message("Found OpenMP, building parallel part.")

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

else(OPENMP_FOUND)
		message(FATAL_ERROR "ERROR: OpenMP not found.")
endif(OPENMP_FOUND)

include(FindMPI)
if(MPI_CXX_FOUND)
	message("Found MPI")
else(MPI_CXX_FOUND)
	message(FATAL_ERROR "Could not find MPI")
endif(MPI_CXX_FOUND)


include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories("${PROJECT_SOURCE_DIR}/utils/include")


add_executable(feddm_serial_cg
		"serial/Serial_CGM.cpp"
		"serial/Serial_CGM.h"
		"serial/main.cpp"
)
target_include_directories(feddm_serial_cg PUBLIC "${PROJECT_SOURCE_DIR}/src")

add_executable(feddm_parallel_cg
		"parallel/main.cpp"
		"parallel/Parallel_CGM.cpp"
		"parallel/Parallel_CGM.h"
		"parallel/Serial_CGM.cpp"
		"parallel/Serial_CGM.h"
)
target_include_directories(feddm_parallel_cg PUBLIC "${PROJECT_SOURCE_DIR}/src")

target_include_directories(feddm_parallel_cg PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(feddm_parallel_cg m ${MPI_CXX_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(feddm_parallel_cg PROPERTIES
  COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(feddm_parallel_cg PROPERTIES
  LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()


