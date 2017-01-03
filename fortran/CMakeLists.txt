project(pyjac)

cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

if(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  set(C_FLAGS "${CMAKE_C_FLAGS}")
	set(C_FLAGS_DEBUG "-g -O0 -Wall ")
  set(C_FLAGS_RELEASE "-O3 -Wall ")
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
  set(C_FLAGS "${CMAKE_C_FLAGS}")
	set(C_FLAGS_DEBUG "-g -O0 -debug all -traceback -Wall -check=conversions,stack,uninit")
	set(C_FLAGS_RELEASE "-O3 -ip -xHost")
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	set(CXX_FLAGS_DEBUG "-g -O0 -Wall ")
  set(CXX_FLAGS_RELEASE "-O3 -Wall ")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	set(CXX_FLAGS_DEBUG "-g -O0 -debug all -traceback -Wall -check=conversions,stack,uninit")
	set(CXX_FLAGS_RELEASE "-O3 -ip -xHost")
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  set(Fortran_FLAGS_DEBUG "-g -O0 -Wall ")
  set(Fortran_FLAGS_RELEASE "-O3 -Wall ")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  set(Fortran_FLAGS_DEBUG "-g -O0 -check all -debug all -traceback -warn all")
  set(Fortran_FLAGS_RELEASE "-O3 -ip -xHost")
endif()

set(CMAKE_C_FLAGS ${C_FLAGS} CACHE STRING "C compiler flags for all builds" FORCE)
set(CMAKE_C_FLAGS_DEBUG ${C_FLAGS_DEBUG} CACHE STRING "Additional C compiler flags for debug builds" FORCE)
set(CMAKE_C_FLAGS_RELEASE ${C_FLAGS_RELEASE} CACHE STRING "Additional C compiler flags for release builds" FORCE)

set(CMAKE_CXX_FLAGS ${CXX_FLAGS} CACHE STRING "CXX compiler flags for all builds" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG ${CXX_FLAGS_DEBUG} CACHE STRING "Additional CXX compiler flags for debug builds" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE ${CXX_FLAGS_RELEASE} CACHE STRING "Additional CXX compiler flags for release builds" FORCE)

set(CMAKE_Fortran_FLAGS ${Fortran_FLAGS} CACHE STRING "Fortran compiler flags for all builds" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG ${Fortran_FLAGS_DEBUG} CACHE STRING "Additional Fortran compiler flags for debug builds" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE ${Fortran_FLAGS_RELEASE} CACHE STRING "Additional Fortran compiler flags for release builds" FORCE)

set(PROJECT_INCLUDE_DIR "${PROJECT_BINARY_DIR}/include")
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_INCLUDE_DIR}")

message(STATUS "Using mechanism: " ${MECH})
#===========================================
#  pyjac sources
#===========================================
include_directories(${MECH})
include_directories(${MECH}/jacobs)
set_source_files_properties("${MECH}/mechanism.c" PROPERTIES LANGUAGE CXX)
FILE(GLOB pyjac_src "${MECH}/*.c")
#===========================================
#  Fortran sources
#===========================================
set(pyjac_src ${pyjac_src} ${MECH}/pyjac_mod.f90)
#===========================================
#  Build library sources
#===========================================
message(STATUS ${pyjac_src})
add_library(pyjac STATIC ${pyjac_src})
add_executable(test_pyjac ${PROJECT_SOURCE_DIR}/test_pyjac.f90)
SET_TARGET_PROPERTIES(test_pyjac PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(test_pyjac pyjac)
install(TARGETS pyjac test_pyjac
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
)
