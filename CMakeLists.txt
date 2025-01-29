cmake_minimum_required(VERSION 3.22)

project(
  MicroMechanics
  VERSION 0.1
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  set(MICROMECH_CPP_COMPILE_OPTIONS /W4 /bigobj)
else()
  set(MICROMECH_CPP_COMPILE_OPTIONS -Wall -Wextra -pedantic)
endif()

include(FetchContent)

# Fetch BioFVM
FetchContent_Declare(
  paraBioFVM
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL file://${CMAKE_CURRENT_SOURCE_DIR}/submodules/paraBioFVM)
FetchContent_MakeAvailable(paraBioFVM)

# Target BioFVMCore
file(GLOB src_files "src/*")
list(FILTER src_files EXCLUDE REGEX ".*main.*")

add_library(MicroMechanicsCore ${src_files})
target_compile_options(
  MicroMechanicsCore
  PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${MICROMECH_CPP_COMPILE_OPTIONS}>)

include_directories(include/MicroMechanics src)

target_include_directories(MicroMechanicsCore
                           PUBLIC ${paraBioFVM_SOURCE_DIR}/include)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(MicroMechanicsCore PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(MicroMechanicsCore PUBLIC BioFVMCore)

# Target MicroMechanics
add_executable(MicroMechanics src/main.cpp)
target_link_libraries(MicroMechanics MicroMechanicsCore)
