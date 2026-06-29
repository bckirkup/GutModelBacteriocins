cmake_minimum_required(VERSION 3.16)
project(GutIBM
  VERSION 0.1.0
  LANGUAGES C CXX
  DESCRIPTION "3D Individual-based Model for Enterobacteriaceae gut dynamics"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ── Options ──────────────────────────────────────────────────────────────
option(GUTIBM_USE_MPI    "Enable MPI parallelism"   ON)
option(GUTIBM_USE_HDF5   "Enable HDF5 output"       ON)
option(GUTIBM_USE_OPENMP "Enable OpenMP parallelism" OFF)
option(GUTIBM_USE_CUDA   "Enable CUDA GPU acceleration" OFF)
option(BUILD_TESTING      "Build unit tests"          ON)

# ── MPI ──────────────────────────────────────────────────────────────────
if(GUTIBM_USE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  add_definitions(-DGUTIBM_MPI)
endif()

# ── OpenMP ─────────────────────────────────────────────────────────────────
if(GUTIBM_USE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  add_definitions(-DGUTIBM_OPENMP)
endif()

# ── CUDA ─────────────────────────────────────────────────────────────────
if(GUTIBM_USE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "60;70;80")
  endif()
  find_package(CUDAToolkit REQUIRED)
  add_definitions(-DGUTIBM_CUDA)
endif()

# ── HDF5 ─────────────────────────────────────────────────────────────────
if(GUTIBM_USE_HDF5)
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 REQUIRED COMPONENTS C)
  add_definitions(-DGUTIBM_HDF5)
endif()

# ── Math ─────────────────────────────────────────────────────────────────
find_library(MATH_LIB m)

# ── Source collection ────────────────────────────────────────────────────
file(GLOB_RECURSE SOURCES
  src/core/*.cpp
  src/fields/*.cpp
  src/diffusion/*.cpp
  src/fixes/*.cpp
  src/genome/*.cpp
  src/io/*.cpp
  src/gpu/*.cpp
)

set(GPU_CUDA_SOURCES "")
if(GUTIBM_USE_CUDA)
  file(GLOB GPU_CUDA_SOURCES src/gpu/*.cu)
endif()

add_executable(gut_ibm src/main.cpp ${SOURCES} ${GPU_CUDA_SOURCES})

target_include_directories(gut_ibm PRIVATE
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/fields
  ${CMAKE_SOURCE_DIR}/src/diffusion
  ${CMAKE_SOURCE_DIR}/src/fixes
  ${CMAKE_SOURCE_DIR}/src/genome
  ${CMAKE_SOURCE_DIR}/src/io
  ${CMAKE_SOURCE_DIR}/src/gpu
)

if(GUTIBM_USE_MPI)
  target_link_libraries(gut_ibm PRIVATE MPI::MPI_CXX)
endif()

if(GUTIBM_USE_HDF5)
  target_include_directories(gut_ibm PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(gut_ibm PRIVATE ${HDF5_LIBRARIES})
endif()

if(GUTIBM_USE_OPENMP)
  target_link_libraries(gut_ibm PRIVATE OpenMP::OpenMP_CXX)
endif()

if(GUTIBM_USE_CUDA)
  target_link_libraries(gut_ibm PRIVATE CUDA::cudart)
  set_target_properties(gut_ibm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

target_link_libraries(gut_ibm PRIVATE ${MATH_LIB})

# ── Tests ────────────────────────────────────────────────────────────────
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

install(TARGETS gut_ibm DESTINATION bin)
