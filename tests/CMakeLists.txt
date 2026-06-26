file(GLOB_RECURSE LIB_SOURCES
  ${CMAKE_SOURCE_DIR}/src/core/*.cpp
  ${CMAKE_SOURCE_DIR}/src/fields/*.cpp
  ${CMAKE_SOURCE_DIR}/src/diffusion/*.cpp
  ${CMAKE_SOURCE_DIR}/src/fixes/*.cpp
  ${CMAKE_SOURCE_DIR}/src/genome/*.cpp
  ${CMAKE_SOURCE_DIR}/src/io/*.cpp
  ${CMAKE_SOURCE_DIR}/src/gpu/*.cpp
)

set(LIB_CUDA_SOURCES "")
if(GUTIBM_USE_CUDA)
  file(GLOB LIB_CUDA_SOURCES ${CMAKE_SOURCE_DIR}/src/gpu/*.cu)
endif()

set(TEST_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/fields
  ${CMAKE_SOURCE_DIR}/src/diffusion
  ${CMAKE_SOURCE_DIR}/src/fixes
  ${CMAKE_SOURCE_DIR}/src/genome
  ${CMAKE_SOURCE_DIR}/src/io
  ${CMAKE_SOURCE_DIR}/src/gpu
)

# Shared static library — compile core once instead of per-test executable.
add_library(gutibm STATIC ${LIB_SOURCES} ${LIB_CUDA_SOURCES})
target_include_directories(gutibm PUBLIC ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(gutibm PUBLIC MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(gutibm PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(gutibm PUBLIC ${HDF5_LIBRARIES})
endif()
if(GUTIBM_USE_OPENMP)
  target_link_libraries(gutibm PUBLIC OpenMP::OpenMP_CXX)
endif()
if(GUTIBM_USE_CUDA)
  target_link_libraries(gutibm PUBLIC CUDA::cudart)
  set_target_properties(gutibm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
target_link_libraries(gutibm PUBLIC ${MATH_LIB})

function(gutibm_add_test name source)
  cmake_parse_arguments(ARG "" "" "LABELS;COMPILE_DEFINITIONS" ${ARGN})
  set(exe test_${name})
  add_executable(${exe} ${source})
  target_link_libraries(${exe} PRIVATE gutibm)
  if(ARG_COMPILE_DEFINITIONS)
    target_compile_definitions(${exe} PRIVATE ${ARG_COMPILE_DEFINITIONS})
  endif()
  add_test(NAME ${name} COMMAND ${exe})
  if(ARG_LABELS)
    set_tests_properties(${name} PROPERTIES LABELS "${ARG_LABELS}")
  endif()
endfunction()

# ── Fast unit tests (default CI gate) ────────────────────────────────────
gutibm_add_test(spatial_hash test_spatial_hash.cpp LABELS "unit")
gutibm_add_test(greens_function test_greens_function.cpp LABELS "unit")
gutibm_add_test(agent test_agent.cpp LABELS "unit")
gutibm_add_test(iron_fallback test_iron_fallback.cpp LABELS "unit")
gutibm_add_test(octree test_octree.cpp LABELS "unit")
gutibm_add_test(fmm test_fmm.cpp LABELS "unit")
gutibm_add_test(conjugation test_conjugation.cpp LABELS "unit")
gutibm_add_test(z_gradient test_z_gradient.cpp LABELS "unit")
gutibm_add_test(domain_decomp test_domain_decomp.cpp LABELS "unit")
gutibm_add_test(acetate_mete test_acetate_mete.cpp LABELS "unit")
gutibm_add_test(advection_peristaltic test_advection_peristaltic.cpp LABELS "unit")
gutibm_add_test(ethanolamine test_ethanolamine.cpp LABELS "unit")
gutibm_add_test(adaptive_dt test_adaptive_dt.cpp LABELS "unit")
gutibm_add_test(agent_transfer test_agent_transfer.cpp LABELS "unit")
gutibm_add_test(fix_registry test_fix_registry.cpp LABELS "unit")
gutibm_add_test(input_parser test_input_parser.cpp
  LABELS "unit"
  COMPILE_DEFINITIONS GUTIBM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
gutibm_add_test(bacteriocin test_bacteriocin.cpp LABELS "unit")
gutibm_add_test(receptor test_receptor.cpp LABELS "unit")
gutibm_add_test(mutation test_mutation.cpp LABELS "unit")

# Scaling benchmark smoke (issue #55) — modest agent counts for CI.
gutibm_add_test(scaling_benchmark test_scaling_benchmark.cpp LABELS "benchmark;slow")

# Slower physics tests — still unit scope but excluded from fast CI shard.
gutibm_add_test(mechanics test_mechanics.cpp LABELS "unit;slow")
gutibm_add_test(immunity_escape test_immunity_escape.cpp LABELS "unit;slow")

# ── Integration tests (simulation / I/O; separate CI job) ────────────────
gutibm_add_test(smoke test_smoke.cpp LABELS "integration")
gutibm_add_test(openmp_parity test_openmp_parity.cpp LABELS "openmp")

if(GUTIBM_USE_HDF5)
  gutibm_add_test(hdf5_roundtrip test_hdf5_roundtrip.cpp LABELS "integration;hdf5")
  gutibm_add_test(hdf5_checkpoint test_hdf5_checkpoint.cpp LABELS "integration;hdf5")

  if(GUTIBM_USE_MPI)
    add_test(
      NAME hdf5_roundtrip_parallel
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_hdf5_roundtrip> ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(hdf5_roundtrip_parallel PROPERTIES
      LABELS "integration;mpi;hdf5"
      PASS_REGULAR_EXPRESSION "All HDF5 round-trip tests passed"
      RESOURCE_LOCK mpi
    )
  endif()
endif()

if(GUTIBM_USE_MPI)
  add_executable(test_mpi_multi_rank test_mpi_multi_rank.cpp)
  target_link_libraries(test_mpi_multi_rank PRIVATE gutibm)
  add_test(
    NAME mpi_multi_rank
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_multi_rank> ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(mpi_multi_rank PROPERTIES
    LABELS "integration;mpi"
    PASS_REGULAR_EXPRESSION "All MPI multi-rank tests passed"
    RESOURCE_LOCK mpi
  )
endif()

# ── GPU tests (CUDA CI job only) ─────────────────────────────────────────
gutibm_add_test(greens_function_gpu test_greens_function_gpu.cpp LABELS "gpu")
gutibm_add_test(gpu_smoke test_gpu_smoke.cpp LABELS "gpu")
