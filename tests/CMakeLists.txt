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
  ${CMAKE_SOURCE_DIR}/tests
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

function(gutibm_enable_test_assertions target)
  if(MSVC)
    set(cxx_undefine_ndebug /UNDEBUG)
  else()
    set(cxx_undefine_ndebug -UNDEBUG)
  endif()
  target_compile_options(${target} PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${cxx_undefine_ndebug}>
    $<$<COMPILE_LANGUAGE:CUDA>:-UNDEBUG>
  )
endfunction()

gutibm_enable_test_assertions(gutibm)

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
  target_compile_definitions(gutibm PUBLIC GUTIBM_CUDA)
  set_target_properties(gutibm PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
target_link_libraries(gutibm PUBLIC ${MATH_LIB})

function(gutibm_add_test name source)
  cmake_parse_arguments(ARG "" "" "LABELS;COMPILE_DEFINITIONS" ${ARGN})
  set(exe test_${name})
  add_executable(${exe} ${source})
  gutibm_enable_test_assertions(${exe})
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
gutibm_add_test(assertions_enabled test_assertions_enabled.cpp LABELS "unit")
gutibm_add_test(spatial_hash test_spatial_hash.cpp LABELS "unit")
gutibm_add_test(greens_function test_greens_function.cpp LABELS "unit")
gutibm_add_test(agent test_agent.cpp LABELS "unit")
gutibm_add_test(iron_fallback test_iron_fallback.cpp LABELS "unit")
gutibm_add_test(octree test_octree.cpp LABELS "unit")
gutibm_add_test(fmm test_fmm.cpp LABELS "unit")
gutibm_add_test(qssa_stoichiometry test_qssa_stoichiometry.cpp LABELS "unit")
gutibm_add_test(conjugation test_conjugation.cpp LABELS "unit")
gutibm_add_test(z_gradient test_z_gradient.cpp LABELS "unit")
gutibm_add_test(nutrient_diffusion test_nutrient_diffusion.cpp LABELS "unit")
gutibm_add_test(domain_decomp test_domain_decomp.cpp LABELS "unit")
gutibm_add_test(acetate_mete test_acetate_mete.cpp LABELS "unit")
gutibm_add_test(protease_decay test_protease_decay.cpp LABELS "unit")
gutibm_add_test(oxygen_gradient test_oxygen_gradient.cpp LABELS "unit")
gutibm_add_test(O2_growth_boost test_O2_growth_boost.cpp LABELS "unit")
gutibm_add_test(mucin_liberation test_mucin_liberation.cpp LABELS "unit")
gutibm_add_test(mechanism_wiring test_mechanism_wiring.cpp LABELS "integration")
gutibm_add_test(advection_peristaltic test_advection_peristaltic.cpp LABELS "unit")
gutibm_add_test(ethanolamine test_ethanolamine.cpp LABELS "unit")
gutibm_add_test(adaptive_dt test_adaptive_dt.cpp LABELS "unit")
gutibm_add_test(agent_transfer test_agent_transfer.cpp LABELS "unit")
gutibm_add_test(fix_registry test_fix_registry.cpp LABELS "unit")
gutibm_add_test(path_utils test_path_utils.cpp LABELS "unit")
gutibm_add_test(input_parser test_input_parser.cpp
  LABELS "unit"
  COMPILE_DEFINITIONS GUTIBM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
gutibm_add_test(config_ingestion test_config_ingestion.cpp
  LABELS "unit"
  COMPILE_DEFINITIONS GUTIBM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
gutibm_add_test(bacteriocin test_bacteriocin.cpp LABELS "unit")
gutibm_add_test(per_receptor_toxin test_per_receptor_toxin.cpp LABELS "unit")
gutibm_add_test(siderophore_depletion test_siderophore_depletion.cpp LABELS "unit")
gutibm_add_test(receptor test_receptor.cpp LABELS "unit")
gutibm_add_test(mutation test_mutation.cpp LABELS "unit")
gutibm_add_test(fur test_fur.cpp LABELS "unit")
gutibm_add_test(cdi test_cdi.cpp LABELS "unit")
gutibm_add_test(motility test_motility.cpp LABELS "unit")

# Scaling benchmark smoke (issue #55) — modest agent counts for CI.
gutibm_add_test(scaling_benchmark test_scaling_benchmark.cpp LABELS "benchmark;slow")

# Slower physics tests — still unit scope but excluded from fast CI shard.
gutibm_add_test(mechanics test_mechanics.cpp LABELS "unit;slow")
gutibm_add_test(immunity_escape test_immunity_escape.cpp LABELS "unit;slow")

# ── Integration tests (simulation / I/O; separate CI job) ────────────────
gutibm_add_test(smoke test_smoke.cpp LABELS "integration")
gutibm_add_test(washout_trap test_washout_trap.cpp LABELS "integration;slow;science")
gutibm_add_test(feature_combinations test_feature_combinations.cpp LABELS "integration")
gutibm_add_test(config_diversity test_config_diversity.cpp
  LABELS "integration"
  COMPILE_DEFINITIONS GUTIBM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
gutibm_add_test(openmp_parity test_openmp_parity.cpp LABELS "openmp")

if(GUTIBM_USE_HDF5)
  gutibm_add_test(hdf5_roundtrip test_hdf5_roundtrip.cpp LABELS "integration;hdf5")
  gutibm_add_test(hdf5_checkpoint test_hdf5_checkpoint.cpp LABELS "integration;hdf5")
  gutibm_add_test(hdf5_schedule test_hdf5_schedule.cpp LABELS "integration;hdf5")
  gutibm_add_test(hdf5_compression test_hdf5_compression.cpp LABELS "integration;hdf5")
  gutibm_add_test(summary_events test_summary_events.cpp LABELS "integration;hdf5")
  gutibm_add_test(hdf5_graceful_shutdown test_hdf5_graceful_shutdown.cpp LABELS "integration;hdf5")

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
  gutibm_enable_test_assertions(test_mpi_multi_rank)
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

  add_executable(test_mpi_four_rank test_mpi_four_rank.cpp)
  gutibm_enable_test_assertions(test_mpi_four_rank)
  target_link_libraries(test_mpi_four_rank PRIVATE gutibm)
  add_test(
    NAME mpi_four_rank
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
            --oversubscribe
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_four_rank> ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(mpi_four_rank PROPERTIES
    LABELS "integration;mpi"
    PASS_REGULAR_EXPRESSION "All MPI four-rank tests passed"
    RESOURCE_LOCK mpi
  )

  add_executable(test_mpi_gpu_multi_rank test_mpi_gpu_multi_rank.cpp)
  gutibm_enable_test_assertions(test_mpi_gpu_multi_rank)
  target_link_libraries(test_mpi_gpu_multi_rank PRIVATE gutibm)
  add_test(
    NAME mpi_gpu_multi_rank
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_gpu_multi_rank> ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(mpi_gpu_multi_rank PROPERTIES
    LABELS "integration;mpi;gpu"
    PASS_REGULAR_EXPRESSION "All MPI GPU multi-rank tests passed"
    RESOURCE_LOCK mpi
  )

  add_executable(test_cuda_aware_mpi_reaction test_cuda_aware_mpi_reaction.cpp)
  gutibm_enable_test_assertions(test_cuda_aware_mpi_reaction)
  target_link_libraries(test_cuda_aware_mpi_reaction PRIVATE gutibm)
  add_test(
    NAME cuda_aware_mpi_reaction
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            --oversubscribe
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_cuda_aware_mpi_reaction> ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(cuda_aware_mpi_reaction PROPERTIES
    LABELS "integration;mpi;gpu"
    PASS_REGULAR_EXPRESSION "All CUDA-aware MPI reaction tests passed"
    RESOURCE_LOCK mpi
  )
endif()

# ── GPU tests (CUDA CI job only) ─────────────────────────────────────────
gutibm_add_test(greens_function_gpu test_greens_function_gpu.cpp LABELS "gpu")
gutibm_add_test(gpu_diffusion test_gpu_diffusion.cpp LABELS "gpu")
gutibm_add_test(gpu_chemical_field test_gpu_chemical_field.cpp LABELS "gpu")
gutibm_add_test(gpu_feature_combinations test_gpu_feature_combinations.cpp LABELS "gpu;integration")
gutibm_add_test(gpu_production_path test_gpu_production_path.cpp LABELS "gpu;integration")
gutibm_add_test(gpu_smoke test_gpu_smoke.cpp LABELS "gpu")
gutibm_add_test(gpu_scaling_benchmark test_gpu_scaling_benchmark.cpp LABELS "gpu;benchmark")
gutibm_add_test(qssa_gpu_parity test_qssa_gpu_parity.cpp LABELS "gpu")
gutibm_add_test(spatial_hash_gpu_csr test_spatial_hash_gpu_csr.cpp LABELS "gpu")
gutibm_add_test(mechanics_gpu_parity test_mechanics_gpu_parity.cpp LABELS "gpu")

add_test(
  NAME rebuild_and_run_script
  COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_rebuild_and_run_script.sh"
)
set_tests_properties(rebuild_and_run_script PROPERTIES LABELS "unit")
