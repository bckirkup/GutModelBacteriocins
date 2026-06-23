file(GLOB_RECURSE LIB_SOURCES
  ${CMAKE_SOURCE_DIR}/src/core/*.cpp
  ${CMAKE_SOURCE_DIR}/src/fields/*.cpp
  ${CMAKE_SOURCE_DIR}/src/diffusion/*.cpp
  ${CMAKE_SOURCE_DIR}/src/fixes/*.cpp
  ${CMAKE_SOURCE_DIR}/src/genome/*.cpp
  ${CMAKE_SOURCE_DIR}/src/io/*.cpp
)

set(TEST_INCLUDE_DIRS
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/fields
  ${CMAKE_SOURCE_DIR}/src/diffusion
  ${CMAKE_SOURCE_DIR}/src/fixes
  ${CMAKE_SOURCE_DIR}/src/genome
  ${CMAKE_SOURCE_DIR}/src/io
)

# ── test_spatial_hash ────────────────────────────────────────────────────
add_executable(test_spatial_hash test_spatial_hash.cpp ${LIB_SOURCES})
target_include_directories(test_spatial_hash PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_spatial_hash PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_spatial_hash PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_spatial_hash PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_spatial_hash PRIVATE ${MATH_LIB})
add_test(NAME spatial_hash COMMAND test_spatial_hash)

# ── test_greens_function ─────────────────────────────────────────────────
add_executable(test_greens_function test_greens_function.cpp ${LIB_SOURCES})
target_include_directories(test_greens_function PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_greens_function PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_greens_function PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_greens_function PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_greens_function PRIVATE ${MATH_LIB})
add_test(NAME greens_function COMMAND test_greens_function)

# ── test_agent ───────────────────────────────────────────────────────────
add_executable(test_agent test_agent.cpp ${LIB_SOURCES})
target_include_directories(test_agent PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_agent PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_agent PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_agent PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_agent PRIVATE ${MATH_LIB})
add_test(NAME agent COMMAND test_agent)

# ── test_iron_fallback ────────────────────────────────────────────────────
add_executable(test_iron_fallback test_iron_fallback.cpp ${LIB_SOURCES})
target_include_directories(test_iron_fallback PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_iron_fallback PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_iron_fallback PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_iron_fallback PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_iron_fallback PRIVATE ${MATH_LIB})
add_test(NAME iron_fallback COMMAND test_iron_fallback)

# ── test_octree ──────────────────────────────────────────────────────────
add_executable(test_octree test_octree.cpp ${LIB_SOURCES})
target_include_directories(test_octree PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_octree PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_octree PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_octree PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_octree PRIVATE ${MATH_LIB})
add_test(NAME octree COMMAND test_octree)

# ── test_smoke (end-to-end integration) ──────────────────────────────────
add_executable(test_smoke test_smoke.cpp ${LIB_SOURCES})
target_include_directories(test_smoke PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_smoke PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_smoke PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_smoke PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_smoke PRIVATE ${MATH_LIB})
add_test(NAME smoke COMMAND test_smoke)

# ── test_z_gradient (z-dependent nutrient gradient) ──────────────────────
add_executable(test_z_gradient test_z_gradient.cpp ${LIB_SOURCES})
target_include_directories(test_z_gradient PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_z_gradient PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_z_gradient PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_z_gradient PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_z_gradient PRIVATE ${MATH_LIB})
add_test(NAME z_gradient COMMAND test_z_gradient)

# ── test_domain_decomp (MPI decomposition unit tests) ────────────────────
add_executable(test_domain_decomp test_domain_decomp.cpp ${LIB_SOURCES})
target_include_directories(test_domain_decomp PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_domain_decomp PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_domain_decomp PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_domain_decomp PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_domain_decomp PRIVATE ${MATH_LIB})
add_test(NAME domain_decomp COMMAND test_domain_decomp)

# ── test_acetate_mete (acetate inhibition of MetE, issue #8) ─────────────
add_executable(test_acetate_mete test_acetate_mete.cpp ${LIB_SOURCES})
target_include_directories(test_acetate_mete PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_acetate_mete PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_acetate_mete PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_acetate_mete PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_acetate_mete PRIVATE ${MATH_LIB})
add_test(NAME acetate_mete COMMAND test_acetate_mete)

# ── test_advection_peristaltic ───────────────────────────────────────────
add_executable(test_advection_peristaltic test_advection_peristaltic.cpp ${LIB_SOURCES})
target_include_directories(test_advection_peristaltic PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_advection_peristaltic PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_advection_peristaltic PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_advection_peristaltic PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_advection_peristaltic PRIVATE ${MATH_LIB})
add_test(NAME advection_peristaltic COMMAND test_advection_peristaltic)

# ── test_immunity_escape (Issue #12: affinity-neutralization matrix) ─────
add_executable(test_immunity_escape test_immunity_escape.cpp ${LIB_SOURCES})
target_include_directories(test_immunity_escape PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_immunity_escape PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_immunity_escape PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_immunity_escape PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_immunity_escape PRIVATE ${MATH_LIB})
add_test(NAME immunity_escape COMMAND test_immunity_escape)

# ── test_mechanics ───────────────────────────────────────────────────────
add_executable(test_mechanics test_mechanics.cpp ${LIB_SOURCES})
target_include_directories(test_mechanics PRIVATE ${TEST_INCLUDE_DIRS})
if(GUTIBM_USE_MPI)
  target_link_libraries(test_mechanics PRIVATE MPI::MPI_CXX)
endif()
if(GUTIBM_USE_HDF5)
  target_include_directories(test_mechanics PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(test_mechanics PRIVATE ${HDF5_LIBRARIES})
endif()
target_link_libraries(test_mechanics PRIVATE ${MATH_LIB})
add_test(NAME mechanics COMMAND test_mechanics)
