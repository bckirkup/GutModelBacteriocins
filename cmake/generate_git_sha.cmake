if(NOT DEFINED GUTIBM_SOURCE_DIR OR NOT DEFINED GUTIBM_OUTPUT)
  message(FATAL_ERROR "GUTIBM_SOURCE_DIR and GUTIBM_OUTPUT are required")
endif()

set(git_sha "")
if(DEFINED GUTIBM_GIT_SHA_OVERRIDE AND NOT
   "${GUTIBM_GIT_SHA_OVERRIDE}" STREQUAL "")
  set(git_sha "${GUTIBM_GIT_SHA_OVERRIDE}")
else()
  find_program(git_executable git)
  if(git_executable AND EXISTS "${GUTIBM_SOURCE_DIR}/.git")
    execute_process(
      COMMAND "${git_executable}" -C "${GUTIBM_SOURCE_DIR}" rev-parse HEAD
      OUTPUT_VARIABLE git_sha
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      RESULT_VARIABLE git_result)
    if(NOT git_result EQUAL 0 OR "${git_sha}" STREQUAL "")
      set(git_sha "unknown-git-unavailable")
    else()
      execute_process(
        COMMAND "${git_executable}" -C "${GUTIBM_SOURCE_DIR}" status
                --porcelain --untracked-files=all
        OUTPUT_VARIABLE git_status
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE status_result)
      if(status_result EQUAL 0 AND NOT "${git_status}" STREQUAL "")
        string(APPEND git_sha "-dirty")
      endif()
    endif()
  else()
    set(git_sha "unknown-git-unavailable")
  endif()
endif()

get_filename_component(output_directory "${GUTIBM_OUTPUT}" DIRECTORY)
file(MAKE_DIRECTORY "${output_directory}")
set(header_content
    "#pragma once\n#define GUTIBM_GIT_SHA \"${git_sha}\"\n")
set(tmp_output "${GUTIBM_OUTPUT}.tmp")
file(WRITE "${tmp_output}" "${header_content}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${tmp_output}" "${GUTIBM_OUTPUT}")
file(REMOVE "${tmp_output}")
