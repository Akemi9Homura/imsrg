cmake_minimum_required(VERSION 3.16)

project(mrimsrg_input_bridge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SIMPLE_NCSM_ROOT "/home/mengziyan/simple-ncsm" CACHE PATH
    "Path to the isolated simple-ncsm validation checkout")

if(NOT EXISTS "${SIMPLE_NCSM_ROOT}/CMakeLists.txt")
    message(FATAL_ERROR "simple-ncsm checkout not found: ${SIMPLE_NCSM_ROOT}")
endif()

# The bridge needs only the validated NN reader and simpleFCI path.  Disabling
# optional 3B readers keeps this small build independent of unused libraries.
set(TBME "NONE" CACHE STRING "simple-ncsm optional 3B reader" FORCE)
set(SIMPLE_FCI_EMAX "2" CACHE STRING "simpleFCI determinant emax")
add_subdirectory("${SIMPLE_NCSM_ROOT}" "${CMAKE_BINARY_DIR}/simple-ncsm" EXCLUDE_FROM_ALL)

execute_process(
    COMMAND git -C "${SIMPLE_NCSM_ROOT}" rev-parse HEAD
    OUTPUT_VARIABLE SIMPLE_NCSM_REVISION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE SIMPLE_NCSM_GIT_RESULT)
if(NOT SIMPLE_NCSM_GIT_RESULT EQUAL 0)
    set(SIMPLE_NCSM_REVISION "unknown")
endif()

add_executable(mrimsrg_prepare cpp/prepare_reference.cpp)
add_executable(mrimsrg_validate cpp/validate_flow.cpp)
add_executable(mrimsrg_export_no2bpack cpp/export_no2bpack.cpp)
add_executable(mrimsrg_minipack_to_j64 cpp/minipack_to_j64.cpp)
add_executable(mrimsrg_check_minipack_subset cpp/check_minipack_subset.cpp)
add_executable(mrimsrg_extract_minipack_subset cpp/extract_minipack_subset.cpp)
foreach(target mrimsrg_prepare mrimsrg_validate mrimsrg_export_no2bpack
               mrimsrg_minipack_to_j64 mrimsrg_check_minipack_subset
               mrimsrg_extract_minipack_subset)
    target_link_libraries(${target} PRIVATE fci nucleus util)
    target_compile_definitions(${target} PRIVATE
        SIMPLE_FCI_EMAX=${SIMPLE_FCI_EMAX}
        # Preserve historical metadata keys while pointing them at the
        # isolated validation checkout.
        SHELL_MODEL_OBS_ROOT="${SIMPLE_NCSM_ROOT}"
        SHELL_MODEL_OBS_REVISION="${SIMPLE_NCSM_REVISION}")
    target_compile_options(${target} PRIVATE
        -Wall -Wextra -Wpedantic
        $<$<CONFIG:Release>:-march=native>)
endforeach()
