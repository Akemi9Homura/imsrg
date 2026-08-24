cmake_minimum_required(VERSION 3.16)

project(mrimsrg_input_bridge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SHELL_MODEL_OBS_ROOT "/home/mengziyan/shell-model-obs" CACHE PATH
    "Path to the validated shell-model-obs checkout")

if(NOT EXISTS "${SHELL_MODEL_OBS_ROOT}/CMakeLists.txt")
    message(FATAL_ERROR "shell-model-obs checkout not found: ${SHELL_MODEL_OBS_ROOT}")
endif()

# The bridge needs only the validated NN reader and simpleFCI path.  Disabling
# optional 3B readers keeps this small build independent of unused libraries.
set(TBME "NONE" CACHE STRING "shell-model-obs optional 3B reader" FORCE)
set(SIMPLE_FCI_EMAX "2" CACHE STRING "simpleFCI determinant emax" FORCE)
add_subdirectory("${SHELL_MODEL_OBS_ROOT}" "${CMAKE_BINARY_DIR}/shell-model-obs" EXCLUDE_FROM_ALL)

execute_process(
    COMMAND git -C "${SHELL_MODEL_OBS_ROOT}" rev-parse HEAD
    OUTPUT_VARIABLE SHELL_MODEL_OBS_REVISION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE SHELL_MODEL_OBS_GIT_RESULT)
if(NOT SHELL_MODEL_OBS_GIT_RESULT EQUAL 0)
    set(SHELL_MODEL_OBS_REVISION "unknown")
endif()

add_executable(mrimsrg_prepare cpp/prepare_reference.cpp)
target_link_libraries(mrimsrg_prepare PRIVATE fci nucleus util)
target_compile_definitions(mrimsrg_prepare PRIVATE
    SIMPLE_FCI_EMAX=2
    SHELL_MODEL_OBS_ROOT="${SHELL_MODEL_OBS_ROOT}"
    SHELL_MODEL_OBS_REVISION="${SHELL_MODEL_OBS_REVISION}")
target_compile_options(mrimsrg_prepare PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-march=native>)
