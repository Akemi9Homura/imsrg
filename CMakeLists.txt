########################################################################################################################
### IMSRG LIBRARY MAIN
########################################################################################################################

cmake_minimum_required(VERSION 3.22)

########################################################################################################################
### PROJECT DESCRIPTION AND META INFO
########################################################################################################################

set(IMSRG_PROJECT_NAME        "imsrg++")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
option(IMSRG_BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(USE_ARMADILLO_EXTERN "Use version of armadillo in extern" ON)
########################################################################################################################
### DECLARE PROJECT
########################################################################################################################

project(imsrg++ VERSION 0.1.0)

set(IMSRG_SRC_DIR ${PROJECT_SOURCE_DIR}/src)


find_package(OpenMP REQUIRED)
set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(GSL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Boost REQUIRED)
if (NOT USE_ARMADILLO_EXTERN)
    find_package(Armadillo)
endif()

# libIMSRG is always a SHARED library built from the OBJECT libraries, so
# they must be position independent even without the Python bindings
# (compilers that do not default to PIE fail to link libIMSRG.so otherwise).
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

add_subdirectory(src)
if(IMSRG_BUILD_PYTHON_BINDINGS)
    #ADD tests to be run directly in the build process
    enable_testing()
    add_executable(UnitTestRecouplingCache test/UnitTestRecouplingCache.cc)
    target_link_libraries(UnitTestRecouplingCache PRIVATE IMSRG OpenMP::OpenMP_CXX)
    target_compile_options(UnitTestRecouplingCache PRIVATE -frounding-math)
    add_test(NAME UnitTestRecouplingCache COMMAND UnitTestRecouplingCache)
    add_test (
        NAME UnitTest
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/UnitTest.py 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestCommutators2
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestCommutators2.py   
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test(
        NAME UnitTestCommutators3
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestCommutators3.py
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}   
        
    )
    add_test (
        NAME UnitTestFactorized
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestFactorized.py 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestTriples
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestTriples.py 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestTensor
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestTensor.py  
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestComm122st
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestComm122st.py
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME UnitTestPandya
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/UnitTestPandya.py
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_test (
        NAME MRPrototype
        COMMAND ${PYTHON_EXECUTABLE} -m unittest discover
                -s ${CMAKE_CURRENT_SOURCE_DIR}/prototype/mrimsrg/tests -v
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()
