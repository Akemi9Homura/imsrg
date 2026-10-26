cmake_minimum_required(VERSION 3.18)
set(INSTDIR $ENV{HOME})
set(BLA_VENDOR OpenBLAS)
# set(CMAKE_SKIP_INSTALL_RPATH ON)
# set(CMAKE_SYSTEM_PREFIX_PATH $ENV{EBROOTGENTOO})
option(IMSRG_BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

# Place to put any settings for quirks of your specific environment
# e.g. maybe you want a different BLAS, or maybe you need to
# explicitly set some additional variables for whatever reason
include(local_settings.cmake OPTIONAL)

# Determine if imsrg++ is built as a subproject (using add_subdirectory)
# or if it is the main project.
if(NOT DEFINED IMSRGPLUSPLUS_MAIN_PROJECT)
  set(IMSRGPLUSPLUS_MAIN_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(IMSRGPLUSPLUS_MAIN_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

# Set the default CMAKE_BUILD_TYPE to Release.
# This should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if(IMSRGPLUSPLUS_MAIN_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(imsrg++ VERSION 0.1.0)

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)
find_package(ZLIB REQUIRED)
# find_package(Armadillo REQUIRED)
# include_directories(${ARMADILLO_INCLUDE_DIRS})


add_compile_options(-frounding-math)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(IMSRGPLUSPLUS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# libIMSRG is always a SHARED library built from the OBJECT libraries, so
# they must be position independent even without the Python bindings
# (compilers that do not default to PIE fail to link libIMSRG.so otherwise).
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

add_subdirectory(extern/armadillo)
add_subdirectory(extern/half)

add_subdirectory(extern/boost_src)

add_subdirectory(Commutator)
add_subdirectory(ReferenceImplementations)
add_subdirectory(ThreeBodyME)
add_subdirectory(imsrg_util)

# add_subdirectory(profiling)


if(IMSRG_BUILD_PYTHON_BINDINGS)
    add_subdirectory(extern/pybind11)
endif()

# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cc
           ${CMAKE_CURRENT_BINARY_DIR}/_version.cc
    COMMAND ${CMAKE_COMMAND} -P
            ${CMAKE_CURRENT_SOURCE_DIR}/git_version.cmake)


add_library(
    AngMom  OBJECT
    AngMom.hh AngMom.cc
    AngMomCache.hh AngMomCache.cc
)
target_link_libraries(
    AngMom
    PRIVATE
    GSL::gsl
)


add_library(
    IMSRGProfiler OBJECT
    IMSRGProfiler.hh IMSRGProfiler.cc
)
target_link_libraries(
    IMSRGProfiler
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_library(
    ModelSpace  OBJECT
    ModelSpace.hh ModelSpace.cc
    TwoBodyChannel.cc
    ThreeBodyChannel.cc
)
target_link_libraries(
    ModelSpace
    PUBLIC
    armadillo
    IMSRGProfiler
    AngMom
    PRIVATE
    OpenMP::OpenMP_CXX
)

add_library(
    PhysicalConstants
    INTERFACE
    PhysicalConstants.hh
)

add_library(
    TwoBodyME  OBJECT
    TwoBodyME.hh TwoBodyME.cc
)
target_link_libraries(
    TwoBodyME
    PUBLIC
    ModelSpace
    PRIVATE
    AngMom
    PhysicalConstants
)

add_library(
    ThreeLegME OBJECT
    ThreeLegME.hh ThreeLegME.cc
)
target_link_libraries(
    ThreeLegME
    PUBLIC
    armadillo
    ModelSpace
    PRIVATE
    # ${ARMADILLO_LIBRARIES}
    PhysicalConstants
)


#add_library(
#    ThreeBodyME OBJECT
#    ThreeBodyME.hh                           ThreeBodyME_dir/ThreeBodyME.cc
#    ThreeBodyStorage.hh      ThreeBodyME_dir/ThreeBodyStorage.cc
#    ThreeBodyME_dir/ThreeBodyStorage_pn.hh   ThreeBodyME_dir/ThreeBodyStorage_pn.cc
#    ThreeBodyME_dir/ThreeBodyStorage_no2b.hh ThreeBodyME_dir/ThreeBodyStorage_no2b.cc
#    ThreeBodyME_dir/ThreeBodyStorage_mono.hh ThreeBodyME_dir/ThreeBodyStorage_mono.cc
#    ThreeBodyME_dir/ThreeBodyStorage_iso.hh  ThreeBodyME_dir/ThreeBodyStorage_iso.cc
#)
#target_include_directories(
#    ThreeBodyME
#    PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}
#    ${CMAKE_CURRENT_SOURCE_DIR}/ThreeBodyME_dir
#)
#target_link_libraries(
#    ThreeBodyME
#    PUBLIC
#    ModelSpace
#    IMSRGBoostZip
#    HalfPrecision
#    PRIVATE
#    IMSRGProfiler
#    AngMom
#    OpenMP::OpenMP_CXX
#)

add_library(
    Parameters
    INTERFACE
)

add_library(
    Operator  OBJECT
    Operator.hh Operator.cc
)
target_include_directories(
    Operator
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ThreeBodyME
)
target_link_libraries(
    Operator
    PUBLIC
    armadillo
    ModelSpace
    TwoBodyME
    ThreeBodyME
    ThreeLegME
    IMSRGProfiler
    PRIVATE
    GSL::gsl
    OpenMP::OpenMP_CXX
    PhysicalConstants
    AngMom
)

add_library(
    MRReference  OBJECT
    MRReference.hh MRReference.cc
)
target_link_libraries(
    MRReference
    PUBLIC
    Operator
    TwoBodyME
    ModelSpace
)

add_library(
    MRCommutator  OBJECT
    MRCommutator.hh MRCommutator.cc
)
target_link_libraries(
    MRCommutator
    PUBLIC
    MRReference
    Operator
    ModelSpace
    PRIVATE
    Commutator
    AngMom
    OpenMP::OpenMP_CXX
)

add_library(
    Version  OBJECT
    version.hh
    ${CMAKE_CURRENT_BINARY_DIR}/version.cc
)
target_include_directories(
    Version
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(
    ReadWrite  OBJECT
    ReadWrite.hh ReadWrite.cc
)
target_compile_definitions(
    ReadWrite
    PUBLIC
    NO_HDF5="true"
)
target_link_libraries(
    ReadWrite
    PUBLIC
    ModelSpace
    Operator
    Jacobi3BME
    ThreeBodyME
    PRIVATE
    OpenMP::OpenMP_CXX
    IMSRGBoostZip
    AngMom
    PhysicalConstants
    Version
)

add_library(
    BCH  OBJECT
    BCH.hh BCH.cc
)
target_link_libraries(
    BCH
    PUBLIC
    Operator
    PRIVATE
    Commutator
    MRCommutator
    MRReference
    OpenMP::OpenMP_CXX
)

add_library(
    IMSRGSolver  OBJECT
    IMSRGSolver.hh IMSRGSolver.cc
    IMSRGSolverPV.hh IMSRGSolverPV.cc
)
target_link_libraries(
    IMSRGSolver
    PUBLIC
    Operator
    Generator
    IMSRGProfiler
    ReadWrite
    PRIVATE
    BCH
    Commutator
    MRCommutator
    MRReference
    Boost::boost
)

#add_library(
#    IMSRGSolverPV  OBJECT
#    IMSRGSolverPV.hh IMSRGSolverPV.cc
#)
#target_link_libraries(
#    IMSRGSolverPV
#    PUBLIC
#    IMSRGSolver
#    GeneratorPV
#    PRIVATE
#    Commutator
#)

#add_library(
#    Commutator   OBJECT
#    Commutator.hh Commutator.cc
#    DaggerCommutators.hh DaggerCommutators.cc
#    TensorCommutators.hh TensorCommutators.cc
#    IMSRG3Commutators.hh IMSRG3Commutators.cc
#    FactorizedDoubleCommutator.hh FactorizedDoubleCommutator.cc
#)
#target_link_libraries(
#    Commutator
#    PUBLIC
#    armadillo
#    Operator
#    TwoBodyME
#    ThreeLegME
#    PRIVATE
#    OpenMP::OpenMP_CXX
#    # ${ARMADILLO_LIBRARIES}
#    PhysicalConstants
#    AngMom
#)


add_library(
    Jacobi3BME OBJECT
    Jacobi3BME.hh Jacobi3BME.cc
)
target_link_libraries(
    Jacobi3BME
    PUBLIC
    ModelSpace
    PRIVATE
    AngMom
    IMSRGProfiler
    OpenMP::OpenMP_CXX
)

add_library(
    HartreeFock OBJECT
    HartreeFock.hh HartreeFock.cc
)
target_link_libraries(
    HartreeFock
    PUBLIC
    armadillo
    ModelSpace
    Operator
    IMSRGProfiler
    Jacobi3BME
    ThreeBodyME
    PRIVATE
    OpenMP::OpenMP_CXX
    GSL::gsl
    AngMom
    PhysicalConstants
)




#add_library(
#    imsrg_util OBJECT
#    imsrg_util.hh
#    imsrg_util_dir/imsrg_util.cc
#    imsrg_util_dir/GaussLaguerre.hh
#    imsrg_util_dir/M0nu.hh
#    imsrg_util_dir/M0nu.cc
#    imsrg_util_dir/DarkMatterNREFT.hh
#    imsrg_util_dir/DarkMatterNREFT.cc
#    imsrg_util_dir/Pwd.hh
#    imsrg_util_dir/Pwd.cc
#    imsrg_util_dir/Atomic.cc
##    Helicity.hh Helicity.cc
#)
#target_include_directories(
#    imsrg_util
#    PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}
#    ${CMAKE_CURRENT_SOURCE_DIR}/imsrg_util_dir
#)
#target_link_libraries(
#    imsrg_util
#    PUBLIC
#    ModelSpace
#    Operator
#    HartreeFock
#    PhysicalConstants
#    GSL::gsl
#    PRIVATE
#    OpenMP::OpenMP_CXX
#)

add_library(
    Generator OBJECT
    Generator.hh Generator.cc
    GeneratorPV.hh GeneratorPV.cc
)
target_link_libraries(
    Generator
    PUBLIC
    ModelSpace
    Operator
    PRIVATE
    GSL::gsl
    imsrg_util
    AngMom
#    Commutator
    PhysicalConstants
    OpenMP::OpenMP_CXX
)

#add_library(
#    GeneratorPV OBJECT
#    GeneratorPV.hh GeneratorPV.cc
#)
#target_link_libraries(
#    GeneratorPV
#    PUBLIC
#    ModelSpace
#    Operator
#    PRIVATE
#    Generator
#    GSL::gsl
#    imsrg_util
#    AngMom
##    Commutator
#    PhysicalConstants
#    OpenMP::OpenMP_CXX
#)

add_library(
    HFMBPT OBJECT
    HFMBPT.hh HFMBPT.cc
)
target_link_libraries(
    HFMBPT
    PUBLIC
    armadillo
    HartreeFock
    PRIVATE
    PhysicalConstants
    OpenMP::OpenMP_CXX
    Commutator
)



add_library(
    RPA OBJECT
    RPA.hh RPA.cc
)
target_link_libraries(
    RPA
    PUBLIC
    armadillo
    ModelSpace
    Operator
    PRIVATE
    PhysicalConstants
    AngMom
)

#add_library(
#    ReferenceImplementations  OBJECT
#    ReferenceImplementations.hh 
#    ReferenceImplementations.cc
#    HFMBPT4.cc
#    IMSRG2Commutators.cc
#    IMSRG3Commutators.cc
#    IMSRG2TensorCommutators.cc
#    IMSRG3TensorCommutators.cc
#    FactorizedNestedCommutators.cc
#)
#
#target_link_libraries(
#    ReferenceImplementations
#    PUBLIC
#    ModelSpace
#    Operator
##    Commutator
#    PRIVATE
#    armadillo
#    PhysicalConstants
#    Version
#    imsrg_util
#    AngMom
#    OpenMP::OpenMP_CXX
#)

add_library(
    UnitTest  OBJECT
    UnitTest.hh UnitTest.cc
)
target_link_libraries(
    UnitTest
    PUBLIC
    ModelSpace
    Operator
    PRIVATE
    Commutator
    armadillo
    PhysicalConstants
    Version
    imsrg_util
    AngMom
    OpenMP::OpenMP_CXX
    ReferenceImplementations
)

add_library(
    IMSRG
    SHARED
)
target_link_libraries(
    IMSRG
#    INTERFACE
    Version
    ReadWrite
    ModelSpace
    TwoBodyME
    ThreeBodyME
    ThreeLegME
    Operator
    MRReference
    MRCommutator
    HartreeFock
    Commutator
    BCH
    IMSRGSolver
#    IMSRGSolverPV
    Generator
#    GeneratorPV
    imsrg_util
#    M0nu
#    Pwd
    AngMom
    IMSRGProfiler
    Jacobi3BME
#    DarkMatterNREFT
    HFMBPT
    UnitTest
    ReferenceImplementations
    PhysicalConstants
    RPA
)
set_target_properties(IMSRG PROPERTIES SUFFIX ".so")
set_target_properties(IMSRG PROPERTIES LINKER_LANGUAGE CXX)


add_executable(
    imsrg++
    imsrg++.cc
)
target_link_libraries(
    imsrg++
    PUBLIC
    IMSRG
    Parameters
    PhysicalConstants
    Version
    OpenMP::OpenMP_CXX
)

add_executable(
    imsrg_operator_validate
    imsrg_operator_validate.cc
)
target_link_libraries(
    imsrg_operator_validate
    PRIVATE
    IMSRG
    OpenMP::OpenMP_CXX
)


#target_link_libraries(
#    UnitTest
#    PUBLIC
#    IMSRG
#    Parameters
#    PhysicalConstants
#    Version
#    OpenMP::OpenMP_CXX
#)

if(IMSRG_BUILD_PYTHON_BINDINGS)
    pybind11_add_module(pyIMSRG pyIMSRG.cc)
    target_link_libraries(
        pyIMSRG
        PRIVATE
        IMSRG
    )
    set_target_properties(pyIMSRG PROPERTIES SUFFIX ".so")
    set_target_properties(pyIMSRG PROPERTIES POSITION_INDEPENDENT_CODE ON)
    add_custom_command(TARGET pyIMSRG PRE_LINK
		       COMMAND ${CMAKE_COMMAND} -E remove -f $<TARGET_FILE:pyIMSRG>
                       COMMENT "Unlinking old shared library to prevent SIGBUS on running jobs" )
endif()

# Install imsrg executable
install(TARGETS imsrg++ 
        DESTINATION ${INSTDIR}/bin 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE) 

# Install library
install(TARGETS IMSRG  
        DESTINATION ${INSTDIR}/lib 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Install the header files
file(GLOB HEADERS *.hh)
install(FILES ${HEADERS} 
        DESTINATION ${INSTDIR}/include
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

include(CTest)
list(APPEND CMAKE_CTEST_ARGUMENTS "--verbose")
# Unit Tests
add_test(NAME Commutators2 COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestCommutators2.py)
add_test(NAME Commutators3 COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestCommutators3.py)
add_test(NAME Factorized COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestFactorized.py)
add_test(NAME Triples COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestTriples.py)
add_test(NAME Tensor COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestTensor.py)
add_test(NAME NormalOrdering COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestNormalOrdering.py)
add_test(NAME MRReference COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRReference.py)
add_test(NAME MRGuards
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRGuards.py)
add_test(NAME MRTensor
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRTensor.py)
add_test(NAME MRReferencePairOrder
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRReferencePairOrder.py)
add_test(NAME MRDriver COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRDriver.py
                               --imsrg-executable $<TARGET_FILE:imsrg++>)
add_test(NAME MRSRDriver COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRSRDriver.py
                                 $<TARGET_FILE:imsrg++>)
add_test(NAME MRCorrelatedDriver
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRCorrelatedDriver.py
                 $<TARGET_FILE:imsrg++>)
add_test(NAME MRDenominators
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRDenominators.py)
add_test(NAME MRMagnus
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRMagnus.py)
add_test(NAME MRJobGenerator
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRJobGenerator.py
                 $<TARGET_FILE:imsrg++>)
add_test(NAME MRTailDiagnostic
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRTailDiagnostic.py)
add_test(NAME MRDownstreamWindow
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRDownstreamWindow.py)
add_test(NAME MRMagnusGate
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/../test/UnitTestMRMagnusGate.py)

# Python itself is not built with the sanitizers.  When an instrumented C++
# extension is loaded later, ASan cannot install its __cxa_throw interceptor
# unless libstdc++ is already present.  Keep libasan first, preload libstdc++ as
# well.  LeakSanitizer is disabled for these Python-hosted tests because the
# non-instrumented CPython interpreter retains process-global allocations;
# address/undefined-behaviour checks remain fatal.  Explicitly
# point every MR Python test at the module from this build tree so a sanitizer
# build cannot silently import the normal Release module.
set(MR_PYTHON_TESTS
    MRReference MRReferencePairOrder MRGuards MRTensor MRDriver MRSRDriver MRCorrelatedDriver MRDenominators MRMagnus
    MRJobGenerator MRTailDiagnostic MRDownstreamWindow MRMagnusGate)
set(MR_PYTHON_TEST_ENVIRONMENT
    "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}"
    "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH}")
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
set(MR_EFFECTIVE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
if(MR_EFFECTIVE_CXX_FLAGS MATCHES "-fsanitize=[^ ]*address")
  execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libasan.so
      OUTPUT_VARIABLE MR_ASAN_RUNTIME
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so.6
      OUTPUT_VARIABLE MR_LIBSTDCXX_RUNTIME
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT EXISTS "${MR_ASAN_RUNTIME}" OR NOT EXISTS "${MR_LIBSTDCXX_RUNTIME}")
    message(FATAL_ERROR "Cannot locate ASan/libstdc++ runtimes for MR Python tests")
  endif()
  list(APPEND MR_PYTHON_TEST_ENVIRONMENT
       "LD_PRELOAD=${MR_ASAN_RUNTIME}:${MR_LIBSTDCXX_RUNTIME}"
       "ASAN_OPTIONS=detect_leaks=0:halt_on_error=1:abort_on_error=1"
       "UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1")
endif()
set_tests_properties(${MR_PYTHON_TESTS}
                     PROPERTIES ENVIRONMENT "${MR_PYTHON_TEST_ENVIRONMENT}")
