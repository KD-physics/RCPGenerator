cmake_minimum_required(VERSION 3.18)

project(rcpgenerator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_PARITY_TESTS "Build legacy parity test executables" OFF)

# OpenMP for parallel force loop (Phase 1 parallelization)
find_package(OpenMP REQUIRED)

add_library(rcpgenerator_core
    src/core/initialize_particles.cpp
    src/core/rcp_generator.cpp
)
set_target_properties(rcpgenerator_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(rcpgenerator_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(rcpgenerator_core PUBLIC OpenMP::OpenMP_CXX)

# Cycle 5c: enable native-target tuning so gcc can emit AVX2 (and any
# other ISA the host supports) for the hot loops. The previous default
# was SSE2-only (x86-64 baseline), which capped forces.loop scaling well
# below the i9-12900K's AVX2 capability.
if(NOT MSVC)
    # Cycle 8 H16: add LTO for cross-translation-unit inlining (range_query
    # is called across the dispatch boundary). -ffast-math is NOT enabled
    # globally because the MIC test must match sort+walk byte-equivalently;
    # selective fast-math is applied via __attribute__((optimize)) on the
    # inner FP loops where reassociation is safe.
    target_compile_options(rcpgenerator_core PRIVATE
        -march=native -funroll-loops -flto)
    target_link_options(rcpgenerator_core PRIVATE -flto)
endif()

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_rcpgenerator
    src/bindings/module.cpp
)
target_link_libraries(_rcpgenerator PRIVATE rcpgenerator_core OpenMP::OpenMP_CXX)
target_include_directories(_rcpgenerator PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

install(TARGETS _rcpgenerator DESTINATION rcpgenerator)
install(
    DIRECTORY python/rcpgenerator/
    DESTINATION rcpgenerator
    FILES_MATCHING
    PATTERN "*.py"
)

if(MSVC)
    target_compile_options(_rcpgenerator PRIVATE /bigobj)
endif()

if(BUILD_PARITY_TESTS)
    enable_testing()

    add_executable(rcpgenerator_initializer_parity
        tests/initializer_parity.cpp
    )
    target_link_libraries(rcpgenerator_initializer_parity PRIVATE rcpgenerator_core)
    target_include_directories(rcpgenerator_initializer_parity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

    add_executable(rcpgenerator_packer_parity
        tests/packer_parity.cpp
    )
    target_link_libraries(rcpgenerator_packer_parity PRIVATE rcpgenerator_core)
    target_include_directories(rcpgenerator_packer_parity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

    add_executable(rcpgenerator_hard_wall_regression
        tests/hard_wall_regression.cpp
    )
    target_link_libraries(rcpgenerator_hard_wall_regression PRIVATE rcpgenerator_core)
    target_include_directories(rcpgenerator_hard_wall_regression PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

    if(MSVC)
        target_compile_options(rcpgenerator_initializer_parity PRIVATE /bigobj)
        target_compile_options(rcpgenerator_packer_parity PRIVATE /bigobj)
        target_compile_options(rcpgenerator_hard_wall_regression PRIVATE /bigobj)
    endif()

    add_test(NAME rcpgenerator_initializer_parity COMMAND rcpgenerator_initializer_parity)
    add_test(NAME rcpgenerator_packer_parity COMMAND rcpgenerator_packer_parity)
    add_test(NAME rcpgenerator_hard_wall_regression COMMAND rcpgenerator_hard_wall_regression)
endif()
