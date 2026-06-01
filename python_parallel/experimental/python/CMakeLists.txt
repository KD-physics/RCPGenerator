cmake_minimum_required(VERSION 3.18)

project(rcpgenerator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_PARITY_TESTS "Build legacy parity test executables" OFF)

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

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_rcpgenerator
    src/bindings/module.cpp
)
target_link_libraries(_rcpgenerator PRIVATE rcpgenerator_core)
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

    add_executable(rcpgenerator_skew_periodic_regression
        tests/skew_periodic_regression.cpp
    )
    target_link_libraries(rcpgenerator_skew_periodic_regression PRIVATE rcpgenerator_core)
    target_include_directories(rcpgenerator_skew_periodic_regression PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

    if(MSVC)
        target_compile_options(rcpgenerator_initializer_parity PRIVATE /bigobj)
        target_compile_options(rcpgenerator_packer_parity PRIVATE /bigobj)
        target_compile_options(rcpgenerator_hard_wall_regression PRIVATE /bigobj)
        target_compile_options(rcpgenerator_skew_periodic_regression PRIVATE /bigobj)
    endif()

    add_test(NAME rcpgenerator_initializer_parity COMMAND rcpgenerator_initializer_parity)
    add_test(NAME rcpgenerator_packer_parity COMMAND rcpgenerator_packer_parity)
    add_test(NAME rcpgenerator_hard_wall_regression COMMAND rcpgenerator_hard_wall_regression)
    add_test(NAME rcpgenerator_skew_periodic_regression COMMAND rcpgenerator_skew_periodic_regression)
endif()
