cmake_minimum_required(VERSION 3.12)
project(Poisson2DSolver VERSION 1.0.0 LANGUAGES CXX)

# set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# compile options
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

# options configuration
option(USE_SPARSE_MATRIX "Use sparse matrix for linear system" ON)
option(USE_QUAD_ELEMENTS "Use quadrilateral elements" ON)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_DOCS "Build documentation" ON)

if(USE_SPARSE_MATRIX)
    add_compile_definitions(USE_SPARSE_MATRIX)
endif()

# find dependencies
find_package(Eigen3 3.3 REQUIRED)
find_package(Boost 1.70 REQUIRED)  # Boost.Math is header-only, no components needed

# Google Test
if(BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    # Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# main library
add_library(Poisson2D_core SHARED
    Src/core/Mesh.cpp
    Src/core/FiniteElement.cpp
    Src/core/Solver.cpp
    Src/core/LinearSystem.cpp
    Src/core/Config.cpp
    Src/elements/TriangleElement.cpp
    Src/elements/QuadElement.cpp
    Src/utils/Parser.cpp
    # Src/utils/Integration.cpp  # Not used - Boost called directly in elements
    Src/utils/VTKWriter.cpp
)

target_include_directories(Poisson2D_core PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Src>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(Poisson2D_core PUBLIC 
    Eigen3::Eigen
    Boost::boost  # Header-only interface target
)

# main executable
add_executable(Poisson2DSolver Src/main.cpp)
target_link_libraries(Poisson2DSolver Poisson2D_core)

# install configuration
install(TARGETS Poisson2DSolver Poisson2D_core
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY Src/ DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)

# tests
if(BUILD_TESTS)
    add_executable(Poisson2D_tests
        tests/test_main.cpp
        tests/test_analytical.cpp
        tests/test_mesh.cpp
        tests/test_elements.cpp
    )
    
    target_link_libraries(Poisson2D_tests 
        Poisson2D_core
        GTest::gtest_main
    )
    
    include(GoogleTest)
    gtest_discover_tests(Poisson2D_tests)
endif()

# documentation generation
if(BUILD_DOCS AND DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()

# performance test script                          
add_custom_target(performance_test
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Poisson2DSolver ${CMAKE_CURRENT_SOURCE_DIR}/examples/performance_test.json
    DEPENDS Poisson2DSolver
    COMMENT "Running performance tests"
)