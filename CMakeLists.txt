cmake_minimum_required(VERSION 3.20)

project(Deep3DRadon
    VERSION 1.0
    DESCRIPTION "A project for 3D Radon transform using deep learning"
    LANGUAGES CXX
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find required packages
find_package(PNG REQUIRED)

# Add source files
set(SOURCES
    source/src/main.cpp
    source/src/orientation.cpp
    source/src/dev_test.cpp
    source/src/detector.cpp
    source/src/microscope.cpp
    source/src/crystal.cpp
    source/src/pattern.cpp
    source/src/radon_transform.cpp
    source/src/k_lines.cpp
    source/src/projector.cpp
)

# Add the executable
add_executable(d3r ${SOURCES})

# Set include directories
target_include_directories(d3r PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/source/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/source/thirdParty/Eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/source/thirdParty/CImg
    ${PNG_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(d3r PRIVATE PNG::PNG)

# Compiler-specific flags for better warnings and optimization
if(MSVC)
    target_compile_options(d3r PRIVATE /W4)
else()
    target_compile_options(d3r PRIVATE -Wall -Wextra -Wpedantic -Dcimg_display=0)
endif()


# Enable testing (optional)
#enable_testing()
