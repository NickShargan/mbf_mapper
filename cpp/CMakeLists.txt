cmake_minimum_required(VERSION 3.16)
project(MbfCli)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(Eigen3 REQUIRED NO_MODULE)

set(SOURCES
    src/main.cpp
    src/utils.cpp
    src/ioDicom.cpp
    src/ioMasks.cpp
    src/mbfPipeline.cpp
    src/MbfDeconvolver.cpp
)

add_executable(mbf_cli ${SOURCES})

target_include_directories(mbf_cli
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(mbf_cli
    PRIVATE
        ${ITK_LIBRARIES}
        Eigen3::Eigen
)

