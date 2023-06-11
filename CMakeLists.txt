cmake_minimum_required( VERSION 3.26 )

project( T-AHH )

set( CMAKE_CXX_STANDARD 17 )

add_executable( ${PROJECT_NAME}
    src/main.cpp
    src/Potential.cpp
)

target_include_directories( ${PROJECT_NAME}
    PUBLIC
        ${CMAKE_SOURCE_DIR}/src
)