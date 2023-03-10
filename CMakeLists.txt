cmake_minimum_required(VERSION 3.15)

include(FetchContent)

# ----------------------------- Eigen library ----------------------------------
FetchContent_Declare(
    Eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    )
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

# ------------------------------ fmt library -----------------------------------
FetchContent_Declare(
    fmt
    URL https://github.com/fmtlib/fmt/archive/refs/tags/9.1.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    )

FetchContent_MakeAvailable(Eigen)
FetchContent_MakeAvailable(fmt)

# ------------------------------------------------------------------------------
project(advection-equation-1d CXX)

add_library(cfd
    INTERFACE
        include/cfd/periodic_boundary.hpp
        include/cfd/problem_parameters.hpp
        include/cfd/riemann_solvers.hpp
        include/cfd/slope_limiters.hpp
        include/cfd/spacial_reconstruction_schemes.hpp
        include/cfd/time_integration_schemes.hpp
        include/cfd/scalar_advection_equation_simulator.hpp
        include/cfd/cfd.hpp
        include/cfd/text_file_writer.hpp
    )

target_include_directories(cfd INTERFACE include/)
target_link_libraries(cfd
    INTERFACE
        Eigen3::Eigen
        fmt::fmt
    )
target_compile_features(cfd INTERFACE cxx_std_17)
target_compile_options(cfd
    INTERFACE
        $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
        $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/W3>
        $<$<CXX_COMPILER_ID:Intel>:$<IF:$<PLATFORM_ID:Windows>,/W3,-w3>>
        $<$<CXX_COMPILER_ID:IntelLLVM>:$<IF:$<PLATFORM_ID:Windows>,/W3,-w3>>
    )
target_compile_definitions(cfd
    INTERFACE
        $<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES NOMINMAX>
    )

function(add_simulator name)
    add_executable(${name}
        src/${name}.cpp
        src/common.cpp
        )
    target_include_directories(${name} PRIVATE src/)
    target_link_libraries(${name} PRIVATE cfd)
endfunction()

add_simulator(first_order_upwind)
add_simulator(lax_wendroff)
add_simulator(beam_warming)
add_simulator(fromm)
add_simulator(tvd_minmod)
add_simulator(tvd_superbee)
add_simulator(tvd_van_leer)
add_simulator(tvd_van_albada)