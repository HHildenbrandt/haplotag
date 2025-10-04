cmake_minimum_required(VERSION 3.10)

set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")
if (WIN32)
    set(VCPKG_TARGET_TRIPLET x64-windows)
endif()

project(haplotag)
find_package(Threads REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

# zlib-ng options
set(ZLIB_COMPAT OFF)
set(WITH_GZFILEOP ON)
set(ZLIB_ENABLE_TESTS OFF)
set(ZLIBNG_ENABLE_TESTS OFF)
set(WITH_GTEST OFF)
add_subdirectory(zlib-ng)


if(MSVC)
    add_definitions(/DNOMINMAX /DWIN32_LEAN_AND_MEAN /D_CRT_SECURE_NO_WARNINGS /DFASTQ_WIN32)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    add_definitions(-DFASTQ_DARWIN)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    add_definitions(-DFASTQ_LINUX)
endif()

# fastq binaries
macro(add_fastq name)
    add_executable(${name} ${CMAKE_SOURCE_DIR}/src/${name}.cpp)
    target_include_directories(${name} 
        PRIVATE ${CMAKE_SOURCE_DIR}
        PRIVATE ${CMAKE_SOURCE_DIR}/zlib-ng
    )
    target_link_libraries(${name} PRIVATE Threads::Threads zlibstatic nlohmann_json::nlohmann_json)
    set_target_properties(${name} PROPERTIES 
        CXX_STANDARD 23
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin$<0:>
    )
endmacro()

add_fastq(fastq_h4)
add_fastq(fastq_cat)
add_fastq(fastq_paste)


# original code with minor changes
add_subdirectory(haplo_demult)
