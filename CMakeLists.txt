cmake_minimum_required(VERSION 3.10)

set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")
if (WIN32)
    set(VCPKG_TARGET_TRIPLET x64-windows)
endif()

project(haplotag)
find_package(Threads)
find_package(nlohmann_json CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

if (CMAKE_CXX_STANDARD LESS 23)
    add_subdirectory(function2)
    set(FU2 function2)
endif()

set(ZLIBPORT zlib-ng)
if (ZLIBPORT STREQUAL "zlib")
    add_subdirectory(zlib)
    include_directories(${CMAKE_SOURCE_DIR}/zlib)
else()
    # zlib-ng options
    set(ZLIB_COMPAT ON)
    set(WITH_NATIVE_INSTRUCTIONS ON)
    set(ZLIB_ENABLE_TESTS OFF)
    set(ZLIBNG_ENABLE_TESTS OFF)
    set(WITH_GTEST OFF)
    add_subdirectory(zlib-ng)
    include_directories(${CMAKE_SOURCE_DIR}/zlib-ng)
endif()

if(MSVC)
    add_definitions(/DNOMINMAX /DWIN32_LEAN_AND_MEAN /D_CRT_SECURE_NO_WARNINGS /DFLUA_USE_WINDOWS /DFASTQ_WIN32)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    add_definitions(-DLUA_USE_MACOSX -DFLUA_USE_MACOS -DFASTQ_DARWIN)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    add_definitions(-DLUA_USE_LINUX -DFLUA_USE_LINUX -DFASTQ_LINUX)
endif()

# cuda
find_package(CUDAToolkit)
if (${CUDAToolkit_FOUND})
    set(FASTQ_HAS_CUDA 1)
endif()


# original binaries
macro(add_haplo_demult name)
    add_executable(${name}
        ${CMAKE_SOURCE_DIR}/haplo_demult/${name}.cpp
        ${CMAKE_SOURCE_DIR}/gzstream/gzstream.cpp
    )
    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/gzstream)
    target_link_libraries(${name} PRIVATE zlibstatic)
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin$<0:>)
endmacro()
add_haplo_demult(demult_fastq)
add_haplo_demult(H4_demult_fastq_with_clipping_8bp-plateBC)
add_haplo_demult(H4_demult_fastq_with_clipping_noPlateBC)
add_haplo_demult(H4_demult_fastq_with_clipping_7bp-plateBC)


# fastq binaries
macro(add_fastq name)
    add_executable(${name} ${CMAKE_SOURCE_DIR}/src/${name}.cpp)
    target_compile_definitions(${name} PRIVATE FASTQ_HAS_CUDA=${FASTQ_HAS_CUDA})
    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR})
    target_link_libraries(${name} PRIVATE Threads::Threads zlibstatic ${FU2} nlohmann_json::nlohmann_json)
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin$<0:>)
endmacro()
add_fastq(fastq_barcode)
add_fastq(fastq_gen)
add_fastq(fastq_write)
add_fastq(fastq_fuzzy)
add_fastq(fastq_read)
add_fastq(fastq_splitter)
add_fastq(fastq_H4)
add_fastq(fastq_cat)

if(NOT WIN32)
    add_executable(pigz
        ${CMAKE_SOURCE_DIR}/pigz/try.c
        ${CMAKE_SOURCE_DIR}/pigz/yarn.c
        ${CMAKE_SOURCE_DIR}/pigz/pigz.c
    )
    target_compile_definitions(pigz PUBLIC -DNOZOPFLI)
    target_link_libraries(pigz PRIVATE Threads::Threads zlibstatic)

    # move to a common destination
    install(TARGETS pigz RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin)
    set_target_properties(pigz PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin$<0:>)
endif()


if (${CUDAToolkit_FOUND})
    add_subdirectory(fastq_cuda)
else()
    message(STATUS "CUDAToolkit not found. Skipping CUDA builds")
endif()
