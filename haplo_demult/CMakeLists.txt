cmake_minimum_required(VERSION 3.10)

set(CMAKE_TOOLCHAIN_FILE "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake")
if (WIN32)
    set(VCPKG_TARGET_TRIPLET x64-windows)
endif()

project(haplo_demult)

find_package(ZLIB REQUIRED)

macro(add_haplo_demult name)
    add_executable(${name}
        ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cpp
        ${CMAKE_SOURCE_DIR}/gzstream/gzstream.cpp
    )
    target_include_directories(${name} 
        PRIVATE ${CMAKE_SOURCE_DIR}/gzstream 
        PRIVATE ${CMAKE_SOURCE_DIR}/zlib
    )
    target_link_libraries(${name} PRIVATE ZLIB::ZLIB)
    set_target_properties(${name} PROPERTIES 
        CXX_STANDARD 11
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin$<0:>
    )
endmacro()

add_haplo_demult(H4_demult_fastq_with_clipping_7bp-plateBC)
add_haplo_demult(H4_demult_fastq_with_clipping_8bp-plateBC)
add_haplo_demult(H4_demult_fastq_with_clipping_noPlateBC)
