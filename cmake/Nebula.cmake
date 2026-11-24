# SPDX-License-Identifier: GPL-3.0-only
set(TN_FINDS_OPENMORESAMPLER_SOURCE "" CACHE PATH "OpenMoresampler source tree for the Nebula F0 backend")
if(NOT TN_FINDS_OPENMORESAMPLER_SOURCE)
    return()
endif()

set(nebula_root "${TN_FINDS_OPENMORESAMPLER_SOURCE}")
if(NOT EXISTS "${nebula_root}/nebula/model/nebula_gmm.inc" OR
   NOT EXISTS "${nebula_root}/nebula/model/nebula_lcal.inc")
    message(FATAL_ERROR "Select the OpenMoresampler source tree containing the Nebula model tables")
endif()
enable_language(C)
if(NOT CMAKE_C_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "The Nebula backend requires the GNU C toolchain")
endif()
find_package(OpenMP REQUIRED COMPONENTS C)
set(nebula_sources
    nebula/nebula_wave.c
    nebula/nebula_f0.c
    nebula/nebula_features.c
    nebula/nebula_model.c
    external/ciglet/ciglet.c
    external/ciglet/external/fftsg_h.c
    external/ciglet/external/fast_median.c
    external/libgvps/gvps_obsrv.c
    external/libgvps/gvps_sampled.c
)
# Keep MinGW object paths short and ASCII for Unicode source locations.
foreach(source IN LISTS nebula_sources)
    get_filename_component(name "${source}" NAME)
    set(unit "${CMAKE_CURRENT_BINARY_DIR}/nebula_units/${name}")
    file(GENERATE OUTPUT "${unit}" CONTENT "#include \"${nebula_root}/${source}\"\n")
    list(APPEND nebula_units "${unit}")
endforeach()
add_library(tn_finds_nebula STATIC ${nebula_units} src/nebula_log.c)
set_target_properties(tn_finds_nebula PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
target_compile_definitions(tn_finds_nebula PRIVATE FP_TYPE=float)
target_compile_options(tn_finds_nebula PRIVATE -Ofast)
target_link_libraries(tn_finds_nebula PRIVATE OpenMP::OpenMP_C)
# MinGW's OpenMP target can omit the driver flag at final executable link.
target_link_options(tn_finds_nebula INTERFACE ${OpenMP_C_FLAGS})

target_sources(tn_finds_engine PRIVATE src/nebula_pitch.cpp)
target_include_directories(tn_finds_engine PRIVATE "${nebula_root}/nebula")
target_compile_definitions(tn_finds_engine PRIVATE TN_FINDS_NEBULA_F0=1)
target_link_libraries(tn_finds_engine PRIVATE tn_finds_nebula)
