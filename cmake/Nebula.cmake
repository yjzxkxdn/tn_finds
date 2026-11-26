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
# Fingerprint the actual source tree, not a Git label: users can supply a
# tarball, a dirty checkout or a different calibration. Normalize CRLF.
file(GLOB_RECURSE nebula_inputs CONFIGURE_DEPENDS LIST_DIRECTORIES false
     RELATIVE "${nebula_root}" "${nebula_root}/nebula/*" "${nebula_root}/external/ciglet/*"
     "${nebula_root}/external/libgvps/*" "${nebula_root}/external/fastapprox/*")
list(FILTER nebula_inputs INCLUDE REGEX "\\.(c|h|hc|inc)$")
list(SORT nebula_inputs)
list(APPEND nebula_inputs logging.h)
set(identity "")
foreach(source IN LISTS nebula_inputs)
    set(path "${nebula_root}/${source}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${path}")
    file(READ "${path}" text)
    string(REPLACE "\r\n" "\n" text "${text}")
    string(SHA256 sha "${text}")
    string(APPEND identity "|${sha}")
endforeach()
string(SHA256 nebula_upstream_id "${identity}")
set(nebula_pinned OFF)
if(nebula_upstream_id STREQUAL "cd3888f56e1309233fd6250a865bcc87306e1894e90e54ff8b5172b71a997b72")
    set(nebula_pinned ON)
endif()
include("${CMAKE_CURRENT_LIST_DIR}/NebulaCold.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/NebulaKernels.cmake")
# Keep MinGW object paths short and ASCII for Unicode source locations.
foreach(source IN LISTS nebula_sources)
    get_filename_component(name "${source}" NAME)
    get_filename_component(stem "${source}" NAME_WE)
    set(unit "${CMAKE_CURRENT_BINARY_DIR}/nebula_units/${name}")
    if(DEFINED ${stem}_override)
        file(GENERATE OUTPUT "${unit}" CONTENT "${${stem}_override}")
    else()
        file(GENERATE OUTPUT "${unit}" CONTENT "#include \"${nebula_root}/${source}\"\n")
    endif()
    list(APPEND nebula_units "${unit}")
endforeach()
add_library(tn_finds_nebula STATIC ${nebula_units} src/nebula_log.c)
target_include_directories(tn_finds_nebula PRIVATE "${nebula_root}/nebula" "${nebula_root}/external/ciglet"
    "${CMAKE_CURRENT_SOURCE_DIR}/src")
set_target_properties(tn_finds_nebula PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
target_compile_definitions(tn_finds_nebula PRIVATE FP_TYPE=float)
target_compile_options(tn_finds_nebula PRIVATE -Ofast)
# The band convolutions and the GMM density are plain float loops, so wider
# vectors double their throughput. x86-64-v3 (AVX2 + FMA, 2013 and later) is
# opt-in: the produced executable then requires such a CPU.
option(TN_FINDS_NEBULA_FAST_ISA "Build the Nebula kernels for x86-64-v3 (AVX2/FMA)" OFF)
if(TN_FINDS_NEBULA_FAST_ISA AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    target_compile_options(tn_finds_nebula PRIVATE -march=x86-64-v3)
    message(STATUS "Nebula: building kernels with -march=x86-64-v3")
endif()
target_link_libraries(tn_finds_nebula PRIVATE OpenMP::OpenMP_C)
# MinGW's OpenMP target can omit the driver flag at final executable link.
target_link_options(tn_finds_nebula INTERFACE ${OpenMP_C_FLAGS})

target_sources(tn_finds_engine PRIVATE src/nebula_pitch.cpp)
target_include_directories(tn_finds_engine PRIVATE "${nebula_root}/nebula")
target_compile_definitions(tn_finds_engine PRIVATE TN_FINDS_NEBULA_F0=1)
target_link_libraries(tn_finds_engine PRIVATE tn_finds_nebula)

if(WIN32)
    # Changing the model, calibration, DSP source, compiler or relevant build
    # options changes the key namespace; OTO/options are not estimator inputs.
    set(cache_inputs "${CMAKE_CURRENT_SOURCE_DIR}/src/nebula_pitch.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/nebula_log.c"
        "${CMAKE_CURRENT_LIST_DIR}/NebulaCold.cmake" "${CMAKE_CURRENT_LIST_DIR}/NebulaKernels.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/nebula_convolution.h" "${CMAKE_CURRENT_LIST_FILE}")
    set(identity "${nebula_upstream_id}|${CMAKE_SYSTEM_NAME}|${CMAKE_SYSTEM_PROCESSOR}|${CMAKE_C_COMPILER_ID}|${CMAKE_C_COMPILER_VERSION}|${CMAKE_CXX_COMPILER_VERSION}|${CMAKE_C_COMPILER_TARGET}|${CMAKE_SIZEOF_VOID_P}|${CMAKE_C_FLAGS}|${CMAKE_CXX_FLAGS}|${OpenMP_C_FLAGS}|${TN_FINDS_NEBULA_FAST_COLD}|${TN_FINDS_NEBULA_FAST_KERNELS}")
    foreach(config IN ITEMS ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE})
        string(TOUPPER "${config}" config)
        string(APPEND identity "|${config}|${CMAKE_C_FLAGS_${config}}|${CMAKE_CXX_FLAGS_${config}}")
    endforeach()
    foreach(source IN LISTS cache_inputs)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${source}")
        file(SHA256 "${source}" sha)
        string(APPEND identity "|${sha}")
    endforeach()
    string(SHA256 cache_id "${identity}")
    target_sources(tn_finds_engine PRIVATE src/nebula_cache.cpp)
    target_compile_definitions(tn_finds_engine PRIVATE "TN_FINDS_NEBULA_CACHE_ID=\"${cache_id}\"")
    target_link_libraries(tn_finds_engine PRIVATE bcrypt)
endif()
