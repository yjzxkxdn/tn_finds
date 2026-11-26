# SPDX-License-Identifier: GPL-3.0-only
if(NOT WIN32)
    return()
endif()

# This cache contains exact DSP intermediates, so include both the source
# analysis implementation and the toolchain/options that produce its values.
# In a Nebula build, cache_id also identifies the upstream data and estimator.
set(analysis_identity "${cache_id}|${CMAKE_SYSTEM_NAME}|${CMAKE_SYSTEM_PROCESSOR}|${CMAKE_CXX_COMPILER_ID}|${CMAKE_CXX_COMPILER_VERSION}|${CMAKE_CXX_COMPILER_TARGET}|${CMAKE_SIZEOF_VOID_P}|${CMAKE_CXX_FLAGS}|${CMAKE_INTERPROCEDURAL_OPTIMIZATION}|${TN_FINDS_DSP_OPENMP}|${OpenMP_CXX_FLAGS}")
foreach(config IN ITEMS ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE})
    string(TOUPPER "${config}" config)
    string(APPEND analysis_identity "|${config}|${CMAKE_CXX_FLAGS_${config}}|${CMAKE_INTERPROCEDURAL_OPTIMIZATION_${config}}")
endforeach()
set(analysis_inputs src/source_analysis.hpp src/source_analysis_cache.cpp src/cache_win.hpp
    src/engine.cpp src/engine.hpp src/pitch.cpp src/pitch.hpp src/harmonic_mask.cpp src/harmonic_mask.hpp
    src/harmonic_ridges.cpp src/harmonic_ridges.hpp src/glottal.cpp src/glottal.hpp src/harmonic_synth.cpp
    src/harmonic_synth.hpp src/harmonic_phase.cpp src/harmonic_phase.hpp src/dsp_fft.hpp src/fractional_sinc.hpp
    src/note_model.hpp CMakeLists.txt cmake/AnalysisCache.cmake)
if(TARGET tn_finds_nebula)
    list(APPEND analysis_inputs src/nebula_cache.cpp)
endif()
foreach(source IN LISTS analysis_inputs)
    set(path "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${path}")
    file(READ "${path}" text)
    string(REPLACE "\r\n" "\n" text "${text}")
    string(SHA256 sha "${text}")
    string(APPEND analysis_identity "|${sha}")
endforeach()
string(SHA256 analysis_id "${analysis_identity}")
target_sources(tn_finds_engine PRIVATE src/source_analysis_cache.cpp)
target_compile_definitions(tn_finds_engine PRIVATE "TN_FINDS_ANALYSIS_CACHE_ID=\"${analysis_id}-$<CONFIG>\"")
target_link_libraries(tn_finds_engine PRIVATE bcrypt)
