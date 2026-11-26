# SPDX-License-Identifier: GPL-3.0-only
cmake_minimum_required(VERSION 3.16)
project(tn_finds VERSION 0.4.9 LANGUAGES CXX)

add_library(tn_finds_engine STATIC src/engine.cpp src/pitch.cpp src/formant.cpp src/glottal.cpp
    src/harmonic_phase.cpp src/harmonic_mask.cpp src/harmonic_ridges.cpp src/harmonic_synth.cpp)
include(cmake/Nebula.cmake)
set(dsp_openmp_default OFF)
if(TARGET tn_finds_nebula)
    set(dsp_openmp_default ON)
endif()
option(TN_FINDS_DSP_OPENMP "Parallelize independent ridge frames of long OTO clips" ${dsp_openmp_default})
if(TN_FINDS_DSP_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_link_libraries(tn_finds_engine PRIVATE OpenMP::OpenMP_CXX)
    if(MINGW)
        target_link_options(tn_finds_engine INTERFACE ${OpenMP_CXX_FLAGS})
    endif()
endif()
include(cmake/AnalysisCache.cmake)
set(cli_sources src/main.cpp src/wav.cpp src/mrq.cpp src/mrq_binding.cpp)
add_executable(tn_finds_resampler ${cli_sources})
add_executable(tn_finds ${cli_sources} src/wavtool.cpp src/bat_host.cpp src/phrase.cpp src/phrase_io.cpp)
if(WIN32)
    target_link_libraries(tn_finds_resampler PRIVATE bcrypt)
    target_link_libraries(tn_finds PRIVATE bcrypt)
endif()
target_compile_definitions(tn_finds PRIVATE TN_FINDS_DUAL_ROLE=1)
set(TN_FINDS_VERSION "${PROJECT_VERSION}" CACHE STRING "Version label reported by the executable")
foreach(target IN ITEMS tn_finds_resampler tn_finds)
    target_compile_definitions(${target} PRIVATE "TN_FINDS_VERSION=\"${TN_FINDS_VERSION}\"")
endforeach()
target_link_libraries(tn_finds_resampler PRIVATE tn_finds_engine)
target_link_libraries(tn_finds PRIVATE tn_finds_engine)
foreach(target IN ITEMS tn_finds_engine tn_finds_resampler tn_finds)
    target_compile_features(${target} PRIVATE cxx_std_17)
    set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /utf-8)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endforeach()
target_include_directories(tn_finds_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
if(MINGW)
    target_link_options(tn_finds PRIVATE -municode -static)
    target_link_options(tn_finds_resampler PRIVATE -municode -static)
endif()
install(TARGETS tn_finds tn_finds_resampler RUNTIME DESTINATION bin)

option(TN_FINDS_TESTS "Register the Python audio regression tests (NumPy and SoundFile required)" OFF)
if(TN_FINDS_TESTS)
    enable_testing()
    add_executable(tn_finds_fft_test tests/test_fft.cpp)
    target_compile_features(tn_finds_fft_test PRIVATE cxx_std_17)
    target_link_libraries(tn_finds_fft_test PRIVATE tn_finds_engine)
    add_test(NAME fft_regression COMMAND tn_finds_fft_test)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_test(NAME audio_regression COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_audio.py --exe $<TARGET_FILE:tn_finds_resampler>)
    add_test(NAME harmonic_flag_regression COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_harmonic_flag.py --exe $<TARGET_FILE:tn_finds_resampler>)
    add_test(NAME voice_flag_regression COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_voice_flags.py --exe $<TARGET_FILE:tn_finds_resampler>)
    add_test(NAME wavtool_regression COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wavtool.py --exe $<TARGET_FILE:tn_finds>
        --resampler $<TARGET_FILE:tn_finds_resampler>)
    add_test(NAME phrase_regression COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_phrase.py --exe $<TARGET_FILE:tn_finds>)
    if(TARGET tn_finds_nebula)
        add_test(NAME nebula_f0_regression COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_nebula_f0.py $<TARGET_FILE:tn_finds_resampler>)
        if(WIN32)
            add_test(NAME nebula_cache_regression COMMAND ${Python3_EXECUTABLE} -B
                ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_nebula_cache.py $<TARGET_FILE:tn_finds_resampler>)
        endif()
    endif()
    if(WIN32)
        # No test uses the user's cache or an earlier test's results.
        get_property(regressions DIRECTORY PROPERTY TESTS)
        set_tests_properties(${regressions} PROPERTIES ENVIRONMENT "TN_FINDS_NEBULA_CACHE=off;TN_FINDS_ANALYSIS_CACHE=off;TN_FINDS_MRQ_BINDINGS=off")
        add_test(NAME source_analysis_cache_regression COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_source_analysis_cache.py
            $<TARGET_FILE:tn_finds_resampler> $<TARGET_FILE:tn_finds>)
    endif()
    add_test(NAME mrq_source_regression COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mrq.py
        $<TARGET_FILE:tn_finds_resampler> $<TARGET_FILE:tn_finds>)
endif()
