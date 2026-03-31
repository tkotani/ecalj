# Make libecaljF.so and all binaries.
#
# Usage:
#   FC=<compiler> cmake .
#   FC=<compiler> cmake . -DCMAKE_BUILD_TYPE=Debug
#   FC=nvfortran cmake . -DBUILD_GPU=ON -DBUILD_MP=ON -DBUILD_MP_GPU=ON
#
# Examples:
#   FC=gfortran cmake .
#   FC=ifort cmake .
#   FC=nvfortran cmake . -DBUILD_GPU=ON
#
# To see compile options:
#   VERBOSE=1 make
#
# To start over:
#   rm -rf CMakeCache.txt CMakeFiles/
#
cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED ENV{FC})
  message(FATAL_ERROR "Fortran compiler must be set via FC environment variable. e.g., FC=gfortran cmake .")
endif()
set(FC $ENV{FC})

project(ecalj LANGUAGES Fortran)

# --- Build Type and Options ---
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type: Debug or Release" FORCE)
endif()
option(BUILD_GPU "Build GPU version (for nvfortran)" OFF)
option(BUILD_MP "Build mixed-precision CPU version" OFF)
option(BUILD_MP_GPU "Build mixed-precision GPU version (for nvfortran)" OFF)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message("BUILD_GPU: ${BUILD_GPU}")
message("BUILD_MP: ${BUILD_MP}")
message("BUILD_MP_GPU: ${BUILD_MP_GPU}")

# --- Enable CUDA if needed ---
if(BUILD_GPU OR BUILD_MP_GPU)
  enable_language(CUDA)
  set(CMAKE_CUDA_FLAGS "-O3")
endif()

# --- Source Files ---
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/../subroutines/*.f90")
# Exclude GPU-only modules from base SOURCES (added back for GPU variants)
set(GPU_ONLY_SOURCES
    "${CMAKE_SOURCE_DIR}/../subroutines/m_batched_diag.f90"
    "${CMAKE_SOURCE_DIR}/../subroutines/m_chefsi.f90"
    "${CMAKE_SOURCE_DIR}/../subroutines/m_ozaki_tridiag.f90"
    "${CMAKE_SOURCE_DIR}/../subroutines/m_qdwh_eigensolver.f90"
    "${CMAKE_SOURCE_DIR}/../subroutines/m_zhev_gpu_handles.f90"
)
list(REMOVE_ITEM SOURCES ${GPU_ONLY_SOURCES})
file(GLOB MAINS "${CMAKE_SOURCE_DIR}/../main/*.f90")
set(MAINS_GPU
    "${CMAKE_SOURCE_DIR}/../main/hsfp0_sc.f90"
    "${CMAKE_SOURCE_DIR}/../main/hrcxq.f90"
    "${CMAKE_SOURCE_DIR}/../main/hx0fp0.f90"
    "${CMAKE_SOURCE_DIR}/../main/hvccfp0.f90"
    "${CMAKE_SOURCE_DIR}/../main/lmf.f90"
)

# --- GEMMul8 paths --- 
set(GEMMUL8_DIR "${CMAKE_BINARY_DIR}/GEMMul8/GEMMul8")
set(GEMMUL8_LIB "${GEMMUL8_DIR}/lib/libgemmul8.so")
set(GEMMUL8_INC "${GEMMUL8_DIR}/include")
set(GEMMUL8_WRAPPER "${CMAKE_SOURCE_DIR}/../subroutines/gemmul8_wrapper.cu")

# --- Compiler Flags ---
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -g")
set(FLAGS -cpp)

if(FC MATCHES "gfortran")
  set(CMAKE_Fortran_COMPILER mpif90)
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fbounds-check")
  set(CMAKE_Fortran_FLAGS "-g -ffree-line-length-none -fimplicit-none -finit-integer=NaN -finit-real=NaN -fpic")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fno-strict-aliasing -fallow-argument-mismatch")
  endif()
  set(LIBM -lmkl_rt)
  set_source_files_properties("../subroutines/mkqp.f90" PROPERTIES COMPILE_FLAGS "-O0")
elseif(FC MATCHES "ifort")
  set(CMAKE_Fortran_COMPILER mpiifort)
  find_program(MPIIFORT_COMPILER mpiifort)
  if(NOT MPIIFORT_COMPILER)
    set(CMAKE_Fortran_COMPILER mpifort)
  endif()
  set(CMAKE_Fortran_FLAGS "-init:snan -fpic -assume byterecl -g -traceback")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  set(LIBM -qmkl=parallel -qopenmp)
  set_source_files_properties("../subroutines/m_qplist.f90" PROPERTIES COMPILE_FLAGS "-O0")
elseif(FC MATCHES "ifx")
  set(CMAKE_Fortran_COMPILER mpiifx)
  set(CMAKE_Fortran_FLAGS "-fpic -assume byterecl -g -traceback")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  set(LIBM -qmkl=parallel -qopenmp)
  set_source_files_properties("../subroutines/m_qplist.f90" PROPERTIES COMPILE_FLAGS "-O0")
elseif(FC MATCHES "nvfortran")
  set(CMAKE_Fortran_COMPILER mpifort)
  set(CMAKE_Fortran_FLAGS "-fpic -Mbackslash -traceback")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  set(LIBM -mp -lmkl_rt)
  set_source_files_properties("../subroutines/pwmat.f90" "../subroutines/rseq.f90" PROPERTIES COMPILE_FLAGS "-O0")
endif()

# --- Dummy executable to aid linker in finding MKL ---
add_executable(dummy ./dummy.f90)
target_link_libraries(dummy ${LIBM})

# --- Main Library and Executables ---
set(LIB ecaljF)
add_library(${LIB} SHARED ${SOURCES})
set_target_properties(${LIB} PROPERTIES
  LINKER_LANGUAGE Fortran
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FC}/mod
)
target_link_options(${LIB} PRIVATE ${LIBM})
target_compile_options(${LIB} PRIVATE ${FLAGS})

set(DELIVER_FILES ${CMAKE_BINARY_DIR}/lib${LIB}.so ${CMAKE_BINARY_DIR}/dummy)

foreach(target IN ITEMS ${MAINS})
  get_filename_component(exename ${target} NAME_WE)
  add_executable(${exename} ${target})
  target_compile_options(${exename} PRIVATE ${FLAGS})
  target_link_libraries(${exename} ${LIB} ${LIBM})
  set_target_properties(${exename} PROPERTIES
    LINKER_LANGUAGE Fortran
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FC}/mod
    BUILD_RPATH "$ORIGIN"
    INSTALL_RPATH "$ORIGIN"
  )
  list(APPEND DELIVER_FILES ${CMAKE_BINARY_DIR}/${exename})
endforeach()

# --- Extended Builds for MP and GPU ---
function(add_extended_build VARIANT)
  set(SOURCES_EXT ${SOURCES})
  # Add GPU-only sources back for GPU variants
  if(VARIANT STREQUAL "gpu" OR VARIANT STREQUAL "mp_gpu")
    list(APPEND SOURCES_EXT ${GPU_ONLY_SOURCES})
  endif()
  # Exclude sources not implemented for GPU or MP
  if(VARIANT STREQUAL "mp" OR VARIANT STREQUAL "gpu" OR VARIANT STREQUAL "mp_gpu")
    list(FILTER SOURCES_EXT EXCLUDE REGEX "main_hmagnon.f90")
    list(FILTER SOURCES_EXT EXCLUDE REGEX "main_mlo_magnon.f90")
    list(FILTER SOURCES_EXT EXCLUDE REGEX "main_hhomogas.f90")
    list(FILTER SOURCES_EXT EXCLUDE REGEX "main_hsfp0.f90")
    list(FILTER SOURCES_EXT EXCLUDE REGEX "sxcf_fal2.f90")
  endif()
  set(LIB ecaljF_${VARIANT})
  if(VARIANT STREQUAL "mp")
    set(CPP_DEFINITIONS -D__MP)
  elseif(VARIANT STREQUAL "gpu")
    set(CPP_DEFINITIONS -D__GPU)
  elseif(VARIANT STREQUAL "mp_gpu")
    set(CPP_DEFINITIONS -D__GPU -D__MP)
  endif()
  set(FLAGS_EXTRA "")
  set(LINK_EXTRA "")
  if (VARIANT STREQUAL "gpu" OR VARIANT STREQUAL "mp_gpu")
    set(FLAGS_EXTRA -Mfree -cuda -acc=gpu -cudalib=cufft,cublas,cusolver)
    set(LINK_EXTRA -acc=gpu -cuda -cudalib=cufft,cublas,cusolver)
  endif()
  message(STATUS "Adding build variant: ${VARIANT}")
  add_library(${LIB} SHARED ${SOURCES_EXT})
  # NVIDIA compiler workaround: remove old library before linking
  add_custom_command(TARGET ${LIB} PRE_LINK
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/lib${LIB}.so)

  target_link_options(${LIB} PRIVATE ${LIBM} ${LINK_EXTRA})
  if (VARIANT STREQUAL "gpu" OR VARIANT STREQUAL "mp_gpu")
    if(EXISTS ${GEMMUL8_LIB})
      set(CPP_DEFINITIONS ${CPP_DEFINITIONS} -D__GEMMUL8)
      target_link_libraries(${LIB} PUBLIC gemmul8wrap)
    endif()
  endif()
  set_target_properties(${LIB} PROPERTIES 
      LINKER_LANGUAGE Fortran
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FC}/mod_${VARIANT}
      BUILD_RPATH "$ORIGIN"
      INSTALL_RPATH "$ORIGIN"
  )
  target_compile_options(${LIB} PRIVATE -cpp ${CPP_DEFINITIONS} ${FLAGS_EXTRA})
  target_compile_definitions(${LIB} PRIVATE ${CPP_DEFINITIONS})
  list(APPEND DELIVER_FILES ${CMAKE_BINARY_DIR}/lib${LIB}.so)

  foreach(target IN ITEMS ${MAINS_GPU})
    get_filename_component(exename ${target} NAME_WE)
    set(exename ${exename}_${VARIANT})
    add_executable(${exename} ${target})
    target_compile_options(${exename} PRIVATE -cpp ${CPP_DEFINITIONS} ${FLAGS_EXTRA})
    target_link_libraries(${exename} ${LIB} ${LIBM} ${LINK_EXTRA})
    set_target_properties(${exename} PROPERTIES
        LINKER_LANGUAGE Fortran
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FC}/mod_${VARIANT}
        BUILD_RPATH "$ORIGIN"
        INSTALL_RPATH "$ORIGIN"
    )
    list(APPEND DELIVER_FILES ${CMAKE_BINARY_DIR}/${exename})
  endforeach()
endfunction()

# --- GEMMul8 Wrapper Library ---
if(BUILD_GPU OR BUILD_MP_GPU)
  if(EXISTS ${GEMMUL8_LIB})
    message(STATUS "GEMMul8 library found at ${GEMMUL8_LIB}. Building gemmul8wrap.")
    add_library(gemmul8wrap SHARED ${GEMMUL8_WRAPPER})
    set_source_files_properties(${GEMMUL8_WRAPPER} PROPERTIES LANGUAGE CUDA)
    target_include_directories(gemmul8wrap PRIVATE ${GEMMUL8_INC})
    target_link_libraries(gemmul8wrap PUBLIC ${GEMMUL8_LIB})
    set_target_properties(gemmul8wrap PROPERTIES
        LINKER_LANGUAGE CUDA
        BUILD_RPATH "$ORIGIN"
        INSTALL_RPATH "$ORIGIN"
    )
  endif()
endif()

if(BUILD_MP)
  add_extended_build("mp")
endif()
if(BUILD_GPU)
  add_extended_build("gpu")
endif()
if(BUILD_MP_GPU)
  add_extended_build("mp_gpu")
endif()

# --- Delivery and Cleanup Targets ---
set(DELIVER_DIR ${CMAKE_SOURCE_DIR}/../exec_${FC})

add_custom_target(deliver ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory ${DELIVER_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${DELIVER_FILES} ${DELIVER_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${DELIVER_FILES} ${CMAKE_SOURCE_DIR}
  DEPENDS ${DELIVER_FILES}
  COMMENT "Copying build products to ${DELIVER_DIR} and ${CMAKE_SOURCE_DIR}"
)

# Create a list of delivered files in the source directory for cleanup
set(DELIVERED_TARGET_FILES "")
foreach(f ${DELIVER_FILES})
  get_filename_component(fname ${f} NAME)
  list(APPEND DELIVERED_TARGET_FILES ${CMAKE_SOURCE_DIR}/${fname})
endforeach()

add_custom_target(clean_delivered
  COMMAND ${CMAKE_COMMAND} -E rm -f ${DELIVERED_TARGET_FILES}
  COMMENT "Cleaning delivered files from ${CMAKE_SOURCE_DIR}"
)
