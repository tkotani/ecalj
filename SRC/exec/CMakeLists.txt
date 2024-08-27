# Make libecaljF.so and all binaries. Run > FC=gfortran cmake . for generating make file
# run >VERBOSE=1 make to see detailed compile options
cmake_minimum_required(VERSION 3.2)

if(DEFINED CACHE{FC} OR DEFINED ENV{FC})
  set(FC $ENV{FC})
else()
  message("--- ")
  message("=== Your choices to run cmake are: ===")
  message(" >FC=ifort cmake .")
  message(" >FC=ifort DEBUG=1 cmake .")
  message(" >FC=gfortran cmake .")
  message(" >FC=gfortran DEBUG=1 cmake .")
  message(" >FC=nvfortran cmake .")
  message(" >FC=nvfortran DEBUG=1 cmake .")
  message("You may run >VERBOSE=1 make -j, for detailed make messages.")
  message("Note that '>rm -rf CMakeCache.txt CMakeFiles' is needed for starting over!")
  message(FATAL_ERROR "==> FC=foobar is needed! e.g, FC=gfortran cmake .")
endif()

set(FCMOD ${FC})

Project(ecalj)
option(BUILD_GPU "Build for GPU" OFF)
option(BUILD_MP "Build for MP" OFF)

message("BUILD_GPU is " ${BUILD_GPU})
message("BUILD_MP is " ${BUILD_MP})

# SET(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

# Note: we always renew show_programinfo.f90 because of Makeshow is always invoked. How to improve this?
add_custom_target(show_prog COMMAND make -f Makeshow WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  DEPENDS "../subroutines/*.f90" "../wanniergw/*.f90" "../main/*.f90"
)

file(GLOB SOURCES "../subroutines/*.f90" "../wanniergw/*.f90" "../subroutines/show_programinfo.f90") # source files
file(GLOB SOURCES1 "../subroutines/*.f90" "../subroutines/show_programinfo.f90") # source files
file(GLOB SOURCES2 "../wanniergw/*.f90") # source files
file(GLOB MAINS "../main/*.f90")
file(GLOB MAINS_BOTTLENECK "../main/hsfp0_sc.f90" "../main/hrcxq.f90" "../main/hx0fp0.f90") 
# list(REMOVE_ITEM MAINS ${MAINS_BOTTLENECK})  #exclude MAINS _BOTTOLECK from MAINS

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # directory for *.so

set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -g -traceback")
# ## compile options ############################################################
if(FC MATCHES "gfortran")
  set(CMAKE_Fortran_COMPILER mpif90) # compiler
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "-g -fimplicit-none -finit-integer=NaN -finit-real=NaN -fpic") # compile flags
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fno-strict-aliasing -fallow-argument-mismatch")
  endif()
  set(LIBM -lmkl_rt) 
elseif(FC MATCHES "ifort")
  set(CMAKE_Fortran_COMPILER mpiifort) # compiler
  set(CMAKE_Fortran_FLAGS "-init:snan -fpic -assume byterecl") # compile flags
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  set(LIBM -qmkl=parallel) 
elseif(FC MATCHES "nvfortran")
  set(CMAKE_Fortran_COMPILER mpifort) # compiler
  set(CMAKE_Fortran_FLAGS "-fpic") # "-Mvect=simd -mcmodel=medium -fpic") # compile flags
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  set(LIBM -lmkl_intel_lp64 -lmkl_sequential -lmkl_core)
endif()

# library compilation
set(FLAGS -cpp)
set(TYPE "")
if(${BUILD_MP})
  set(TYPE ${TYPE}_mp)
  set(FLAGS ${FLAGS} -D__MP)
endif()
if(${BUILD_GPU})
  set(TYPE ${TYPE}_gpu)
  set(LIBM ${LIBM} -acc=gpu -cuda -cudalib=cufft,cublas,cusolver)
  set(FLAGS ${FLAGS} -D__GPU -Minfo -Mfree -cuda -acc=gpu -cudalib=cufft,cublas,cusolver)
endif()
set(LIB ecaljF${TYPE})
set(MOD mod${TYPE})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FCMOD}/${MOD})
add_library(${LIB} SHARED ${SOURCES})
target_link_options(${LIB} PRIVATE ${LIBM})
target_compile_options(${LIB} PRIVATE ${FLAGS})

# emacs tags
add_custom_target(etags ALL COMMAND etags ../subroutines/*.f90 ../wanniergw/*.f90 ../main/*.f90
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# For setting mkl location
add_executable(dummy ./dummy.f90) # exe
target_link_libraries(dummy ${BLAS} ${LIBM} ${MPI}) # linked
# add_custom_target(fcheck ALL COMMAND ldd dummy|grep mkl >mklloc.txt  WORKING_DIRECTORY ${CMAKE_BINARY_DIR} DEPENDS dummy)
# Generate library and binaries ###############################
#add_library(${LIB1} SHARED ${SOURCES1}) # Library
#add_library(${LIB2} SHARED ${SOURCES2})

# ## source-dependent compile options
if(FC MATCHES "gfortran")
  message(STATUS "===== gfortran release  =======================")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/mkqp.f90" OR ${target} MATCHES "../wannier/iopen.f90")
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    endif()
  endforeach()
elseif(FC MATCHES "ifort")
  message(STATUS "====== ifort release =======================")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/m_qplist.f90")# OR ${target} MATCHES "../subroutines/m_bzintegration2.f90" )
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    endif()
  endforeach()
elseif(FC MATCHES "nvfortran")
  message(STATUS "====== nvfortran release =======================")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/pwmat.f90" OR ${target} MATCHES "../subroutines/rseq.f90")
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    endif()
  endforeach()
endif()

if(TYPE) #if TYPE is not empty
  foreach(target IN ITEMS ${MAINS_BOTTLENECK})
    get_filename_component(exename ${target} NAME_WE)
    set(exename ${exename}${TYPE})
    add_executable(${exename} ${target})
    add_dependencies(${exename} ${LIB})
    target_compile_options(${exename} PRIVATE ${FLAGS})
    target_link_libraries(${exename} ${LIB} ${LIBM})
  endforeach()
else()
  foreach(target IN ITEMS ${MAINS})
    get_filename_component(exename ${target} NAME_WE)
    add_executable(${exename} ${target})
    add_dependencies(${exename} ${LIB})
    target_link_libraries(${exename} ${LIB} ${LIBM})
  endforeach()
endif()

file(MAKE_DIRECTORY   ${CMAKE_BINARY_DIR}/../exec${FCMOD})
set(CMAKE_Fortran_BINS ${CMAKE_BINARY_DIR}/../exec${FCMOD}) # directory for *.mod
add_custom_target(copyexe ALL COMMAND cp `find . -mindepth 1 -maxdepth 1 -executable -type f` ${CMAKE_Fortran_BINS} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
