# Make libecaljF.so and all binaries. Run > FC=gfortran cmake . for generating make file
# run >VERBOSE=1 make to see detailed compile options
cmake_minimum_required(VERSION 3.2)

if(DEFINED CACHE{FC} OR DEFINED ENV{FC})
  set(FC $ENV{FC})
else()
  message("--- ")
  message("=== Your choices to run cmake are: ===")
  message(" >FC=ifort cmake .")
  message(" >FC=ifort cmake . -DCMAKE_BUILD_TYPE=Debug")
  message(" >FC=ifx cmake .")
  message(" >FC=ifx cmake . -DCMAKE_BUILD_TYPE=Debug")
  message(" >FC=gfortran cmake .")
  message(" >FC=gfortran cmake . -DCMAKE_BUILD_TYPE=Debug")
  message(" >FC=nvfortran cmake .")
  message(" >FC=nvfortran cmake . -DCMAKE_BUILD_TYPE=Debug")
  message(" >FC=nvfortran cmake . -DBUILD_GPU=ON -DBUILD_MP=ON -DBUILD_MP_GPU=ON")
  message("You may run >VERBOSE=1 make -j, for detailed make messages.")
  message("Note that '>rm -rf CMakeCache.txt CMakeFiles' is needed for starting over!")
  message(FATAL_ERROR "==> FC=foobar is needed! e.g, FC=gfortran cmake .")
endif()

set(FCMOD ${FC})

Project(ecalj)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
option(BUILD_GPU "Build GPU version" OFF)
option(BUILD_MP "Build mixed precision CPU version" OFF)
option(BUILD_MP_GPU "Build mixed precision GPU version" OFF)
message("CMAKE_BUILD_TYPE is " ${CMAKE_BUILD_TYPE})
message("BUILD_GPU is " ${BUILD_GPU})
message("BUILD_MP is " ${BUILD_MP})
message("BUILD_MP_GPU is " ${BUILD_MP_GPU})
# SET(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

file(GLOB SOURCES "../subroutines/*.f90") # existing source files 
file(GLOB MAINS "../main/*.f90")
file(GLOB MAINS_GPU "../main/hsfp0_sc.f90" "../main/hrcxq.f90" "../main/hx0fp0.f90" "../main/hvccfp0.f90" "../main/lmf.f90")  #main files for GPU versions
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # directory for *.so
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -g") # -traceback")
set(FLAGS -cpp)
# ## compile options ############################################################
if(FC MATCHES "gfortran")
  set(CMAKE_Fortran_COMPILER mpif90) # compiler
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -fbounds-check")
  set(CMAKE_Fortran_FLAGS "-g -ffree-line-length-none -fimplicit-none -finit-integer=NaN -finit-real=NaN -fpic") # compile flags
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fno-strict-aliasing -fallow-argument-mismatch")
  endif()
  set(LIBM -lmkl_rt) 
elseif(FC MATCHES "ifort")
  set(CMAKE_Fortran_COMPILER mpiifort)
  find_program(MPIIFORT_COMPILER mpiifort)
  if(NOT MPIIFORT_COMPILER) 
    set(CMAKE_Fortran_COMPILER mpifort)
    # set(FLAGS ${FLAGS} -D__Nomm3m)               #use zgemm instaed of zgemm3m
    # set(FLAGS ${FLAGS} -D__MKL_ZMINV_SEQUENTIAL) #use sequental version of zgetri 
  endif()
  set(CMAKE_Fortran_FLAGS "-init:snan -fpic -assume byterecl -g -traceback") # compile flags
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  set(LIBM -qmkl=parallel -lmkl_rt -qopenmp)
elseif(FC MATCHES "ifx")
  set(CMAKE_Fortran_COMPILER mpiifx)
  set(CMAKE_Fortran_FLAGS "-fpic -assume byterecl -g -traceback") # ifx: -init:snan option causes segmentation fault in zhegvx
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  set(LIBM -qmkl=parallel -qopenmp)
elseif(FC MATCHES "nvfortran")
  set(CMAKE_Fortran_COMPILER mpifort) # compiler
  #-Mbackslash is necessary to write backslash without escape character
  set(CMAKE_Fortran_FLAGS "-fpic -Mbackslash -traceback") # "-Mvect=simd -mcmodel=medium -fpic") # compile flags
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback")
  # set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -Mbounds -Mchkptr -Mchkstk")
  set(LIBM -mp -lmkl_intel_lp64 -lmkl_pgi_thread -lmkl_core)
endif()

## source-dependent compile options
if(FC MATCHES "gfortran")
  message(STATUS "===== gfortran release  =======================")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/mkqp.f90")
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    endif()
  endforeach()
elseif(FC MATCHES "ifort" OR FC MATCHES "ifx")
  message(STATUS "====== ifort/ifx release =======================")
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

# makeshow which is always invoked to save the comilation time
# execute_process(COMMAND uname -orin  OUTPUT_VARIABLE OUTPUT_UNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND bash -c "LANG=us; date"  OUTPUT_VARIABLE OUTPUT_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND bash -c "${FC} --version |sed '/^$/d' |head -1" OUTPUT_VARIABLE OUTPUT_FCVERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND bash -c "git log |head -1" OUTPUT_VARIABLE OUTPUT_GITCOMMIT1 OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND bash -c "git log |head -2|tail -n 1" OUTPUT_VARIABLE OUTPUT_GITCOMMIT2 OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND bash -c "git log |head -3|tail -n 1" OUTPUT_VARIABLE OUTPUT_GITCOMMIT3 OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND bash -c "grep -v \"^\$\" /etc/issue | head -1 " OUTPUT_VARIABLE OUTPUT_ISSUE OUTPUT_STRIP_TRAILING_WHITESPACE)
# if(CMAKE_BUILD_TYPE STREQUAL "Release")
#   string(JOIN " "  FLAGS_ALL ${FLAGS} ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}) #conbine lists to single variable
# elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
#   string(JOIN " "  FLAGS_ALL ${FLAGS} ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
# endif()
# add_custom_command(
#   OUTPUT ../subroutines/show_programinfo.f90
#   COMMAND cpp -P -D___UNAME_A___=\"${OUTPUT_UNAME}\"
#               -D___ETC_ISSUE___=\"${OUTPUT_ISSUE}\" 
#               -D___FFLAGS___=\"${FLAGS_ALL}\"
#               -D___PLATFORM___=\"${FC}\" 
#               -D___FC_VERSION___=\"${OUTPUT_FCVERSION}\" 
#               -D___LIBLOC___=\"${LIB}\" 
#               -D___LINK_TIME___=\"${OUTPUT_DATE}\" 
#               -D___GIT_COMMIT1___=\"${OUTPUT_GITCOMMIT1}\" 
#               -D___GIT_COMMIT2___=\"${OUTPUT_GITCOMMIT2}\" 
#               -D___GIT_COMMIT3___=\"${OUTPUT_GITCOMMIT3}\"
#               ../subroutines/show_programinfo.template > ../subroutines/show_programinfo.f90
#   COMMENT "Generating show_programinfo.f90!!"
#   VERBATIM
# )
# add_custom_target(show_prog ALL
#   DEPENDS ../subroutines/show_programinfo.f90 #which is generated by custom_command
# )

# emacs tags
add_custom_target(etags ALL COMMAND etags ../subroutines/*.f90 ../main/*.f90
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# copy executables etc to ../exec${FCMOD}
file(MAKE_DIRECTORY   ${CMAKE_BINARY_DIR}/../exec${FCMOD})
set(CMAKE_Fortran_BINS ${CMAKE_BINARY_DIR}/../exec${FCMOD}) # directory for *.mod
add_custom_target(copyexe ALL 
  COMMAND cp `find . -mindepth 1 -maxdepth 1 -executable -type f` ${CMAKE_Fortran_BINS} 
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# For setting mkl location
add_executable(dummy ./dummy.f90) # exe
target_link_libraries(dummy ${BLAS} ${LIBM} ${MPI}) # linked
# add_custom_target(fcheck ALL COMMAND ldd dummy|grep mkl >mklloc.txt  WORKING_DIRECTORY ${CMAKE_BINARY_DIR} DEPENDS dummy)

# library compilation
set(LIB ecaljF) # ecaljF, ecaljF_gpu ecaljF_mp, or ecaljF_mp_gpu
add_library(${LIB} SHARED ${SOURCES})
set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FCMOD}/mod)
#add_dependencies(${LIB} show_prog)
target_link_options(${LIB} PRIVATE ${LIBM})
target_compile_options(${LIB} PRIVATE ${FLAGS})

foreach(target IN ITEMS ${MAINS})
  get_filename_component(exename ${target} NAME_WE)
  add_executable(${exename} ${target})
  target_compile_options(${exename} PRIVATE ${FLAGS})
  target_link_libraries(${exename} ${LIB} ${LIBM})
  #library's path is set to $ORIGIN same with binary's directory
  set_target_properties(${exename} PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FCMOD}/mod
      BUILD_RPATH "\$ORIGIN"  
      INSTALL_RPATH "\$ORIGIN"
  ) 
  add_dependencies(copyexe ${exename})
endforeach()

### GPU or MP version ####
# 1:MP, 2:GPU, 3:GPU+MP
# Some of sources are excluded because GPU or MP are not implemented.
if(BUILD_MP OR BUILD_GPU OR BUILD_MP_GPU)
  list(FILTER SOURCES EXCLUDE REGEX "main_hmagnon.f90")
  list(FILTER SOURCES EXCLUDE REGEX "main_hhomogas.f90")
  list(FILTER SOURCES EXCLUDE REGEX "main_hsfp0.f90")
  list(FILTER SOURCES EXCLUDE REGEX "sxcf_fal2.f90")
endif()
foreach(i RANGE 1 3)
  if(BUILD_MP AND i EQUAL 1)
    set(TYPE _mp)
    set(LIBM2 )
    set(FLAGS2 -D__MP)
  elseif(BUILD_GPU AND i EQUAL 2)
    set(TYPE _gpu)
    set(LIBM2 -acc=gpu -cuda -cudalib=cufft,cublas,cusolver)
    set(FLAGS2 -D__GPU -Minfo -Mfree -cuda -acc=gpu -cudalib=cufft,cublas,cusolver)
  elseif(BUILD_MP_GPU AND i EQUAL 3)
    set(TYPE _mp_gpu)
    set(LIBM2 -acc=gpu -cuda -cudalib=cufft,cublas,cusolver)
    set(FLAGS2 -D__MP -D__GPU -Minfo -Mfree -cuda -acc=gpu -cudalib=cufft,cublas,cusolver)
  else()
    continue()
  endif()
  set(LIB2 ecaljF${TYPE})
  message(STATUS "LIB is " ${LIB2})
  add_library(${LIB2} SHARED ${SOURCES})
  set_target_properties(${LIB2} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FCMOD}/mod${TYPE})
#  add_dependencies(${LIB2} show_prog)
  target_link_options(${LIB2} PRIVATE ${LIBM} ${LIBM2})
  target_compile_options(${LIB2} PRIVATE ${FLAGS} ${FLAGS2})

  foreach(target IN ITEMS ${MAINS_GPU})
    get_filename_component(exename ${target} NAME_WE)
    set(exename ${exename}${TYPE})
    add_executable(${exename} ${target})
    target_compile_options(${exename} PRIVATE ${FLAGS} ${FLAGS2})
    target_link_libraries(${exename} ${LIB2} ${LIBM} ${LIBM2})
    set_target_properties(${exename} PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FCMOD}/mod${TYPE}
        BUILD_RPATH "\$ORIGIN"
        INSTALL_RPATH "\$ORIGIN"
    ) 
    add_dependencies(copyexe ${exename})
  endforeach()
endforeach()
