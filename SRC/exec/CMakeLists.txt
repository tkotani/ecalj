#Make libecaljF.so and all binaries.
#Run >  PLATFORM=gfortran cmake .
Project (ecalj)
cmake_minimum_required(VERSION 2.8)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)
file(GLOB SOURCES  "../subroutines/*.f90" "../wanniergw/*.f90") # source files
file(GLOB MAINS    "../main/*.f90") # main source files
set(LIBRARY_NAME ecaljF)            # library name is libecaljF.so
if(DEFINED CACHE{PLATFORM} OR DEFINED ENV{PLATFORM})
  set(PLATFORM $ENV{PLATFORM})
else()
  set(PLATFORM "gfortran") #default
endif() 
message(STATUS "PLATFORM = ${PLATFORM}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # directory for *.so
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${PLATFORM}/mod) # directory for *.mod

### show_programinfo.f90 get date and so on for binary info. This should be newer than others ###
### but I don't know how to do it. So I use add_custom_target with ALLch
add_custom_target(show_programinfo.f90 ALL COMMAND make -f Makeshow WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#find_package(BLAS)
#find_package(MPI)
#find_package(LAPACK)
#message(STATUS ${BLAS_LIBARARIES})
#message(STATUS ${LAPACK})

### compile options ############################################################
if(PLATFORM MATCHES "gfortran" OR PLATFORM MATCHES "gfortran_dbg")
  set(CMAKE_Fortran_COMPILER mpif90)    # compiler
  set(CMAKE_Fortran_FLAGS "-g -fimplicit-none -finit-integer=NaN -finit-real=NaN -fPIC") # compile flags
  set(LIB ${LIBRARY_NAME} -lmkl_rt)  #set LIB variable
  add_link_options(-lmkl_rt)
elseif(PLATFORM MATCHES "ifort" OR PLATFORM MATCHES "ifort_dbg")
  set(CMAKE_Fortran_COMPILER mpiifort)    # compiler
  set(CMAKE_Fortran_FLAGS "-init:snan -fPIC") # compile flags
  set(LIB ${LIBRARY_NAME} -lmkl)  #set LIB variable
  add_link_options(-lmkl)
endif()

### source-dependent compile options ##########################################
if(PLATFORM MATCHES "gfortran")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/mkqp.f90" OR ${target} MATCHES "../wannier/iopen.f90" )
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    else()
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O2")
    endif()
  endforeach()
elseif(PLATFORM MATCHES "gfortran_dbg")
  foreach(target IN ITEMS ${SOURCES})
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
  endforeach()
elseif(PLATFORM MATCHES "ifort")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/m_qplist.f90") # OR ${target} MATCHES "../subroutines/x0kf_v4h.f90" )
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    else()
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O2")
    endif()
  endforeach()
elseif(PLATFORM MATCHES "ifort_dbg")
  foreach(target IN ITEMS ${SOURCES})
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0  -g -traceback")
  endforeach()
endif()

#not working. add_custom_target(dummys ALL COMMAND echo "call sgemm();end" >dummy.f90)
add_executable(dummy ./dummy.f90 ) # exe
target_link_libraries(dummy ${LIB} ${BLAS} ${MPI}) #linked 
add_custom_target(fcheck ALL COMMAND ldd dummy|grep mkl >mklloc.txt  WORKING_DIRECTORY ${CMAKE_BINARY_DIR} DEPENDS dummy) #check mkl location

#### Generate library and binaries ###############################
add_library(${LIBRARY_NAME} SHARED ${SOURCES}) # Library
foreach(target IN ITEMS ${MAINS}) #main routines
  get_filename_component(exename ${target} NAME_WE ) #set DIR variable
  add_executable(${exename} ${target} ) # exe
  target_link_libraries(${exename} ${LIB} ${BLAS} ${MPI}) #linked 
endforeach()



#set(BINDIR $ENV{HOME}/bin)
#set(BINDIR $ENV{BINDIR}) 
#message(STATUS "BINDIR = ${BINDIR}")
#add_custom_target(inst cp `find . -mindepth 1 -maxdepth 1 -executable -type f` ${BINDIR})
