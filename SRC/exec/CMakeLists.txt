# Make libecaljF.so and all binaries. Run > FC=gfortran cmake . for generating make file
# run >VERBOSE=1 make to see detailed compile options
cmake_minimum_required(VERSION 3.2)
if(DEFINED CACHE{FC} OR DEFINED ENV{FC})
  set(FC $ENV{FC})
else()
    message("")
    message("Your choices to run cmake are")
    message(" >FC=ifort cmake .")
    message(" >FC=ifort DEBUG=1 cmake .")
    message(" >FC=gfortran cmake .")
    message(" >FC=gfortran DEBUG=1 cmake .")
    message("You may run >VERBOSE=1 make -j, for detailed make messages.")
    message("Note that '>rm CMakeCache.txt' is needed for starting over!")
    message(FATAL_ERROR "Set FC=foobar explained above!")
endif() 
set(FCMOD ${FC})

Project (ecalj)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)
# Note: we always renew show_programinfo.f90 because of Makeshow is always invoked. How to improve this?
add_custom_target(show_programinfo.f90 ALL COMMAND make -f Makeshow WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
file(GLOB SOURCES  "../subroutines/*.f90" "../wanniergw/*.f90" "../subroutines/show_programinfo.f90") # source files
file(GLOB MAINS    "../main/*.f90") # main source files
set(LIBRARY_NAME ecaljF)            # library name is libecaljF.so

#if(DEFINED CACHE{FC} OR DEFINED ENV{FC})
#  set(FC $ENV{FC})
#else()
#  set(FC "gfortran") #default
#endif() 

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # directory for *.so
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${FCMOD}/mod) # directory for *.mod

### show_programinfo.f90 get date and so on for binary info. This should be newer than others ###
### but I don't know how to do it. So I use add_custom_target with ALLch
#find_package(BLAS)
#find_package(MPI)
#find_package(LAPACK)
#message(STATUS ${BLAS_LIBARARIES})
#message(STATUS ${LAPACK})

if(DEFINED CACHE{DEBUG} OR DEFINED ENV{DEBUG})
  set(DEBUG $ENV{DEBUG})
endif() 
message(STATUS "DEBUG= ${DEBUG}")

### compile options ############################################################
if(FC MATCHES "gfortran")
  set(CMAKE_Fortran_COMPILER mpif90)    # compiler
  set(CMAKE_Fortran_FLAGS "-g -fimplicit-none -finit-integer=NaN -finit-real=NaN -fPIC") # compile flags
  set(LIB ${LIBRARY_NAME} -lmkl_rt)  #set LIB variable
  add_link_options(-lmkl_rt)
elseif(FC MATCHES "ifort")
  set(CMAKE_Fortran_COMPILER mpiifort)    # compiler
  set(CMAKE_Fortran_FLAGS "-init:snan -fPIC -assume byterecl") # compile flags
  set(LIB ${LIBRARY_NAME} -mkl)  #set LIB variable
  add_link_options(-mkl)
endif()

### source-dependent compile options 
if(FC MATCHES "gfortran" AND DEFINED DEBUG)
  message(STATUS "===== gfortran debug  =======================")
  foreach(target IN ITEMS ${SOURCES})
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
  endforeach()
elseif(FC MATCHES "gfortran")
  message(STATUS "===== gfortran release  =======================")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/mkqp.f90" OR ${target} MATCHES "../wannier/iopen.f90" )
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    else()
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O2")
    endif()
  endforeach()
elseif(FC MATCHES "ifort" AND DEFINED DEBUG)
  message(STATUS "===== ifort debug  =======================")
  foreach(target IN ITEMS ${SOURCES})
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0  -g -traceback")
  endforeach()
elseif(FC MATCHES "ifort")
  message(STATUS "====== ifort release =======================")
  foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/m_qplist.f90") # OR ${target} MATCHES "../subroutines/x0kf_v4h.f90" )
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    else()
      set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O2")
    endif()
  endforeach()
endif()

if(FC MATCHES "gfortran")   #Set mkl location
  add_executable(dummy ./dummy.f90 ) # exe
  target_link_libraries(dummy ${LIB} ${BLAS} ${MPI}) #linked 
  add_custom_target(fcheck ALL COMMAND ldd dummy|egrep mkl >mklloc.txt  WORKING_DIRECTORY ${CMAKE_BINARY_DIR} DEPENDS dummy) 
else()
  add_custom_target(fcheck ALL COMMAND echo -n >mklloc.txt WORKING_DIRECTORY ${CMAKE_BINARY_DIR}) 
endif()

# Generate library and binaries ###############################
add_library(${LIBRARY_NAME} SHARED ${SOURCES}) # Library
foreach(target IN ITEMS ${MAINS}) #main routines
  get_filename_component(exename ${target} NAME_WE ) #set DIR variable
  add_executable(${exename} ${target} ) # exe
  target_link_libraries(${exename} ${LIB} ${BLAS} ${MPI}) #linked 
endforeach()

#set(BINDIR $ENV{HOME}/bin)
#message(STATUS "BINDIR = ${BINDIR}")
#add_custom_target(inst cp `find . -mindepth 1 -maxdepth 1 -executable -type f` ${BINDIR})
