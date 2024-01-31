#Make libecaljF.so and all binaries.
#Run >  PLATFORM=gfortran cmake .
Project (ecalj)
cmake_minimum_required(VERSION 3.8)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)
file(GLOB SOURCES  "../subroutines/*.f90" "../wanniergw/*.f90" "./show_programinfo.f90") # source files
file(GLOB MAINS    "../main/*.f90") # main source files
set(PLATFORM $ENV{PLATFORM}) 
message(STATUS "PLATFORM = ${PLATFORM}")
#set(PLATFORM "gfortran")
set(LIBRARY_NAME ecaljF) # library name
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PLATFORM}/lib) # directory for *.so
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${PLATFORM}/mod) # directory *.mod
# compile options ##########################################
add_custom_target(show_programinfo.f90  ALL COMMAND make -f Makeshow WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
if(PLATFORM MATCHES "gfortran")
    set(CMAKE_Fortran_COMPILER mpif90)    # compiler
    set(CMAKE_Fortran_FLAGS "-g -fimplicit-none -finit-integer=NaN -finit-real=NaN -fPIC") # compile flags
    foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/mkqp.f90" OR
       ${target} MATCHES "../wannier/iopen.f90" )
        set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    else()
        set_source_files_properties(${target} P ROPERTIES COMPILE_FLAGS "-O2")
    endif()
    endforeach()
    set(LIB ${LIBRARY_NAME} mkl_rt)  #set LIB variable
elseif(PLATFORM MATCHES "gfortran_dbg")
    set(CMAKE_Fortran_COMPILER mpif90)    # compiler
    set(CMAKE_Fortran_FLAGS "-g -fimplicit-none -finit-integer=NaN -finit-real=NaN -fPIC") # compile flags
    foreach(target IN ITEMS ${SOURCES})
        set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    endforeach()
    set(LIB ${LIBRARY_NAME} mkl_rt)  #set LIB variable
elseif(PLATFORM MATCHES "ifort")
    set(CMAKE_Fortran_COMPILER mpiifort)    # compiler
    set(CMAKE_Fortran_FLAGS "-init:snan -fPIC") # compile flags
    foreach(target IN ITEMS ${SOURCES})
    if(${target} MATCHES "../subroutines/m_qplist.f90") # OR ${target} MATCHES "../subroutines/x0kf_v4h.f90" )
        set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    else()
        set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O2")
    endif()
    endforeach()
    set(LIB ${LIBRARY_NAME} mkl)  #set LIB variable
elseif(PLATFORM MATCHES "ifort_dbg")
    set(CMAKE_Fortran_COMPILER mpiifort)    # compiler
    set(CMAKE_Fortran_FLAGS "-traceback -init:snan -fPIC") # compile flags
    foreach(target IN ITEMS ${SOURCES})
        set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS "-O0")
    endforeach()
    set(LIB ${LIBRARY_NAME} mkl)  #set LIB variable
endif()

#### Generate library and binaries ###############################
set_source_files_properties(${CMAKE_BINARY_DIR}/show_programinfo.f90 PROPERTIES COMPILE_FLAGS "-O0")
add_library(${LIBRARY_NAME} SHARED ${SOURCES}) # Library
foreach(target IN ITEMS ${MAINS})
    get_filename_component(exename ${target} NAME_WE ) #set DIR variable
    add_executable(${exename} ${target} ) # exe
    target_link_libraries(${exename} ${LIB})
endforeach()

#set(BINDIR $ENV{HOME}/bin)
#set(BINDIR $ENV{BINDIR}) 
#message(STATUS "BINDIR = ${BINDIR}")
#add_custom_target(inst cp `find . -mindepth 1 -maxdepth 1 -executable -type f` ${BINDIR})
