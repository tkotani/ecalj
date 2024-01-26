#Make libecaljF.so and all binaries.
Project (ecalj)
cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)

set(LIBRARY_NAME ecaljF) # library name
set(CMAKE_Fortran_COMPILER mpif90)    # compiler
set(CMAKE_Fortran_FLAGS "-O2 -fPIC ") # compile flags
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) # directory for *.so 
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod) # directory *.mod

# librcaljF.so Library
file(GLOB SOURCES  "../subroutines/*.f90"  "../wanniergw/*.f90")
add_library(${LIBRARY_NAME} SHARED ${SOURCES}) 
# executables
set(LIB ${LIBRARY_NAME} ${MPI_Fortran_LIBRARIES} mkl_rt)  #set LIB variable
set(DIR ${CMAKE_BINARY_DIR}/../main)                      #set DIR variable
foreach(target IN ITEMS lmfa lmchk lmf-MPIK lmfham1 lmfham2 gwinit hbasfp0 heftet hmagnon 
    hmaxloc hpsig_MPI hqpe hqpe_sc hrcxq hsfp0 hsfp0_sc hsocmat huumat_MPI 
    hvccfp0 hwmatK_MPI hx0fp0 hx0init wanplot)
    add_executable(${target} ${DIR}/${target}.f90)
    target_link_libraries(${target} ${LIB})
endforeach()