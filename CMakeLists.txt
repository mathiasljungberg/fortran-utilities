cmake_minimum_required(VERSION 2.4)

PROJECT(fortran-utilities Fortran)

IF(CMAKE_Fortran_COMPILER MATCHES ifort)

  MESSAGE(STATUS  "ifort compiler found")

  SET(FLAGS "-O0 -warn all -g -CB")
  SET(INCDIR "-I/usr/include")
  #SET(OMP_FLAGS "")
  #SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")
  #SET(MKLPATH "/opt/intelcompiler/mkl/lib/intel64")
  SET(LAPACK_LIB "-Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_sequential.a ${MKLPATH}/libmkl_core.a -Wl,--end-group")

ELSE(CMAKE_Fortran_COMPILER MATCHES gfortran)

  MESSAGE(STATUS  "gfortran compiler found")

  IF(USER_OS STREQUAL mac)

    MESSAGE(STATUS  "option -DUSER_OS=mac, gfortran for mac")

    # gfortran/mac
    #SET(FLAGS "-O0  -fbounds-check -ffree-line-length-none -framework vecLib -DUSE_OWN_CZDOT")
    SET(FLAGS "-O0  -fbounds-check -framework vecLib -DUSE_OWN_CZDOT")  
    SET(INCDIR "-I/usr/local/include -I/usr/include")
    #SET(OMP_FLAGS "")
    #SET(FFTW_LIB "-L/usr/local/lib -lfftw3 -lm")
    SET(LAPACK_LIB "")

  ELSE(USER_OS STREQUAL mac)

    MESSAGE(STATUS  "gfortran for Linux")

    # gfortran/linux
    SET(FLAGS "-O0  -fbounds-check -ffree-line-length-none -fmax-errors=5")
    SET(INCDIR "-I/usr/include")
    SET(OMP_FLAGS "")
    SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")
    SET(MKLPATH "")
    SET(LAPACK_LIB "-L/usr/lib -lblas -llapack")

  ENDIF(USER_OS STREQUAL mac)

ENDIF(CMAKE_Fortran_COMPILER MATCHES ifort)


SET(ALL_MODULES_F90 
m_FFT.f90               hist_2d_class.f90       spline_easy2d.f90
delta.f90               hist_class.f90          match_wild.f90          spline_easy.f90         spline_m.f90
fcn3params.f90          linspace.f90            parameters.f90          qsort_c_module.f90      spline_easy2.f90
m_splines.F90
)

SET(ALL_MODULES_F 
bubblesort.f    crossprod.f     factorials.f    rnd.f           sffteu.f        splib.f         splines.f
)


ADD_LIBRARY(fortran-utilities STATIC ${ALL_MODULES_F90} ${ALL_MODULES_F})

add_subdirectory(morsefitting morsefitting)
add_subdirectory(plotting plotting)


#FOREACH(EXEC ${EXECUTABLES})
#  ADD_EXECUTABLE(${EXEC} ${EXEC}.F90 ${ALL_MODULES})
#  TARGET_LINK_LIBRARIES(${EXEC} ${LAPACK_LIBRARIES} 
#    ${LINK_ADDONS})
#ENDFOREACH(EXEC)
