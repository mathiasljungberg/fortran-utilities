#F90 =/opt/intel_fc_80/bin/ifort # kompilator på qcllx07
#F77=g77
#F90 =/opt/intel/fc/9.0/bin/ifort # kompilator på qcllx18
#F90 =/opt/intel/fc/10.1.015/bin/ifort # kompilator p qcllx21
#F90=pgf90 # qcllx22
F90=gfortran
#F90FLAGS = -O3
#F77FLAGS = -O3
#LFLAGS = -static-libcxa
#LFLAGS =
#LOBJS = /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/bubblesort.o \
#         /usr/libexec/CERNLIB/2004/lib/libblas.a
#LOBJS2 = /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/bubblesort.o \
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/crossprod.o \
#	/usr/libexec/CERNLIB/2004/lib/libblas.a
#LOBJS3 = /usr/libexec/CERNLIB/2004/lib/libblas.a

filesF90 = linspace.f90 spline_easy.f90 parameters.f90 hist_class.f90 spline_m.f90 FFT_m.f90
filesF77 = splines.f factorials.f splib.f sffteu.f


allF90: $(filesF90)
	$(F90) $(F90FLAGS) -c $(filesF90)

allF77: $(filesF77)
	$(F90) $(F90FLAGS) -c $(filesF77)

parameters: parameters.f90  
	$(F90) $(F90FLAGS) -c parameters.f90

hist_class: hist_class.f90
	 $(F90) $(F90FLAGS) -c hist_class.f90
 
freq_mukamel: freq_mukamel.f90
	 $(F90) $(F90FLAGS) -c freq_mukamel.f90 

qsort_c_module: qsort_c_module.f90
	 $(F90) $(F90FLAGS) -c qsort_c_module.f90 

hist_2d_class: hist_2d_class.f90
	 $(F90) $(F90FLAGS) -c hist_2d_class.f90

all: @.f90
	$(F90) $(F90FLAGS) -c @.f90


#%.o: %.f90
#	$(F90) $(F90FLAGS) -c $*.f90
#%.mod: %.f90
#	$(F90) $(F90FLAGS) -c $*.f90

#%.o: %.f
#	$(F77) $(F77FLAGS) -c $*.f

#clean: 
#	rm -f *.o timing_scalarprod.out

#depend:
#	makedepend -- $(FLAGS) -- *.cpp
