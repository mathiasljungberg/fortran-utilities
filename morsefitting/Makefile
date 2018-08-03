#F90 =/opt/intel_fc_80/bin/ifort # kompilator på qcllx07
#F90 =/opt/intel/fc/9.1.036/bin/ifort # ny kompilator på qcllx07
F90 =/opt/intel/fc/9.0/bin/ifort # kompilator på qcllx18
F90FLAGS = -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB
LFLAGS = -static-libcxa 
#LFLAGS = -static -static-libcxa
#F77 =f77
#F77FLAGS = 
#LDIR =  /afs/physto.se/home/m/mathiasl/fortranfiler/LIB
#LOBJS = bubblesort.o 
LOBJS = /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/bubblesort.o \
         /usr/libexec/CERNLIB/2004/lib/libblas.a
LOBJS2 = /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/bubblesort.o \
	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/crossprod.o \
	/usr/libexec/CERNLIB/2004/lib/libblas.a
LOBJS3 = /usr/libexec/CERNLIB/2004/lib/libblas.a
LOBJS4=
LOBJS5=/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/rnd.o
LOBJS6 = /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/bubblesort.o \
	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/crossprod.o \
	/usr/libexec/CERNLIB/2004/lib/libblas.a \
	/afs/physto.se/home/m/mathiasl/fortranfiler/test/hist3_class.o
LOBJS7 = /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/bubblesort.o \
	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/crossprod.o \
	/usr/libexec/CERNLIB/2004/lib/libblas.a \
	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/freq_mukamel.o \
	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/hist_class.o

LOBJS8 = /usr/libexec/CERNLIB/2004/lib/libblas.a 

morselmdr_fixdist: morselmdr_fixdist.o lmder.o fcn2params.o fitting_routines.o
	$(F90) $(LFLAGS) -o morselmdr_fixdist.out morselmdr_fixdist.o lmder.o\
	fcn2params.o fitting_routines.o $(LOBJS8)
polyfit: polyfit.o lmder.o fcn2params.o fitting_routines.o
	$(F90) $(LFLAGS) -o polyfit.out polyfit.o lmder.o\
	fcn2params.o fitting_routines.o $(LOBJS8)
morse_fit: morse_fit.o lmder.o fcn2params.o fitting_routines.o
	$(F90) $(LFLAGS) -o morse_fit.x morse_fit.o lmder.o\
	fcn2params.o fitting_routines.o $(LOBJS8)
morse_fit_fixfreq: morse_fit_fixfreq.o lmder.o fcn2params.o fcn8params.o \
	fitting_routines.o
	$(F90) $(LFLAGS) -o morse_fit_fixfreq morse_fit_fixfreq.o lmder.o\
	fcn2params.o fcn8params.o fitting_routines.o $(LOBJS8)
spline_energy: spline_energy.o
	$(F90) $(LFLAGS) -o spline_energy spline_energy.o ../LIB/linspace.o \
	../LIB/spline_easy.o ../LIB/splines.o

#calculate_e-field.o: calculate_e-field.f90
#	$(F90) $(F90FLAGS) -c calculate_e-field.f90
#watersort.o: watersort.f90
#	$(F90) $(F90FLAGS) -c watersort.f90
#replaceSPC.o: replaceSPC.f90
#	$(F90) $(F90FLAGS) -c replaceSPC.f90
#makeperiodic.o: makeperiodic.f90
#	$(F90) $(F90FLAGS) -c makeperiodic.f90
#cutoutnearest.o: cutoutnearest.f90
#	$(F90) $(F90FLAGS) -c cutoutnearest.f90
#makeinput.o: makeinput.f90
#	$(F90) $(F90FLAGS) -c makeinput.f90	
#frequencies.o: frequencies.f90
#	$(F90) $(F90FLAGS) -c frequencies.f90
#mdtoxyz.o: mdtoxyz.f90
#	$(F90) $(F90FLAGS) -c mdtoxyz.f90
#cutoutneighbors.o:cutoutneighbors.f90
#	$(F90) $(F90FLAGS) -c cutoutneighbors.f90 
#distribution.o: distribution.f90
#	$(F90) $(F90FLAGS) -c distribution.f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $*.f90
%.mod: %.f90
	$(F90) $(F90FLAGS) -c $*.f90

#morselmdr_fixdist.o: morselmdr_fixdist.f90
#	$(F90) $(F90FLAGS) -c $*.f90
#lmder.o: morselmdr_fixdist.f
#	$(F90) $(F90FLAGS) -c $*.f90


#depend:
#	makedepend -- $(FLAGS) -- *.cpp
