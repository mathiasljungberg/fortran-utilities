
#F90 =/opt/intel/fc/9.1.036/bin/ifort # ny kompilator på qcllx07
#F90 =/opt/intel/fc/9.0/bin/ifort # kompilator på qcllx18
#F90 =/opt/intel/fc/10.1.015/bin/ifort # kompilator p qcllx21
#F90FLAGS = -O3 -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/
#LFLAGS = -static -static-libcxa #qcllx07
#LFLAGS = -static-libcxa #qcllx18 

F90=gfortran
F90FLAGS = -I/Users/mathiasl/work_bcn/fortranfiler/LIB -J/Users/mathiasl/work_stockholm/fortranfiler/plotting -fbounds-check #-ffixed-line-length-none

#LOBJS8= -L/opt/intel/fc/9.1.036/lib -L/opt/intel/fc/9.1.036/lib -L/usr/libexec/CERNLIB/2004/lib/ -L/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/ -llapack3 -lblas 
#LOBJS8= /usr/libexec/CERNLIB/2004/lib/liblapack3.a  /usr/libexec/CERNLIB/2004/lib/libblas.a
#LOPT= -L/usr/lib/  -L/usr/libexec/CERNLIB/2004/lib  -lifport -lifcoremt -limf -lm \
#-lcxa -lpthread -lirc -lunwind -lc -lirc_s -llapack3 -lblas -lg2c 

#LOPT= -L/usr/lib/  -L/usr/libexec/CERNLIB/2004/lib  -llapack3 -lblas  -lg2c 

# order important in the following sequence: -llapack3 -lblas -lg2c 

#-L/opt/intel/fc/9.1.036/lib -L/usr/libexec/CERNLIB/2004/lib \
#-lifport -lifcoremt -limf -lm -lcxa -lpthread -lirc -lunwind -lc -lirc_s \
#-llapack3 -lblas -lg2c

#LFLAGS = -static -static-libcxa
#F77 =f77
#F77FLAGS = 
#LDIR =  /afs/physto.se/home/m/mathiasl/fortranfiler/LIB
#LOBJS = bubblesort.o 

LOBJS = /Users/mathiasl/work_bcn/fortranfiler/LIB/bubblesort.o \
         /usr/libexec/CERNLIB/2004/lib/libblas.a
LOBJS2 = /Users/mathiasl/work_bcn/fortranfiler/LIB/bubblesort.o \
	/Users/mathiasl/work_bcn/fortranfiler/LIB/crossprod.o \
	/usr/libexec/CERNLIB/2004/lib/libblas.a
LOBJS3 = /usr/libexec/CERNLIB/2004/lib/libblas.a
LOBJS4=/Users/mathiasl/work_bcn/fortranfiler/LIB/splines.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/spline_easy.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/linspace.o

LOBJS5=/Users/mathiasl/work_bcn/fortranfiler/LIB/splines.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/spline_easy.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/linspace.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/hist_class.o


LOBJS6=/Users/mathiasl/work_bcn/fortranfiler/LIB/splines.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/spline_easy.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/linspace.o\
	/Users/mathiasl/work_bcn/fortranfiler/LIB/hist_2d_class.o


addmanygraphs:addmanygraphs.o
	$(F90) $(LFLAGS) -o addmanygraphs.out addmanygraphs.o $(LOBJS4)
addmanygraphs_new:addmanygraphs_new.o
	$(F90) $(LFLAGS) -o addmanygraphs_new addmanygraphs_new.o $(LOBJS4)
addmanygraphs_new_henrik:addmanygraphs_new_henrik.o
	$(F90) $(LFLAGS) -o addmanygraphs_new_henrik addmanygraphs_new_henrik.o $(LOBJS4)
broaden:broaden.o
	$(F90) $(F90FLAGS) -o broaden broaden.o $(LOBJS5)
broaden_2:broaden_2.o
	$(F90) $(LFLAGS) -o broaden_2 broaden_2.o $(LOBJS5)
broaden_3:broaden_3.o
	$(F90) $(F90FLAGS) -o broaden_3 broaden_3.o $(LOBJS5)
broaden_2d:broaden_2d.o
	$(F90) $(LFLAGS) -o broaden_2d broaden_2d.o $(LOBJS6)



#	$(F90) $(F90FLAGS) -c distribution.f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $*.f90



clean: 
	rm -f *.o 

#depend:
#	makedepend -- $(FLAGS) -- *.cpp
