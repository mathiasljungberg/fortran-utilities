program broaden
use parameters
use hist_class

implicit none 

!input/output
integer::ninput, noutput,broadening
character(75):: inputfile,outputfile
real(kind=wp),dimension(:),allocatable::x,y
real(kind=wp):: fwhm

!loop variables
integer :: i

!functions
!real(kind=wp)::dnrm2,ddot

!other variables
type(hist)::hist0
real(kind=wp):: xmin, xmax


!read from standard input
read(5,*) inputfile
read(5,*) outputfile
read(5,*) ninput, noutput
read(5,*) xmin,xmax
read(5,*) broadening, fwhm

allocate(x(ninput),y(ninput)) 
open(10,file=inputfile,status='old')
!read file
do i=1,ninput
   read(10,*) x(i),y(i)
end do

!xmin=minval(x)
!xmax=maxval(x)


!initialize histogram
call hist_init(hist0, noutput, xmin, xmax)


! put in histogram
do i=1,ninput
   call hist_add(hist0,x(i), y(i))
end do

!broaden spectra
if (broadening .eq. 1) then
   call hist_broadening(hist0, fwhm)
end if

!write histogram
call hist_write(hist0, outputfile )

end program broaden

