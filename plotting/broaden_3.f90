program broaden
use parameters
use hist_class

implicit none 

!input/output
integer::ninput, noutput,broadening,ncols
character(75):: inputfile,outputfile
real(kind=wp),dimension(:),allocatable::x
real(kind=wp),dimension(:,:),allocatable::y
real(kind=wp):: fwhm

!loop variables
integer :: i,j

!functions
!real(kind=wp)::dnrm2,ddot

!other variables
type(hist)::hist0
real(kind=wp):: xmin, xmax


!read from standard input, only use data from last column!
read(5,*) inputfile, ncols 
read(5,*) outputfile
read(5,*) ninput, noutput
read(5,*) xmin,xmax
read(5,*) broadening, fwhm

allocate(x(ninput),y(ninput,ncols)) 
open(10,file=inputfile,status='old')
!read file
do i=1,ninput
   read(10,*) x(i),(y(i,j),j=1,ncols)
end do

!xmin=minval(x)
!xmax=maxval(x)


!initialize histograms
call hist_init(hist0, noutput, xmin, xmax)

! put in histogram
do i=1,ninput
   call hist_add(hist0,x(i), y(i,ncols))
end do

!broaden spectra
if (broadening .eq. 1) then
   call hist_broadening(hist0, fwhm)
end if

!write histogram
call hist_write(hist0, outputfile )

end program broaden

