program broaden
use parameters
use hist_2d_class

implicit none 

!input/output
integer::ninput1,ninput2, noutput1,noutput2, broadening
character(75):: inputfile,outputfile
real(kind=wp),dimension(:),allocatable::x1,x2
real(kind=wp),dimension(:,:),allocatable::y
real(kind=wp):: fwhm

!loop variables
integer :: i, i1, i2

!functions
!real(kind=wp)::dnrm2,ddot

!other variables
type(hist_2d)::hist0
real(kind=wp):: xmin1, xmax1,xmin2, xmax2


!read from standard input
read(5,*) inputfile
read(5,*) outputfile
read(5,*) ninput1, noutput1
read(5,*) ninput2, noutput2
read(5,*) xmin1,xmax1
read(5,*) xmin2,xmax2
read(5,*) broadening, fwhm

allocate(x1(ninput1), x2(ninput1) , y(ninput1,ninput2)) 
open(10,file=inputfile,status='old')

!read file
do i1=1,ninput1
   do i2=1,ninput2
      read(10,*) x1(i1),x2(i2), y(i1,i2)
   end do
end do

!xmin=minval(x)
!xmax=maxval(x)


!initialize histogram
call hist_2d_init(hist0, noutput1, xmin1, xmax1,noutput2, xmin2, xmax2)


! put in histogram
do i1=1,ninput1
   do i2=1,ninput2

      call hist_2d_add(hist0, x1(i1),x2(i2), y(i1,i2))

   end do
end do

!broaden spectra
if (broadening .eq. 1) then
   call hist_2d_broadening(hist0, fwhm)
end if

!write histogram
call hist_2d_write(hist0, outputfile )

end program broaden

