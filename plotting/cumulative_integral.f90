program cumulative_integral
use parameters
!use hist_class

implicit none 

!input/output
integer::ninput
character(75):: inputfile,outputfile
real(kind=wp),dimension(:),allocatable::x,y
real(kind=wp):: y_out, dx

!loop variables
integer :: i

!other variables
!type(hist)::hist0
real(kind=wp):: xmin, xmax


!read from standard input
read(5,*) inputfile
read(5,*) outputfile
read(5,*) ninput

allocate(x(ninput),y(ninput)) 
open(10,file=inputfile,status='old')
!read file
do i=1,ninput
   read(10,*) x(i),y(i)
end do

close(10)

open(10,file=outputfile,status='unknown')
! write cumulative integral
y_out =0.0_wp

write(10,*) x(1), 0.0_wp
do i=2, ninput
  dx = x(i)-x(i-1)
  y_out = y_out + y(i-1) * dx
  write(10,*) x(i), y_out
end do

close(10)

end program cumulative_integral

