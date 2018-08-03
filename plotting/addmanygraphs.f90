program addmanygraphs
use parameters
implicit none 

!parameters
integer,parameter::maxsize=5000 ! max, if changed, must change in spine functions

!input/output
!character(75),dimension(maxinfiles)::inputfile
character(75),dimension(:),allocatable::inputfile
character(75)::outfile
integer:: ninputfiles, noutput
!real(kind=wp),dimension(maxinfiles):: amount
!real(kind=wp),dimension(maxsize,maxinfiles)::x,y
!real(kind=wp),dimension(maxsize):: ynew,xnew
real(kind=wp),dimension(:),allocatable:: amount,shift,ninput
real(kind=wp),dimension(:,:),allocatable::x,y
real(kind=wp),dimension(:),allocatable:: ynew,xnew

!loop variables
integer :: i,j

!other variables
real(kind=wp):: dx,dx2,xmin,xmin2,xmax,xmax2
real(kind=wp):: xrange,xrange2
integer :: rx
!real(kind=wp),dimension(maxsize,maxinfiles):: ynew2
real(kind=wp),dimension(:,:),allocatable:: ynew2

! sums together an arbitrary number of graphs with different weights

!read from standard input
read(5,*) ninputfiles

allocate(inputfile(ninputfiles), ninput(i), amount(ninputfiles),shift(ninputfiles),&
     x(maxsize,ninputfiles),y(maxsize,ninputfiles))

do i=1,ninputfiles
   read(5,*) inputfile(i), ninput(i), amount(i), shift(i)
end do

read(5,*) outfile, noutput

allocate(xnew(noutput), xnew(noutput), ynew2(noutput,ninputfiles))

read(5,*) range, start, end
! done reading

x=0
y=0
xnew=0
ynew=0
ynew2=0

! read inputfiles
do i=1,ninputfiles
   open(10,file=inputfile(i),status='old')
   do j=1,ninput(i)
      read(10,*) x(j,i), y(j,i)
   end do
   close(10)
end do

! read inputfile2
!open(10,file=inputfile2,status='old')
!do i=1,ninput
!   read(10,*) x2(i), y2(i)
!end do
!close(10)


!
if(range.eq.1) then 
   xmin = start
   xmax = end
   
end if

xmin = minval(x(1:ninput,1:ninputfiles))
xmax = maxval(x(1:ninput,1:ninputfiles))

print *, xmin, xmax

xrange= xmax-xmin
dx = xrange/noutput



! create new x-values
do i=1,noutput
   xnew(i) = xmin + (i-1)*dx !-dx/2.
end do

ynew2=0
! put graphs in new range and add them (do not care about eventual spillover)
do i=1,ninputfiles
   
   xmin2 = minval(x(1:ninput,i))
   xmax2 = maxval(x(1:ninput,i))
   xrange2 = xmax2-xmin2
   dx2 = xrange2/ninput

   print *, xmin2, xmax2, xrange2,dx2

   do j=1,noutput

      !graph 1
      
      rx =int((xnew(j)-xmin2)/dx2 )
      
      ! force them inte a bin
      !if(rx.lt.1) rx=1
      !if(rx.ge.ninput) rx=ninput
      
      !no contribution if out of range 
      if(rx.lt.1.or.rx.gt.(ninput-1)) then
         ynew2(j,i) = 0
    
      else         
         ! linear inperpolation between rx 
         ynew2(j,i)= y(rx,i)+(xnew(j)-x(rx,i))*y(rx+1,i)
         !print *, rx,y(rx,i),ynew2(j,i)
      end if

   end do
end do

ynew=0
!add alla graphs 
do i=1,ninputfiles
   do j=1,noutput
      ynew(j) = ynew(j) + amount(i)*ynew2(j,i)
   end do
end do


! write to file
open(10,file=outfile,status='unknown')
do i=1,noutput
   write(10,*) xnew(i), ynew(i)
end do
close(10)

! write to file2
!open(10,file=outfile2,status='unknown')
!do i=1,ninput
!   write(10,*) xnew(i), amount1*ynew1(i)
!end do
!close(10)


! write to file3
!open(10,file=outfile3,status='unknown')
!do i=1,ninput
!   write(10,*) xnew(i),  amount2*ynew2(i)
!end do
!close(10)


end program addmanygraphs
