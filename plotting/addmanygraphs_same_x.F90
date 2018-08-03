program addmanygraphs
!use parameters
implicit none 

!parameters
integer,parameter::maxsize=5000 ! max, if changed, must change in spine functions

!input/output
!character(75),dimension(maxinfiles)::inputfile
character(75),dimension(:),allocatable::inputfile
character(75)::outfile
integer:: ninputfiles, noutput
real(kind=8),dimension(:),allocatable:: amount,shift
integer,dimension(:),allocatable::ninput
real(kind=8),dimension(:,:),allocatable::x,y
real(kind=8),dimension(:),allocatable:: ynew,xnew,ytemp

!loop variables
integer :: i,j,ind

!other variables
real(kind=8):: dx,dx2,xmin,xmin2,xmax,xmax2
real(kind=8):: xrange,xrange2,xmax_tmp,xmin_tmp
integer :: rx,flag
real(kind=8),dimension(:,:),allocatable:: ynew2
real(kind=8),dimension(maxsize):: x_tmp,y_tmp

! sums together an arbitrary number of graphs with different weights

!read from standard input
read(5,*) ninputfiles 

allocate(inputfile(ninputfiles), ninput(ninputfiles), amount(ninputfiles),&
     shift(ninputfiles))

write(6,*) ninputfiles
 
do i=1,ninputfiles
   read(5,*) inputfile(i), ninput(i), amount(i)
   write(6,*) inputfile(i), ninput(i), amount(i)
end do

allocate( x(maxval(ninput),ninputfiles),y(maxval(ninput),ninputfiles))

read(5,*) outfile, noutput
write(6,*) outfile, noutput

allocate(xnew(noutput), ynew(noutput))!, ynew2(noutput,ninputfiles))
!allocate(ytemp(noutput))

! done reading

x=0
y=0

! read inputfiles
do i=1,ninputfiles
   open(10,file=inputfile(i),status='old')
   do j=1,ninput(i)
      read(10,*) x(j,i), y(j,i)
   end do
   close(10)
end do

ynew=0
!add alla graphs 
do i=1,ninputfiles
   do j=1,noutput
      ynew(j) = ynew(j) + amount(i)*y(j,i)
   end do
end do

! write to file
open(10,file=outfile,status='unknown')
do i=1,noutput
   write(10,*) x(i,1), ynew(i)
end do
close(10)

end program addmanygraphs
