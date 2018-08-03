program addmanygraphs
use parameters
implicit none 

!parameters
integer,parameter::maxsize=5000 ! max, if changed, must change in spine functions

!input/output
!character(75),dimension(maxinfiles)::inputfile
character(75),dimension(:),allocatable::inputfile
character(75)::outfile, dummy
integer:: ninputfiles, noutput,range, normalize
!real(kind=wp),dimension(maxinfiles):: amount
!real(kind=wp),dimension(maxsize,maxinfiles)::x,y
!real(kind=wp),dimension(maxsize):: ynew,xnew
real(kind=wp):: start,end
real(kind=wp),dimension(:),allocatable:: amount,shift
integer,dimension(:),allocatable::ninput,column
real(kind=wp),dimension(:,:),allocatable::x,y
real(kind=wp),dimension(:),allocatable:: ynew,xnew,ytemp

!loop variables
integer :: i,j,ind

!other variables
real(kind=wp):: dx,dx2,xmin,xmin2,xmax,xmax2
real(kind=wp):: xrange,xrange2,xmax_tmp,xmin_tmp
integer :: rx,flag
!real(kind=wp),dimension(maxsize,maxinfiles):: ynew2
real(kind=wp),dimension(:,:),allocatable:: ynew2
real(kind=wp),dimension(maxsize):: x_tmp,y_tmp

! sums together an arbitrary number of graphs with different weights
! column is the column to be summed, 2 should be default

!read from standard input
read(5,*) ninputfiles 

allocate(inputfile(ninputfiles), ninput(ninputfiles), amount(ninputfiles),&
     shift(ninputfiles), x(maxsize,ninputfiles),y(maxsize,ninputfiles),column(ninputfiles))

write(6,*) ninputfiles
 
do i=1,ninputfiles
   read(5,*) inputfile(i), column(i),ninput(i), amount(i), shift(i)
   write(6,*) inputfile(i), column(i), ninput(i), amount(i), shift(i)
end do


read(5,*) outfile, noutput
write(6,*) outfile, noutput

allocate(xnew(noutput), ynew(noutput), ynew2(noutput,ninputfiles))
allocate(ytemp(noutput))

read(5,*) range, start, end 
 write(6,*) range, start, end

read(5,*) normalize


! done reading

x=0
y=0
xnew=0
ynew=0
ynew2=0

! read inputfiles

do i=1,ninputfiles
   open(10,file=inputfile(i),status='old')

   if(column(i).eq.2) then
      do j=1,ninput(i)
         read(10,*) x(j,i), y(j,i)
      end do
   end if

   if(column(i).eq.3) then
      do j=1,ninput(i)
         read(10,*) x(j,i), dummy, y(j,i)
      end do
   end if

   close(10)
end do

! first shift
do i=1,ninputfiles
   write(6,*) "so far"
      x(1:ninput(i),i) = x(1:ninput(i),i) - shift(i)
      write(6,*) "so far"
end do

write(6,*) "so far"

!find xmin,xmax
if(range.eq.1) then 
   xmin = start
   xmax = end
else 
   flag=0
   do i=1,ninputfiles
      xmin_tmp = minval(x(1:ninput(i),i))
      xmax_tmp = maxval(x(1:ninput(i),i))
      if (flag.eq.0) then
         xmin =xmin_tmp
         xmax =xmax_tmp
         flag=1
      end if
      if(xmin_tmp.lt.xmin) then
         xmin=xmin_tmp
      end if
      if(xmax_tmp.gt.xmax) then
         xmax=xmax_tmp
      end if
   end do
end if

call linspace(xnew,xmin,xmax,noutput)

print *, xmin, xmax

! take away points outside of new x_range
if(range .eq. 1) then
   do i=1,ninputfiles
      x_tmp = 0
      y_tmp = 0
      ind =0
      do j=1,ninput(i)
         if(x(j,i) .le. xmax .and. x(j,i) .ge. xmin) then
            ind = ind + 1            
            x_tmp(ind) = x(j,i)
            y_tmp(ind) = x(j,i)
         end if
      end do
      x(:,i) = x_tmp
      y(:,i) = y_tmp
      ninput(i) = ind
   end do
end if


ynew2=0
! put graphs in new range and add them
do i=1,ninputfiles
   write(6,*) ninput(i)
   call spline_easy(x(1:ninput(i),i),y(1:ninput(i),i),ninput(i), xnew, &
          ynew2(:,i) ,noutput)
end do

ynew=0
!add alla graphs 
do i=1,ninputfiles
   do j=1,noutput
      ynew(j) = ynew(j) + amount(i)*ynew2(j,i)
   end do
end do


! if normalize=1 area normalize spectrum to 1
if (normalize .eq. 1) then
   write(6,*) "normalizing spectrum"
   ynew = ynew / (sum(ynew)*(xnew(2)-xnew(1))) 
end if
     



! write to file
open(10,file=outfile,status='unknown')
do i=1,noutput
   write(10,*) xnew(i), ynew(i)
end do
close(10)

end program addmanygraphs
