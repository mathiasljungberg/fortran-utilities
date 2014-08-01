      subroutine bubblesort(d,xin,xout,perm)
      implicit none
c variables
      integer d,perm(d) 
      real*8 xin(d),xout(d)
c local variables
      integer i,j,temp2
      real*8 temp
c xin är vektorn som ska sorteras, ut kommer xout, perm är 
c ordningen på xout relativt xin

      do i=1,d
         xout(i)=xin(i)
         perm(i)=i
      enddo
c bubblesort      
      do i = 1,d-1
         do j = i+1,d
            if (xout(j).lt.xout(i)) then
               temp = xout(i)
               temp2 =perm(i)

               xout(i) = xout(j)
               xout(j) = temp

               perm(i) = perm(j)
               perm(j) = temp2

            end if
         end do
      end do

      end subroutine bubblesort
