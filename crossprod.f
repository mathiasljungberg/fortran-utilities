c
c     mina egna vektor-/matrisfunktioner
c
      subroutine crossprod(A,B,C)
c A är kryssprodukten mellan B och C 
      implicit none
      real*8 A(*),B(*),C(*)
      
      A(1)= B(2)*C(3)-B(3)*C(2) 
      A(2)= -( B(1)*C(3)-B(3)*C(1) ) 
      A(3)= B(1)*C(2)-B(2)*C(1)
      end subroutine crossprod

      subroutine multmm(A,B,C) 
c  A är matrisprodukten mellan B och C     
      implicit none
      real*8 A(3,3),B(3,3),C(3,3)
      integer i,j,k
   

      do i=1,3
         do j=1,3
            A(i,j)=0
            do k=1,3
                 A(i,j)=A(i,j)+B(i,k)*C(k,j)
            enddo
         enddo
      enddo
      end subroutine multmm
      
      subroutine multmv(A,B,C) 
c  A är produkten mellan B och C (vektor)     
      implicit none
      real*8 A(3),B(3,3),C(3)
      integer i,j
     
      do i=1,3
         A(i)=0
         do j=1,3
            A(i)=A(i)+B(i,j)*C(j)
         enddo
      enddo
      
      end subroutine multmv
      

      subroutine unitmat(A)
c enhetsmatrisen 
      implicit none
      real*8 A(3,3)
      integer i,j
      do i=1,3
         do j=1,3
            if(i==j) A(i,j)=1
            if(i.ne.j) A(i,j)=0
         enddo
      enddo
      end subroutine unitmat

      subroutine multscm(A,B, s) 
c A är B multiplicerat med en skalär
      implicit none
      real*8 A(3,3),B(3,3),s
      integer i,j
      do i=1,3
         do j=1,3
            A(i,j)=B(i,j)*s
         enddo
      enddo
      end subroutine multscm

      subroutine multscv(A,B, s) 
c A är B multiplicerat med en skalär
      implicit none
      real*8 A(3),B(3),s
      integer i
      do i=1,3
         A(i)=B(i)*s
      enddo
      end subroutine multscv


      subroutine addmm(A,B,C) 
c A är B+C
      implicit none
      real*8 A(3,3),B(3,3),C(3,3)
      integer i,j
      do i=1,3
         do j=1,3
            A(i,j)=B(i,j) +C(i,j)
         enddo
      enddo
      end subroutine addmm
      
      subroutine subtrmm(A,B,C) 
c A är B+C
      implicit none
      real*8 A(3,3),B(3,3),C(3,3)
      integer i,j
      do i=1,3
         do j=1,3
            A(i,j)=B(i,j) -C(i,j)
         enddo
      enddo
      end subroutine subtrmm

      real*8  function norm(A) 
c normen av vektor A
      implicit none
      real*8 A(*)
      
      norm=sqrt(A(1)**2 +A(2)**2 +A(3)**2)
      return
      end

c      real(3) function norm2(A) 
c normen av vektor A
c      implicit none
c      real A(*)
c      
c      norm2(1)=sqrt(A(1)**2 +A(2)**2 +A(3)**2)
c      norm2(2)=0
c c     norm2(3)=0
ccc     A(1)=10000.
c      return
c      end

