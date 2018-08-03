program morselmdr
  use fcn2params
  implicit none
  integer:: m,n,ldfjac,info,lwa
!  parameter (m=8,n=3,ldfjac=m,lwa=30)
  integer,allocatable:: ipvt(:)
  real(8):: tol
  real(8),allocatable:: x(:),fvec(:),fjac(:,:),wa(:)
  external fcn2,fcn3,fcn4,fcn5,fcn6,fcn7

!  real(8),allocatable :: y1(:),y2(:)
  real(8)::dummy,my,freq,p2,p3,p4,omega,corr3,corr4,startx,starty,my_in
  integer:: i,alt,alt2,alt3
  
  
  ! read the input file
  !read(5,*) alt
  
  ! read the input file
  read(5,*) alt
  read(5,*) alt2, startx, starty
  read(5,*) alt3, my_in  
  read(5,*) m

  ! alt=0, fit all parameters, alt=1, don't fit dist, alt=2 don't fit energy
  !if(alt.eq.1.or.alt.eq.2) then
  !   read(5,*) dist
  !end if
  !if(alt.eq.2) then
  !   read(5,*) energy
  !end if
  
  !read(5,*) m
 
  ! set up parameters
  !  m=8
  if(alt.eq.0) n=4  ! full morse
  if(alt.eq.1) n=3  ! don't fit dist
  if(alt.eq.2) n=2  ! dont fit dist or energy
  if(alt.eq.3) n=6  ! 4-deg polynomial in morse var, all parameters fitted 
  if(alt.eq.4) n=5  ! 3-deg polynomial in morse var, all parameters fitted 
  if(alt.eq.5) n=5  ! 3-deg polynomial in morse var, all parameters fitted, only 2,4 
  
 
  ldfjac=m
  lwa= 5*n+m
  
  ! allocate
  allocate( ipvt(n),x(n),fvec(m),fjac(ldfjac,n),wa(lwa))
  
  ! allocate in module fcn2params
  allocate( y1(m),y2(m))
  
  ! continue to read the input file
  do i=1,m        
     read(5,*) dummy,y1(i),y2(i) 
  enddo
 
!     the program calculates frequencies for both OH and OD stretches
!     lwa>5*n+m, ldfjac =m

      tol=1.0e-6


!     start guesses

      x(1)=0.10
      x(2)=2.5

      if(alt.eq.0) then 
         
         if(alt2.eq.0) then
            x(3) = y1(1)
            x(4) = y2(1)  
         else if (alt2.eq.1) then
            x(3) = startx
            x(4) = starty
         end if
         
         call lmder1(fcn2,m,n,x,fvec,fjac,ldfjac,tol,info, &
              ipvt,wa,lwa)
         
         ! not supported anymore
      else if(alt.eq.1) then
         
         x(3)=y2(1)
         
         call lmder1(fcn3,m,n,x,fvec,fjac,ldfjac,tol,info, &
              ipvt,wa,lwa)
         
         ! not supported
      else if(alt.eq.2) then
         
         call lmder1(fcn4,m,n,x,fvec,fjac,ldfjac,tol,info, &
              ipvt,wa,lwa)
         
      
      ! poly4 in morse var
   else if(alt.eq.3) then
      
      if(alt2.eq.0) then
         x(3) = y1(1)
         x(6) = y2(1)  
      else if (alt2.eq.1) then
         x(3) = startx
         x(6) = starty
      end if
      
      x(4) = 0.1
      x(5) = 0.1

      call lmder1(fcn5,m,n,x,fvec,fjac,ldfjac,tol,info, &
           ipvt,wa,lwa)
      
      !poly3 i morse var
   else if(alt.eq.4) then
      
      if(alt2.eq.0) then
         x(3) = y1(1)
         x(5) = y2(1)  
      else if (alt2.eq.1) then
         x(3) = startx
         x(5) = starty
      end if
      
      x(4) = 0.1
      
      call lmder1(fcn6,m,n,x,fvec,fjac,ldfjac,tol,info, &
           ipvt,wa,lwa)
    
      !poly4 i morse var, 2,4
   else if(alt.eq.5) then
      
      if(alt2.eq.0) then
         x(3) = y1(1)
         x(5) = y2(1)  
      else if (alt2.eq.1) then
         x(3) = startx
         x(5) = starty
      end if
      
      x(4) = 0.1
      
      call lmder1(fcn7,m,n,x,fvec,fjac,ldfjac,tol,info, &
           ipvt,wa,lwa)
      
      
   endif
   


   ! get polynomial coefficients

   if(alt.eq.3) then
      ! for polynomial 4 in morse variables
      call coeff4(x(1),x(2),x(4),x(5),p2,p3,p4)   
   
   else if(alt.eq.4) then
      ! for polynomial 3 in morse variables
      call coeff4(x(1),x(2),x(4),0,p2,p3,p4)   
      
   else if(alt.eq.5) then
      ! for polynomial 3 in morse variables 2,4 
      call coeff4(x(1),x(2),0,x(4),p2,p3,p4)   
      
   else

      call coeff(x(1),x(2),p2,p3,p4)
      
   end if


   !     f r equency for OH stretch
   !      m y=18./19. 
   !    my=1.0672
   
   my=18./19.
   
   ! alt3=1, use my_inp as reduced mass
   if(alt3.eq.1) my = my_in

   

       !call morsefreq(x(1),x(2),my,freq)
      call perturb(p2,p3,p4,my,freq,omega,corr3,corr4)
  
      write(6,*) 'x',x      
      write(6,'(A24,3ES20.10)') 'Polynomial coefficients',p2, p3, p4
      write(6,*)
      write(6,*) 'Frequency, OH',freq, 'cm-1'

      call harmfreq(x(1),x(2),my,freq)
      write(6,*) 'Frequency, OH, harm',freq, 'cm-1'
      write(6,'(A21,3ES20.10)') 'Contribution 2,3,4 OH',omega, corr3,corr4

!     frequency for OD stretch
      my=34./19. 
      !call morsefreq(x(1),x(2),my,freq)
      call perturb(p2,p3,p4,my,freq,omega,corr3,corr4)
      
      write(6,*)     
      write(6,*) 'Frequency, OD',freq, 'cm-1'
      
      call harmfreq(x(1),x(2),my,freq)
      write(6,*) 'Frequency, OD, harm',freq, 'cm-1'
      write(6,'(A21,3ES20.10)') 'Contribution 2,3,4 OD',omega, corr3,corr4
   
     end program morselmdr



      subroutine fcn2(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn2params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision y1(m),y2(m)
      integer:: i

!     This function fits all parameters

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 + x(4) &
             - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
         fjac(i,1)=(1-exp(-x(2)*(y1(i)-x(3) )))**2

         fjac(i,2)=2*x(1)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
             exp(-x(2)*(y1(i)-x(3))))


         fjac(i,3)=-2*x(1)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3))))

         fjac(i,4)=1
         enddo
      endif

      return
    end subroutine fcn2

   subroutine fcn3(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn2params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision y1(m),y2(m)
      integer:: i

!     this function doesn't fit the equilibrium distance

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-dist )))**2 + x(3) &
             - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
         fjac(i,1)=(1-exp(-x(2)*(y1(i)-dist )))**2

         fjac(i,2)=2*x(1)*exp(-x(2)*(y1(i)-dist))*(y1(i)-dist)*(1- &
             exp(-x(2)*(y1(i)-dist)))


!         fjac(i,3)=-2*x(1)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1-
!     +        exp(-x(2)*(y1(i)-x(3))))

         fjac(i,3)=1
         enddo
      endif

      return
    end subroutine fcn3


   subroutine fcn4(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn2params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision y1(m),y2(m)
      integer:: i

!     this function doesn't fit the equilibrium distance nor the
!     equilibrium energy


!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-dist )))**2 + energy &
             - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
         fjac(i,1)=(1-exp(-x(2)*(y1(i)-dist )))**2

         fjac(i,2)=2*x(1)*exp(-x(2)*(y1(i)-dist))*(y1(i)-dist)*(1- &
             exp(-x(2)*(y1(i)-dist)))


!         fjac(i,3)=-2*x(1)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1-
!     +        exp(-x(2)*(y1(i)-x(3))))

!         fjac(i,3)=1
         enddo
      endif

      return

    end subroutine fcn4


 subroutine fcn5(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn2params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision y1(m),y2(m)
      integer:: i

!     This function fits all parameters, a polynomial of degree 4 in the morse varibles

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 &
              + x(4)*(1-exp(-x(2)*(y1(i)-x(3) )))**3 &
              + x(5)*(1-exp(-x(2)*(y1(i)-x(3) )))**4 + x(6) &
             - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
         fjac(i,1)=(1-exp(-x(2)*(y1(i)-x(3) )))**2
         
         fjac(i,2)=2*x(1)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
              exp(-x(2)*(y1(i)-x(3)))) &
              + 3*x(4)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
              exp(-x(2)*(y1(i)-x(3))))**2  &
              + 4*x(5)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
              exp(-x(2)*(y1(i)-x(3))))**3 

         fjac(i,3)= - 2*x(1)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3)))) &
              -  3*x(4)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3))))**2  &
              -  4*x(5)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3))))**3 
         
         fjac(i,4)=(1-exp(-x(2)*(y1(i)-x(3) )))**3
         fjac(i,5)=(1-exp(-x(2)*(y1(i)-x(3) )))**4
        
         fjac(i,6)=1

         enddo
      endif

      return
    end subroutine fcn5

 subroutine fcn6(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn2params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision y1(m),y2(m)
      integer:: i

!     This function fits all parameters, a polynomial of degree 3 in the morse varibles

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 &
              + x(4)*(1-exp(-x(2)*(y1(i)-x(3) )))**3 &
              + x(5) &
             - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
         fjac(i,1)=(1-exp(-x(2)*(y1(i)-x(3) )))**2
         
         fjac(i,2)=2*x(1)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
              exp(-x(2)*(y1(i)-x(3)))) &
              + 3*x(4)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
              exp(-x(2)*(y1(i)-x(3))))**2 

         fjac(i,3)= - 2*x(1)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3)))) &
              -  3*x(4)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3))))**2 
         
         fjac(i,4)=(1-exp(-x(2)*(y1(i)-x(3) )))**3
         fjac(i,5)=1

         enddo
      endif

      return
    end subroutine fcn6



 subroutine fcn7(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn2params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision y1(m),y2(m)
      integer:: i

!     This function fits all parameters, a polynomial of degree 3 in the morse varibles
!     only 2, 4 
!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 &
              + x(4)*(1-exp(-x(2)*(y1(i)-x(3) )))**4 &
              + x(5) &
             - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
         fjac(i,1)=(1-exp(-x(2)*(y1(i)-x(3) )))**2
         
         fjac(i,2)=2*x(1)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
              exp(-x(2)*(y1(i)-x(3)))) &
              + 4*x(4)*exp(-x(2)*(y1(i)-x(3)))*(y1(i)-x(3))*(1- &
              exp(-x(2)*(y1(i)-x(3))))**3 

         fjac(i,3)= - 2*x(1)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3)))) &
              -  4*x(4)*exp(-x(2)*(y1(i)-x(3)))*x(2)*(1- &
              exp(-x(2)*(y1(i)-x(3))))**3 
         
         fjac(i,4)=(1-exp(-x(2)*(y1(i)-x(3) )))**4
         fjac(i,5)=1

         enddo
      endif

      return
    end subroutine fcn7
