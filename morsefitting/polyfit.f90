program polyfit
  use fcn2params
  implicit none
  integer:: m,n,ldfjac,info,lwa
!  parameter (m=8,n=3,ldfjac=m,lwa=30)
  integer,allocatable:: ipvt(:)
  real(8):: tol
  real(8),allocatable:: x(:),fvec(:),fjac(:,:),wa(:)
  external poly4,poly6,poly8
!  real(8),allocatable :: y1(:),y2(:)
  real(8)::dummy,my,freq,omega,p2,p3,p4,corr3,corr4,startx,starty
  integer:: i,alt,alt2
  

  ! read the input file
  read(5,*) alt
  read(5,*) alt2, startx, starty
  read(5,*) m
 
  ! alt =0, poly4, alt=1 ploy6, alt2=0, startvalues y1(1),y2(1)
  ! alt2 =1 startvalues startx, starty

  ! set up parameters, n= number of unknowns
  !  m=8
  if(alt.eq.0) n=5 ! fjärdegrads
  if(alt.eq.1) n=7 ! sjättegrads
  if(alt.eq.2) n=9 ! åttondegrads
  !if(alt.eq.2) n=2

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

      if(alt.eq.0) then
         
!     start guesses

         x(2)=0.1
         x(3)=0.1
         x(4)=0.1
         
         if(alt2.eq.0) then
            x(1) = y1(1) !0.962
            x(5) = y2(1) !-76.458
         else if (alt2.eq.1) then
            x(1) = startx
            x(5) = starty
         end if
         
         
         call lmder1(poly4,m,n,x,fvec,fjac,ldfjac,tol,info, &
              ipvt,wa,lwa)

      else if(alt.eq.1) then

!        start guesses

         x(2)=0.1
         x(3)=0.1
         x(4)=0.1
         x(5)=0.1
         x(6)=0.1

         if(alt2.eq.0) then
            x(1) = y1(1) !0.962
            x(7) = y2(1) !-76.458
         else if (alt2.eq.1) then
            x(1) = startx
            x(7) = starty
         end if
         
         call lmder1(poly6,m,n,x,fvec,fjac,ldfjac,tol,info, &
              ipvt,wa,lwa)


      else if(alt.eq.2) then
         
!        start guesses
         
         x(2)=0.1
         x(3)=0.1
         x(4)=0.1
         x(5)=0.1
         x(6)=0.1
         x(7)=0.1
         x(8)=0.1

         if(alt2.eq.0) then
            x(1) = y1(1) !0.962
            x(9) = y2(1) !-76.458
         else if (alt2.eq.1) then
            x(1) = startx
            x(9) = starty
         end if
         
         call lmder1(poly8,m,n,x,fvec,fjac,ldfjac,tol,info, &
              ipvt,wa,lwa)

      end if

!     frequency for OH stretch
      my=18./19. 
!      call morsefreq(x(1),x(2),my,freq)            
      write(6,*) 'x',x   
      write(6,'(A24,3ES20.10)') 'Polynomial coefficients',x(2), x(3), x(4)
      write(6,*)

      call perturb(x(2),x(3),x(4),my,freq,omega,corr3,corr4)
      write(6,*) 'Frequency, OH',freq, 'cm-1'
      write(6,*) 'Frequency, OH, harm',omega, 'cm-1'
      write(6,'(A21,3ES20.10)') 'Contribution 2,3,4 OH', omega, corr3, corr4
      write(6,*)
!     frequency for OD stretch
      my=34./19. 
      
      call perturb(x(2),x(3),x(4),my,freq,omega,corr3,corr4)
      write(6,*) 'Frequency, OD',freq, 'cm-1'
      write(6,*) 'Frequency, OD, harm',omega, 'cm-1'

      write(6,'(A21,3ES20.10)') 'Contribution 2,3,4 OD', omega, corr3, corr4


    end program polyfit




    subroutine poly4(m,n,x,fvec,fjac,ldfjac,iflag)
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

!         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 + x(4) &
!             - y2(i)

         fvec(i)= x(2)*(y1(i)-x(1))**2  + x(3)*(y1(i)-x(1))**3 &
              + x(4)*(y1(i)-x(1))**4  + x(5) &
              - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(y1(i)-x(1)) - 3.*x(3)*(y1(i)-x(1))**2 &
                 - 4.*x(4)*(y1(i)-x(1))**3
            fjac(i,2)= (y1(i)-x(1))**2
            fjac(i,3)= (y1(i)-x(1))**3
            fjac(i,4)= (y1(i)-x(1))**4
            fjac(i,5)=1.

         enddo

     endif

      return
    end subroutine poly4

   subroutine poly6(m,n,x,fvec,fjac,ldfjac,iflag)
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

!         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 + x(4) &
!             - y2(i)

         fvec(i)= x(2)*(y1(i)-x(1))**2  + x(3)*(y1(i)-x(1))**3 &
              + x(4)*(y1(i)-x(1))**4  + x(5)*(y1(i)-x(1))**5   &
              + x(6)*(y1(i)-x(1))**6 + x(7)                    &
              - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(y1(i)-x(1)) - 3.*x(3)*(y1(i)-x(1))**2 &
                 - 4.*x(4)*(y1(i)-x(1))**3 - 5.*x(5)*(y1(i)-x(1))**4 &
                 -6.*x(6)*(y1(i)-x(1))**5
            fjac(i,2)= (y1(i)-x(1))**2
            fjac(i,3)= (y1(i)-x(1))**3
            fjac(i,4)= (y1(i)-x(1))**4
            fjac(i,5)= (y1(i)-x(1))**5
            fjac(i,6)= (y1(i)-x(1))**6
            fjac(i,7)=1.

         enddo

     endif

      return
    end subroutine poly6




 subroutine poly8(m,n,x,fvec,fjac,ldfjac,iflag)
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

!         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 + x(4) &
!             - y2(i)

         fvec(i)= x(2)*(y1(i)-x(1))**2  + x(3)*(y1(i)-x(1))**3 &
              + x(4)*(y1(i)-x(1))**4  + x(5)*(y1(i)-x(1))**5   &
              + x(6)*(y1(i)-x(1))**6 + x(7)*(y1(i)-x(1))**7    &
              + x(8)*(y1(i)-x(1))**8  +  x(9)                  &
              - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(y1(i)-x(1)) - 3.*x(3)*(y1(i)-x(1))**2 &
                 - 4.*x(4)*(y1(i)-x(1))**3 - 5.*x(5)*(y1(i)-x(1))**4 &
                 -6.*x(6)*(y1(i)-x(1))**5 -7.*x(7)*(y1(i)-x(1))**6 &
                 -8.*x(8)*(y1(i)-x(1))**7
            fjac(i,2)= (y1(i)-x(1))**2
            fjac(i,3)= (y1(i)-x(1))**3
            fjac(i,4)= (y1(i)-x(1))**4
            fjac(i,5)= (y1(i)-x(1))**5
            fjac(i,6)= (y1(i)-x(1))**6
            fjac(i,7)= (y1(i)-x(1))**7
            fjac(i,8)= (y1(i)-x(1))**8
            fjac(i,9)=1.

         enddo

     endif

      return
    end subroutine poly8
