subroutine poly_n_fit(degree, x_inp, y_inp, m, start, coeff)
  use fcn3params
  implicit none
  
  !passed variables
  integer, intent(in)::m, degree
  real(8),dimension(m),intent(in)::x_inp, y_inp
  real(8),dimension(degree+1),intent(in):: start
  real(8),dimension(degree+1),intent(out):: coeff
  !real(8),intent(in):: startx, starty
  
  
  ! local variables
  integer:: n,ldfjac,info,lwa
  integer,allocatable:: ipvt(:)
  real(8):: tol
  real(8),allocatable:: x(:),fvec(:),fjac(:,:),wa(:)
  external poly2,poly4,poly6,poly8,poly10,poly12

  !loop variables
  integer::i

  ! real(8)::dummy,my,freq,startx,starty
 
  ! degree: integer, degree of polynomial fitted, 2,4,6, or8
  ! x_inp: real(8) array of dimension m, the x values of the points to be fitted
  ! y_inp: real(8) array of dimension m, the x values of the points to be fitted
  ! m: integer, the number of points to be fitted
  ! start: real(8), dimension degree+1, the starting values of the coefficients of
  !        the polynomial, x0 =x(1), c=x(degree+1),
  ! coeff: real(8) dimension degree +1, output of the fitting coeficients
  !
  ! x_in, y_in is located in the module fcn3params
  ! 

  n = degree +1 ! andragrads
  tol=1.0e-6
  ldfjac=m
  lwa= 5*n+m

 ! allocate
  allocate( ipvt(n),x(n),fvec(m),fjac(ldfjac,n),wa(lwa))
  
! allocate in module fcn2params
  allocate( x_in(m),y_in(m))
  
  x_in = x_inp
  y_in = y_inp

  !     start guesses
  !x(1) = start(1) !y1(1) !0.962
  !x(2) = start(2)    
  !x(3) = start(3) !y2(1) !-76.458
     

  do i=1,n
     x(i)=start(i)
  end do

  if(degree.eq.2) then
     call lmder1(poly2,m,n,x,fvec,fjac,ldfjac,tol,info, &
          ipvt,wa,lwa)
     
  else if(degree.eq.4) then
     call lmder1(poly4,m,n,x,fvec,fjac,ldfjac,tol,info, &
          ipvt,wa,lwa)
  
  else if(degree.eq.6) then
     call lmder1(poly6,m,n,x,fvec,fjac,ldfjac,tol,info, &
          ipvt,wa,lwa)
  
  else if(degree.eq.8) then
     call lmder1(poly8,m,n,x,fvec,fjac,ldfjac,tol,info, &
          ipvt,wa,lwa)

  else if(degree.eq.10) then
     call lmder1(poly10,m,n,x,fvec,fjac,ldfjac,tol,info, &
          ipvt,wa,lwa)

  else if(degree.eq.12) then
     call lmder1(poly12,m,n,x,fvec,fjac,ldfjac,tol,info, &
          ipvt,wa,lwa)
  
  end if
  
  ! return values
  do i=1,n
     coeff(i)=x(i)
  end do

! deallocate
  deallocate( ipvt,x,fvec,fjac,wa)
  
! deallocate in module fcn2params
  deallocate( x_in,y_in)

end subroutine poly_n_fit


   subroutine poly2(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn3params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision y1(m),y2(m)
      integer:: i

!     This function fits all parameters
!     x is the coefficients to be fitted, x_in is the array of ingoing 
!     x-values and y_in the corresponding y-values, they are allocated in the 
!     module fcn3params

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
      if(iflag == 1) then
         do i=1,m

!         fvec(i)=x(1)*(1-exp(-x(2)*(y1(i)-x(3) )))**2 + x(4) &
!             - y2(i)

         fvec(i)= x(2)*(x_in(i)-x(1))**2 &
              + x(3) &
              - y_in(i)

         enddo

      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then
         
         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(x_in(i)-x(1)) 
            fjac(i,2)= (x_in(i)-x(1))**2
            fjac(i,3)= 1.
            
         enddo
         
      endif

      return
    end subroutine poly2


 subroutine poly4(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn3params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision x_in(m),y_in(m)
      integer:: i

!     This function fits all parameters

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

!         fvec(i)=x(1)*(1-exp(-x(2)*(x_in(i)-x(3) )))**2 + x(4) &
!             - y_in(i)

         fvec(i)= x(2)*(x_in(i)-x(1))**2  + x(3)*(x_in(i)-x(1))**3 &
              + x(4)*(x_in(i)-x(1))**4  + x(5) &
              - y_in(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(x_in(i)-x(1)) - 3.*x(3)*(x_in(i)-x(1))**2 &
                 - 4.*x(4)*(x_in(i)-x(1))**3
            fjac(i,2)= (x_in(i)-x(1))**2
            fjac(i,3)= (x_in(i)-x(1))**3
            fjac(i,4)= (x_in(i)-x(1))**4
            fjac(i,5)=1.

         enddo

     endif

      return
    end subroutine poly4

   subroutine poly6(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn3params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision x_in(m),y_in(m)
      integer:: i

!     This function fits all parameters

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

!         fvec(i)=x(1)*(1-exp(-x(2)*(x_in(i)-x(3) )))**2 + x(4) &
!             - y_in(i)

         fvec(i)= x(2)*(x_in(i)-x(1))**2  + x(3)*(x_in(i)-x(1))**3 &
              + x(4)*(x_in(i)-x(1))**4  + x(5)*(x_in(i)-x(1))**5   &
              + x(6)*(x_in(i)-x(1))**6 + x(7)                    &
              - y_in(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(x_in(i)-x(1)) - 3.*x(3)*(x_in(i)-x(1))**2 &
                 - 4.*x(4)*(x_in(i)-x(1))**3 - 5.*x(5)*(x_in(i)-x(1))**4 &
                 -6.*x(6)*(x_in(i)-x(1))**5
            fjac(i,2)= (x_in(i)-x(1))**2
            fjac(i,3)= (x_in(i)-x(1))**3
            fjac(i,4)= (x_in(i)-x(1))**4
            fjac(i,5)= (x_in(i)-x(1))**5
            fjac(i,6)= (x_in(i)-x(1))**6
            fjac(i,7)=1.

         enddo

     endif

      return
    end subroutine poly6




 subroutine poly8(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn3params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision x_in(m),y_in(m)
      integer:: i

!     This function fits all parameters

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

!         fvec(i)=x(1)*(1-exp(-x(2)*(x_in(i)-x(3) )))**2 + x(4) &
!             - y_in(i)

         fvec(i)= x(2)*(x_in(i)-x(1))**2  + x(3)*(x_in(i)-x(1))**3 &
              + x(4)*(x_in(i)-x(1))**4  + x(5)*(x_in(i)-x(1))**5   &
              + x(6)*(x_in(i)-x(1))**6 + x(7)*(x_in(i)-x(1))**7    &
              + x(8)*(x_in(i)-x(1))**8  +  x(9)                  &
              - y_in(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(x_in(i)-x(1)) - 3.*x(3)*(x_in(i)-x(1))**2 &
                 - 4.*x(4)*(x_in(i)-x(1))**3 - 5.*x(5)*(x_in(i)-x(1))**4 &
                 -6.*x(6)*(x_in(i)-x(1))**5 -7.*x(7)*(x_in(i)-x(1))**6 &
                 -8.*x(8)*(x_in(i)-x(1))**7
            fjac(i,2)= (x_in(i)-x(1))**2
            fjac(i,3)= (x_in(i)-x(1))**3
            fjac(i,4)= (x_in(i)-x(1))**4
            fjac(i,5)= (x_in(i)-x(1))**5
            fjac(i,6)= (x_in(i)-x(1))**6
            fjac(i,7)= (x_in(i)-x(1))**7
            fjac(i,8)= (x_in(i)-x(1))**8
            fjac(i,9)=1.

         enddo

     endif

      return
    end subroutine poly8 

 subroutine poly10(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn3params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision x_in(m),y_in(m)
      integer:: i

!     This function fits all parameters

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

!         fvec(i)=x(1)*(1-exp(-x(2)*(x_in(i)-x(3) )))**2 + x(4) &
!             - y_in(i)

         fvec(i)= x(2)*(x_in(i)-x(1))**2  + x(3)*(x_in(i)-x(1))**3 &
              + x(4)*(x_in(i)-x(1))**4  + x(5)*(x_in(i)-x(1))**5   &
              + x(6)*(x_in(i)-x(1))**6 + x(7)*(x_in(i)-x(1))**7    &
              + x(8)*(x_in(i)-x(1))**8  +  x(9) *(x_in(i)-x(1))**9  &
              + x(10)*(x_in(i)-x(1))**10  +  x(11)                  &
              - y_in(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(x_in(i)-x(1)) - 3.*x(3)*(x_in(i)-x(1))**2 &
                 - 4.*x(4)*(x_in(i)-x(1))**3 - 5.*x(5)*(x_in(i)-x(1))**4 &
                 -6.*x(6)*(x_in(i)-x(1))**5 -7.*x(7)*(x_in(i)-x(1))**6 &
                 -8.*x(8)*(x_in(i)-x(1))**7  -9.*x(9)*(x_in(i)-x(1))**8 &
                 -10.*x(10)*(x_in(i)-x(1))**9  

            fjac(i,2)= (x_in(i)-x(1))**2
            fjac(i,3)= (x_in(i)-x(1))**3
            fjac(i,4)= (x_in(i)-x(1))**4
            fjac(i,5)= (x_in(i)-x(1))**5
            fjac(i,6)= (x_in(i)-x(1))**6
            fjac(i,7)= (x_in(i)-x(1))**7
            fjac(i,8)= (x_in(i)-x(1))**8
            fjac(i,9)= (x_in(i)-x(1))**9
            fjac(i,10)= (x_in(i)-x(1))**10
            fjac(i,11)=1.

         enddo

     endif

      return
    end subroutine poly10 

subroutine poly12(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn3params
      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n)
!
!      double precision x_in(m),y_in(m)
      integer:: i

!     This function fits all parameters

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,alfa,r0,c)
      if(iflag == 1) then
         do i=1,m

!         fvec(i)=x(1)*(1-exp(-x(2)*(x_in(i)-x(3) )))**2 + x(4) &
!             - y_in(i)

         fvec(i)= x(2)*(x_in(i)-x(1))**2  + x(3)*(x_in(i)-x(1))**3 &
              + x(4)*(x_in(i)-x(1))**4  + x(5)*(x_in(i)-x(1))**5   &
              + x(6)*(x_in(i)-x(1))**6 + x(7)*(x_in(i)-x(1))**7    &
              + x(8)*(x_in(i)-x(1))**8  +  x(9) *(x_in(i)-x(1))**9  &
              + x(10)*(x_in(i)-x(1))**10  +  x(11) *(x_in(i)-x(1))**11  &
              + x(12)*(x_in(i)-x(1))**12  +  x(13)                  &
              - y_in(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         do i=1,m
            
            fjac(i,1)= - 2.*x(2)*(x_in(i)-x(1)) - 3.*x(3)*(x_in(i)-x(1))**2 &
                 - 4.*x(4)*(x_in(i)-x(1))**3 - 5.*x(5)*(x_in(i)-x(1))**4 &
                 -6.*x(6)*(x_in(i)-x(1))**5 -7.*x(7)*(x_in(i)-x(1))**6 &
                 -8.*x(8)*(x_in(i)-x(1))**7  -9.*x(9)*(x_in(i)-x(1))**8 &
                 -10.*x(8)*(x_in(i)-x(1))**9  -11.*x(9)*(x_in(i)-x(1))**10 &
                 -12.*x(10)*(x_in(i)-x(1))**11  

            fjac(i,2)= (x_in(i)-x(1))**2
            fjac(i,3)= (x_in(i)-x(1))**3
            fjac(i,4)= (x_in(i)-x(1))**4
            fjac(i,5)= (x_in(i)-x(1))**5
            fjac(i,6)= (x_in(i)-x(1))**6
            fjac(i,7)= (x_in(i)-x(1))**7
            fjac(i,8)= (x_in(i)-x(1))**8
            fjac(i,9)= (x_in(i)-x(1))**9
            fjac(i,10)= (x_in(i)-x(1))**10
            fjac(i,11)= (x_in(i)-x(1))**11
            fjac(i,12)= (x_in(i)-x(1))**12
            fjac(i,13)=1.

         enddo

     endif

      return
    end subroutine poly12 
