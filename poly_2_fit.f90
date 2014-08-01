subroutine poly_2_fit(x_inp, y_inp, m, start, coeff)
  use fcn3params
  implicit none
  
  !passed variables
  integer, intent(in)::m
  real(8),dimension(m),intent(in)::x_inp, y_inp
  real(8),dimension(3),intent(in):: start
  real(8),dimension(3),intent(out):: coeff
  !real(8),intent(in):: startx, starty
  
  
  ! local variables
  integer:: n,ldfjac,info,lwa
  integer,allocatable:: ipvt(:)
  real(8):: tol
  real(8),allocatable:: x(:),fvec(:),fjac(:,:),wa(:)
  external poly2
 ! real(8)::dummy,my,freq,startx,starty
 

  ! x_inp is the input array, x_in is located in the module fcn3params

  n=3 ! andragrads
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
  x(1) = start(1) !y1(1) !0.962
  x(2) = start(2)    
  x(3) = start(3) !y2(1) !-76.458
     
    
  call lmder1(poly2,m,n,x,fvec,fjac,ldfjac,tol,info, &
       ipvt,wa,lwa)
  
  ! return values
  coeff(1) = x(1)
  coeff(2) = x(2)
  coeff(3) = x(3)

! deallocate
  deallocate( ipvt,x,fvec,fjac,wa)
  
! deallocate in module fcn2params
  deallocate( x_in,y_in)

end subroutine poly_2_fit


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
