subroutine spline_easy2d(x1,x2, y, M,N, xx1, xx2, y2,M2,N2)
  implicit none
  
  !passed variables
  integer,intent(in):: M,N,M2,N2
  real(8), intent(in), dimension(M):: x1
  real(8), intent(in), dimension(N):: x2
  real(8), intent(in), dimension(M2):: xx1
  real(8), intent(in), dimension(N2):: xx2
  real(8), intent(in), dimension(M,N) :: y
  real(8), intent(out),dimension(M2,N2) :: y2
  
  !local variables
  real(8),dimension(M,N):: yder
  integer::i,j
  
  ! This routine works like octaves spline function. 
  !
  !   x1 :  input array of tabulated x-values in direction 1
  !   x2 :  input array of tabulated x-values in direction 2
  !   y :  input 2d array of tabulated y-values
  !   M,N : size of x,x1 and of y 
  !   xx1 :  input array of new x-values in direction 1
  !   xx2 :  input array of new x-values in direction 2
  !   y2 :   2d array of new interpolated y-values
  !   M2,N2 : size of xx1,xx2 and of y2 
  ! 
  !   The subroutine uses splie2.f and splin2.f from numerical recepies
  !   modified for real(8) 
  !   
  !   Mathias Ljungberg, 2010-01-05
  !
  
  call splie2(x1,x2,y,M,N,yder)


  ! call splin2
  do i=1,M2
     do j=1,N2
        call splin2(x1,x2, y, yder,M,N ,xx1(i),xx2(j), y2(i,j))
     end do
  end do

end subroutine spline_easy2d
