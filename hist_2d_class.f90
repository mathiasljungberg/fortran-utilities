! last changed 2009-08-10
! Warning, not fully tested

module hist_2d_class
  use parameters
  implicit none
  type hist_2d
     integer:: nbins1, nbins2
     real(kind=wp),dimension(:),allocatable :: x1,x2
     real(kind=wp),dimension(:,:),allocatable :: y     
     real(kind=wp):: start1,end1, interval1, dx1, start2,end2, interval2, dx2
  end type hist_2d

contains

!constructor
  subroutine hist_2d_init(a, nb1, start1, end1, nb2, start2, end2)
    !passed variables
    type (hist_2d),intent(out):: a
    integer,intent(in):: nb1,nb2
    real(kind=wp):: start1, end1, start2, end2

    !local variables
    integer:: i1,i2

    !set class variables
    a%nbins1=nb1
    a%nbins2=nb2
    
    allocate( a%x1(nb1), a%x2(nb2), a%y(nb1,nb2) )
  
    a%y=0
    a%start1 = start1
    a%end1 = end1
    a%interval1= end1-start1
    a%dx1 = (end1 - start1)/dfloat(nb1-1)

    a%start2 = start2
    a%end2 = end2
    a%interval2= end2-start2
    a%dx2 = (end2 - start2)/dfloat(nb2-1)

    ! points defined at left of interval
    do i1=1,nb1
       a%x1(i1)= start1 + (i1-1)*a%dx1
    end do

    do i2=1,nb2
       a%x2(i2)= start2 + (i2-1)*a%dx2
    end do


  end subroutine hist_2d_init

!destructor
  subroutine hist_2d_final(a)
    !passed variables
    type (hist_2d),intent(out):: a
 
    deallocate(a%x1)
    deallocate(a%x2)
    deallocate(a%y)

  end subroutine hist_2d_final

! add point to hist_2dogram
  subroutine hist_2d_add(a, xin1,xin2, yin)
    !passed variables
    type (hist_2d),intent(inout):: a
    real(kind=wp):: xin1,xin2, yin
    
    !local variables
    integer:: rx1, rx2
    
    rx1 =int( (xin1 -a%start1 )/ a%dx1) + 1 
    rx2 =int( (xin2 -a%start2 )/ a%dx1) + 1 
    
    ! force them inte a bin
    if(rx1.lt.1) rx1=1
    if(rx1.ge.a%nbins1) rx1=a%nbins1

    if(rx2.lt.1) rx2=1
    if(rx2.ge.a%nbins2) rx2=a%nbins2

    
    !add intensities
    a%y(rx1,rx2) = a%y(rx1,rx2)+ yin 
   
  end subroutine hist_2d_add
  
  ! write hist_2dogram to file
  subroutine hist_2d_write(a, outfile)
    !passed variables
    type (hist_2d),intent(inout):: a
    character(LEN=*):: outfile
    
    !local variables
    integer:: i1,i2
    
    open(10,file=outfile,status='unknown')
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          write(10,'((2ES16.6E3) (ES16.6E3))') a%x1(i1), a%x2(i2), a%y(i1,i2)
       end do
       write(10,*)
    end do

    close(10)

  end subroutine hist_2d_write

  subroutine hist_2d_broadening(a,fwhm)
    !passed variables
    type (hist_2d),intent(inout):: a
    real(kind=wp)::fwhm
    !local variables
    real(kind=wp), dimension(:,:), allocatable::intenshi
    real(kind=wp):: alpha, sum1, sum2
    integer::i1,i2,j1,j2

    ! calculate sum of initial hist_2dogram
    sum1=0
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          sum1 = sum1 + a%y(i1,i2) 
       end do
    end do

    allocate( intenshi(a%nbins1,a%nbins2) )
    intenshi=0
    alpha=4.0_wp*log(2.0_wp)/(fwhm**2.0_wp)

    !broadening with 2d gaussian
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          do j1=1,a%nbins1
             do j2=1,a%nbins2
                intenshi(i1,i2) = intenshi(i1,i2) &
               + a%y(j1,j2) * (alpha / pi)**0.5_wp * exp(-alpha*(a%x1(i1)-a%x1(j1))**2.0_wp -alpha*(a%x2(i2)-a%x2(j2))**2.0_wp)
             end do
          end do
       end do
    end do


    ! calculate sum of broadened hist_2dogram
    sum2=0
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          sum2 = sum2 + a%y(i1,i2) 
       end do
    end do

    ! set hist_2dogram to normalized broadened hist_2dogram
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          a%y(i1,i2) = intenshi(i1,i2) *(sum1/sum2)
       end do
    end do

    deallocate(intenshi)
    

  end subroutine hist_2d_broadening


  subroutine normalize_integral(a, norm)
    !passed variables
    type (hist_2d),intent(inout):: a
    real(kind=wp):: norm
    
    !local variables
    integer:: i1,i2
    real(kind=wp):: sum2
    
    !calculate integral      
    sum2=0
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          sum2 = sum2 + a%y(i1,i2) 
       end do
    end do

    sum2=sum2 * a%dx1 * a%dx2 
    
    !normalize
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          a%y(i1,i2) = a%y(i1,i2)*(norm/sum2) 
       end do
    end do

  end subroutine normalize_integral
  
  subroutine normalize_sum(a, norm)
    !passed variables
    type (hist_2d),intent(inout):: a
    real(kind=wp):: norm
    
    !local variables
    integer:: i1,i2
    real(kind=wp):: sum2
    
    !calculate integral      
    sum2=0
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          sum2 = sum2 + a%y(i1,i2) 
       end do
    end do
    
    !normalize
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          a%y(i1,i2) = a%y(i1,i2)*(norm/sum2) 
       end do
    end do
    
  end subroutine normalize_sum
  
  real(kind=wp) function get_sum(a)
    !passed variables
    type (hist_2d),intent(inout):: a
    
    !local variables
    integer:: i1,i2
    real(kind=wp):: sum
    
    !calculate integral      
    sum=0
    do i1=1,a%nbins1
       do i2=1,a%nbins2
          sum = sum + a%y(i1,i2) 
       end do
    end do    

    get_sum=sum
    
  end function get_sum
  
  
  

end module hist_2d_class
