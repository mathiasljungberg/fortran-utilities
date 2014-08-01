module m_histogram_array
  use m_globals, only: PREC
  use m_constants 
  use m_common_io, only: get_unit
  implicit none
  type histogram_array
     integer:: nbins, dim
     real(kind=PREC),dimension(:),allocatable :: x     
     real(kind=PREC),dimension(:,:),allocatable :: y     
     real(kind=PREC):: start,end, interval, dx
  end type histogram_array

contains

  !constructor
  subroutine histogram_array_init(a, nb, start, end, dim)
    !passed variables
    type (histogram_array),intent(out):: a
    integer,intent(in):: nb, dim
    real(kind=PREC), intent(in):: start, end

    !local variables
    integer:: i

    !set class variables
    a%nbins=nb
    a%dim =dim

    allocate( a%x(nb), a%y(nb,dim) )

    a%y=0.0_PREC
    a%start = start
    a%end = end
    a%interval= end-start
    a%dx = (end - start)/dfloat(nb-1)

    ! points defined at left of interval
    do i=1,nb
      a%x(i)= start + (i-1)*a%dx
    end do

  end subroutine histogram_array_init

  !destructor
  subroutine histogram_array_final(a)
    !passed variables
    type (histogram_array),intent(out):: a

    deallocate(a%x)
    deallocate(a%y)

  end subroutine histogram_array_final

  ! add point to histogramogram
  subroutine histogram_array_add(a, xin, yin)
    !passed variables
    type (histogram_array),intent(inout):: a
    real(kind=PREC), intent(in), dimension(:):: xin
    real(kind=PREC), intent(in):: yin


    !local variables
    integer:: rx, i

    do i=1, a%dim
      rx =int( (xin(i) -a%start )/ a%dx) + 1 
      
      ! force them inte a bin
      if(rx.lt.1) rx=1
      if(rx.ge.a%nbins) rx=a%nbins
      
      !add intensities
      a%y(rx,i) = a%y(rx,i)+ yin 
    end do

  end subroutine histogram_array_add

  ! write histogramogram to file
  subroutine histogram_array_write(a, outfile)
    !passed variables
    type (histogram_array),intent(inout):: a
    character(LEN=*), intent(in):: outfile
    
    !local variables
    integer:: i, j, funit, ios
    character(80):: string


    write(string,*) a%dim
    string = '((ES25.16) (' // trim(adjustl(string)) // 'ES25.16))'
    string = trim(adjustl(string))

    ios = 0
    call get_unit(funit,ios)
    open(funit,file=outfile,status='unknown')
    do i=1,a%nbins
      write(funit, string) a%x(i), (a%y(i,j), j = 1, a%dim)
    end do
    close(funit)

  end subroutine histogram_array_write

  subroutine histogram_array_broadening(a,fwhm)
    !passed variables
    type (histogram_array),intent(inout):: a
    real(kind=PREC), intent(in)::fwhm
    !local variables
    real(kind=PREC), dimension(:), allocatable::intenshi
    real(kind=PREC):: alpha, sum1, sum2
    integer::i,j,k

    ! calculate sum of initial histogramogram
    do k  =1, a % dim
      
      sum1=0.0_PREC
      do i=1,a%nbins
        sum1 = sum1 + a%y(i,k) 
      end do
      
      allocate( intenshi(a%nbins) )
      intenshi=0.0_PREC
      alpha=4.0_PREC*log(2.0_PREC)/(fwhm**2.0_PREC)
      
      !broadening with gaussian
      do i=1,a%nbins
        do j=1,a%nbins
          intenshi(i) = intenshi(i) &
               + a%y(j,k) * (alpha / const_pi)**0.5_PREC * exp(-alpha*(a%x(i)-a%x(j))**2.0_PREC)
        end do
      end do
      
      ! calculate sum of broadened histogramogram
      sum2=0.0_PREC
      do i=1,a%nbins
        sum2 = sum2 + intenshi(i) 
      end do
      
      ! set histogramogram to normalized broadened histogramogram
      do i=1,a%nbins
        a%y(i,k) = intenshi(i) *(sum1/sum2)
      end do
      deallocate(intenshi)
  
    end do

  end subroutine histogram_array_broadening


  subroutine histogram_array_normalize_integral(a, norm)
    !passed variables
    type (histogram_array),intent(inout):: a
    real(kind=PREC), intent(in):: norm

    !local variables
    integer:: i,k
    real(kind=PREC):: sum2

    do k =1, a%dim
      !calculate integral      
      sum2=0.0_PREC
      do i=1,a%nbins
        sum2 = sum2 + a%y(i,k) 
      end do
      sum2=sum2*a%dx

      !normalize
      do i=1,a%nbins
        a%y(i,k) = a%y(i,k)*(norm/sum2) 
      end do
    end do

  end subroutine histogram_array_normalize_integral


  subroutine histogram_array_normalize_sum(a, norm)
    !passed variables
    type (histogram_array),intent(inout):: a
    real(kind=PREC), intent(in):: norm

    !local variables
    integer:: i,k 
    real(kind=PREC):: sum2

    do k=1, a%dim

      !calculate integral      
      sum2=0.0_PREC
      do i=1,a%nbins
        sum2 = sum2 + a%y(i,k) 
      end do
      
      !normalize
      do i=1,a%nbins
        a%y(i,k) = a%y(i,k)*(norm/sum2) 
      end do
    
    end do

  end subroutine histogram_array_normalize_sum


  real(kind=PREC) function histogram_array_get_sum(a,dim)
    !passed variables
    type (histogram_array),intent(inout):: a
    integer:: dim

    !local variables
    integer:: i
    real(kind=PREC):: sum

    !calculate integral      
    sum=0.0_PREC
    do i=1,a%nbins
      sum = sum + a%y(i,dim) 
    end do
    histogram_array_get_sum = sum

  end function histogram_array_get_sum


end module m_histogram_array
