program morselmdr
  use fcn2params
  use fcn8params
  use parameters
  implicit none
  integer:: m,n,ldfjac,info,lwa
!  parameter (m=8,n=3,ldfjac=m,lwa=30)
  integer,allocatable:: ipvt(:)
  real(8):: tol
  real(8),allocatable:: x(:),fvec(:),fjac(:,:),wa(:)
  external fcn8 !2,fcn3,fcn4,fcn5,fcn6,fcn7
  character(40)::inputfile

!  real(8),allocatable :: y1(:),y2(:)
  real(8)::dummy,my,freq,p2,p3,p4,omega,corr3,corr4,startx,starty,my_in
  real(8):: infreq,a
  integer:: i,alt,alt2,alt3
    
  ! read the input file
  read(5,*)  inputfile
  read(5,*)  infreq
  read(5,*)  alt2, startx, starty
  read(5,*)  alt3, my_in  
  read(5,*)  m

  n=3
  ldfjac=m
  lwa= 5*n+m


  ! allocate
  allocate( ipvt(n),x(n),fvec(m),fjac(ldfjac,n),wa(lwa))
  
  ! allocate in module fcn2params
  allocate( y1(m),y2(m))
  
  ! read the inputfile
  open(10,file=inputfile,status='old')
  do i=1,m        
     read(10,*) y1(i),y2(i) 
  enddo
  close(10)

  ! set parameters in fcn8params
  my=18./19.
  if(alt3.eq.1) my = my_in

  fcn8_my = my 
  fcn8_omega = infreq / (hbar * cm)
  fcn8_k = sqrt(amu / hartree)*1.0d-10



 
!     fit only D and r0, fix a to omega*sqrt(m/2D)
!     lwa>5*n+m, ldfjac =m

      tol=1.0e-6
!     start guesses

      x(1)=0.10   ! D
      if(alt2.eq.0) then
         x(2) = y1(1)
         x(3) = y2(1)  
      else if (alt2.eq.1) then
         x(2) = startx
         x(3) = starty
      end if
         
      call lmder1(fcn8,m,n,x,fvec,fjac,ldfjac,tol,info, &
           ipvt,wa,lwa)
        

      a=fcn8_omega*sqrt(fcn8_my/(2 * x(1)))*fcn8_k
      ! get polynomial coefficients
      call coeff(x(1),a,p2,p3,p4)


   !     f r equency for OH stretch
   !      m y=18./19. 
   !    my=1.0672
   

      call morsefreq(x(1), a,my,freq)       
      write(6,*) freq

      call perturb(p2,p3,p4,my,freq,omega,corr3,corr4)
  
      write(6,*) 'x',x      
      write(6,*) 'a',a
      write(6,'(A24,3ES20.10)') 'Polynomial coefficients',p2, p3, p4
      write(6,*)
      write(6,*) 'Frequency, OH',freq, 'cm-1'

      call harmfreq(x(1),a,my,freq)
      write(6,*) 'Frequency, OH, harm',freq, 'cm-1'
      write(6,'(A21,3ES20.10)') 'Contribution 2,3,4 OH',omega, corr3,corr4

!     frequency for OD stretch
      my=34./19. 
      !call morsefreq(x(1),x(2),my,freq)
      call perturb(p2,p3,p4,my,freq,omega,corr3,corr4)
      
      write(6,*)     
      write(6,*) 'Frequency, OD',freq, 'cm-1'
      
      call harmfreq(x(1),a,my,freq)
      write(6,*) 'Frequency, OD, harm',freq, 'cm-1'
      write(6,'(A21,3ES20.10)') 'Contribution 2,3,4 OD',omega, corr3,corr4
   
     end program morselmdr


    subroutine fcn8(m,n,x,fvec,fjac,ldfjac,iflag)
      use fcn2params
      use fcn8params

      implicit none
      integer:: m,n,ldfjac,iflag
      real(8):: x(n),fvec(m),fjac(ldfjac,n),a,test
!
!      double precision y1(m),y2(m)
      integer:: i

!     This function fixes the frequency 

!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     parametrarna x= (D,r0,c)



      if(iflag == 1) then

         a = fcn8_omega*sqrt(fcn8_my/(2 * x(1)))*fcn8_k   
         do i=1,m

         fvec(i)=x(1)*(1-exp(- a * (y1(i)-x(2) )))**2 + x(3) &
             - y2(i)

         enddo
   
      endif

!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.      

      if(iflag == 2) then

         a = fcn8_omega*sqrt(fcn8_my/(2 * x(1)))*fcn8_k   
         do i=1,m
            
            fjac(i,1)=(1-exp(-a *(y1(i)-x(2) )))**2 - &
                 exp(- a * (y1(i)-x(2))) * (1- &
                 exp(-a * (y1(i)-x(2)))) * a * (y1(i)-x(2)) 

            fjac(i,2)=-2 * x(1) * exp(- a * (y1(i)-x(2))) * a * (1- &
                 exp(-a * (y1(i)-x(2))))

            fjac(i,3)=1

         enddo
      endif

      return
    end subroutine fcn8
