module m_poisson_1d_finite_diff
  implicit none
contains

  ! this routine solves [ d^2/dx + eps'(x)*d/dx + d_ext(x) ] v(x) = c * \rho(x) / eps(x) 
  ! where eps and eps' is the dielectric function and its derivative, c is a constant
  ! and left and right boundary conditions should be specified quite freely.
  ! d_ext(x) is some diagonal term coming from integrating out the other dimensions.
  ! For now, we use fixed bc, that is v(1) = v_1 and v(N) = v_N
  ! this can also treat dipole correction boundary conditions
  ! 
  subroutine poisson_1d_fd_eps(x, rho, d_ext, c, flag_use_eps, eps, boundary_cond, v_1, v_N, v)
    real(8), intent(in):: x(:)
    real(8), intent(in):: rho(:)
    real(8), intent(in):: d_ext(:)
    real(8), intent(in):: c
    logical, intent(in):: flag_use_eps
    real(8), intent(in):: eps(:)
    character(*), intent(in):: boundary_cond
    real(8), intent(in):: v_1, v_N
    real(8), intent(out):: v(:)

    real(8):: dx
    integer:: N, i

    real(8), allocatable:: diag(:), subdiag(:), superdiag(:), eps_p(:)
    real(8):: m_0, m_1
    
    ! for solver
    real(8), allocatable:: DLF(:), DF(:), DUF(:), DU2(:), FERR(:), &
         BERR(:), WORK(:)
    real(8), allocatable:: B(:,:), v_sol(:,:) 
    integer, allocatable:: IPIV(:), IWORK(:)
    real(8):: RCOND
    integer:: NRHS, LDB, LDX, INFO

    
    N = size(x)
    dx = x(2)-x(1)
    allocate(diag(N))
    allocate(subdiag(N-1))
    allocate(superdiag(N-1))
    allocate(eps_p(N))

    ! create derivative of epsilon
    eps_p(1) = (eps(2)-eps(1)) / dx
    eps_p(N) = (eps(N)-eps(N-1)) / dx
    do i=2, N-1
      eps_p(i) = (eps(i+1) -eps(i-1)) / (2.0d0 * dx)
    end do

    ! set up splver variables
    NRHS = 1
    LDB = N 
    LDX = N
    allocate(DLF(N-1), DF(N), DUF(N-1), DU2(N-2), FERR(NRHS), &
         BERR(NRHS), WORK(3*N))
    allocate( B(LDB,NRHS), v_sol(LDX,NRHS)) 
    allocate(IPIV(N), IWORK(N))
        
    ! diagonal
    
    ! contribition from second derivative
    diag(:) = diag(:) -2.0d0 / dx**2

    ! contribution from external diagonal terms
    diag(:) = diag(:) + d_ext(:) 
        
    ! sub/super-diagonalen
    subdiag= 0.0d0
    subdiag(:) = subdiag(:) + 1.0d0 / dx**2
    subdiag(:) = subdiag(:) - eps_p(2:N) / ( eps(2:N) * 2.0d0 * dx)
    
    superdiag= 0.0d0
    superdiag(:) = superdiag(:) + 1.0d0 / dx**2
    superdiag(:) = superdiag(:) + eps_p(1:N-1) / ( eps(1:N-1) * 2.0d0 * dx)
    
    ! right hand side
    B(:,1) = c * rho(:) / eps(:)  

    ! impose boundary conditions
    if (boundary_cond .eq. "explicit") then 
      B(1,1) =B(1,1) - v_1 * (1.0d0 / dx**2 + eps_p(1) / ( eps(1) * 2.0d0 * dx))
      B(N,1) =B(N,1) - v_N * (1.0d0 / dx**2 - eps_p(N) / ( eps(N) * 2.0d0 * dx))
    else if (boundary_cond .eq. "periodic") then
      ! can we wrap around and still have tridiagonal system? 
      write(6,*) "boundary_cond 'periodic' not implementend"
      stop
    else if (boundary_cond .eq. "dipole_corr") then
      ! find zeroth moment of distibution, no eps
      m_0 = sum(c * rho(:)) * dx
      write(6,*) "m_0", m_0

      ! find first moment of distibution, centered at m_0,  no eps
      m_1 = sum(x(:) * c * rho(:)) * dx
      write(6,*) "m_1", m_1
      write(6,*) "m_1/m_0", m_1/m_0

      ! set boundary condition
      B(1,1) = B(1,1) - m_1 / (2.0d0  *dx**2) 
      B(N,1) = B(N,1) + m_1 / (2.0d0  *dx**2) 
    else if (boundary_cond .eq. "zero") then
      ! do nothing
    else
      ! do nothing
    end if
    
    ! solve eigenvalue problem
    !abstol=2d0*dlamch('s')
    !call dstevx('v','i',npoints,diag,subdiag, 0.d0,1.d0, 1, nstates, abstol, &
    !     n_eigfound, eigenval , eigenvec, npoints,work,iwork,ifail,info)
    !
    !eigvec(:,1:nstates) = eigenvec(:,1:nstates)
    !eigval(1:nstates) = eigenval(1:nstates)
    
    ! solve tridiagonal system of equations
    ! SUBROUTINE DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,
    !$                   DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,
    !$                   WORK, IWORK, INFO )

    call DGTSVX('N', 'N', N, 1, subdiag, diag, superdiag, DLF, DF, DUF,&
         DU2, IPIV, B, LDB, v_sol, LDX,  RCOND, FERR, BERR, &
         WORK, IWORK, INFO )

    if (INFO > 0) then
      write(6,*) "WARNING, INFO > 0 in poisson_1d_fd_eps", INFO
    end if

    v = v_sol(:,1)
    
  end subroutine poisson_1d_fd_eps

  subroutine test_poisson_1d_fd_eps()

    integer:: N, i
    real(8):: x_start, x_end, x_mid, c, alpha, v_1, v_N, dx, pi
    real(8), allocatable:: x(:), rho(:), eps(:), d(:), v(:), eps_p(:) 
    
    N=1000

    allocate( x(N), rho(N), eps(N), d(N), v(N), eps_p(N)) 
    
    x_start= 0.0d0
    x_end = 1000.0d0
    x_mid = (x_end+x_start)/2
    alpha= 1.0d0
    pi = 3.14159265359
    !c=-1.0d0 / (4.0d0 * pi)
    c=-1.0
    
    do i=1, N
      x(i) = x_start + (i-1)*(x_end-x_start)/N
      rho(i) = (x(i)-x_mid) * exp(-alpha * (x(i) -x_mid )**2)
      !rho(i) =  (1.0d0/(2.0d0 * sqrt(alpha*pi))) * exp(-alpha * (x(i) -x_mid )**2)
      rho(i) = rho(i) + 3.0d0 * (x(i)-x_mid)**3 * exp(-alpha * (x(i) -x_mid )**2)
      rho(i) = rho(i) +  0.2d0 * exp(-alpha * (x(i) -x_mid +10)**2)
      rho(i) = rho(i) -  0.2d0 * exp(-alpha * (x(i) -x_mid -10)**2)
    end do

    ! normalisera rho
    dx = x(2)-x(1)
    !rho(N/2-2) = 0.1d0 / dx 
    !rho(N/2-1) = -0.5d0 / dx 
    !rho(N/2) = 0.5d0 / dx 
    !rho(N/2+1) = -0.5d0 / dx
    !rho(N/2+2) = 0.5d0 / dx
    !rho(N/2+3) = -0.1d0 / dx
    !rho = rho * dx / sum(rho) 
        
    eps= 1.0d0
    !eps(1:N/3) = 3.0d0
    
    v_1 = -0.00d0
    v_N = 0.00d0
    d = -0.0001d0
    
    !call poisson_1d_fd_eps(x, rho, d, c, .true., eps, "dipole_corr", v_1, v_N, v)
    call poisson_1d_fd_eps(x, rho, d, c, .true., eps, "explicit", v_1, v_N, v)

   write(6,*) "integral rho", sum(rho)* dx
   write(6,*) "integral v", sum(v)* dx
    
    ! create derivative of epsilon
    eps_p(1) = (eps(2)-eps(1)) / dx
    eps_p(N) = (eps(N)-eps(N-1)) / dx
    do i=2, N-1
      eps_p(i) = (eps(i+1) -eps(i-1)) / (2.0d0 * dx)
    end do

    
    open(10, file="test_poisson_1d_fd_eps.txt")
    do i=1, N
      write(10,*) x(i), rho(i), v(i), eps(i), eps_p(i)
    end do
    close(10)

  end subroutine test_poisson_1d_fd_eps

  
end module m_poisson_1d_finite_diff
