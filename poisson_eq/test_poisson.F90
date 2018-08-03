program test_poisson
  use m_poisson_1d_finite_diff, only: test_poisson_1d_fd_eps
  implicit none

  call test_poisson_1d_fd_eps()
  
end program test_poisson
