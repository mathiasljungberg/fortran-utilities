program spline_energy
  use parameters
  implicit none

  ! input/output
  character(80)::infile,outfile
  integer::nin, nout
  real(kind=wp):: first,last
  real(kind=wp),dimension(:),allocatable::x_in,V_in
  ! loop variables
  integer::i,j

  ! other variables
  !real(kind=wp):: 
  real(kind=wp),dimension(:),allocatable::x_out,V_out

  read(5,*) infile
  read(5,*) outfile
  read(5,*) nin
  read(5,*) nout ,first , last
 
  allocate(x_in(nin),V_in(nin))
  allocate(x_out(nout),V_out(nout))

  open(10,file=infile,status='unknown')
  do i=1,nin
     read(10,*) x_in(i), V_in(i)
  end do
  close(10)

  if(first.lt.x_in(1)) then
     write(6,*) "point lower than supplied potential range!"
     stop
  end if
  if(last.gt.x_in(nin)) then
     write(6,*) "point higher than supplied potential range!"
     stop
  end if

  call linspace(x_out,first,last,nout)

  call spline_easy(x_in, V_in, nin, x_out, V_out, nout)

  open(10,file=outfile,status='unknown')
  do i=1,nout
     write(10,*) x_out(i), V_out(i)
  end do
  close(10)

end program spline_energy

