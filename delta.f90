integer function delta(i,j)
  implicit none
  
  integer,intent(in):: i,j


  if(i.eq.j) then
     delta = 1
  else
     delta = 0
  end if


end function delta

