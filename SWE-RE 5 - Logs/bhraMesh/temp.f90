program temp
implicit none

  integer(kind=4)::i,j,k
  real(kind=8)::tr1,tr2,tr3

  open(11,file='Tecplot_bnd_points.dat')

  do i=1,566
    read(11,*)tr1,tr2,tr3
    write(12,'(3F20.10)')tr1,tr2,tr3
  enddo



end program temp