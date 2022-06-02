subroutine dpscop(afrom,ato,nel,n1,n2,fac)
  !- shift and copy.
  !i nel number of elements
  !i n1  offset in afrom
  !i n2  ofset in ato
  !     implicit none
  integer :: n1,n2,i,iadd,ntop,nel
  double precision :: afrom(1),ato(1),fac

  if (fac /= 1d0) goto 100
  call dcopy(nel,afrom(n1),1,ato(n2),1)
  return

  ! --- fac not unity ---
100 continue
  iadd = n2-n1
  ntop = n1+nel-1
  do    i = n1, ntop
     ato(i+iadd) = fac*afrom(i)
  enddo
end subroutine dpscop
