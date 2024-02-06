subroutine hello() bind(C)
  write(6,*)'hello world'
endsubroutine 
subroutine hello2(iii) bind(C)
  integer::iii
  write(6,*)'hello world2! iii=',iii
endsubroutine
module m_test2
  contains
 subroutine hello3(lll,intx,r1,r2) bind(C)
   use m_comm,only: comm,size,rank
   implicit none
   include "mpif.h"
   integer :: ierr,intx
   logical:: lll
   real(8):: r1,r2
   call MPI_barrier(comm,ierr)
   if(mod(rank,2)==0) print *,'rank/size=',rank,size, "Hello World333 even", lll,intx,'sum=',r1+r2
   if(mod(rank,2)==1) print *,'rank/size=',rank,size, "Hello World333 odd ", lll,intx,'sum=',r1+r2
 end subroutine hello3
end module m_test2
subroutine hello4(cname,nsize) bind(C)
   implicit none
   integer:: nsize,i
   character(1):: cname(nsize)
   character(512):: string=''
   forall(i=1:nsize) string(i:i)=cname(i)
   write(6,*)' ',string(1:nsize),' ##########'
end subroutine hello4
