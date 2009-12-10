      subroutine mkiaxd(npr,lmx,ips,iax)
C- Makes basis-related parts of iax, and re
C     implicit none
      integer npr,niax,lmx(*),ips(*)
      parameter (niax=10)
      integer iax(niax,npr)
      integer nli,is,i,jbas
c     ibas0 = 0
      do  10  i = 1, npr
c        ibas = iax(1,i)        
        jbas = iax(2,i)
        is = ips(jbas)
        nli = (lmx(is)+1)**2
        iax(9,i) = nli
C        if (ibas .ne. ibas0) then
C          ibas0 = ibas
CC         iax(10,i) = nli
C        else
CC         iax(10,i) = iax(10,i-1) + nli 
C        endif
   10 continue
C      print *, isum(50,iax(9,1),niax)
C      print *, isum(npr,iax(9,1),niax)
      end

