      subroutine pgbasp(nbas,npadl,npadr,bas,plat,platl,platr,basp)
C- Padds" basis to accommodate extra bulk layers on LHS and RHS.
C ----------------------------------------------------------------
Co basp: padded basis
C ----------------------------------------------------------------
C     implicit none
      integer nbas,npadl,npadr
      double precision plat(3,3),bas(3,nbas),basp(3,*),
     .  platl(3,3),platr(3,3)
      double precision wk(3),pgfn(3),xx,ddot
      integer ib,nbasp

C     call prmx('original basis',bas,3,3,nbas)

      call dpcopy(bas,basp,1,3*nbas,1d0)

C --- make pgfn ---
      call cross(plat,plat(1,2),pgfn)
      xx = ddot(3,pgfn,1,plat(1,3),1)
      if (xx .lt. 0) call dscal(3,-1d0,pgfn,1)

C --- Padding on LHS ---
      nbasp = nbas
      xx = ddot(3,pgfn,1,platl(1,3),1)
      call dpcopy(platl(1,3),wk,1,3,dsign(1d0,xx))
      do  30  ib = 1, npadl
        call dpcopy(bas(1,ib),basp(1,nbasp+1),1,3,1d0)
        call dpadd(basp(1,nbasp+1),wk,1,3,-1d0)
        nbasp = nbasp+1
   30 continue

C --- Padding on RHS ---
      xx = ddot(3,pgfn,1,platr(1,3),1)
      call dpcopy(platr(1,3),wk,1,3,dsign(1d0,xx))
      do  40  ib = nbas-npadr+1, nbas
        call dpcopy(bas(1,ib),basp(1,nbasp+1),1,3,1d0)
        call dpadd(basp(1,nbasp+1),wk,1,3,1d0)
        nbasp = nbasp+1
   40 continue

C      call prmx('padded basis',basp,3,3,nbasp)
C      call fopna('out',66,0)
C      call dfdump(basp,3*nbasp,-66)
C      stop
C      do  50  ib = 1, nbasp
C        print 333, ib, basp(1,ib), basp(2,ib), basp(3,ib)
C  333   format(i4,3f12.6)
C   50 continue

      end

