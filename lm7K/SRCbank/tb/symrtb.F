      subroutine symrtb(nl,nsp,nbas,ng,g,gd,istab,rwk,rhol)
C- Symmetrise rhol
C ----------------------------------------------------------------------
Ci Inputs:
Ci   rhol (unsymmetrised, from mkrhol in tbfrce)
Co Outputs:
Co   rhol (symmetrised)
Cr Remarks
Cr   Symmetrise rhol, in the usual way.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nl,nsp,nbas,ng,istab(nbas,ng)
      double precision g(3,3,ng),gd(5,5,ng),
     .                  rwk(0:nl**2-1,0:nl**2-1,nsp,1),
     .                 rhol(0:nl**2-1,0:nl**2-1,nsp,1)
C Local Variables
      integer ic,jc,ig,m1,m2,m3,m4,ib,i,j,i1mach,iprint

      call tcn('symrtb')

      if (iprint() .gt. 60) then
        print *, 'SYMRTB: Before symmetrisation Rho_LL''='
        do  ib = 1, nbas
          call awrit1('Atom %i',' ',60,i1mach(2),ib)
          do  i = 0, nl**2-1
            write (*,300) (rhol(i,j,1,ib),j=0,nl**2-1)
          enddo
        enddo
      endif

      call dpcopy(rhol,rwk,1,nl**4*nsp*nbas,1d0/ng)
      call dpzero(rhol,nl**4*nsp*nbas)

      do  100  ig = 1, ng
        do  90  ic = 1, nbas
          jc = istab(ic,ig)
C ... s block
          rhol(0,0,1,ic) = rhol(0,0,1,ic) + rwk(0,0,1,jc) 
C ... p block
          if (nl .gt. 1) then
            do  4  m1 = 1, 3
              do  3  m2 = 1, 3
                do  2  m3 = 1, 3
                  do  1  m4 = 1, 3
                    rhol(m1,m2,1,ic) = rhol(m1,m2,1,ic) +
     .                         g(m1,m3,ig)*g(m2,m4,ig)*rwk(m3,m4,1,jc)
    1             continue
    2           continue
    3         continue
    4       continue
          endif
C ... d block
          if (nl .gt. 2) then
            do  8  m1 = 4, 8
              do  7  m2 = 4, 8
                do  6  m3 = 4, 8
                  do  5  m4 = 4, 8
                    rhol(m1,m2,1,ic) = rhol(m1,m2,1,ic) +
     .              gd(m1-3,m3-3,ig)*gd(m2-3,m4-3,ig)*rwk(m3,m4,1,jc)
    5             continue
    6           continue
    7         continue
    8       continue
          endif
C ... sp block
          if (nl .gt. 1) then
            do  10  m1 = 1, 3
              do  9  m3 = 1, 3
                rhol(0,m1,1,ic) = rhol(0,m1,1,ic) +
     .                            g(m1,m3,ig)*rwk(0,m3,1,jc)
                rhol(m1,0,1,ic) = rhol(m1,0,1,ic) +
     .                            g(m1,m3,ig)*rwk(m3,0,1,jc)
    9         continue
   10       continue 
          endif
C ... sd block
          if (nl .gt. 2) then
            do  12  m1 = 4, 8
              do  11   m3 = 4, 8
                rhol(0,m1,1,ic) = rhol(0,m1,1,ic) +
     .                            gd(m1-3,m3-3,ig)*rwk(0,m3,1,jc)
                rhol(m1,0,1,ic) = rhol(m1,0,1,ic) +
     .                            gd(m1-3,m3-3,ig)*rwk(m3,0,1,jc)
   11         continue
   12       continue
          endif
C ... pd block
          if (nl .gt. 2) then
            do  16  m1 = 1, 3
              do  15  m2 = 4, 8
                do  14  m3 = 1, 3
                  do  13  m4 = 4, 8
                    rhol(m1,m2,1,ic) = rhol(m1,m2,1,ic) +
     .                     g(m1,m3,ig)*gd(m2-3,m4-3,ig)*rwk(m3,m4,1,jc)
                    rhol(m2,m1,1,ic) = rhol(m2,m1,1,ic) +
     .                     gd(m2-3,m4-3,ig)*g(m1,m3,ig)*rwk(m4,m3,1,jc)
   13             continue
   14           continue
   15         continue 
   16       continue 
          endif
   90   continue
  100 continue

      if (iprint() .gt. 60) then
        print *, 'SYMRTB: After symmetrisation Rho_LL''='
        do  ib = 1, nbas
          call awrit1('Atom %i',' ',60,i1mach(2),ib)
          do  i = 0, nl**2-1
            write (*,300) (rhol(i,j,1,ib),j=0,nl**2-1)
          enddo
        enddo
      endif
  300 format (5x,9f10.6)
      call tcx('symrtb')

      end
