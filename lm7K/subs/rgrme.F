      subroutine rgrme(nsp,lmxs,lmx,nr,ri,wi,phi,dphi,gradrm)

C- Matrix elements of a gradients of a radial wave function.
C ----------------------------------------------------------------
Ci Inputs
Ci   nsp       2 for spin polarized
Ci   phi,dphi  radial wave functions and energy derivatives
Ci   lmxs      fourth dimension of gradrm is 0:lmxs
Ci   lmx       make integrals from 0 to lmx
Ci   ri,wi,nr  radial points, weights, and number
Co Outputs
Co   gradrm(i,j,ll,l,nsp):  <g2_l grad_r g1_l'> with l' = l +/- 1
Co         i = 1 for g1=phi,    g2=phidot
Co             2 for g1=phidot, g2=phi
Co             3 for g1=phi,    g2=phidot
Co             4 for g1=phidot, g2=phidot
Co         j = 1 for <g2 | grad_r g1> - (l+1) < g2 | 1/r g1 >
Co             2 for <g2 | grad_r g1> +    l  < g2 | 1/r g1 >
Co        ll = 1 for l' = l-1  for l' = l+1
Cr Remarks
Cr   Adapted from an optics package by V. Antropov
Cu Updates
Cu   05 Jan 04 dimensions of gradrm distinct from lmx
C ----------------------------------------------------------------
C     implicit none
      integer nr,nsp,lmxs,lmx
      double precision ri(nr),wi(nr)
      double precision phi(nr,0:lmx,nsp),dphi(nr,0:lmx,nsp)
C Outputs:
      double precision gradrm(4,2,2,0:lmxs,nsp)
C Local variables
      integer nrmx,lun,isp,loc,ir,ll
      parameter (nrmx=1501)
      double precision rq1,rq2,wir(nrmx),wbr(nrmx)
      
C     call prrmsh('phi',ri,phi,nr,nr,lmx+1)

C ... Extra weights for radial mesh
      wir(1) = 0
      wbr(1) = 0
      do  5  ir = 2, nr
      wir(ir) = wi(ir)*ri(ir)
    5 wbr(ir) = wi(ir)/ri(ir)
      call dpzero(gradrm,16*(lmxs+1)*nsp)
      
C --- Loop over spins ---
      do  10  isp = 1, nsp

C --- Matrix elements <g1 grad g2> using Wigner-Eckhart theorem ---
      do  20  lun = 0, lmx

C   ... Save radial integrals g1 grad g2 for each l
        do  22  loc = 0, lmx
          if (abs(loc-lun) .ne. 1) goto 22
          ll = 1
          if (loc-lun .eq. 1) ll = 2

          call prgrme(phi(1,loc,isp),phi(1,lun,isp),ri,wi,wir,wbr,
     .     nr,loc,gradrm(1,1,ll,loc,isp),gradrm(1,2,ll,loc,isp),rq1,rq2)
          
          call prgrme(dphi(1,loc,isp),phi(1,lun,isp),ri,wi,wir,wbr,
     .     nr,loc,gradrm(2,1,ll,loc,isp),gradrm(2,2,ll,loc,isp),rq1,rq2)
        
          call prgrme(phi(1,loc,isp),dphi(1,lun,isp),ri,wi,wir,wbr,
     .     nr,loc,gradrm(3,1,ll,loc,isp),gradrm(3,2,ll,loc,isp),rq1,rq2)
        
          call prgrme(dphi(1,loc,isp),dphi(1,lun,isp),ri,wi,wir,wbr,
     .     nr,loc,gradrm(4,1,ll,loc,isp),gradrm(4,2,ll,loc,isp),rq1,rq2)
   22   continue

   20 continue
C     call gradme(isp,nsp,1,1,lmx,nl*nl,gradrm,aamt,abmt,bamt,bbmt)
   10 continue

C      write(99,1004) aamt
C      write(99,1004) abmt
C      write(99,1004) bamt
C      write(99,1004) bbmt
C 1004 format(1p5d17.9)

Cgetarg       end
       end subroutine rgrme 


      subroutine prgrme(g1,g2,ri,wi,wir,wbr,nr,lc,r1,r2,rq1,rq2)

C- Make radial integrals of gradients
C ----------------------------------------------------------------
Ci Inputs
Co    r1  = <g2 | grad g1>   - (lc+1) < g2 | 1/r g1 >
Co    r2  = <g2 | grad g1>   +     lc < g2 | 1/r g1 >
Co    rq1 = <g2 | r*grad g1> - (lc+1) < g2 | g1 >
Co    rq2 = <g2 | r*grad g1> +     lc < g2 | g1 >
C ----------------------------------------------------------------
C     implicit none
      integer nr,lc
      double precision g1(nr),g2(nr),ri(nr),wi(nr),wir(nr),wbr(nr)
      integer nrmx,lerr
      parameter (nrmx=1501)
      double precision dg(nrmx)
      double precision r1,r2,rq1,rq2
      double precision sum1,sum2,dot3

C     call radgra(a,b,nr,ri,g1,dg)
      call poldvm(ri,g1,nr,6,.false.,1d-12,lerr,dg)
C     call prrmsh('dg2',ri,dg,nr,nr,1)
      sum1 = dot3(nr,g2,dg,wi)
      sum2 = dot3(nr,g2,g1,wbr)
      r1 = sum1 - sum2*(lc+1)
      r2 = sum1 + sum2*lc
      sum1 = dot3(nr,g2,dg,wir)
      sum2 = dot3(nr,g2,g1,wi)
      rq1 = sum1 - sum2*(lc+1)
      rq2 = sum1 + sum2*lc

Cgetarg       end
       end subroutine prgrme 


