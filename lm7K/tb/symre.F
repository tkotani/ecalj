      subroutine symre(lov,trh,nl,nsp,nbas,ng,istab,rwk,rwklm,wk,rho,
     .                 rholm,e)
C- Symmetrise e,rho,rholm
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lov, trh
Ci   e, rho, rholm (unsymmetrised, from tbfrce)
Co Outputs:
Co   e, rho, rholm (symmetrised)
Cr Remarks
Cr   Symmetrise rho, in the usual way.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nl,nsp,nbas,ng,istab(nbas,ng)
      double precision rho(0:nl-1,2,nbas),rwk(0:nl-1,2,nbas),
     .                 rholm(nl**2,2,nbas),rwklm(nl**2,2,nbas),
     .                 e(nbas,nsp),wk(nbas,nsp)
      logical trh,lov
C Local Variables
      integer isp,ib,jb,ig,l,ilm

      call dpcopy(rho,rwk,1,nl*2*nbas,1d0/ng)
      call dpzero(rho,nl*2*nbas)
      if (trh) then
        call dpcopy(e,wk,1,nbas*nsp,1d0/ng)
        call dpzero(e,nbas*nsp)
      endif
      if (lov) then
        call dpcopy(rholm,rwklm,1,nl**2*2*nbas,1d0/ng)
        call dpzero(rholm,nl**2*nsp*nbas)
      endif

      do   isp = 1, nsp
        do   ig = 1, ng
          do   ib = 1, nbas
            jb = istab(ib,ig)
            if (trh) e(ib,isp) = e(ib,isp) + wk(jb,isp)
            do  l = 0, nl-1
              rho(l,isp,ib) = rho(l,isp,ib) + rwk(l,isp,jb)
            enddo
            if (lov) then
              do  ilm = 1, nl**2
                rholm(ilm,isp,ib) = rholm(ilm,isp,ib)
     .                            + rwklm(ilm,isp,jb)
              enddo
            endif
          enddo
        enddo
      enddo
      end
