      subroutine mkrho0(lr,nl,nbas,nclass,ipc,dclabl,nsp,qnu,rho0,rhol0,
     .                  rholm0,mmom)
C- Construct rho_in
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lov,nl,nbas,nclass,ipc,npsu,qnu
Ci   lr: if 1 make also full on site density matrix, otherwise just
Ci       Mulliken charges
Co Outputs:
Co   rho0: density matrix, rhol0,rholm0: Mulliken charges
Co   if nsp .eq. 2, make also total moment at each site.
Cr Remarks
Cr   Makes the "input" density matrix from the Q's. Also if lov=T makes
Cr   input Mulliken charges.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer lr,nl,nbas,nsp,nclass,ipc(nbas)
      double precision qnu(3,0:nl-1,nsp,nclass),mmom(nbas),
     .  rho0(nl**2,nl**2,nbas,nsp),rhol0(0:nl-1,2,nbas),
     .  rholm0(nl**2,2,nbas),dclabl(1)
C Local Variables
      integer ilm,l,m,ib,isp,ic,iprint,i1mach
      double precision q,dsum
      character*8 clabl

      if (lr .eq. 1) call dcopy(nl**4*nbas*nsp,0d0,0,rho0,1)

      do  ib = 1, nbas
        do  isp = 1, nsp
          ilm = 0
          do  l = 0, nl-1
            ic = ipc(ib)
            q = qnu(1,l,isp,ic) / (2*l + 1)
            rhol0(l,isp,ib) = qnu(1,l,isp,ic)
            do  m = -l, l
              ilm = ilm + 1
              if (lr .eq. 1) then
                rho0(ilm,ilm,ib,isp) = q
              endif
              rholm0(ilm,isp,ib) = q
            enddo
          enddo
        enddo
        if (nsp .eq. 2) then
          mmom(ib) = dsum(nl,qnu(1,0,1,ic),3) - dsum(nl,qnu(1,0,2,ic),3)
        endif
      enddo
      if (iprint() .lt. 30) then
        return
      endif
C --- Printout ---
      call awrit0('%N MKRHO0: input Mulliken charges ..',' ',120,
     .            i1mach(2))
      do  ib = 1, nbas
        ic = ipc(ib)
        call r8tos8(dclabl(ic),clabl)
        call awrit1('Atom %i '//clabl,' ',60,i1mach(2),ib)
        if (nsp .eq. 1) then
           call awrit4('   by l: %n:2d%N   by L: %n:2d',' ',120,
     .                 i1mach(2),nl,rhol0(0,1,ib),nl**2,rholm0(1,1,ib))
        else
          do  isp = 1, nsp
            call awrit1('    spin %i',' ',120,i1mach(2),isp)
            call awrit4('   by l: %n:2d%N   by L: %n:2d',' ',120,
     .             i1mach(2),nl,rhol0(0,isp,ib),nl**2,rholm0(1,isp,ib))
          enddo
         call awrit1('   moment:  %d%N',' ',120,i1mach(2),mmom(ib))
        endif
      enddo
      end
