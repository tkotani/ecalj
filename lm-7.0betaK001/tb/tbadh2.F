      subroutine tbadh2(nl,nbas,pos,plat,nspu,nsites,ipc,iax,npr,
     .                  dh,vm,vu,h0,hrs,ors)
C- Add mixed increments to H to real-space H
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nl,nbas,nspu,nsites,iax,ipc,npr
Ci   dh, dho: increments to on-site and off-site hamiltonian
Ci   h0 : input hamiltonian, H_in
Ci   ors : real space overlap matrix
Co Outputs:
Co   hrs : real space hamiltonian: H_0 + H'
Cr Remarks
Cr   This routine for density mixing branch. The density matrix, rhoc
Cr   or rhon, is mixed after bndtb and before tbesel. The components
Cr   of H' returned from tbesel are added here to H_in
C ----------------------------------------------------------------------
C     implicit none
      integer niax
      parameter (niax=10)
C Passed Parameters
      integer nl,nbas,nspu,nsites,npr(0:1,nbas),
     .        ipc(1),iax(niax,nsites)
      double precision dh(nl**2,nl**2,nbas,nspu),
     .                 vm(25,nbas),vu(0:nl-1,nbas),pos(3,1),plat(3,3),
     .                 h0(nl**2,nl**2,nsites,nspu),
     .                 hrs(nl**2,nl**2,nsites,nspu),
     .                 ors(nl**2,nl**2,nsites)

C Local Variables
      integer ispu,ib,j,ilm,ilmp,iprint

      if (iprint() .gt. 40) print 100
      do  ib = 1, nbas
        do ispu = 1, nspu
          if (iprint() .gt. 40 .and. nspu .eq. 2) write (*,150) ispu
          j = npr(1,ib) + 1
          do  ilm = 1, nl**2
            do  ilmp = 1, nl**2

C --- verbose output ---
              if ((ilm .eq. ilmp .and. iprint() .gt. 40)
     .          .or. iprint() .gt. 50) then
                write (*,200)
     .          ib,ilm,ilmp,hrs(ilm,ilmp,j,ispu),
     .          h0(ilm,ilmp,j,ispu)+dh(ilm,ilmp,ib,ispu)
              endif
C ----------------------

              hrs(ilm,ilmp,j,ispu) = h0(ilm,ilmp,j,ispu)
     .                             + dh(ilm,ilmp,ib,ispu)
            enddo
          enddo
        enddo
      enddo
      
  100 format(' TBADH2: site ilm ilm''   H_old   H_new')
  150 format ('        spin ',i1)
  200 format (7x,3i4,2x,2f10.6)


      end
