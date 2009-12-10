      subroutine tbmpol(nbas,nsp,nl,qnu,ipc,ltb,gaunt,qpol,rho,rhoc,
     .                  qmpol,mmom)
C- Make multipole (and magnetic) moments from TB "density matrix"
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas,nl,nsp,qnu,ipc,ltb,gaunt
Ci   qpol: polarisation parameters (delta_ll'l'')
Ci   rho, rhoc (see tbfrce); if mixrho then these are from rhomix.
Co Outputs:
Co   qmpol
Cr Remarks
Cr   Our multipoles are defined such that the monopole moment is the
Cr   total Mulliken charge, while the higher multipoles are approximated
Cr   as an onsite sum without reference to overlap (Finnis, p. 216)
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer ltb,nbas,nl,nsp,ipc(1)
      double precision qpol(10,1),qnu(3,0:nl-1,nsp,1),gaunt(9,9,25),
     .                 rho(nl,2,nbas),rhoc(nl**2,nl**2,nbas),
     .                 qmpol(9,nbas),mmom(nbas)
C Local Variables
      integer ib,ic,ilm,ilmp,ilmpp,iprint
      integer ll
      double precision q0,M,dsum

      call dcopy(9*nbas,0d0,0,qmpol,1)
      call dcopy(nbas,0d0,0,mmom,1)

C --- get multipole moments ---
      do  ib = 1, nbas
        ic = ipc(ib)
        if (nsp .eq. 1) then
          q0 = dsum(nl,qnu(1,0,1,ic),3)
          qmpol(1,ib) = dsum(nl,rho(1,1,ib),1) - q0
        else
          q0 = dsum(nl,qnu(1,0,1,ic),3) + dsum(nl,qnu(1,0,2,ic),3)
          qmpol(1,ib) = dsum(nl,rho(1,1,ib),1)
     .                + dsum(nl,rho(1,2,ib),1) - q0
        endif
        if (iprint() .gt. 60) write (*,100) 
        do  ilm = 2, 9
          do  ilmp  = 1, nl**2
            do  ilmpp = 1, nl**2
              if (ilm .gt. 1) then
                call getM(ilm,ilmp,ilmpp,qpol(1,ic),M)
                qmpol(ilm,ib) = qmpol(ilm,ib) +
     .               rhoc(ilmp,ilmpp,ib) * M * gaunt(ilmp,ilmpp,ilm)
              endif
                
C --- verbose output ---
              if (iprint() .gt. 60 .and. M .ne. 0d0
     .            .and. gaunt(ilmp,ilmpp,ilm) .ne. 0d0) then
                write (*,200) 
     .              ilmp,ilmpp,ilm,ll(ilmp),ll(ilmpp),ll(ilm),
     .              M,gaunt(ilmp,ilmpp,ilm),rhoc(ilmp,ilmpp,ib)
              endif
C ---------------------

            enddo
          enddo
        enddo
C        if (iprint() .gt. 30) then
C          call awrit3('        Q^e/e=%d, %3:1d, %5:1d',' ',180,
C     .                i1mach(2),qmpol(1,ib),qmpol(2,ib),qmpol(5,ib))
C        endif

C --- make magnetic moments ---
        if (nsp .eq. 2) then
          mmom(ib) = dsum(nl,rho(1,1,ib),1) - dsum(nl,rho(1,2,ib),1)
        endif
      enddo
  100 format ('  L''   L''''  L    l''   l''''  l      M         CG
     .      rho_L''L''''')
  200 format (6(i3,2x),2(2x,f6.2,2x),2f10.6)
      end
