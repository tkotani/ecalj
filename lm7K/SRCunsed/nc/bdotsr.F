      subroutine bdotsr(nbas,nl,ipc,rhos,nrhos,qnu,bdots,
     .  lihdim,indxsh,mode,bsigr)
C- Double-counting term <B.sigma.rho>
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   rhos  :spin density matrix (for mode=0), by class
Ci         :rhos should be hermitian in spin space, but may not be owing
Ci         :to energy integration errors in the complex plane.
Ci         :bdotsr uses a symmetrized form to minimize errors.
Ci   qnu   :moments (for mode=1)
Ci   bdots :(magnetic field . Pauli matrix), local coordinate system,
Ci         :downfolding order
Ci   mode  :0, use spin density matrix to make moments along T
Ci         :1, use qnus to make moments along qnu
Co Outputs
Cr Remarks
Cr   Definition of rho in terms of M: (standard definition of sigma)
Cr      rho = M . sigma/2
Cr   Pauli matrices sigma: (NB: normal def is twice this)
Cr
Cr              (0  1)             (0 -i)           (1  0)
Cr     sigmax = (    )    sigmay = (    )  sigmaz = (    )
Cr              (1  0)             (i  0)           (0 -1)
Cr
Cr   Given rho, M can be obtain from:
Cr     M_x =  2 Re(rho21) = Re (rho12+rho21)
Cr     M_y =  2 Im(rho21) = Im (rho21-rho12)
Cr     M_z =  (rho11)-(rho22)
Cr
Cr   Second (symmetrized) form is used because for numerical reasons,
Cr   rhos may not be quite hermitian; e.g. when rhos is generated
Cr   by a Green's function technique.
Cr
Cr   Double counting term is
Cr     Tr <(B.sigma)(rho)>
Cr
Cr   Input B.sigma is
Cr                1 (Bz   Bx-iBy)
Cr    B.sigma =   - (           )
Cr                2 (Bx+iBy  -Bz)
Cr
Cr   Then (CHECK wrong factors of 2 in both b and sigma)
Cr      Tr <(B.sigma)(rho)>
Cr      = 1/2 (Bz(rho11-rho22) + (Bx-iBy) rho21 + (Bx+iBy) rho12)
Cr      = Bz(rho11-rho22)/2 + (Bx-iBy)(Mx+iMy)/4 + (Bx+iBy)(Mx-iMy)/4
Cr      = Bz(rho11-rho22)/2 + Bx Mx / 2 + By My /2
Cr      = 1/2 (B . M)
Cr   This formula can be computed either with the moments qnu
Cr   or from the spin-density matrix.
Cr Bugs
Cr   bdots is orbital-resolved; rho is not.
Cu Updates
Cu   21 Apr 04 Revised to accomodate an m-dependent density-matrix
Cu   15 Feb 03 Created from amagnc.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer nbas,nl,nrhos,ipc(*),mode,lihdim,indxsh(lihdim)
      double precision rhos(2,0:2,nrhos,2,2,1),
     .  qnu(3,nl,2,1),bdots(2,2,2,*),bsigr
C Local variables
      logical lrhol
      integer i,ib,ic,lgunit,ipr,k,stdo,lmr,l,m,ilmr,ilm
      double precision aml(3),samg(3),bm(3),bloc(3),facB
C     logical pass1,lprt
C     facB covers convention for B:
C     facB= 1 => +B induces positive -M
C     facB=-1 => +B induces positive +M
      parameter (facB=-1d0)

C --- Setup ---
      call getpr(ipr)
      stdo = lgunit(1)
      call isanrg(mode,0,1,' bdotsr','mode',.true.)
      lrhol = nrhos .eq. nl .or. mode .eq. 1
      if (ipr .gt. 50) then
        if (mode.eq.0) write(stdo,332) 'density matrix'
        if (mode.eq.1) write(stdo,332) 'sphere charges'
        write(stdo,335)
      endif
  332 format(/' BDOTSR: magnetic moments from ',a,':',:,
     .  ' <B.sig.rho>=',f12.7)
  335 format(3x,'ib  l     rhox*bx     rhoy*by     rhoz*bz      Sum')

C --- Accumulate double counting orbital by orbital ---
      bsigr = 0
      lmr = 0
      do  ib = 1, nbas
        ic = ipc(ib)

        call dpzero(samg,3)
        bm(1) = 0
        bm(2) = 0
        bm(3) = 0
        ilm = 0
        do  l = 0, nl-1
          k = l+1
          if (indxsh(lmr+1) .gt. lihdim) then
            lmr = lmr + 2*l+1
            ilm = ilm + 2*l+1
            goto 2
          endif
          do   m = -l, l
            ilm = ilm + 1
            lmr = lmr+1
            ilmr = indxsh(lmr)

            bloc(1) = 2*bdots(1,2,1,ilmr) *facB
            bloc(2) = 2*bdots(2,2,1,ilmr) *facB
            bloc(3) = 2*bdots(1,1,1,ilmr) *facB

C       ... local Mx,My,Mz from density-matrix or from qnu
            if (mode .eq. 0) then
              k = l+1
              if (.not. lrhol) k = ilm
              aml(1) = rhos(1,0,k,1,2,ic) + rhos(1,0,k,2,1,ic)
              aml(2) = rhos(2,0,k,2,1,ic) - rhos(2,0,k,1,2,ic)
              aml(3) = rhos(1,0,k,1,1,ic) - rhos(1,0,k,2,2,ic)
              if (lrhol) call dscal(3,1/dble(2*l+1),aml,1)
            else
              aml(1) = 0
              aml(2) = 0
              aml(3) = qnu(1,k,1,ic) - qnu(1,k,2,ic)
              call dscal(3,1/dble(2*l+1),aml,1)
            endif

            do  i = 1, 3
              bm(i) = bloc(i) * aml(i) / 2    * facB
            enddo

C       ... Add to double-counting term
            bsigr = bsigr + bm(1) + bm(2) + bm(3)

C       ... Printout
            if (ipr .gt. 50 .and. .not. lrhol)
     .      write(stdo,'(i5,i3,3f12.7,f14.7)') ib,l,bm,bm(1)+bm(2)+bm(3)
          enddo
          if (ipr .gt. 50 .and. lrhol)
     .      write(stdo,'(i5,i3,3f12.7,f14.7)') ib,l,bm,bm(1)+bm(2)+bm(3)

        enddo
    2   continue
      enddo

      if (ipr .gt. 30 .and. ipr .le. 50) then
        if (mode.eq.0) write(stdo,332) 'density matrix',bsigr
        if (mode.eq.1) write(stdo,332) 'sphere charges',bsigr
      elseif (ipr .ge. 30) then
        write(stdo,'('' double-counting <B.sigma.rho>='',f12.7)') bsigr
      endif

      end

