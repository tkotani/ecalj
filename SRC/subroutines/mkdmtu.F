      subroutine mkdmtu(isp,iq,qp,nev,evec,dmatu)
      use m_lmfinit,only: ssite=>v_ssite,sspec=>v_sspec,nbas,nlmax,nsp,nspc,nl,n0,nppn,nlibu,lmaxu,nlibu,lldau
      use m_mkpot,only: ppnl=>ppnl_rv 
      use m_subzi, only: wtkb=>rv_a_owtkb
      use m_igv2x,only: ndimh
      use m_suham,only:  ndham=>ham_ndham,ndhamx=>ham_ndhamx 
C- Calculate density matrix for LDA+U channels
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored: *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idu
Ci     Stored: *
Ci     Passed to: *
Ci   wtkb  :eigenvalue weights for BZ integration of occupied states
Ci   isp   :current spin channel (1 or 2)
Ci   iq    :qp index, used only to address element in wtkb
Ci         :NB: aus is stored only for current qp
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   ndham :dimensions wtkb,aus
Ci   nlmax :1st dimension of aus (maximum nlma over all sites)
Ci   nbas  :size of basis
Ci   nev   :actual number of eigenvectors generated
Ci   ppnl  :nmto-like pot pars
Ci   aus   :coefficients to phi and phidot made previously by makusqldau
Ci  lldau  :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat beginning at dmats(*,lldau(ib))
Co Outputs
Co   dmatu :density matrix for specified LDA+U channels
Cb Bugs
Cb   Never checked for noncollinear case
Cu Updates
Cu   09 Nov 05 Convert dmat to complex form
Cu   28 Jun 05 bug fix for nspc=2
Cu   09 Jun 05 (MvS) extended to local orbitals
Cu   30 Apr 05 (WRL) first created
C ----------------------------------------------------------------------
      implicit none
C     ... Passed parameters
      integer:: i_copy_size
      integer isp,iq,nev
c      real(8):: wtkb(ndham,nsp,iq)
c      double complex aus(nlmax,ndham*nspc,3,nsp,nbas)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
c      double precision ppnl(nppn,n0,nsp,nbas)
C ... Local parameters
      double complex add,au,as,az,ap1,ap2
      double precision dlphi,rmt,dlphip,phi,phip,dphi,dphip,r(2,2),det,
     .phz,dphz
      integer lmxa,ilm1,ilm2,l,iv,m1,m2,ib,is,idu(4),igetss,iblu,ispc,
     .     ksp
      complex(8),allocatable ::aus(:,:,:,:,:)
      double complex evec(ndimh,nsp,nev)
      real(8)::qp(3)
      allocate(aus(nlmax,ndham*nspc,3,nsp,nbas))! aus(2*nlmax*ndhamx*3*nsp*nbas))
      aus=0d0
      call makusq(0, nbas,0 , nev,  isp, 1, qp, evec, aus )

      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) .eq. 0) goto 10
        is = int(ssite(ib)%spec)

        lmxa=sspec(is)%lmxa
        i_copy_size=size(sspec(is)%idu)
        call icopy(i_copy_size,sspec(is)%idu,1,idu,1)
        rmt=sspec(is)%rmt

        do  l = 0, min(lmxa,3)

          if (idu(l+1) .ne. 0) then
            iblu = iblu+1

C           In noncollinear case, isp=1 always => need internal ispc=1..2
C           ksp is the current spin index in both cases:
C           ksp = isp  in the collinear case
C               = ispc in the noncollinear case
C           ispc=1 for independent spins, and spin index when nspc=2

            do  ispc = 1, nspc
              ksp = max(ispc,isp)

C             For rotation (u,s) to (phi,phidot)
              dlphi  = ppnl(3,l+1,ksp,ib)/rmt
              dlphip = ppnl(4,l+1,ksp,ib)/rmt
              phi    = ppnl(5,l+1,ksp,ib)
              phip   = ppnl(6,l+1,ksp,ib)
              dphi   = phi*dlphi/rmt
              dphip  = dlphip/rmt*phip
              det    = phi*dphip - dphi*phip
              r(1,1) = dphip/det
              r(1,2) = -dphi/det
              r(2,1) = -phip/det
              r(2,2) = phi/det
C             For projection from loc. orb. onto (u,s)
              phz    = ppnl(11,l+1,ksp,ib)
              dphz   = ppnl(12,l+1,ksp,ib)

              ilm1 = l*l
              do  m1 = -l, l
                ilm1 = ilm1+1
                ilm2 = l*l
                do  m2 = -l, l
                  ilm2 = ilm2+1
                  add = (0d0,0d0)
C                 If (au,as) are coefficients to (phi,phidot), use this
C                 do  iv = 1, nev
C                   ap1 = aus(ilm1,iv,1,ksp,ib)
C                   ap2 = aus(ilm2,iv,1,ksp,ib)
C                   add = add + ap1*dconjg(ap2)*wtkb(iv,ksp,iq)
C                 enddo

C                 If (au,as) are coefficients to (u,s), use this
C                 Local orbital contribution adds to u,s
C                 deltau = -phi(rmax) * az   deltas = -dphi(rmax) * az
                  do  iv = 1, nev
                    az = aus(ilm1,iv,3,ksp,ib)
                    au = aus(ilm1,iv,1,ksp,ib) - phz*az
                    as = aus(ilm1,iv,2,ksp,ib) - dphz*az
                    ap1 = au*r(1,1) + as*r(2,1)
C                   ad1 = au*r(1,2) + as*r(2,2)

                    az = aus(ilm2,iv,3,ksp,ib)
                    au = aus(ilm2,iv,1,ksp,ib) - phz*az
                    as = aus(ilm2,iv,2,ksp,ib) - dphz*az
                    ap2 = au*r(1,1) + as*r(2,1)
C                   ad2 = au*r(1,2) + as*r(2,2)
                    add = add + ap1*dconjg(ap2)*wtkb(iv,isp,iq)
                  enddo
                  dmatu(m1,m2,ksp,iblu) = dmatu(m1,m2,ksp,iblu) + add
C                  if (m1 .eq. m2 .and. l .eq. 3) call info5(0,0,0,
C     .              'l=%i m=%i iq=%i isp=%i add=%d',l,m1,iq,ksp,add)

                enddo
              enddo
            enddo
          endif
        enddo
   10   continue
      enddo
      deallocate(aus)
      end subroutine mkdmtu


