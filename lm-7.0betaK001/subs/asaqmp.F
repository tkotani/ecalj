      subroutine asaqmp(iopt,sctrl,spot,slat,sarray,lmxa,lmxf,nlmf,qmp)
C- Multipole moments from 2nd gen ASA wave function products
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt  :not used: must be zero
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nl nspin
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opmpol oqpp
Ci      qpp(*,1,ib) phi-phi   products for (L<=L')
Ci      qpp(*,2,ib) phi-dphi  products for all (L,L')
Ci      qpp(*,4,ib) dphi-dphi products for (L<=L')
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat qlat nsgrp oistab oag osymgr ocg ojcg 
Ci                oidxcg ocy
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp oclabl oipc opos oics
Ci   lmxa  :array of augmentation l-cutoffs
Ci   lmxf  :array of l-cutoffs for qpp
Ci   nlmf  :leading dimension and global L-cutoff for qmp
Co Outputs
Ci   qmp   :multipole moments
Cr Remarks
Cr   For nonspherical multipole moments, let
Cr        I^(pp)_l'l''m = int (phi_l' phi_l'' r^m)
Cr        I^(pd)_l'l''m = int (phi_l' phidot_l'' r^m)
Cr        I^(dd)_l'l''m = int (phidot_l' phidot_l'' r^m)
Cr   Then the multipole moments inside sphere R are
Cr     q_M = sum_L',L'' CG_ML'L''
Cr           [ (1+oh)z+_RL' I^(pp)_l'l''m (1+oh)z_RL''   +
Cr             (1+oh)z+_RL' I^(pd)_l'l''m (hz)_RL'' + h.c. +
Cr                (hz)+_RL' I^(dd)_l'l''m (hz)_RL'' ]
Cr         = sum_L',L'' CG_ML'L''
Cr           [ q_pp(L',L'') I^(pp)_l'l''m +
Cr           2*q_pd(L',L'') I^(pd)_l'l''m
Cr             q_dd(L',L'') I^(dd)_l'l''m ]
Cu Updates
Cu   08 Nov 07 (J. Xu) qpp is complex
Cu   23 Aug 01 adapted from makqmp.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer iopt,nlmf,lmxa(1),lmxf(1)
      double precision sctrl(1),spot(1),slat(1),sarray(1),qmp(*)
C ... Local parameters
      double precision alat,plat(9),qlat(9)
      integer ngrp,nrclas,nclasp,ic,nlml,nl,nn,ipr,iprint,
     .  nlmx,lgunit,i,nbas,nqpp,nsp,oistab,oag,og,oipa,oipc,oics,
     .  opos,oposc,oqwk,opmpol,oqpp,osym,ocy,ocg,ojcg,
     .  oidxcg,oclabl
      double precision dval
      character*8 clabl
C ... Heap
      integer w(1)
      common /w/ w

C ... Setup
      call upack('array nclasp oclabl oipc opos oics',sarray,nclasp,
     .  oclabl,oipc,opos,oics)
      call upack2('pot opmpol oqpp',spot,opmpol,oqpp)
      call upack('ctrl nbas nl nspin',sctrl,nbas,nl,nsp,0,0)
      call upack('lat alat plat qlat',slat,alat,plat,qlat,0,0)
      call upack('lat nsgrp oistab oag osymgr',slat,ngrp,oistab,oag,
     .  og,0)
      call upack('lat ocg ojcg oidxcg',slat,ocg,ojcg,oidxcg,0,0)
      call upack1('lat ocy',slat,ocy)
      ipr = iprint()
      if (iopt .ne. 0) call rx('asaqmp: bad  iopt')

C --- Multipole moments at all sites from the qpp ---
      i = nl**2
      nqpp = (i*(i+1))/2
      call qpp2mp(nqpp,nl,nsp,nbas,nlmf,w(oipc),lmxa,
     .  w(ojcg),w(oidxcg),w(ocg),w(oqpp),w(opmpol),qmp)
c     call psymqp(nlmf,nlmf,1,nbas,qmp)

      if (ipr .ge. 20)
     .  call awrit1('%N ASAQMP: Make and symmetrize multipole moments'
     .  //' for %i classes',' ',80,lgunit(1),nclasp)

      call defi (oipa,    nbas)
      call defrr(oposc, 3*nbas)
      call defrr(oqwk, nlmf)

C --- For each class, do ---
      do  10  ic = 1, nclasp
        nlml = (lmxf(ic)+1)**2
        nlmx = max0(nlml,4)

C   ... Make nrclas,ipa,posc
        call psymr0(-2,ic,nbas,w(oipc),w(opos),w(oposc),w(oipa),nrclas)
        do  i  = 1, nrclas
          w(oipa+i-1) = (w(oipa+i-1)-1)*nlmf
        enddo

C   ... Symmetrize qmp for members of this class
        call defrr (osym,   nlmx*nlmx*nrclas)
        call symqmp(nrclas,nlml,nlmx,plat,w(oposc),ngrp,w(og),w(oag),
     .    w(oqwk),w(oipa),w(osym),qmp,nn)
        call rlse (osym)

        call r8tos8(dval(w(oclabl),ic),clabl)
        if (ipr .ge. 30)
     .    call awrit3('%N Class '//clabl//'%a: %16p'//
     .    'nrc = %i,  %i nonspherical elements of %i ',' ',80,lgunit(1),
     .    nrclas,nn,nlml)

        if (ipr .ge. 30 .and. ipr .lt. 50)
     .    call psymqp(nlmf,nlmf,w(oipa),w(oipa),qmp)

   10 continue

      call rlse (oipa)

      if (ipr .ge. 50) call psymqp(nlmf,nlmf,1,nbas,qmp)

      end

      subroutine psymqp(nlmf,nlml,ib1,ib2,qmp)
C- Printout
C     implicit none
      integer nlmf,nlml,ib1,ib2
      double precision qmp(nlmf,ib2),fpi,y0
      integer j1,ib,ilm

      print 221
      fpi = 16d0*datan(1d0)
      y0 = 1d0/dsqrt(fpi)
      j1 = 1
      do  10  ib = ib1, ib2
        print 220, ib,1,qmp(j1,ib),qmp(j1,ib)/y0
        do  20  ilm = 2, nlml
   20   if (dabs(qmp(ilm,ib)).gt.1d-6) print 222, ilm,qmp(ilm,ib)
   10 continue

  220 format(i4,i6,f12.6,f12.6,2f9.2)
  222 format(4x,i6,f12.6)
  221 format('  ib   ilm      qmom',8x,'Qval')

      end


      subroutine qpp2mp(nqpp,nl,nsp,nbas,nlmf,ipc,lmxa,
     .  jcg,indxcg,cg,qpp,pmpol,qmp)
C- Make multipole moments from the qpp and pmpol
C     implicit none
      integer nqpp,nbas,nsp,nl,nlmf,lmxa(1),ipc(nbas),
     .  jcg(1),indxcg(1)
      double precision qmp(nlmf,nbas),
     .  pmpol(nl,nl,2*nl-1,3,nsp,1),cg(1)
      double complex qpp(nqpp,4,nsp,nbas)
C Local
      integer ib,nlm,iqpp,iqpd,ilm1,l1p1,ll,ilm2,l2p1,ix,icg,mlm,lm,
     .  isp,ic
      double precision wij

      call dpzero(qmp,nlmf*nbas)

C --- For each site, make qmp ---
      do  10  isp = 1, nsp
      do  10  ib = 1, nbas
        ic = ipc(ib)
        nlm = (lmxa(ic)+1)**2

C   ... debugging
C        iqpp = 0
C        snot = 0
C        do  ilm1 = 1, nlm
C        do  ilm2 = 1, nlm
C          if (ilm2 .le. ilm1) then
C            wij = 2
C            if (ilm1 .eq. ilm2) wij = 1
C            iqpp = iqpp+1
C            print *, ilm1,ilm2,iqpp,qpp(iqpp,1,isp,ib)
C            snot(ilm1,ilm2) = qpp(iqpp,1,isp,ib)
C          endif
C        enddo
C        enddo
C        print *, 'ib,isp=',ib,isp
C        call prmx('qpp',snot,16,16,16)

        iqpp = 0
        iqpd = 0
        do  20  ilm1 = 1, nlm
        l1p1 = ll(ilm1)+1
        do  20  ilm2 = 1, nlm
        l2p1 = ll(ilm2)+1

          ix = max0(ilm1,ilm2)
          ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)

C     ... phi-phi and dot-dot terms
          if (ilm2 .le. ilm1) then
            wij = 2
            if (ilm1 .eq. ilm2) wij = 1
            iqpp = iqpp+1
            do  32  icg = indxcg(ix), indxcg(ix+1)-1
              mlm = jcg(icg)
              if (mlm .le. nlmf) then
                lm = ll(mlm)
                qmp(mlm,ib) = qmp(mlm,ib) + cg(icg)*wij*
     .      (dble(qpp(iqpp,1,isp,ib))*pmpol(l1p1,l2p1,lm+1,1,isp,ic) +
     .      dble(qpp(iqpp,4,isp,ib))*pmpol(l1p1,l2p1,lm+1,3,isp,ic))
              endif
   32       continue
          endif

C     ... phi-dot terms
          iqpd = iqpd+1
          do  34  icg = indxcg(ix), indxcg(ix+1)-1
            mlm = jcg(icg)
            if (mlm .le. nlmf) then
              lm = ll(mlm)+1
C              if (mlm .eq. 2 .and. ib .eq. 1) then
C                call awrit8(
C     .            'l1,m1= %i %i l2,m2= %i %i iqpd=%i %32p%g %g %g',
C     .            ' ',80,6,l1p1-1,ilm1,l2p1-1,ilm2,iqpd,
C     .            qpp(iqpd,2,isp,ib),cg(icg),
C     .            pmpol(l1p1,l2p1,lm,2,isp,ic))
C              endif
              qmp(mlm,ib) = qmp(mlm,ib) + cg(icg)*2*
     .    (dble(qpp(iqpd,2,isp,ib))*pmpol(l1p1,l2p1,lm,2,isp,ic))

            endif
   34     continue
   20   continue

   10 continue

      end

