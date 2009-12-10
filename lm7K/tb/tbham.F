      subroutine tbham(nsp,nspc,sctrl,sarray,slat,
     .  fitpar,nvar,ip1,ivar,nlmesp,memode,
     .  decay,deccf,decov,dcocf,iam,npm,nterm,nset,
     .  tabme,tabcf,tabov,tbocf,nsites,npr,oiax,
     .  oh,oh0,odh,oov,odov,odhcf,odovcf)
C- Real-space Slater-Koster tight-binding Hamiltonian for all sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp: TB+U or TB-L spin polarised
Ci   nspc: coupled spins (JEK's S-O coupling)
Ci   nlmesp
Ci   fitpar: if true, get derivatives for fitting, not normal ones
Ci   nvar: total number of TB parameters to vary (only used if fitpar=T)
Ci   ip1: pointer to locations in full list of variables (fitpar=T)
Ci   ivar(1,i): points to position of ith variable (fitpar=T)
Ci   ivar(2,i): paramter type of ith variable (fitpar=T)
Ci   liax: use iax, nsites as passed rather than making them (order N)
Ci   memode,decay,deccf,decov,dcocf,iam,npm,nterm,tabme,tabcf,tabov,
Ci     tbocf: see rdtbh
Co Outputs
Co   nsites: accumulated number of neighbors in all clusters
Co   npr(0,i): number of neighbors within rmax for ith atom
Co   npr(1,i): offset in accumulated list of all neighbors to the
Co             ith cluster associated with the ith atom
Co   oiax (and contents iax): neighbor lists, for the jth atom pair:
Co     iax(1,j) = cental atom of cluster
Co     iax(2,j) = second atom of the pair
Co     connecting vector: dr(i) = plat(i,1)*iax(3,j)
Co                              + plat(i,2)*iax(4,j) + plat(i,3)*iax(5,j)
Co   oh,oh0,odh (and contents s,dh): Real-space Ham and its derivative
Co   oov,odov (and contents ov,dov): Real-space overlap and deriv.
Co   odhcf,odovcf (and contents dhcf,dovcf): Crystal field derivatives
Cr Remarks
Cr   In TB+U the spin up and down hamiltonians are identical to
Cr   start with. Hence make one spin only and copy into second spin in
Cr   TBDIAG. Spin polarisation then happens during self-consistency
Cr   and the terms are added in tbaddh along with the electrostatic
Cr   potential.
Cr   oh0 is allocated here to hold a copy of the real space H before
Cr   the diagonal terms are added. This is then H^in for use in TB-L
Cr   and TB+U. (See tbdiag)
Cu Updates
Cu    8 Jun 07 (MvS) Merged Klepeis's additions, esp bands fitting
Cu   23 May 07 (MvS) bug fix, GSP mode
Cu   28 Sep 98 (JEK) updated to handle fitting of TB parameters (TBFIT)
Cu             and arbitrary number of MEs (prelude to including f's)
Cu   1  Dec 05 (ATP) Added oh0
Cu   18 May 05 (ATP) Add the spin index nsp for TB+U
Cu   19 May 96 (MvS) require now that liax true, since it is passed
Cu             iax(6) data, and can handle padding cases.
Cu             NB: should replace npr with ntab also.
Cu   11 Sep 97 (MvS) adapted for data contained in structures.
Cu   11 Aug 98 (ATP) updated for new addtos (ov and crf needs testing)
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      double precision sctrl(1),sarray(1),slat(1)
      integer nvar,nlmesp,memode,nterm,nset,nsites
      integer oiax,oh,oh0,odh,oov,odov,odhcf,odovcf
      integer ip1(1),ivar(1),iam(3,1),npm(2,1),npr(0:1,1)
      double precision decay(nlmesp,*),deccf(1),decov(1),dcocf(1),
     .  tabme(nterm,nlmesp,*),tabcf(1),tabov(1),tbocf(1)
C     character*8 clabl
      logical fitpar
C heap:
      integer w(1)
      common /w/ w
C Local parameters
      integer nbas,nl,nspc,nsp,ltb
      integer iat,ndimL,nl2,nsts,nsited,ilme,nlme,pointr,
     .  i1mach,isp1,isp2,pntcf,pntov,pntocf,niax
      parameter (niax=10)
      integer oiaxc,ohj,odhxj,odhyj,odhzj,odhrj,owk1,owk2,owk3,owk4,
     .  odv,ohcf,oipc,opos
      double precision alat,plat(3,3)
      logical hdiffs,cryf,lov,ocryf,bittst
C     double precision avw,avwsr

      call upack('ctrl nbas nl',sctrl,nbas,nl,0,0,0)
      call upack('ctrl ltb',sctrl,ltb,0,0,0,0)
      call upack('lat alat plat',slat,alat,plat,0,0,0)
      call upack('array oipc opos',sarray,oipc,opos,0,0,0)

      nl2 = nl**2
      lov = bittst(ltb,1)
      cryf = bittst(ltb,2)
      ocryf = bittst(ltb,4)
C     hdiffs = bitand(ltb,16+128) .ne. 0
      hdiffs = bittst(ltb,16) .or. bittst(ltb,128)
      if ((ocryf .or. cryf) .and. nsp .eq. 2) then
        call rx('TBHAM: cannot have CRYF and TBU')
      endif
      if (nspc .eq. 2 .and. nsp .eq. 2) then
        call rx('TBHAM: cannot have S-O and TBU')
      endif

      call rxx(lov .and. memode .eq. 1,
     .  'TBHAM: No overlap matrix for universal Hamiltonian')
      call rxx((cryf .or. ocryf) .and. memode .eq. 1,
     .  'TBHAM: No crystal field terms for universal Hamiltonian')
      call rxx(fitpar .and. memode .eq. 1,
     .  'TBHAM: Universal Hamiltonian not allowed for fitting')

C --- Get neighbor table iax for each atom ---
C oiaxc is the offset to iax in the current cluster
C nsites is the accumulated number of neighbors in all clusters
C      if (.not. liax) then
C        call rx('this branch eliminated; missing iax(6)')
C        if (mxnbr0 .eq. 0) then
C          avw = avwsr(plat,1d0,avw,nbas)
C          mxnbr = 4*(rmax/avw)**3*nbas
C        else
C          mxnbr = mxnbr0*nbas
C        endif
C        nsites = 0
C        call defi(oiax,niax*mxnbr)
C        call defdr(owk,mxnbr)
C        if (iprint() .ge. 30) print *
C        do  20  iat = 1, nbas
C          oiaxc = oiax + niax*nsites
C          call nghbor(nbas,plat,w(opos),rmax,rmax,iat,
C     .      mxnbr-nsites,npr(0,iat),w(oiaxc),w(owk))
C          nsites = nsites + npr(0,iat)
C          if (.not. fitpar .and. iprint() .ge. 30 .or. iprint() .gt. 40)
C     .      write(*,10) iat, clabl(ipc(iat)), npr(0,iat)-1, rmax
C   10     format(' TBHAM:  Site',i4,' (',a,')',':',i4,
C     .      ' neighbours within ',f8.4,'a')
C   20   continue
CC ... Reallocate table for iax ...
C        call rlse(oiax)
C        call defi(oiax,niax*nsites)
C      endif

C --- Allocate memory for H, H_0, dH, O, dO ---
      call defdr(oh,-nl2**2*nsites*min(nspc**2*nsp,4))
      if (bittst(ltb,2**15) .or. bittst(ltb,2**13)) then
        call defdr(oh0,-nl2**2*nsites*min(nspc**2*nsp,4))
      else
        oh0 = oh
      endif
      if (lov) call defdr(oov,-nl2**2*nsites*nspc**2)
      if (fitpar) then
        call defdr(odh,-nvar*nl2**2*nsites*nspc**2)
      elseif (hdiffs) then
        call defdr(odh,-4*nl2**2*nsites*nspc**2)
        if (cryf) call defdr(odhcf,-4*nl2**2*nsites*nspc**2)
        if (lov) then
          call defdr(odov,-4*nl2**2*nsites*nspc**2)
          if (ocryf) call defdr(odovcf,-4*nl2**2*nsites*nspc**2)
        endif
      endif

C --- Make the Hamiltonian and overlap, one site at a time ---
      nsts = nsites*nspc**2
      nsited = 0
      ohcf = 1
      ilme = nlme(nl)
      call defdr(owk1,   ilme)
      call defdr(owk2,   ilme)
      call defdr(owk3, 9*ilme)
      call defdr(owk4, nterm*ilme)
      if (hdiffs .or. fitpar) call defdr(odv, ilme)
      if (cryf .or. ocryf) call defdr(ohcf,nl2**2)

C --- Loop over spin combinations (spin-orbit) ---
      do  30  isp1 = 1, nspc
      do  30  isp2 = 1, nspc
      nsites = 0

C --- Loop over sites (clusters) ---
C oiaxc is the offset to iax in the current cluster
      do  30  iat = 1, nbas
        oiaxc = oiax + niax*nsites
        ndimL = npr(0,iat)*nl2
        call defdr(ohj,-ndimL*nl2)
        if (fitpar) then
          call defdr(odhrj,-ndimL*nl2*nvar)
        elseif (hdiffs) then
          call defdr(odhxj,-ndimL*nl2)
          call defdr(odhyj,-ndimL*nl2)
          call defdr(odhzj,-ndimL*nl2)
          call defdr(odhrj,-ndimL*nl2)
        endif

C --- Add to Hamiltonian ---
        call tbham1(nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,ndimL,
     .    memode,decay,iam,npm,nterm,tabme,npr(0,iat),isp1,isp2,
     .    w(oiaxc),cryf,deccf,tabcf,w(ohcf),w(ohj))
        call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(ohj),nsited,
     .              w(oh))

C --- Add to overlap matrix ---
        if (lov) then
          call dpzero(w(ohj),ndimL*nl2)
          call tbham1(nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,
     .      ndimL,memode,decov,iam,npm,nterm,tabov,npr(0,iat),isp1,isp2,
     .      w(oiaxc),ocryf,dcocf,tbocf,w(ohcf),w(ohj))
          call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(ohj),nsited,
     .                w(oov))
        endif

C --- Add to derivatives of Hamiltonian ---
        if (fitpar) then
C --- Parameter derivatives of Hamiltonian ---
          call tbdhm2(1,nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,
     .      ndimL,memode,decay,iam,npm,nterm,nset,tabme,npr(0,iat),isp1,
     .      isp2,w(oiaxc),ip1,w(owk2),w(odv),w(owk1),w(owk3),w(owk4),
     .      w(odhrj))

C --- Parameter derivatives of crystal field matrix ---
          if (cryf) call tbdhm2(2,nbas,w(oipc),alat,plat,w(opos),nl,
     .      nlmesp,ilme,ndimL,memode,deccf,iam,npm,nterm,nset,tabcf,
     .      npr(0,iat),isp1,isp2,w(oiaxc),ip1,w(owk2),w(odv),w(owk1),
     .      w(owk3),w(owk4),w(odhrj))

C --- Parameter derivatives of overlap matrix ---
          if (lov) then
            call tbdhm2(3,nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,
     .        ndimL,memode,decov,iam,npm,nterm,nset,tabov,npr(0,iat),
     .        isp1,isp2,w(oiaxc),ip1,w(owk2),w(odv),w(owk1),w(owk3),
     .        w(owk4),w(odhrj))

C --- Parameter derivatives of overlap crystal field matrix ---
            if (ocryf) call tbdhm2(4,nbas,w(oipc),alat,plat,w(opos),nl,
     .        nlmesp,ilme,ndimL,memode,dcocf,iam,npm,nterm,nset,
     .        tbocf,npr(0,iat),isp1,isp2,w(oiaxc),ip1,w(owk2),
     .        w(odv),w(owk1),w(owk3),w(owk4),w(odhrj))
          endif

C --- Add to parameter derivatives ---
          call xaddts(ndimL,nl2,nvar,nsts,npr(0,iat),w(oiaxc),nsited,
     .      w(odhrj),w(odh))

        elseif (hdiffs) then
C --- Add to derivatives of Hamiltonian ---
          call tbdhm1(nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,
     .      ndimL,memode,decay,iam,npm,nterm,tabme,npr(0,iat),isp1,isp2,
     .      w(oiaxc),.false.,w(odhxj),w(odhyj),w(odhzj),w(odhrj))
          pointr = odh
          call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhxj),
     .                nsited,w(pointr))
          pointr = pointr + nl2**2*nsts*i1mach(18)
          call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhyj),
     .                nsited,w(pointr))
          pointr = pointr + nl2**2*nsts*i1mach(18)
          call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhzj),
     .                nsited,w(pointr))
          pointr = pointr + nl2**2*nsts*i1mach(18)
          call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhrj),
     .                nsited,w(pointr))

C --- Add to derivatives of crystal field matrix ---
          if (cryf) then
            call dpzero(odhxj,-ndimL*nl2)
            call dpzero(odhyj,-ndimL*nl2)
            call dpzero(odhzj,-ndimL*nl2)
            call dpzero(odhrj,-ndimL*nl2)
            call tbdhm1(nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,
     .        ndimL,memode,deccf,iam,npm,nterm,tabcf,npr(0,iat),isp1,
     .        isp2,w(oiaxc),cryf,w(odhxj),w(odhyj),w(odhzj),w(odhrj))
            pntcf = odhcf
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhxj),
     .                  nsited,w(pntcf))
            pntcf = pntcf + nl2**2*nsts*i1mach(18)
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhyj),
     .                  nsited,w(pntcf))
            pntcf = pntcf + nl2**2*nsts*i1mach(18)
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhzj),
     .                  nsited,w(pntcf))
            pntcf = pntcf + nl2**2*nsts*i1mach(18)
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhrj),
     .                  nsited,w(pntcf))
          endif

C --- Add to derivatives of overlap matrix ---
          if (lov) then
            call dpzero(odhxj,-ndimL*nl2)
            call dpzero(odhyj,-ndimL*nl2)
            call dpzero(odhzj,-ndimL*nl2)
            call dpzero(odhrj,-ndimL*nl2)
            call tbdhm1(nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,
     .        ndimL,memode,decov,iam,npm,nterm,tabov,npr(0,iat),isp1,
     .        isp2,w(oiaxc),.false.,w(odhxj),w(odhyj),w(odhzj),w(odhrj))
            pntov = odov
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhxj),
     .                  nsited,w(pntov))
            pntov = pntov + nl2**2*nsts*i1mach(18)
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhyj),
     .                  nsited,w(pntov))
            pntov = pntov + nl2**2*nsts*i1mach(18)
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhzj),
     .                  nsited,w(pntov))
            pntov = pntov + nl2**2*nsts*i1mach(18)
            call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhrj),
     .                  nsited,w(pntov))

C --- Add to derivatives of overlap crystal field matrix ---
            if (ocryf) then
              call dpzero(odhxj,-ndimL*nl2)
              call dpzero(odhyj,-ndimL*nl2)
              call dpzero(odhzj,-ndimL*nl2)
              call dpzero(odhrj,-ndimL*nl2)
              call tbdhm1(nbas,w(oipc),alat,plat,w(opos),nl,nlmesp,ilme,
     .          ndimL,memode,dcocf,iam,npm,nterm,tbocf,npr(0,iat),isp1,
     .          isp2,w(oiaxc),ocryf,w(odhxj),w(odhyj),w(odhzj),w(odhrj))
              pntocf = odovcf
              call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhxj),
     .                    nsited,w(pntocf))
              pntocf = pntocf + nl2**2*nsts*i1mach(18)
              call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhyj),
     .                    nsited,w(pntocf))
              pntocf = pntocf + nl2**2*nsts*i1mach(18)
              call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhzj),
     .                    nsited,w(pntocf))
              pntocf = pntocf + nl2**2*nsts*i1mach(18)
              call addtos(nl2,ndimL,w(oiaxc),nl2,1,npr(0,iat),w(odhrj),
     .                    nsited,w(pntocf))
            endif
          endif
        endif

C --- Book keeping ---
        npr(1,iat) = nsites
        nsited = nsited + npr(0,iat)
        nsites = nsites + npr(0,iat)
        call rlse(ohj)
   30 continue
      if (cryf .or. ocryf) call rlse(ohcf)
      call rlse(owk1)

C --- Fix up inequivalent pairs: (sp,ps), (sd,ds), (pd,dp) ---
      if (nl2 .lt. 2) return
      call swapV(nl2,nsites,nspc,w(oiax),w(oh))
      if (lov) call swapV(nl2,nsites,nspc,w(oiax),w(oov))
      if (fitpar) then
        call xswap(nl2,nsites,nspc,nvar,ivar,w(oiax),w(odh))
      elseif (hdiffs) then
        call swapdV(nl2,nsites,nspc,w(oiax),w(odh))
        call swapV(nl2,nsites,nspc,w(oiax),w(pointr))
        if (lov) then
          call swapdV(nl2,nsites,nspc,w(oiax),w(odov))
          call swapV(nl2,nsites,nspc,w(oiax),w(pntov))
        endif
      endif

C --- iax(niax,i) = i ---
C      oiaxc = oiax-1
C      do  40  i = 1, nsites
C   40 w(oiaxc+niax*i) = i

      end

      subroutine tbham1(nbas,ipc,alat,plat,bas,nl,nlmesp,ilme,ndimL,
     .  memode,decay,iam,npm,nterm,tabme,npr,isp1,isp2,iax,cryf,deccf,
     .  tabcf,hcf,hj)
C- Real-space Slater-Koster tight-binding Hamiltonian for 1 site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas,ipc,alat,plat,bas,nl
Ci   nlmesp: =nlme(nl)*(2*nspc-1), number of TB MEs incl. spin
Ci   ilme: =nlme(nl), number of TB MEs
Ci   ndimL: first dimension of hj
Ci   memode,decay,iam,npm,nterm,tabme,deccf,tabcf: see rdtbh
Ci   npr: number of neighbors for this cluster
Ci   isp1,isp2: spin combination (spin orbit)
Ci   iax: neighbor list for this cluster
Ci   cryf: T, include crystal field terms; F, do not
Ci   hcf: work array for crystal field terms
Cl Local variables
Cl   wk1,wk2,wk3,wk4: work arrays
Co Outputs
Co   hj(ndimL,nl**2): Hamiltonian for one site
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nl,nlmesp,ilme,ndimL,memode,nterm,npr,isp1,isp2,niax
      parameter (niax=10)
      integer ipc(1),iam(3,1),npm(0:1,1),iax(niax,1)
      double precision alat
      double precision plat(3,3),bas(3,nbas),decay(nlmesp,2),
     .  tabme(nterm,nlmesp,1),deccf(nlmesp,2),
     .  tabcf(nterm,nlmesp,1),hcf(nl**2,nl**2),hj(ndimL,nl**2)
      logical cryf
C Local parameters
      integer i,j,k,l,m,nl2,ic1,ic2,idum
      integer itab(2,2)
      double precision wk(0:3),wk1(ilme),wk2(ilme),wk3(3,3,ilme),
     .  wk4(nterm,ilme)
      double precision n,nc,r0,rc,r,r1,A,c(0:5),p,val,slo,curv
      double precision d1mach
      data itab /0,2,2,1/
      p(r) = c(0)+c(1)*r+c(2)*r**2+c(3)*r**3+c(4)*r**4+c(5)*r**5 

      nl2 = nl**2
      idum = 1 + itab(isp1,isp2)*ilme
      j = 1
C --- Loop over neighbors ---
      do  110  i = 2, npr
        if (cryf) call dpzero(hcf,nl2**2)
        j = j + nl2
        call dlmn(nbas,plat,bas,iax(1,i),wk)
        ic1 = ipc(iax(1,i))
        ic2 = ipc(iax(2,i))
        if (memode .ne. 1) then
          call meptr(ic1,ic2,iam,npm,k)
          if (k .eq. 0) goto 110
        endif

C --- Fixed MEs ---
        if (memode .eq. 0 .or. memode .eq. 6) then
          call skham(wk(1),wk(2),wk(3),tabme(1,idum,k),nl,ndimL,hj(j,1))
          if (memode .eq. 6) call
     .      gskham(wk(1),wk(2),wk(3),tabme(1,idum,k),nl,ndimL,hj(j,1))
C --- Crystal field terms ---
          if (cryf) then
            call skham(wk(1),wk(2),wk(3),tabcf(1,idum,k),nl,nl2,hcf)
            if (memode .eq. 6) call
     .        gskham(wk(1),wk(2),wk(3),tabcf(1,idum,k),nl,nl2,hcf)
          endif

C --- Harrison universal MEs ---
        elseif (memode .eq. 1) then
          call rxx(nl .gt. 2,'TBHAM1: No universal Hamiltonian for l>1')
          call ume(ilme,alat*wk(0),wk1)
          call skham(wk(1),wk(2),wk(3),wk1,nl,ndimL,hj(j,1))

C --- Exponential or power decay ---
        elseif (memode .eq. 2 .or. memode .eq. 3) then
          call dcopy(ilme,tabme(1,idum,k),1,wk1,1)
          call dcopy(ilme,decay(idum,k),1,wk2,1)
          do  10  l = 1, ilme
            if (memode .eq. 2) then
              wk1(l) = wk1(l)*dexp(-alat*wk(0)*wk2(l))
            elseif (memode .eq. 3) then
              wk1(l) = wk1(l)*(alat*wk(0))**(-wk2(l))
            endif
   10     continue
          call skham(wk(1),wk(2),wk(3),wk1,nl,ndimL,hj(j,1))
C --- Crystal field terms ---
          if (cryf) then
            call dcopy(ilme,tabcf(1,idum,k),1,wk1,1)
            call dcopy(ilme,deccf(idum,k),1,wk2,1)
            do  20  l = 1, ilme
              if (memode .eq. 2) then
                wk1(l) = wk1(l)*dexp(-alat*wk(0)*wk2(l))
              elseif (memode .eq. 3) then
                wk1(l) = wk1(l)*(alat*wk(0))**(-wk2(l))
              endif
   20       continue
            call skham(wk(1),wk(2),wk(3),wk1,nl,nl2,hcf)
          endif

C --- ME = \sum_i=1,3 a_i d^b_i exp(-c_i d) ---
        elseif (memode .eq. 4) then
          call dcopy(9*ilme,tabme(1,idum,k),1,wk3,1)
          call dpzero(wk1,ilme)
          do  40  l = 1, ilme
            do  30  m = 1, 3
              wk1(l) = wk1(l) + wk3(1,m,l)*(alat*wk(0))**wk3(2,m,l)
     .                                 *dexp(-alat*wk(0)*wk3(3,m,l))
   30       continue
   40     continue
          call skham(wk(1),wk(2),wk(3),wk1,nl,ndimL,hj(j,1))
C --- Crystal field terms ---
          if (cryf) then
            call dcopy(9*ilme,tabcf(1,idum,k),1,wk3,1)
            call dpzero(wk1,ilme)
            do  60  l = 1, ilme
              do  50  m = 1, 3
                wk1(l) = wk1(l) + wk3(1,m,l)*(alat*wk(0))**wk3(2,m,l)
     .                                   *dexp(-alat*wk(0)*wk3(3,m,l))
   50         continue
   60       continue
            call skham(wk(1),wk(2),wk(3),wk1,nl,nl2,hcf)
          endif

C --- Goodwin-Skinner-Pettifor: V (r0/d)^n exp[n(-{d/rc}^nc+{r0/rc}^nc)]
        elseif (memode .eq. 5) then
          call dcopy(5*ilme,tabme(1,idum,k),1,wk4,1)
          do  70  l = 1, ilme
            n  = wk4(2,l)
            nc = wk4(3,l)
            r0 = wk4(4,l)
            rc = wk4(5,l)
            r  = alat*wk(0)
            if (nc .gt. -d1mach(1)) then
              if (wk4(1,l) .eq. 0) then
                wk1(l) = 0
              else
                A  = wk4(1,l)*r0**n*exp(n*(r0/rc)**nc)
                wk1(l) = A*exp(-n*(r/rc)**nc)/r**n
              endif
            else
C --- Znam et al. augmented power law
              r1 = -nc
              if (r .gt. rc .or. wk4(1,l) .eq. 0) then
                wk1(l) = 0d0
              elseif (r .lt. r1) then
                wk1(l) = wk4(1,l)*r0**n / r**n
              else
                val = wk4(1,l)*r0**n / r1**n
                slo = -n*wk4(1,l)*r0**n / r1**(n+1)
                curv = n*(n+1)*wk4(1,l)*r0**n / r1**(n+2)
                call rcut(val,slo,curv,r1,rc,c)
                wk1(l) = p(r)
              endif
            endif
   70     continue
          call skham(wk(1),wk(2),wk(3),wk1,nl,ndimL,hj(j,1))
C --- Crystal field terms ---
          if (cryf) then
            call dcopy(5*ilme,tabcf(1,idum,k),1,wk4,1)
            do  80  l = 1, ilme
              if (wk4(1,l) .eq. 0) then
                wk1(l) = 0
              else
                n  = wk4(2,l)
                nc = wk4(3,l)
                r0 = wk4(4,l)
                rc = wk4(5,l)
                r  = alat*wk(0)
                A  = wk4(1,l)*r0**n*exp(n*(r0/rc)**nc)
                wk1(l) = A*exp(-n*(r/rc)**nc)/r**n
              endif
   80       continue
            call skham(wk(1),wk(2),wk(3),wk1,nl,nl2,hcf)
          endif

C --- ME = a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc) ---
        elseif (memode .eq. 7) then
          call dcopy(4*ilme,tabme(1,idum,k),1,wk4,1)
          do  90  l = 1, ilme
            wk1(l) = wk4(1,l)*(alat*wk(0))**(-wk4(2,l))
     .             / (1 + dexp(wk4(3,l)*(alat*wk(0) - wk4(4,l))))
   90     continue
          call skham(wk(1),wk(2),wk(3),wk1,nl,ndimL,hj(j,1))
C --- Crystal field terms ---
          if (cryf) then
            call dcopy(4*ilme,tabcf(1,idum,k),1,wk4,1)
            do  100  l = 1, ilme
              wk1(l) = wk4(1,l)*(alat*wk(0))**(-wk4(2,l))
     .               / (1 + dexp(wk4(3,l)*(alat*wk(0) - wk4(4,l))))
  100       continue
            call skham(wk(1),wk(2),wk(3),wk1,nl,nl2,hcf)
          endif

        else
          call rx('TBHAM1: bad memode')
        endif

C --- Add crystal field terms to hj ---
        if (cryf) call dmadd(hcf,nl2,1,1d0,hj(1,1),ndimL,1,1d0,
     .    hj(1,1),ndimL,1,nl2,nl2)
  110 continue
C      call shoblk(npr,1,ndimL,nl2,nl2,iax,hj,0)
      end

      subroutine tbdhm1(nbas,ipc,alat,plat,bas,nl,nlmesp,ilme,ndimL,
     .  memode,decay,iam,npm,nterm,tabme,npr,isp1,isp2,iax,cryf,dhxj,
     .  dhyj,dhzj,dhrj)
C- Real-space Slater-Koster tight-binding Hamiltonian deriv. for 1 site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas,ipc,alat,plat,bas,nl,nspc
Ci   ndimL: first dimension of dh
Ci   memode,decay,iam,npm,nterm,tabme: see rdtbh
Ci   npr: number of neighbors for this cluster
Ci   isp1,isp2: spin combination (spin orbit)
Ci   iax: neighbor list for this cluster
Ci   cryf: if true then crystal field derivatives
Co Outputs
Co   dhxj,dhyj,dhzj,dhrj: x,y,z, and radial deriv. of Ham for one site
Cl Local variables
Cl   V,dV,wk1,wk3,wk4: work arrays
Cr Remarks
Cu Updates
Cu   23 May 07 (MvS) bug fix, GSP mode
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nl,nlmesp,ilme,ndimL,memode,nterm,npr,isp1,isp2,niax
      parameter (niax=10)
      integer ipc(1),iam(3,1),npm(0:1,1),iax(niax,1)
      double precision alat
      double precision plat(3,3),bas(3,nbas),decay(nlmesp,2),
     .  tabme(nterm,nlmesp,1),dhrj(ndimL,nl**2),
     .  dhxj(ndimL,nl**2),dhyj(ndimL,nl**2),dhzj(ndimL,nl**2)
      logical cryf
C Local parameters
      integer i,j,k,l,m,nl2,ic1,ic2,idum
      integer itab(2,2)
      double precision wk(0:3),wk1(ilme),wk3(3,3,ilme),wk4(nterm,ilme),
     .  V(ilme),dV(ilme)
      double precision n,nc,r0,rc,r,r1,A,c(0:5),p,pp,val,slo,curv
      double precision d1mach
      data itab /0,2,2,1/
C 5th order polynomial and its derivatives
      p(r) = c(0)+c(1)*r+c(2)*r**2+c(3)*r**3+c(4)*r**4+c(5)*r**5 
      pp(r) = c(1)+2d0*c(2)*r+3d0*c(3)*r**2+4d0*c(4)*r**3+5d0*c(5)*r**4

      nl2 = nl**2
      idum = 1 + itab(isp1,isp2)*ilme
      j = 1
C --- Loop over neighbors ---
      do  70  i = 2, npr
        j = j + nl2
        call dlmn(nbas,plat,bas,iax(1,i),wk)
        ic1 = ipc(iax(1,i))
        ic2 = ipc(iax(2,i))
        if (memode .ge. 2) then
          call meptr(ic1,ic2,iam,npm,k)
          if (k .eq. 0) goto 70
        endif

C --- Fixed MEs ---
        if (memode .eq. 0 .or. memode .eq. 6) then
          call rx('TBDHM1: No Hamiltonian derivatives for modes 0 or 6')

C --- Harrison unversal MEs ---
        elseif (memode .eq. 1) then
          call rxx(nl .gt. 2,'TBDHM1: No universal Hamiltonian for l>1')
          call ume(ilme,alat*wk(0),V)
          call dcopy(ilme,V,1,dV,1)
          call dscal(ilme,-2d0/(wk(0)*alat),dV,1)

C --- Exponential decay ---
        elseif (memode .eq. 2) then
          call dcopy(ilme,tabme(1,idum,k),1,V,1)
          call dcopy(ilme,decay(idum,k),1,wk1,1)
          do  10  l = 1, ilme
            V(l)  = V(l)*dexp(-wk1(l)*alat*wk(0))
            dV(l) = -wk1(l)*V(l)
   10     continue

C --- Power decay ---
        elseif (memode .eq. 3) then
          call dcopy(ilme,tabme(1,idum,k),1,V,1)
          call dcopy(ilme,decay(idum,k),1,wk1,1)
          do  20  l = 1, ilme
            V(l)  = V(l)*(alat*wk(0))**(-wk1(l))
            dV(l) = -wk1(l)/(alat*wk(0))*V(l)
   20     continue

C --- ME = \sum_i=1,3 a_i d^b_i exp(-c_i d) ---
        elseif (memode .eq. 4) then
          call dcopy(9*ilme,tabme(1,idum,k),1,wk3,1)
          call dpzero(V,ilme)
          call dpzero(dV,ilme)
          do  40  l = 1, ilme
            do  30  m = 1, 3
              wk1(l) = wk3(1,m,l)*(alat*wk(0))**wk3(2,m,l)
     .                        *dexp(-alat*wk(0)*wk3(3,m,l))
              V(l) = V(l) + wk1(l)
              dV(l) = dV(l) + wk3(2,m,l)/(alat*wk(0))*wk1(l)
     .                      - wk3(3,m,l)*wk1(l)
   30       continue
   40     continue

C --- Goodwin-Skinner-Pettifor: V (r0/d)^n exp[n(-{d/rc}^nc+{r0/rc}^nc)]
        elseif (memode .eq. 5) then
          call dcopy(5*ilme,tabme(1,idum,k),1,wk4,1)
          do  50  l = 1, ilme
            n  = wk4(2,l)
            nc = wk4(3,l)
            r0 = wk4(4,l)
            rc = wk4(5,l)
            r  = alat*wk(0)
            if (nc .gt. -d1mach(1)) then
              if (wk4(1,l) .eq. 0) then
                V(l) = 0
                dV(l) = 0
              else
                A  = wk4(1,l)*r0**n*exp(n*(r0/rc)**nc)
                V(l) = A*exp(-n*(r/rc)**nc)/r**n
                dV(l) = -n*V(l)*(1d0 + nc*(r/rc)**nc)/r
              endif
            else
C --- Znam et al. augmented power law
              r1 = -nc
              if (r .gt. rc) then
                V(l) = 0d0
                dV(l) = 0d0
              elseif (r .lt. r1 .or. wk4(1,l) .eq. 0) then
                V(l) = wk4(1,l)*r0**n / r**n
                dV(l) = -n*V(l) / r
              else
                val = wk4(1,l)*r0**n / r1**n
                slo = -n*wk4(1,l)*r0**n / r1**(n+1)
                curv = n*(n+1)*wk4(1,l)*r0**n / r1**(n+2)
                call rcut(val,slo,curv,r1,rc,c)
                V(l) = p(r)
                dV(l) = pp(r)
              endif
            endif
   50     continue

C --- ME = a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc) ---
        elseif (memode .eq. 7) then
          call dcopy(4*ilme,tabme(1,idum,k),1,wk4,1)
          do  60  l = 1, ilme
            A = dexp(wk4(3,l)*(alat*wk(0) - wk4(4,l)))
            V(l) = wk4(1,l)*(alat*wk(0))**(-wk4(2,l)) / (1 + A)
            dV(l) = -wk4(2,l)/(alat*wk(0))*V(l) 
     .            - wk4(3,l)*A*V(l) / (1 + A)
   60     continue

        else
          call rx('TBDHM1: bad memode')
        endif

C --- Get Hamiltonian derivatives ---
        call dskham(wk(1),wk(2),wk(3),wk(0)*alat,V,dV,nl,ndimL,cryf,
     .    dhxj(j,1),dhyj(j,1),dhzj(j,1),dhrj(j,1))

   70 continue

C      call shoblk(npr,1,ndimL,nl2,nl2,iax,dhrj,0)
C      call shoblk(npr,1,ndimL,nl2,nl2,iax,dhxj,0)
C      call shoblk(npr,1,ndimL,nl2,nl2,iax,dhyj,0)
C      call shoblk(npr,1,ndimL,nl2,nl2,iax,dhzj,0)
C      pause

      end

      subroutine tbdhm2(ityp,nbas,ipc,alat,plat,bas,nl,nlmesp,ilme,
     .  ndimL,memode,decay,iam,npm,nterm,nset,tabme,npr,isp1,isp2,iax,
     .  ip1,V,dV,wk1,wk3,wk4,dhj)
C- Real-space Slater-Koster Hamiltonian deriv. wrt parameters for 1 site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ityp: ME type (1=Ham, 2=CF, 3=Overlap, 4=Overlap CF)
Ci   nbas,ipc,alat,plat,bas,nl
Ci   nlmesp: =nlme(nl)*(2*nsp-1), number of TB MEs incl. spin
Ci   ilme: =nlme(nl), number of TB MEs
Ci   ndimL: first dimension of dhj
Ci   memode,decay,iam,npm,nterm,nset,tabme: see rdtbh
Ci   npr: number of neighbors for this cluster
Ci   isp1,isp2: spin combination (spin orbit)
Ci   iax: neighbor list for this cluster
Ci   ip1: pointer to locations in full list of variables
Ci   V,dV,wk1,wk3,wk4: work arrays
Co Outputs
Co   dhj: deriv. wrt TB parameters for one site
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ityp,nbas,nl,nlmesp,ilme,ndimL,memode,nterm,nset,npr,isp1,
     .  isp2
      integer ipc(1),iam(3,1),npm(0:1,1),iax(10,1),
     .  ip1(nterm+1,nlmesp,nset,4)
      double precision alat
      double precision plat(3,3),bas(3,nbas),decay(nlmesp,2),
     .  tabme(nterm,nlmesp,1),V(ilme),dV(ilme),wk1(ilme),wk3(3,3,ilme),
     .  wk4(nterm,ilme),dhj(ndimL,nl**2,2)
C Local parameters
      integer i,j,k,l,m,nl2,ic1,ic2,idum,iv,ii,kold
      integer itab(2,2)
      double precision wk(0:3)
      double precision n,nc,r0,rc,r,A
      data itab /0,2,2,1/

      nl2 = nl**2
      idum = 1 + itab(isp1,isp2)*ilme
      j = 1
C --- Loop over neighbors ---
      do  110  i = 2, npr
        j = j + nl2
        call dlmn(nbas,plat,bas,iax(1,i),wk)
        ic1 = ipc(iax(1,i))
        ic2 = ipc(iax(2,i))
        do  100  ii = 1, 2
          if (memode .ne. 1) then
            if (ii .eq. 1) then
              call meptr(ic1,ic2,iam,npm,k)
              if (k .eq. 0) goto 100
              kold = k
            else
              call meptr(ic2,ic1,iam,npm,k)
              if (k .eq. 0 .or. k .eq. kold) goto 100
              call dscal(3,-1d0,wk(1),1)
            endif
          endif

C --- Fixed MEs ---
          if (memode .eq. 0) then
            do  10  l = 1, ilme
              iv = ip1(1,l+idum-1,k,ityp)
              if (iv .eq. 0) goto 10
              call dpzero(dV,ilme)
              dV(l) = 1d0
              call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
   10       continue

C --- Fixed MEs + extension ---
          elseif (memode .eq. 6) then
            do  20  l = 1, ilme
              iv = ip1(1,l+idum-1,k,ityp)
              if (iv .eq. 0) goto 20
              call dpzero(wk4,2*ilme)
              wk4(l,1) = 1d0
              call skham(wk(1),wk(2),wk(3),wk4,nl,ndimL,dhj(j,1,iv))
   20       continue
            do  30  l = ilme+1, 2*ilme
              iv = ip1(2,l+idum-1,k,ityp)
              if (iv .eq. 0) goto 30
              call dpzero(wk4,2*ilme)
              wk4(l,1) = 1d0
              call gskham(wk(1),wk(2),wk(3),wk4,nl,ndimL,dhj(j,1,iv))
   30       continue

C --- Exponential decay ---
          elseif (memode .eq. 2) then
            call dcopy(ilme,tabme(1,idum,k),1,V,1)
            call dcopy(ilme,decay(idum,k),1,wk1,1)
            do  40  l = 1, ilme
              call dpzero(dV,ilme)
              dV(l) = dexp(-wk1(l)*alat*wk(0))
              iv = ip1(1,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              iv = ip1(2,l+idum-1,k,ityp)
              if (iv .eq. 0) goto 40
              dV(l) = -alat*wk(0)*V(l)*dV(l)
              call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
   40       continue

C --- Power decay ---
          elseif (memode .eq. 3) then
            call dcopy(ilme,tabme(1,idum,k),1,V,1)
            call dcopy(ilme,decay(idum,k),1,wk1,1)
            do  50  l = 1, ilme
              call dpzero(dV,ilme)
              dV(l) = (alat*wk(0))**(-wk1(l))
              iv = ip1(1,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              iv = ip1(2,l+idum-1,k,ityp)
              if (iv .eq. 0) goto 50
              dV(l) = -dlog(alat*wk(0))*V(l)*dV(l)
              call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
   50       continue

C --- ME = \sum_i=1,3 a_i d^b_i exp(-c_i d) ---
          elseif (memode .eq. 4) then
            call dcopy(9*ilme,tabme(1,idum,k),1,wk3,1)
            do  70  l = 1, ilme
              call dpzero(dV,ilme)
              do  60  m = 1, 3
                dV(l) = (alat*wk(0))**wk3(2,m,l)
     .            *dexp(-alat*wk(0)*wk3(3,m,l))
                iv = ip1(1+3*(m-1),l+idum-1,k,ityp)
                if (iv .gt. 0)
     .            call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
                V(l) = wk3(1,m,l)*dV(l)
                dV(l) = dlog(alat*wk(0))*V(l)
                iv = ip1(2+3*(m-1),l+idum-1,k,ityp)
                if (iv .gt. 0)
     .            call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
                iv = ip1(3+3*(m-1),l+idum-1,k,ityp)
                if (iv .eq. 0) goto 60
                dV(l) = -alat*wk(0)*V(l)
                call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
   60         continue
   70       continue

C --- Goodwin-Skinner-Pettifor: V (r0/d)^n exp[n(-{d/rc}^nc+{r0/rc}^nc)]
          elseif (memode .eq. 5) then
            call dcopy(5*ilme,tabme(1,idum,k),1,wk4,1)
            do  80  l = 1, ilme
              call dpzero(dV,ilme)
              n  = wk4(2,l)
              nc = wk4(3,l)
              r0 = wk4(4,l)
              rc = wk4(5,l)
              r  = alat*wk(0)
              A  = r0**n*exp(n*(r0/rc)**nc)
              dV(l) = A*exp(-n*(r/rc)**nc)/r**n
              iv = ip1(1,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              V(l) = wk4(1,l)*dV(l)
              dV(l) = (dlog(r0/r) - (r/rc)**nc + (r0/rc)**nc)*V(l)
              iv = ip1(2,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              dV(l) = (-dlog(r/rc)*(r/rc)**nc
     .              + dlog(r0/rc)*(r0/rc)**nc)*n*V(l)
              iv = ip1(3,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              dV(l) = (1 + nc*(r0/rc)**nc)*n*V(l)/r0
              iv = ip1(4,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              iv = ip1(5,l+idum-1,k,ityp)
              if (iv .eq. 0) goto 80
              dV(l) = ((r/rc)**nc - (r0/rc)**nc)*n*nc*V(l)/rc
              call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
   80       continue

C --- ME = a d^-b / {1 + exp[c(d - d0)]} (Sawada, Kohyama, etc) ---
          elseif (memode .eq. 7) then
            call dcopy(4*ilme,tabme(1,idum,k),1,wk4,1)
            do  90  l = 1, ilme
              call dpzero(dV,ilme)
              A = dexp(wk4(3,l)*(alat*wk(0) - wk4(4,l)))
              dV(l) = (alat*wk(0))**(-wk4(2,l)) / (1 + A)
              iv = ip1(1,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              V(l) = wk4(1,l)*dV(l)
              dV(l) = -dlog(alat*wk(0))*V(l)
              iv = ip1(2,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              dV(l) = -(alat*wk(0) - wk4(4,l))*A*V(l) / (1 + A)
              iv = ip1(3,l+idum-1,k,ityp)
              if (iv .gt. 0)
     .          call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
              iv = ip1(4,l+idum-1,k,ityp)
              if (iv .eq. 0) goto 90
              dV(l) = wk4(3,l)*A*V(l) / (1 + A)
              call skham(wk(1),wk(2),wk(3),dV,nl,ndimL,dhj(j,1,iv))
   90       continue

          else
            call rx('TBDHM2: bad memode')
          endif

  100   continue
  110 continue

      end
      subroutine xaddts(ndimL,nl2,nvar,nsts,npr,iaxc,nsited,dhj,dh)
C- Help routine to add to matrix of parameter derivatives
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl2,nsts,iaxc
Ci   ndimL: first dimension of dhj
Ci   nvar: total number of TB parameters to vary
Ci   npr: number of neighbors for this cluster
Ci   iaxc: iax for the current cluster
Ci   nsited: total number of sites, updated
Ci   dhj: deriv. wrt TB parameters for one site
Co Outputs
Co   dh: deriv. wrt TB parameters for all sites, updated
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndimL,nl2,nvar,nsts,npr(0:1),iaxc(1),nsited
      double precision dhj(ndimL*nl2,nvar),dh(nl2**2*nsts,nvar)
C Local parameters
      integer iv

      do  10  iv = 1, nvar
      call addtos(nl2,ndimL,iaxc,nl2,1,npr(0),dhj(1,iv),nsited,dh(1,iv))
   10 continue

      end
      subroutine xswap(nl2,nsites,nsp,nvar,ivar,iax,dh)
C- Help routine to swap inequivalent pairs of parameter derivatives
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl2,nsites,nsp,iax
Ci   nvar: total number of TB parameters to vary
Ci   ivar(1,i): points to position of ith variable
Ci   ivar(2,i): paramter type of ith variable (1=Ham, 3=Overlap)
Co Outputs
Co   dh: deriv. wrt TB parameters for all sites, MEs swapped
Cr Remarks
Cr   Only Hamiltonian and overlap (if any) MEs are swapped.  Crystal
Cr   field MEs do not need to be swapped.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl2,nsites,nsp,nvar
      integer ivar(2,nvar),iax(1)
      double precision dh(nl2**2*nsites*nsp**2,nvar)
C Local parameters
      integer iv

C --- Swap Hamiltonian and overlap derivative MEs ---
      do  10  iv = 1, nvar
        if (ivar(2,iv) .eq. 1 .or. ivar(2,iv) .eq. 3)
     .    call swapV(nl2,nsites,nsp,iax,dh(1,iv))
   10 continue

      end
