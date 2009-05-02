      subroutine suham(sctrl,slat,sspec,ssite,sarray,sham,spot,sstrn)
C- Hamiltonian setup
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Cio  sctrl :struct for program flow parameters; see routine uctrl
Cio    Elts read: nbasp nspec nl lham lcgf lgen3 lfp lncol nspin lpgf
Cio    Stored:    *
Cio    Passed to: *
Cio  slat  :struct for lattice information; see routine ulat
Cio    Elts read: ng ogv alat tolft plat qlat
Cio    Stored:    *
Cio    Passed to: *
Cio  sspec :struct for species-specific information; see routine uspec
Cio    Elts read: lmxa lmxl p pz
Cio    Stored:    *
Cio    Passed to: atfold makidx showbs sugcut
Cio  ssite :struct for site-specific information; see routine usite
Cio    Elts read: spec
Cio    Stored:    pnu pz
Cio    Passed to: showbs pvioeu
Cio  sarray:structure containing offsets to various arrays
Cio    Elts read: oips oipc nclasp opos
Cio    Stored:    *
Cio    Passed to: *
Cio  sham  :struct for parameters defining hamiltonian; see routine uham
Cio    Elts read: pwmode pwemax lncol oeula neula qss omagf nbf
Cio    Stored:    ndham ndofH ooffH oindxo ldham lmxax hord obdots
Cio    Passed to: *
Cio  spot  :struct for information about the potential; see routine upot
Cio    Elts read: opp vmtz
Cio    Stored:    nlma nlml opti
Cio    Passed to: *
Cio  sstrn :struct for global strings
Cio    Elts read: gfopt
Cio    Stored:    *
Cio    Passed to: *
Ci   gfopts:string containing switches for GF program
Cl Local variables
Cl   ndim  :total number of lmto orbitals = nl**2 * nbas
Cl   npwmin:lower limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   npwmax:upper limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   nqdiv: loop over mesh of q-points to estimate npwmax
Cl        : nqdiv is fineness of q-mesh.
Cl  npwpad: a 'safety' padding to npwmax in case npwmax
Cl        : underestimates actual upper limit
Cr Remarks
Cr   This routine generates energy-independent hamiltonian setup.
Cr  *It generates and packs a table of hamiltonian offsets offH,
Cr   orbital permutation indices oindxo.
Cr
Cr  *For the ASA 2nd generation LMTO:
Cr   Extract order of potential function from gfopts
Cr   Transform pp's to alpha representation
Cu Updates
Cu   07 Jul 08 Make sham->ndham = estimate for upper dimension of
Cu             hamiltonian, including possible PW part
Cu             Make sham->lmxax = largest lmax in basis
Cu   18 Apr 05 Force small parameter p -> 0 in 2-center turned on (ASA)
Cu   14 Feb 03 Makes and packs magnetic field B.sigma
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision sctrl(1),slat(1),sspec(1),ssite(1),sarray(1),
     .  sham(1),spot(1)
      character*(*) sstrn
C ... Local parameters
      integer nkap0,n0H
      parameter (nkap0=3,n0H=5)
      integer hord,i,i1,i2,iprint,lcgf,ldham(16),ldim,lgen3,
     .  lgunit,lham,lidim,lihdim,lncol,lpgf,nbasp,nclasp,ndim,neul,nl,
     .  nlspcp,nsp,nspc,nspec,nspx,nttab,partok,igets,lfp,nvi,nvl,
     .  ib,is,lmxa,lmxl,stdo,nglob,nkaph,isw,nbf,lmxax,pwmode,
     .  j1,j2,j3,m,npw,npwmin,npwmax,ndham
      double precision pwemin,pwemax,plat(3,3),qlat(3,3),q(3),Gmin,Gmax,
     .  xx
      integer n0,nqdiv,npwpad
      parameter (n0=10,nqdiv=12)
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))

C ... Needed for iostr
      double precision ckbas,cksumf,kap2(10)
      integer oalph,oiax,ontab,os,ng,ogv,okv,oigv2
      integer obas,offH,oidxsh,oipc,oips,oo,opp,oeula,opti,obs,omagf
      logical ltmp,iostr,bittst,adnf
      character*80 outs,gfopts
      double precision qss(4),vmtz,pnu(n0,2),pz(n0,2),alat,tolgv
C ... For NMTO
C     integer olmx,ohcr,otral,otrad,oalpha,oadot,nmto
C     double precision kmto(10),avw,dglob
C ... Heap
      integer w(1)
      common /w/ w

C --- Setup ---
      call upack('array oips oipc',sarray,oips,oipc,0,0,0)
      call upack('ctrl nbasp nspec nl',sctrl,nbasp,nspec,nl,0,0)
      call upack('ctrl lham lcgf lgen3 lfp',sctrl,lham,lcgf,lgen3,
     .  lfp,0)
      call upack('ctrl lncol nspin',sctrl,lncol,nsp,0,0,0)
      lpgf  = igets('ctrl lpgf',sctrl)
      nkaph = nglob('nkaph')
      call upack('ham pwmode pwemin pwemax npwpad',sham,pwmode,pwemin,
     .  pwemax,npwpad,0)

      nspc  = 1
      nspx  = nsp
      if (lncol .ne. 0) then
        nspc = 2
        nspx = 1
      endif
      ndim = nbasp * nl**2 * nkaph
      call defi(offH, -n0H*nkap0*(nbasp+1))
      call defi(oidxsh,ndim+3)
      stdo = lgunit(1)

C --- Automatic downfolding (2nd generation only?) ---
      if (lgen3 .eq. 0 .and. lfp .eq. 0) then
        adnf = bittst(lham,4)
C       Suppress if in gamma repsn (not needed for GF mode)
C       if (lgors('ctrl lham,128',sctrl)) adnf = 0
        call upack2('pot opp vmtz',spot,opp,vmtz)
        call atfold(0,adnf,nl,nsp,vmtz,w(opp),sspec)
      endif

C --- Hamiltonian offsets, orbital permutation table ---
      call iinit(ldham,16)
      if (mod(pwmode,10) .eq. 2) then
        call info0(20,0,0,' suham: LMTO basis will be excluded')
      else
        call makidx(nl,nkaph,1,nbasp,0,sspec,w(oips),w(offH),w(oidxsh),
     .    ldham)
      endif
C     Default value for hamiltonian dimension
      call pack1('ham ndham',sham,ldham)
C     print *, '!!'; call pack1('ham ndham',sham,ldham+3)

C ... Make rest of ldham
      ldham(4) = nspc
      ldham(5) = ldham(1) * nspc
      ldham(6) = ldham(2) * nspc
      ldham(7) = ldham(3) * nspc
      ldham(8) = nspx
      call pack5('ham ndofH ooffH oindxo ldham',sham,nbasp+1,offH,
     .  oidxsh,ldham,0)

C ... Printout orbital positions in hamiltonian, resolved by l
      if (iprint() .gt. 50) then
        call showbs(ssite,sspec,nkaph,0,w(oidxsh),ldham)
      endif

C ------- Potential- and implementation- specific setup -------

C --- FP setup ---
      if (lfp .ne. 0) then
        nvi = 0
        nvl = 0
        lmxax = -1
        do  20  ib = 1, nbasp
          call upack2('site spec',ssite,ib,is)
          call upack('spec lmxa lmxl p pz',sspec,is,lmxa,lmxl,pnu,pz)
cccccccccccccc
c        print *,'suham 111 pnu=',i,pnu(1:10,1:2)
cccccccccccccc
          call dcopy(n0,pz,1,pz(1,2),1)
C         Augmentation dimensioning parameters
          nvi = nvi + (lmxa+1)**2
          nvl = nvl + (lmxl+1)**2
C         Poke spec starting pnu to site
          call pack5('site pnu pz',ssite,ib,pnu,pz,0,0)
C         Find largest lmxa
          lmxax = max(lmxax,lmxa)
   20   continue
        call pack1('ham lmxax',sham,lmxax)
        call pack2('pot nlma nlml',spot,nvi,nvl)
        call info5(30,0,0,' suham :  %i augmentation'//
     .    ' channels, %i local potential channels  Maximum lmxa=%i',
     .    nvi,nvl,lmxax,0,0)

        call upack('lat ng ogv alat tolft',slat,ng,ogv,alat,tolgv,0)
        call sugcut(1,nspec,sspec,alat,ng,w(ogv),tolgv)

C   ... PW setup : estimate upper bound to number of G vectors
C       to set up upper bound to hamiltonian dimension
        call upack2('ham pwmode pwemax',sham,pwmode,pwemax)
        if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then

          call upack3('lat alat plat qlat',slat,alat,plat,qlat)
          Gmin = dsqrt(pwemin)
          Gmax = dsqrt(pwemax)
          if (mod(pwmode/10,10) .eq. 1) then
            call info0(70,1,0,' Estimate max size of PW basis from'//
     .        'combinations of recip. lattice vectors ...')
            npwmax = -1
            npwmin = 99999
            do  j1 = 0, nqdiv
            do  j2 = 0, nqdiv
            do  j3 = 0, nqdiv
              do   m = 1, 3
                q(m) = (qlat(m,1)/nqdiv)*j1 +
     .                 (qlat(m,2)/nqdiv)*j2 +
     .                 (qlat(m,3)/nqdiv)*j3
              enddo
              call pshpr(iprint()-40)
ctakao nqdiv is large but fixed... So this is safer?
c              call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
c     .          xx,xx,xx)
c             write(6,"('GVLST2:goto q=',3f13.5)")q
c             call shorbz(q,q,qlat,plat)
             call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
     &          xx,xx,xx)
              call poppr
              npwmin = min(npwmin,npw)
              npwmax = max(npwmax,npw)
            enddo
            enddo
            enddo
            if (npwpad .lt. 0) then
              npwpad = max(nint((npwmax-npwmin)*0.2d0),3)
            endif
          else
            call dpzero(q,3)
            call pshpr(iprint()-40)
            call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
     .        xx,xx,xx)
            call poppr
            npwmin = npw
            npwmax = npw
            npwpad = 0
          endif
          ndham = npwmax + npwpad
          if (mod(pwmode,10) .ne. 2) ndham = ldham(1) + npwmax + npwpad
          call pack3('ham npwmin npwpad ndham',sham,npwmin,npwpad,ndham)
          if (mod(pwmode/10,10) .eq. 1) then
            call info2(20,1,0,' suham:  q-dependent PW basis with'//
     .        '  Emin = %d < E < %d.',pwemin,pwemax)
            call info5(30,0,0,'%9fEst. min,max PW dimension = %i,%i.'//
     .        '  Use npwpad = %i => ndham = %i',
     .        npwmin,npwmax,npwpad,ndham,0)
          else
            call info5(20,0,0,' suham:  PW basis with  %d < E < '//
     .        '%d  =>  npw = %i,  ndham = %i',
     .        pwemin,pwemax,npw,ndham,0)
          endif

C    ...  Printout
          if (iprint() .ge. 40) then
            call pshpr(0)
            call dpzero(q,3)
            call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
     .        xx,xx,xx)
            call poppr
            call info0(40,1,-1,' G vectors at the Gamma point:')
            call defi(oigv2,3*npw)
            call defi(okv,3*npw)
            call defrr(ogv,3*npw)
            call pshpr(iprint())
            if (iprint() .ge. 50) call setpr(100)
            call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,8+2,npw,npw,
     .        w(okv),w(ogv),xx,w(oigv2))
            call poppr
            call rlse(oigv2)
          endif
        endif
      endif

C --- Third-generation NMTO ---
      if (lgen3 .ne. 0) then

C        call upack('ham nmto kmto',sham,nmto,kmto,0,0,0)
C        call isanrg(nmto,2,5,'suham','nmto',.true.)
C        call defdr(otral,-4*nl**2*nbasp*nmto)
C        call defdr(otrad,-4*nl**2*nbasp*nmto)
C
C        call sp2cls('spec lmxb',sspec,0,1,1,nspec,olmx)
C        call sp2cls('spec hcr',sspec,0,nl,1,nspec,ohcr)
C        call defdr(oalpha,-nl**2*nbasp*nmto)
C        call defdr(oadot,-nl**2*nbasp*nmto)
C        avw = dglob('avw',0d0,0)
C        call dscal(nmto,avw**2,kmto,1)
C        call dscal(nl*nspec,1/avw,w(ohcr),1)
C        call pshpr(50)
C        call mktra2(1,loka,nbasp,w(oips),nl,w(olmx),avw,4,kmto,nmto,
C     .    w(ohcr),1,w(otral),w(otrad),w(oalpha),w(oadot),w,w)
C        call poppr
C
C        call dscal(nmto,1/avw**2,kmto,1)
C        call dscal(nl*nspec,avw,w(ohcr),1)
C        call rlse(olmx)
      endif

C --- Second-generation LMTO, ASA ---
C ... Green's function-specific initialization
      if (lcgf .ne. 0 .or. lpgf .ne. 0) then
        call upacks('strn gfopt',i1,i2)
        gfopts = ' '
        if (i2 .ge. i1) then
          gfopts = sstrn(i1:i2)
          call partk0(0,len(gfopts),1,-1,0,len(gfopts),-1,31,.false.)
          hord = 2
          ltmp = .false.
          i = partok(gfopts,'pz',  ' ;',ltmp,' ',0,0,0,0)
          if (ltmp) hord = 4
          i = partok(gfopts,'p3',  ' ;',ltmp,' ',0,0,0,0)
          if (ltmp) hord = 3
          i = partok(gfopts,'p1',  ' ;',ltmp,' ',0,0,0,0)
          if (ltmp) hord = 1
          call pack1('ham hord',sham,hord)
          if (iprint() .ge. 30) call awrit1(' SUHAM : 2nd generation'//
     .      ' ASA; potential functions P%?#n<4#%-1j%i#(z)#',' ',
     .      80,stdo,hord)
        endif
      endif

C ... 2nd gen LMTO-specific initialization
C     Make pti here since pp in gamma rep
      if (lcgf .eq. 0 .and. lfp .eq. 0 .and. lgen3 .eq. 0) then

        call defdr(opti,-ndim*nsp)
        call makipt(nl,nbasp,nsp,w(oipc),lihdim,w(oidxsh),w(opp),ndim,
     .    w(opti))
        call pack1('pot opti',spot,opti)

      endif

C ... Transform pp's to alpha representation
      if (lfp .eq. 0 .and. lgen3 .eq. 0) then
C       Get screening alpha's
        call upack('array nclasp opos',sarray,nclasp,obas,0,0,0)
        nlspcp = nl * nsp * nclasp
        ckbas = cksumf(w(obas),3*nbasp)
        ltmp = iostr(8,'STR',nl,nbasp,1,kap2,0,ckbas,-1,nttab,oalph,
     .    oiax,ontab,os)
        i = 0
        if (bittst(lham,1) .or. bittst(lham,2)) i = 10
C       Do the transformation
        call defdr(oo,nlspcp)
        call pptrns(i,nl,w(oipc),nclasp,nsp,w(oalph),1,w(opp),w(oo))
        call rlse(oalph)
      endif

C ... Printout Euler angles
      call upack('ham lncol oeula neula qss',sham,lncol,oeula,neul,qss,
     .  0)
      if (iprint() .ge. 30 .and. lncol .ne. 0) then
        call awrit4('%x          %?#n#Qss%3:1;6d angle %;6d.  #%2j#'//
     .    '%?#n#Euler angles:#noncollinear hamiltonian',outs,80,
     .    -stdo,bittst(lncol,2),qss,qss(4),isw(bittst(lncol,1).and.
     .    iprint().gt.30))
        if (iprint().gt.30 .and. bittst(lncol,1))
     .    call pvioeu(1,ssite,w(oeula),nbasp,neul)
      endif

C ... Make B.sigma for magnetic field
      obs = 1
      if (lfp.eq.0 .and. lgen3.eq.0 .and. bittst(lncol,8)) then
        call upack('ham oeula omagf nbf',sham,oeula,omagf,nbf,0,0)
        call defcc(obs,-4*lihdim)
        call mkbfld(nl,nbasp,lihdim,w(oidxsh),w(omagf),nbf,w(oeula),
     .    neul,w(obs))
      endif
      call pack1('ham obdots',sham,obs)

      end
      subroutine showbs(ssite,sspec,nkaph,iprma,iprmb,ldham)
C- Makes hamiltonian offsets and permutation indices
C ----------------------------------------------------------------------
Ci Inputs
Ci   nkaph :number of types of one l-quantum number in the basis
Ci   iprma :if supplied, iprma holds a permutation table of site
Ci         :indices.  Orbital of the hamiltonian are ordered by
Ci         :sites iprma(1),iprma(2),iprma(3),...
Ci         :iprma(1)=0=> iprma is not supplied.  Then the orbitals
Ci         :are ordered by sites 1,2,3,...
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:   spid
Ci     Stored:
Co Outputs
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      double precision sspec(1),ssite(1)
      integer nkaph,iprma(*)
      integer iprmb(*),ldham(4)
C Local parameters
      integer nkap0,n0
      parameter (nkap0=3,n0=10)
      integer l,ndim,ipr,nglob,nbas,ldim,off,offs,
     .  specw,fieldw,iorb,offsi,ib,is,norb,
     .  ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      logical lprmib
      character sout*120,spid*8,fmt*20,spdlab*7,lsym*1

C --- Setup ---
      nbas = nglob('nbas')
      ldim = ldham(1)
      if (nkaph .gt. nkap0) call rx('showbs: increase nkap0')
C     stdo = lgunit(1)
C     stdl = lgunit(2)
      lprmib = iprma(1) .gt. 0
      if (lprmib) call rx('showbs not ready for iprma')
      call getpr(ipr)
      if (ipr .lt. 10) return
      spdlab = 'spdfghi'

C     Field width = number of digits each offset uses * 2 + padding
      fieldw = int(dlog(dble(ldim))/dlog(10d0)+1d-14)+1

C     Get longest spec name + site
      specw = 0
      do  ib = 1, nbas
        call upack('site spec',ssite,ib,is,0,0,0)
        call spacks(0,'spec name',sspec,spid,is,is)
        specw = max(specw,len(trim(spid)))
      enddo
      call awrit2('(i4,3x,a%i,i%i)',fmt,len(fmt),0,specw+2,fieldw+1)
      offs = 4+2 + specw+2 + 2*fieldw + 6

      call info0(0,1,0,
     .  ' Orbital positions in hamiltonian, resolved by l:%N'//
     .  ' Site  Spec  Total    By l ...')

C     Loop over sites
      do  ib = 1, nbas
        call upack('site spec',ssite,ib,is,0,0,0)
        call spacks(0,'spec name',sspec,spid,is,is)
        call orbl(ib,0,ldim,iprmb,norb,ltab,ktab,off,offl,ndim)

        write(sout,fmt) ib, spid, offl(1)+1
        call awrit1('%a:%i',sout,len(sout),0,offl(1)+ndim)

C       Loop over orbitals
        do  iorb = 1, norb
          l   = ltab(iorb)
C           = ktab(iorb)  <- not needed for single-kappa hamiltonians
          off = offl(iorb)

          offsi = (iorb-1)*(2*fieldw + 5)
          lsym = spdlab(l+1:l+1)
          call awrit3('%np%i:%i('//lsym//')',sout,len(sout),0,
     .      offs+offsi,offl(iorb)+1,offl(iorb)+2*l+1)
        enddo
        call info0(0,0,0,sout)

      enddo

      end

