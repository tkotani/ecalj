      subroutine secmtn(sctrl,sspec,slat,sham,spot,nl,nsp,nbas,ips,
     .  indxsh,qss,eula,neul,ikp,nkp,ldim,lidim,lihdim,qp,ppn,sop,isp,
     .  nevmx,efmax,nev,z,eb)
C- Hamiltonian and Overlap, NMTO
C ----------------------------------------------------------------
Ci Inputs:
Ci   nl,nsp,nbas,pp
Ci   ldim: dimension of l - wave block of hamiltonian matrix
Ci   idim: dimension of i - wave block of hamiltonian matrix
Co   pph : vector of pot. par's. in alpha rep'n (suham.f)
Ci   eula: Euler angles for spin rotations
Ci   nevmx max. no. evecs to generate.
Ci         -1 suppresses generation of z
Ci         -2 Do not diagonalize, but return overlap in z,
Ci            allocate oc for hamiltonian and place there
Ci   z:    used as a work array, whether or not evecs generated
Co Outputs:
Ci   ccd:  diagonal matrices for 1- 2- & 3-centre CCOR integrals
Co   eigenvalues and eigenvectors are returned in eb, z
Co   nev:  number of evecs generated.
Cr Remarks
Cr   Downfolding automatically turns on the combined correction.
Cr   bittst(lham,8) can be used to transform
Cr   structure constants to an arbitrary representation.
Cr   Hybridisation is turned off when bittst(lham,16) is set (see remhyb.f).
Cr   Dimensions of pph,eb,z are doubled when spins are coupled.
Cu Updates
Cu   08 Feb 03 altered dmensions of sop
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nl,nsp,nbas,neul,isp,ldim,lidim,lihdim,nevmx,nev,
     .  ips(nbas),ikp,nkp,indxsh(1)
      integer nppn
      parameter (nppn=12)
      double precision sctrl(1),sspec(1),eb(ldim*2),ppn(nppn,lihdim,*),
     .  qp(3),z(ldim,ldim*2),eula(nbas,neul,3),efmax,qss(4),
     .  sop(0:nl-1,nsp,nsp,9,*)
      double precision slat(1),sham(1),spot(1)
C ... Local parameters
      logical lewald,lgors,lnc,lx,ltmp,iostr
      integer hdim,i,ii,ipr,j,l2,ldimx,lgunit,linv,lncol,nglob,nspec,
     .  nmto,nkaps,n0,nsite,oalph,oadot,oiax,ontab,itral
      integer ohcr,ohk,olmx,ook,osdll,oshl,osll,otrad,otral,owk,os,osd,
     .  lmaxw,lbloch,obas,nl2,idim,nitab,ndoti
      parameter (n0=10)
      double precision vmtz,dgets,avw,dglob,kaps(n0),kmto(n0),ckbas0,
     .  plat(3,3),ckbas,cksumf
      character*80 outs
C ... Heap
      integer w(1)
      common /w/ w

      call upack('ctrl lncol',sctrl,lncol,0,0,0,0)
      call upack('ham nmto kmto',sham,nmto,kmto,0,0,0)

      call getpr(ipr)
      lewald = lgors('ctrl lgen3,2',sctrl)
      nspec = nglob('nspec')
      avw = dglob('avw',0d0,0)
      vmtz = dgets('pot vmtz',spot)
      lnc = lncol .ne. 0
      nl2 = nl**2
      hdim = lihdim - lidim
C     idim = lidim - ldim
      ldimx = ldim
C     ld2 = ldim**2
      if (lnc) ldimx = 2*ldim
      l2 = ldimx**2
      lbloch = 10
      call upack('lat plat opos',slat,plat,obas,0,0,0)
      ckbas = cksumf(w(obas),3*nbas)

C     Space for nmto strux
      call defdc(osll,lidim**2*nmto)
      call defdc(osdll,lidim**2*nmto)
      if (hdim .ne. 0) call defdc(oshl,hdim*lidim)

C --- Screened strux by Ewald summation ---
      if (lewald) then

C   ... Make the tral matrix on the fly
        call defdr(otral,-4*nl**2*nbas*nmto)
        call defdr(otrad,-4*nl**2*nbas*nmto)
C       call defdr(oalpha,-nl**2*nbas*nmto)
C       call defdr(oadot,-nl**2*nbas*nmto)
        call sp2cls('spec lmxb',sspec,0,1,1,nspec,olmx)
        call sp2cls('spec hcr',sspec,0,nl,1,nspec,ohcr)
        call dscal(nmto,avw**2,kmto,1)
        call dscal(nl*nspec,1/avw,w(ohcr),1)
        call pshpr(0)
        if (ikp.eq.1 .and. isp.eq.1 .and. ipr.ge.50) call setpr(50)
        call mktra2(1,1,nbas,ips,nl,w(olmx),avw,4,kmto,nmto,w(ohcr),1,
     .    w(otral),w(otrad),w,w,w,w)
        call poppr
        call dscal(nmto,1/avw**2,kmto,1)
        call dscal(nl*nspec,avw,w(ohcr),1)

        call sstrxq(1,slat,sham,qp,nmto,kmto,w(otral),w(otrad),w(osll),
     .    w(osdll),w(oshl),w)

        call rlse(otral)

C --- Screened strux by Bloch sum of R.S. strux ---
      else
C       Read header data
        ltmp = iostr(10006,'STR',nl,nbas,nkaps,kaps,itral,ckbas0,lmaxw,
     .    nitab,oalph,oiax,ontab,osd)
C       Some sanity checks
        if (ckbas .ne. ckbas0)
     .    call rx('secmtn: basis mismatch in strux file')
        call isanrg(itral,4,4,'secmtn','tral',.true.)
C       Read the strux and energy derivatives
        ltmp = iostr(10008,'SDOT',nl,nbas,nkaps,kaps,itral,ckbas,lmaxw,
     .    nitab,oalph,oiax,ontab,osd)
        ltmp = iostr(10008,'STR',nl,nbas,nkaps,kaps,itral,ckbas,lmaxw,
     .    nitab,oadot,oiax,ontab,os)
        nsite = w(ontab+nbas)
C       Energy-interpolate and Bloch sum strux
        ndoti = min(nkaps,2)
        call dscal(nmto,avw**2,kmto,1)
        call blochi(lbloch,qp,nl,plat,indxsh,1,nsite,w(oiax),nkaps,
     .    ndoti,kaps,nmto,kmto,w(os),w(osd),nl2,nitab,0,lidim,0,lidim,
     .    lidim,lidim,0,w(osll))
        call blochi(10000+lbloch,qp,nl,plat,indxsh,1,nsite,w(oiax),
     .    nkaps,ndoti,kaps,nmto,kmto,w(os),w(osd),nl2,nitab,0,lidim,0,
     .    lidim,lidim,lidim,0,w(osdll))
        call dscal(nmto,1/avw**2,kmto,1)
        call rlse(oalph)
      endif

C --- Kink matrix from structure constants ---
      call dscal(2*lidim**2*nmto,avw**2,w(osdll),1)
      call sstr2k(ldim,lidim,lihdim,nmto,isp,ppn,w(osll),w(osdll))

C --- NMTO Hamiltonian and overlap ---
      call defdc(ohk,-l2)
      call defdc(ook,-l2)
      call defdc(owk, l2)
      call nmham(1,lidim,ldim,nmto,kmto,vmtz,w(osll),w(osdll),w(owk),
     .  w(ohk),w(ook))
      call rlse(owk)

C --- Diagonalize ---
#if BLAS3
      lx = .true.
#else
      lx = .false.
#endif
C     Diagonalize by inverse iteration, or not
      linv = 0
      if (nevmx .gt. 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1
      call defdr(owk,ldimx*11)
      call zhev(ldimx,w(ohk),w(ook),.true.,lx,nevmx,efmax,nev,w(owk),
     .  linv,0,eb,z)
      call rlse(owk)

C --- Printout ---
      if (ipr .ge. 30) then
        j = min(9,ldimx)
        if (ipr .ge. 35) j = ldimx
#if LINUX_PGI
        do  18  ii = 1, 1
#else
        do  18  ii = 1, 2
#endif
        call awrit3(' SECMTN:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(ii),ikp,nkp,qp)
   18   write(lgunit(ii),'(255(9f8.4:/))') (eb(i), i=1,j)
        if (ipr .ge. 36 .and. nev .gt. 0) call awrit5(
     .    ' nev, nevmx, ldim=  %i  %i  %i  ev(nev) = %1;5d  efmax '//
     .    '= %1;5d',' ',80,lgunit(1),nev,nevmx,ldimx,eb(nev),efmax)
        call ftflsh(lgunit(1))
      endif

      if (ipr .ge. 110) then
        outs = 'evec'
        call yprm(outs,2,z,ldimx*nev,ldimx,ldimx,nev)
        call zprm('eigenvectors',2,z,ldimx,ldimx,ldimx)
        call yprm('eval',1,eb,ldimx*1,ldimx,nev,1)
        call query('V<110 to skip matrix printing',-1,0)
      endif

      call rlse(osll)
      end

      subroutine sstr2k(ldim,lidim,lihdim,nmto,isp,ppn,sll,sdll)
C- Overwrite screened structure constants with Kink matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   lidim :number of lower+intermediate orbitals
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   nmto :number of types of one l-quantum number in the basis
Ci   isp   :current spin channel (1 or 2)
Ci   ppn  :nmto potential parameters, in downfolding order
Cio Inputs/Outputs
Cio  sll   :On input, screened structure constants
Cio        :On output, K = aD - S
Cr Remarks
Cr  For downfolding, the KKR-equations have the simple form:
Cr
Cr    K^a *  u = 0   with  K^a = S^a + aD
Cr
Cr and where K^a is the kink matrix in the a-representation.
Cr Divide orbitals into lower and intermediate sets:
Cr Then u_i can be expressed by u_l:
Cr
Cr    u_i = - (K^a_ii)^-1 * K^a_il * u_l
Cr
Cr  and the equation for u_l is:
Cr
Cr    (K^a_ll - K^a_li * (K^a_ii)^-1 * K^a_il) u_l \equiv K^b_ll u_l = 0
Cr
Cr   with beta = inverse potential function P^0.
Cr
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ldim,lidim,lihdim,nppn,nmto,isp
      parameter (nppn=12)
      double precision ppn(nppn,lihdim,nmto,isp)
      double complex sll(lidim,lidim,nmto),sdll(lidim,lidim,nmto)
C ... Local parameters
      integer inc,ik,i,j,idim,ierr,owk,iprint
      double precision ad,add,phia,phida
      double complex alpha,beta
      character*40 outs
C ... Heap
      integer w(1)
      common /w/ w

C --- Downfolded kink matrix and its energy derivative ---
      idim = lidim-ldim
      if (idim .ne. 0) then

      do  ik = 1, nmto

        alpha = 1
        beta = 0

C   ... Scale the SSW ii block to match the back extrapolated phi(a)
        do  i = ldim+1, lidim
          do  j = ldim+1, lidim
            sdll(i,j,ik) = sll(i,j,ik) * ppn(6,j,ik,isp) +
     .                    sdll(i,j,ik) * ppn(5,j,ik,isp)
            sll(i,j,ik) =  sll(i,j,ik) * ppn(5,j,ik,isp)
          enddo
        enddo

C   ... Add diagonal D to make iwaves K, using SSW rescaled to phi0(a)
        do  i = ldim+1, lidim
             aD = ppn(3,i,ik,isp)
            aDd = ppn(4,i,ik,isp)
           phia = ppn(5,i,ik,isp)
          phida = ppn(6,i,ik,isp)

          sll(i,i,ik)  =  sll(i,i,ik) + aD*phia
          sdll(i,i,ik) = sdll(i,i,ik) + aDd*phia + aD*phida
        enddo

C   ... Kii -> Kii^-1; Overwrite Sil with Kii^-1 Sil = Kii^-1 Sli+
        call defrr(owk,idim*max(idim+1,3))
        i = ldim+1
        call zqinv('g',sll(i,i,ik),lidim,0,idim,w(owk),idim,ierr)
        if (ierr .ne. 0) call rx('sstr2k:  Kii is singular')
        call rlse(owk)
        call zgemm('N','C',idim,ldim,idim,alpha,sll(i,i,ik),lidim,
     .    sll(1,i,ik),lidim,beta,sll(i,1,ik),lidim)

C   ... Store in wk : Sdil - Kdii Kii^-1 Sil
        call defcc(owk,idim*ldim)
        call zgemm('N','N',idim,ldim,idim,alpha,sdll(i,i,ik),lidim,
     .    sll(i,1,ik),lidim,beta,w(owk),idim)
        call dscal(idim*ldim*2,-1d0,w(owk),1)
        call zmsadd(idim,ldim,lidim,idim,0,0,0,0,alpha,alpha,
     .    sdll(i,1,ik),w(owk))

C   ... Scale the SSW li block to match the back extrapolated phi(a)
C       Also scale by -1 to prepare for the final multiplication steps
        do  i = ldim+1, lidim
          phia  = ppn(5,i,ik,isp)
          phida = ppn(6,i,ik,isp)
          do  j = 1, ldim
            sdll(j,i,ik) =-sll(j,i,ik) * phida - sdll(j,i,ik) * phia
            sll(j,i,ik) = -sll(j,i,ik) * phia
          enddo
        enddo

C   ... Overwrite Sdil with : Kii^-1 (Sdil - Kdii Kii^-1 Sil)
        i = ldim+1
        beta = 0
        call zgemm('N','N',idim,ldim,idim,alpha,sll(i,i,ik),lidim,
     .    w(owk),idim,beta,sdll(i,1,ik),lidim)
        call rlse(owk)

C   ... Kll(b) = Kll(a) - Kli * (Kii^-1 * Kil)
        i = ldim+1
        beta = 1
C        call zprm('-Kli',2,sll(1,i,ik),lidim,ldim,idim)
C        call zprm('Kii^-1 Kil',2,sll(i,1,ik),lidim,idim,ldim)
        call zgemm('N','N',ldim,ldim,idim,alpha,sll(1,i,ik),lidim,
     .    sll(i,1,ik),lidim,beta,sll(1,1,ik),lidim)

C       call zprm('K^a - Kli Kii^-1 Kil',2,sll(1,1,ik),lidim,ldim,ldim)

C   ... Kdll(b) = Kdll(a) - Kdli * (Kii^-1 Kil)dot - Kdli * (Kii^-1 Kil)
        call zgemm('N','N',ldim,ldim,idim,alpha,sll(1,i,ik),lidim,
     .    sdll(i,1,ik),lidim,beta,sdll(1,1,ik),lidim)
        call zgemm('N','N',ldim,ldim,idim,alpha,sdll(1,i,ik),lidim,
     .    sll(i,1,ik),lidim,beta,sdll(1,1,ik),lidim)

C       call zprm('Kd^b',2,sdll(1,1,ik),lidim,ldim,ldim)

      enddo
      endif

C --- Kink matrix of ll block ---
      inc = 2*(1+lidim)
      do  ik = 1, nmto
        call daxpy(ldim,1d0,ppn(3,1,ik,isp),nppn,sll(1,1,ik),inc)
        call daxpy(ldim,1d0,ppn(4,1,ik,isp),nppn,sdll(1,1,ik),inc)

        if (iprint() .gt. 110) then
          call awrit2('%xK^b%?#n#-dot##, ik=%i',outs,len(outs),0,0,ik)
          call zprm(outs,2,sll(1,1,ik),lidim,ldim,ldim)
          call awrit2('%xK^b%?#n#-dot##, ik=%i',outs,len(outs),0,1,ik)
          call zprm(outs,2,sdll(1,1,ik),lidim,ldim,ldim)
        endif

      enddo


      end

