Cgetarg...info...           structure ['mcasim', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['mcasim', 'lat', 'slat']
Cgetarg...info...           structure ['mcasim', 'array', 'sarray']
Cgetarg...info...           structure ['mcasim', 'spec', 'sspec']
Cgetarg...info...           structure ['mcasm2', 'move', 'smove']
Cgetarg...info...           structure ['mcmet', 'move', 'smove']
Cgetarg...info...           use_to_add ['mcmet', 'm_struc_def']
Cgetarg...info...           use_to_add ['mcasim', 'm_struc_def']
Cgetarg...info...           use_to_add ['mcasm2', 'm_struc_def']
      subroutine mcasim(sctrl,sarray,sspec,slat,smove,sstrn)
       
       use m_struc_def  !Cgetarg

C- Monte-Carlo simulations for alloy configurations
C ----------------------------------------------------------------------
Ci   sctrl :struct containing parameters governing program flow
Ci   sarray:structure containing offsets to various arrays
Ci   sspec :struct containing species-specific information
Ci   slat  :struct containing information about the lattice
Ci   sstrn :struct containing global strings
Co Outputs
Cr Remarks
Cr   This is a special-purpose for monte-carlo simulations of an
Cr   alloy.
Cr
Cr   Check for special constraints added directly into the code in mcmet!
Cr   For example, paper with Alan Wright,
Cr        if (ips(ib) .gt. 2 .or. ips(jb) .gt. 2) goto 12
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mxcsiz
      character*(*) sstrn
Cgetarg       double precision sctrl(1),sspec(1),slat(1),sarray(1),smove(1)
       type(s_ctrl)::sctrl
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_array)::sarray
       type(s_move)::smove

C Local
      logical cmdopt
      integer i,i1,i2,igets,iosiga,j,lpbc,ltype(2),nbas,nbasp,
     .  nspec,nttab,ntype(2),oepair,oesite,oiax,oipair,oiprm,oiprmi,
     .  oips,oisite,ontab,opgfsl,opos,orham,ortab,osigb,osigl,osigma
      character htype(2)*5, outs*80
      double precision alat,plat(3,3),avw
C ... Heap
      integer w(1)
      common /w/ w

      data htype /'h1:','h2:'/

Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nbasp , integer(8) :: sctrl%nbasp 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg       call upack('ctrl nbas nbasp nspec',sctrl,nbas,nbasp,nspec,0,0)
       
       nbas=sctrl%nbas
       nbasp=sctrl%nbasp
       nspec=sctrl%nspec

Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: avw , real(8) :: slat%avw 1
Cgetarg       call upack('lat alat plat avw',slat,alat,plat,avw,0,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       avw=slat%avw

Cgetarg...info...           integer :: opgfsl , integer(8) :: sarray%opgfsl 1
Cgetarg...info...           integer :: oips , integer(8) :: sarray%oips 1
Cgetarg...info...           integer :: opos , integer(8) :: sarray%opos 1
Cgetarg       call upack('array opgfsl oips opos',sarray,opgfsl,oips,opos,0,0)
       
       opgfsl=sarray%opgfsl
       oips=sarray%oips
       opos=sarray%opos

      call defrr(orham,nspec)
Cgetarg       call spackv(10,'spec rham',sspec,1,nspec,w(orham))
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_r8_r8('u',sspec(i_spackv)%rham,i_copy_size,i_spackv+1-1,w(orham))
       enddo

Cgetarg       lpgf = igets('ctrl lpgf',sctrl)
c       lpgf = int(sctrl%lpgf(1)) 

      lpbc = 0
c      if (lpgf .gt. 0) lpbc = 11

C --- Make iax and rtab ---
      if (lpbc .eq. 0) then
        i = 3
        j = -1
      elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
        i = 2
        j = 1
      else
        call rx('MMHAM: not implemented for lpbc>1')
      endif
      mxcsiz = 0
      call pairs(nbas,nbasp,alat,plat,w(orham),w(opos),
     .  w(oips),i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
      call defrr(ortab,3*nttab)
      call mkrtab(000,alat,plat,w(opos),w(oiax),nttab,w(opos),w(ortab))

C --- Make the Hamiltonian ---
      call defrr(oesite, nbasp)
      call defi(oisite,  nbasp)
      call defrr(oepair, nttab)
      call defi(oipair,  nttab)
      call upacks('strn mmham',i1,i2)
      ntype(1) = 1
      ntype(2) = 1
      ltype(1) = 1
      ltype(2) = 1
      call suemph(sspec,sstrn(i1:i2),w(ontab),w(oiax),w(ortab),nbasp,
     .  nttab,w(oips),ntype,ltype,htype,w(oesite),w(oisite),w(oepair),
     .  w(oipair))

C --- sigmas, permutation table from file; if missing create and exit---
      call defi(oiprm, nbasp)
      call defi(oiprmi,nbasp)
      call defi(osigb,  nbasp)
      call defi(osigl,  nbasp)
      call defi(osigma,  nbasp)
      if  (iosiga(0,'sigma',nbasp,w(oips),w(osigl),w(osigb),w(oiprm),
     .  w(opos)).lt. 0) then
C       Initial permutation table
        do  10  i = 1, nbasp
   10   w(oiprm+i-1) = i
C       Both sigb and sigl have same convention for now
        call mcasms(0,1,nbas,w(oips),w(osigb))
        call mcasms(0,1,nbas,w(oips),w(osigl))
        i = iosiga(1,'sigma',nbasp,w(oips),w(osigl),w(osigb),w(oiprm),
     .    w(opos))
        call rx0(' wrote new sigma file')
      endif

      call mcasm2(smove,w(ontab),w(oiax),nbasp,nttab,w(oips),w(oesite),
     .  w(oepair),w(osigma),w(oiprm),w(oiprmi),w(osigb),w(osigl))

      if (.not. cmdopt('-wsigma=',8,0,outs)) return
      i = iosiga(1,outs(9:),nbasp,w(oips),w(osigl),w(osigb),w(oiprm),
     .  w(opos))

Cgetarg       end
       end subroutine mcasim 

      subroutine mcasm2(smove,ntab,iax,nbas,nttab,ips,
     .  esite,epair,sigma,iprm,iprmi,sigb,sigl)
       
       use m_struc_def  !Cgetarg

C- Monte-carlo simulations of an an alloy
C ----------------------------------------------------------------------
Ci Inputs
Ci   smove
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nbas  :size of basis
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   iprm  :site permutation table: {ib} are permuted to {iprm(ib)}
Ci   iprmi :inverse of iprm
Ci   esite :one-center hamiltonian
Ci   epair :pairwise energy
Co Outputs
Co   sig   :p
Cl Local variables
Cl         :
Cr Remarks
Cr One-center Hamiltonian is  sum_i e_i sigl(i) sigb(iprm(i))
Cr Two-center Hamiltonian is  1/2 sum_ib,jb e_ib,jb sigma_ib sigma_jb
Cr Sigma is defined in mcasms.
Cr
Cr To update site list efficiently, hold a permutation table
Cr that maps site ib into iprm(ib).
Cr iax table is NOT updated; indices (1..2,*) refer to UNPERMUTED ib,bj
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
Cgetarg       double precision smove(1)
       type(s_move)::smove

      integer nbas,ntab(1),niax,nttab,ips(1),sigma(1),iprm(1),iprmi(1)
      parameter (niax=10)
      integer iax(niax,1),sigb(nbas),sigl(nbas)
      double precision esite(nbas),epair(nttab)
C Local
      logical cmdopt,a2bin
      integer lgunit,iprint,ntot,nequ,i,ib,iseed
      double precision kt,ts,ttot,tequ
      character*80 outs
C ... heap:
      integer w(1)
      common /w/ w

C --- Setup and initial printout ---
Cgetarg...info...           real(8) :: kt , real(8) :: smove%kt 1
Cgetarg...info...           real(8) :: ts , real(8) :: smove%ts 1
Cgetarg...info...           real(8) :: ttot , real(8) :: smove%tstot 1
Cgetarg...info...           real(8) :: tequ , real(8) :: smove%tsequ 1
Cgetarg       call upack('move kt ts tstot tsequ',smove,kt,ts,ttot,tequ,0)
       
       kt=smove%kt
       ts=smove%ts
       ttot=smove%tstot
       tequ=smove%tsequ

C ... Initial printout
      ntot = ttot/ts
      nequ = tequ/ts
      call awrit5(' mcasim:   %i mc steps (ttot=%;4d, ts=%;4d)'//
     .  ' incl. %i equil. steps,  kT=%;4d',
     .  outs,80,-lgunit(2),ntot,ttot,ts,nequ,kt)
      if (iprint() .ge. 20) then
        print *
        call awrit0('%a',outs,80,-lgunit(1))
      endif
C     Create inverse of iprm and initial sigma
      do  10  i = 1, nbas
        ib = iprm(i)
        iprmi(ib) = i
        sigma(ib) = sigb(ib)*sigl(i)
   10 continue

      i = 7
      if (cmdopt('--seed=',i,0,outs)) then
        if (.not. a2bin(outs,iseed,2,0,' ',i,-1))
     .    call rxs('mcasim failed to parse argument: ',outs)
        call ran1in(iseed)
      endif

      call mcmet(smove,ntab,iax,nbas,nttab,ips,esite,epair,sigma,
     .  iprm,iprmi,sigb,sigl)

Cgetarg       end
       end subroutine mcasm2 


      subroutine mcmet(smove,ntab,iax,nbas,nttab,ips,esite,epair,
     .  sigma,iprm,iprmi,sigb,sigl)
       
       use m_struc_def  !Cgetarg

C- Metropolis algorithm
C     implicit none
C Passed parameters
Cgetarg       double precision smove(1)
       type(s_move)::smove

      integer nbas,ntab(nbas+1),niax,nttab,ips(1),sigma(1),
     .  iprm(nbas),iprmi(nbas)
      parameter (niax=10)
      integer iax(niax,1),sigb(nbas),sigl(nbas)
      double precision esite(nbas),epair(nttab)
C Local
      logical lundo
      integer ntry,i1,ib,i2,jb,nclus,ip,ipr,nswap,nequ,ntot,stdo,lgunit
      real ran1
      double precision kt,etot,e1c,e2c,edc,eold,de,ebar,esqr
      double precision ttot,ts,tequ

      call tcn('mcmet')

      if (ran1() .eq. 0.)
     .  call rx('mcmet: random number generator not seeded')

C     Inititializations
      call getpr(ipr)
Cgetarg...info...           real(8) :: kt , real(8) :: smove%kt 1
Cgetarg...info...           real(8) :: ts , real(8) :: smove%ts 1
Cgetarg...info...           real(8) :: ttot , real(8) :: smove%tstot 1
Cgetarg...info...           real(8) :: tequ , real(8) :: smove%tsequ 1
Cgetarg       call upack('move kt ts tstot tsequ',smove,kt,ts,ttot,tequ,0)
       
       kt=smove%kt
       ts=smove%ts
       ttot=smove%tstot
       tequ=smove%tsequ

      ntot = ttot/ts
      nequ = tequ/ts
      call mcasan(1,nbas,iprm,iprmi,ntab,iax,ips,esite,epair,sigb,sigl,
     .  sigma,eold)
      call mcasho(1,nbas,iprm,iprmi,ntab,iax,esite,epair,sigma,sigb,
     .  eold)
      etot = eold
      ebar = 0
      esqr = 0
      nswap = 0
      stdo = lgunit(1)

C --- For ntot time steps, do ---
      ntry = 0
   10 ntry = ntry+1

C   ... Decide on first of pair of atoms to swap
   12   continue
        ib = nbas*ran1() + 1
C   ... If no sigma, no contribution to hamiltonian ... skip
C       if (sigb(ib) .eq. 0) goto 12
        i1 = iprmi(ib)

C   ... Decide on second of pair of atoms to swap
C       Require that it belong to cluster i1
        nclus = ntab(i1+1)-ntab(i1)
        ip = nclus*ran1() + 1 + ntab(i1)
        i2 = iax(2,ip)
        jb = iprm(i2)
C   ... If no sigma, no contribution to hamiltonian ... skip
C        if (sigb(jb) .eq. 0) goto 12
C   ... If equivalent sigma, no change in  hamiltonian ... skip
        if (sigb(jb) .eq. sigb(ib)) goto 12
C   ... special constraints on movement
        if (ips(ib) .gt. 2 .or. ips(jb) .gt. 2) goto 12

C       Energy connected with the two sites before swap
        call mcasme(ib,ib,ntab,iax,esite,epair,iprm,iprmi,
     .    sigma,sigb,e1c,e2c,edc)
        etot = etot - e1c - (e2c-edc)*2 - edc
        call mcasme(jb,jb,ntab,iax,esite,epair,iprm,iprmi,
     .    sigma,sigb,e1c,e2c,edc)
        etot = etot - e1c - (e2c-edc)*2 - edc

C       Swap the atoms
        call mcasw(ib,jb,iprm,iprmi,sigb,sigl,sigma)

C       Energy connected with the two sites after swap
        call mcasme(ib,ib,ntab,iax,esite,epair,iprm,iprmi,
     .    sigma,sigb,e1c,e2c,edc)
        etot = etot + e1c + (e2c-edc)*2 + edc
        call mcasme(jb,jb,ntab,iax,esite,epair,iprm,iprmi,
     .    sigma,sigb,e1c,e2c,edc)
        etot = etot + e1c + (e2c-edc)*2 + edc

        de = etot - eold

C  ...  Metropolis step : retain swap with probability exp(-de/kt)
        lundo = .false.
        if (de .gt. 0) then
          if (ran1() .gt. dexp(-de/kt)) then
C           Undo atom swap
            call mcasw(jb,ib,iprm,iprmi,sigb,sigl,sigma)
            etot = eold
            lundo = .true.
          endif
        endif
        if (.not. lundo) nswap = nswap+1

C   ... Keep statistics
        if (.not. lundo .and. nswap .gt. nequ) then
          ebar = ebar + etot
          esqr = esqr + etot**2
        endif

#if DEBUG
C       Calculate the energy, brute force
C        call mcasme(1,nbas,ntab,iax,esite,epair,iprm,iprmi,
C     .    sigma,sigb,e1c,e2c,edc)
C        if (abs(etot-(e1c+e2c)).gt.1d-8) call rxi('bug in mcmet',ntry)
C        call awrit7('  mcmet swap i1(ib1)=%i(%i) with i2(ib2)=%i(%i)'//
C     .  ' etot=%;4d de=%;4d undo=%l',' ',80,6,i1,ib,i2,jb,etot,de,lundo)
C        call awrit2('  iprm=%n:1i',' ',80,6,nbas,iprm)
#endif
        if (ipr .le. 50 .and. ipr .gt. 30) then
          if (mod(ntry,1000) .eq. 0)
     .      call awrit4(' mcmet completed %i+%i swaps in %i attempts'
     .      //' e=%d',' ',80,lgunit(1),nswap-nequ,nequ,ntry,eold)
        endif
        if (ipr .gt. 50 .or. ipr .ge. 50 .and. .not. lundo)  then
          print 357, ib,jb,etot,de,lundo
  357     format(' swap',i4,i4,' enow = ',f12.6,' de = ',f12.6,' undo=',
     .      L1)
        endif

        eold = etot

C ... end of big loop
      if (nswap .lt. ntot) goto 10

      ebar = ebar/(nswap-nequ)
      esqr = dsqrt(esqr/(nswap-nequ) - ebar**2)

      if (ipr .ge. 10) then
        do  90  ip = 1, 2
          call awrit3('%N mcmet completed %i+%i swaps in %i attempts:',
     .      ' ',80,lgunit(ip),nswap-nequ,nequ,ntry)
          call awrit3(' kt = %;4d  ebar = %;4d  rms(e-ebar) = %;4d ',
     .      ' ',80,lgunit(ip),kt,ebar,esqr)
   90   continue
      endif

      call mcasan(1,nbas,iprm,iprmi,ntab,iax,ips,esite,epair,sigb,sigl,
     .  sigma,eold)
      call mcasho(1,nbas,iprm,iprmi,ntab,iax,esite,epair,sigma,sigb,
     .  eold)

      call tcx('mcmet')
Cgetarg       end
       end subroutine mcmet 

      subroutine mcasw(ib1,ib2,iprm,iprmi,sigb,sigl,sigma)

C- swap two entries in permutation table
C ----------------------------------------------------------------------
Ci Inputs
Ci  ib1,ib2:pairs of sites to swap
Ci   iprm  :site permutation table: atom ib=iprm(i) is moved to site i
Ci   iprmi :inverse of iprm: atom ib is moved to site iprmi(ib)
Co Outputs
Co   iprm  :entries iprm(i1) and iprm(i2) are swapped, with
Co          ib1 = iprm(i1), ib2=iprm(ib2)
Co   iprmi :entries iprmi(ib1) and iprmi(ib2) are swapped
Co   sigma :sigma(ib1) = sigb(ib1)*sigl(i1) and
Co         :sigma(ib2) = sigb(ib2)*sigl(i2) are updated.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ib1,ib2,iprm(ib2),iprmi(ib2),
     .  sigb(ib2),sigl(ib2),sigma(ib2)
C Local
      integer i1,i2

      i1 = iprmi(ib1)
      i2 = iprmi(ib2)
      iprmi(ib1) = i2
      iprmi(ib2) = i1
#if DEBUG
      call awrit4('  ... swap i1(ib1)=%i(%i) with i2(ib2)=%i(%i)',
     .  ' ',80,6,i1,ib1,i2,ib2)
      if (iprm(i1) .ne. ib1) stop 'bug in mcasw'
      if (iprm(i2) .ne. ib2) stop 'bug in mcasw'
#endif
      iprm(i1) = ib2
      iprm(i2) = ib1
C     sigma(ib1) = sigb(ib1)*sigl(iprmi(ib1))
      sigma(ib1) = sigb(ib1)*sigl(i2)
C     sigma(ib2) = sigb(ib2)*sigl(iprmi(ib2))
      sigma(ib2) = sigb(ib2)*sigl(i1)

Cgetarg       end
       end subroutine mcasw 


      subroutine mcasms(mode,i1,i2,ips,sigma)

C- Generate sigma
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :determines how sigma is made
Ci  i1,i2 :range of sites for which sigma is calculated
Ci   ips   :species table: site  belongs to species ips()
Co Outputs
Co   sigma :sigma is generated for lattice sites i=i1..i2
Cr Remarks
Cr   mode=0:  sigma(i) = 1 if ips(i) = 1, sigma = -1 if ips(i) = 2
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer mode,i1,i2,ips(i2),sigma(i2)
C Local
      integer i

      if (mode .ne. 0) call rxi('mcasms: mode not recognized, ',mode)
      do  10  i = i1, i2
        if (ips(i) .eq. 1) then
          sigma(i) = 1
        elseif (ips(i) .eq. 2) then
          sigma(i) = -1
        else
          sigma(i) = 0
        endif
   10 continue

Cgetarg       end
       end subroutine mcasms 

      subroutine mcasme(ib1,ib2,ntab,iax,esite,epair,iprm,iprmi,
     .  sig1,sig2,e1c,e2c,edc)

C- Hamiltonian total energy
C     implicit none
C Passed parameters
      integer ib1,ib2,sig1(ib2),sig2(ib2),niax
      parameter (niax=10)
      integer iax(niax,1),ntab(1),iprm(ib2),iprmi(ib2)
      double precision esite(ib2),epair(ib2),e1c,e2c,edc
C Local
      integer ib,ibp,ip,jbp,jb

      call tcn('mcasme')

C --- One-center hamiltonian ---
      e1c = 0
      do  10  ib = ib1, ib2
        e1c = e1c + sig1(ib) * esite(ib)
C       print *, 'e1c after ib=',ib,e1c,sig1(ib)*esite(ib)
   10 continue

C --- Pair hamiltonian ---
      e2c = 0
      edc = 0
C     For each site do
      do  20  ib = ib1, ib2
        ibp = iprmi(ib)
C       Loop over each pair in the corresponding pair table
C       e2ci = 0
        do  22  ip = ntab(ibp)+1, ntab(ibp+1)
          jbp = iax(2,ip)
          jb = iprm(jbp)
          e2c = e2c + sig2(ib)*sig2(jb)*(epair(ip)/2)
          if (ib .eq. jb) edc = edc + sig2(ib)*sig2(jb)*(epair(ip)/2)
   22   continue
   20 continue

C          e2ci = e2ci + sig2(ib)*sig2(jb)*(epair(ip)/2)
C          if (sig2(ib)*sig2(jb)*(epair(ip)/2) .ne. 0)
C     .     print 333, ip,iax(6,ip),ib,jb,sig2(ib)*sig2(jb)*epair(ip)/2
C  333     format(4i4,f10.4)
C   22   continue
C        print *, 'e2c after ib=',ib,e2c,e2ci,edc
C   20 continue

      call tcx('mcasme')

Cgetarg       end
       end subroutine mcasme 

      integer function iosiga(lread,fnam,nbas,ips,sigl,sigb,iprm,pos)

Cr Remarks
Cr   Loaded into variables table for each of the ib=1..nbas lines read:
Cr     is = ips(ib)
C     implicit none
      integer lread,nbas,sigl(nbas),sigb(nbas),iprm(nbas),ips(nbas)
      character*(*) fnam
      double precision pos(3,nbas)
      integer ifi,ib,fopna,j1,j2,is,ival
      double precision px,py,pz
      character*1 ch
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),aa*(recl),
     .  vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev),a2bin
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)

      iosiga = 0
      ifi = fopna(fnam,-1,0)
      rewind ifi
      if (lread .eq. 0) then
        read(ifi,'(a)',end=99,err=99) ch
        nr = 0
        do  10  ib = 1, nbas
          call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .                ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
          if (nr .lt. 0) goto 99
C         read(a,*,end=99,err=99) sigl(ib),sigb(ib),iprm(ib),px,py,pz
          is = ips(ib)
          call lodsyv('is',1,dble(is),ival)
          j1 = 0
          if  (.not. a2bin(a,sigl(ib),2,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,sigb(ib),2,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,iprm(ib),2,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,px,4,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,py,4,0,' ',j1,-1)) goto 99
          if (.not. a2bin(a,pz,4,0,' ',j1,-1)) goto 99
          if (abs(px-pos(1,ib))+abs(py-pos(2,ib))+abs(pz-pos(3,ib))
     .      .gt. 1d-6) goto 99
   10   continue
        print *, 'iosiga:  read sigl,sigb,iprm from file ',fnam
      else
        write(ifi,
     .    '(''#   sigl sigb  iprm              pos'',21x,''spec'')')
        do  20  ib = 1, nbas
          write(ifi,333) sigl(ib),sigb(ib),iprm(ib),
     .      pos(1,ib),pos(2,ib),pos(3,ib),ips(iprm(ib))
  333     format(3i6,3f12.7,i6)
   20   continue
        call word(fnam,1,j1,j2)
        print *, 'iosiga:  wrote sigl,sigb,iprm to file ',fnam(j1:j2)
      endif
      call fclose(ifi)

      return
C --- Error exit ---
   99 continue
      print *, 'iosiga: file read failed or file mismatch'
      iosiga = -1
      call fclose(ifi)
Cgetarg       end
       end function iosiga 


      subroutine mcasho(ib1,ib2,iprm,iprmi,ntab,iax,esite,epair,sig1,
     .  sig2,etot)

C- Printout
C ----------------------------------------------------------------------
Ci Inputs
Ci  ib1,ib2:range of sites for which sig1 is calculated
Ci   iprm  :site permutation table: atom ib=iprm(i) is moved to site i
Ci   iprmi :inverse of iprm: atom ib is moved to site iprmi(ib)
Ci   sig1  :sigma for 1-center part of hamiltonian, sig1_i e1_i
Ci   sig2  :sigma for 2-center part of hamiltonian, sig2_i sig2_j e2_ij
Co Outputs
Co   etot  :
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ib1,ib2,iprm(ib2),iprmi(ib2),ntab(ib2),niax
      parameter (niax=10)
      integer sig2(ib2),sig1(ib2),iax(niax,1)
      double precision esite(ib2),epair(1),etot,edc
C Local
      integer i,ncol,nskip
      double precision e1c,e2c

      ncol = min(4,ib2-ib1+1)
      nskip = (ib2-ib1+ncol)/ncol

      call mcasme(ib1,ib2,ntab,iax,esite,epair,iprm,iprmi,
     .  sig1,sig2,e1c,e2c,edc)
      etot = e1c + e2c
      call awrit5('%N  mcasho :  sites %i..%i  e1c=%;4d  epair=%;4d'//
     .  '  etot=%;4d',' ',80,6,ib1,ib2,e1c,e2c,e1c+e2c)
      call arrprt('  ib site  s1  s2','%,4i%,4i%,4i%,4i','Iiii',
     .  ib2-ib1+1,nskip,ncol,nskip,'  | ',i,iprmi,sig1,sig2,i,i,i,i)

Cgetarg       end
       end subroutine mcasho 

      subroutine mcasan(ib1,ib2,iprm,iprmi,ntab,iax,ips,esite,epair,
     .  sigb,sigl,sig1,etot)

C- Sanity checking
C     implicit none
C Passed parameters
      integer ib1,ib2,iprm(ib2),iprmi(ib2),ntab(ib2),niax
      parameter (niax=10)
      integer sigb(ib2),sigl(ib2),sig1(ib2),ips(ib2),iax(niax,1)
      double precision esite(ib2),epair(1),etot
C Local
      integer i,ib


C     Check that sigma tables are sensible
      do  10  i = ib1, ib2
        ib = iprm(i)
        if (iprmi(ib) .ne. i) call rx('bug in mcasim')
        if (sig1(ib) .ne. sigb(ib)*sigl(i)) then
          call rx('bug in mcasim')
        endif
   10 continue

Cgetarg       end
       end subroutine mcasan 


