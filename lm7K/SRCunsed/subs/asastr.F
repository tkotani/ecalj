Cgetarg...info...           structure ['asastr', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['asastr', 'array', 'sarray']
Cgetarg...info...           structure ['asastr', 'lat', 'slat']
Cgetarg...info...           structure ['asastr', 'str', 'sstr']
Cgetarg...info...           use_to_add ['asastr', 'm_struc_def']
      subroutine asastr(prgnam,sctrl,slat,sspec,sstr,sarray,sham)
       
       use m_struc_def  !Cgetarg

C- Screened, real-space structure constants
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of main program
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nspec nl nspin lpgf
Ci     Stored:
Ci     Passed to: chkstr
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy alat plat avw
Ci     Stored:
Ci     Passed to: chkstr
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxb hcr
Ci     Duplicate: lmxb
Ci     Stored:
Ci     Passed to: hcr2a maadot chkstr
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci     Elts read: nkaps kaps rmax rfit amode lmaxw drwats iinv loka
Ci                mxnbr
Ci     Stored:
Ci     Passed to: lgors chkstr
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: ohave oics oipc oclabl oips opos nclasp npadl npadr
Ci                opgfsl
Ci     Stored:
Ci     Passed to: chkstr
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read:
Ci     Stored:
Ci     Passed to:
Co Outputs
Co   strx are written to disk
Cl Local variables
Cl   amode :determines kind of strux to make:
Cl         :  0 2nd generation kappa=0
Cl         :  1 SSSW
Cl         :  2 NMTO
Cl   rmaxs :sphere radius for structure constants, a.u.
Cl   nkapn :Number of families of functions to make
Cl         (1 for now, nkap for NMTO)
Cl   ncplx :1 if strux are real; 2 if complex
Cr Remarks
Cu Updates
Cu   04 Jul 08 (S. Lozovoi) val-lap structure constants
Cu   06 Aug 06 New 2-kappa strux
Cu   19 May 04 make strux work with Methfessel conventions for H,J
Cu   31 Aug 00 some adaptation to fit with NMTO
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) prgnam*8
Cgetarg       double precision sarray(1),sctrl(1),slat(1),sspec(1),sstr(1),
Cgetarg      .  sham(1)
       type(s_array)::sarray
       type(s_ctrl)::sctrl
       type(s_lat)::slat
       type(s_spec)::sspec(*)
       type(s_str)::sstr
       type(s_ham)::sham

C ... Local parameters
      character*120 outs
      logical cmdopt,lgors,strdif,ltmp,lsdot
      integer oadot,oalpha,ocg,oclabl,ocy,oiax,oics,oidxcg,oinitc,oipc,
     .  oips,ojcg,olmx,ontab,opgfsl,opos,os,osd,ohcr,otral,otrad,orwats,
     .  o1,o2,ontabg,okmx,orsmh,orsma,osg,oiax1,ontab1,ontabg1,osg1,
     .  oehvl
      integer amode,awrite,fopn,fopno,i,igets,iprint,isw,
     .  itral,j,nglob,lmaxw,lpbc,lpgf,mordrn,mxcsiz,n0,nbas,
     .  nbasp,nbaspp,nclasp,niax,nitab,nkaps,nl,nla,npadl,npadr,
     .  nttabg,liog,ivl,
     .  nsp,nspec,nttab,stdo,streqv,loka,nkapn,nds,lio,ncplx,nmto
      parameter (niax=10,n0=10)
      double precision avw,kap2(n0),rmaxs,rfit,plat(3,3),alat,
     .  drwats,xx,ckbas,cksumf,siinv(5),rmxg0,rmxg
C ... Heap
      integer w(1)
      common /w/ w

C --- Unpack global variables ---
      stdo = nglob('stdo')
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg       call upack('ctrl nbas nspec nl nspin',sctrl,nbas,nspec,nl,nsp,0)
       
       nbas=sctrl%nbas
       nspec=sctrl%nspec
       nl=sctrl%nl
       nsp=sctrl%nspin

Cgetarg       lpgf  = igets('ctrl lpgf',sctrl)
       lpgf = int(sctrl%lpgf(1)) 

C     lgen3 = igets('ctrl lgen3',sctrl)
C     lpbc  = 0 for pbc in 3 dimensions, 11 pgf padded geometry
      lpbc  = 0
      if (lpgf .gt. 0) lpbc = 11
Cgetarg...info...           integer :: oinitc , integer(8) :: sarray%ohave 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg...info...           integer :: oipc , integer(8) :: sarray%oipc 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: oips , integer(8) :: sarray%oips 1
Cgetarg       call upack('array ohave oics oipc oclabl oips',sarray,
Cgetarg      .  oinitc,oics,oipc,oclabl,oips)
       
       oinitc=sarray%ohave
       oics=sarray%oics
       oipc=sarray%oipc
       oclabl=sarray%oclabl
       oips=sarray%oips

Cgetarg...info...           integer :: opos , integer(8) :: sarray%opos 1
Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: npadl , integer(8) :: sarray%npadl 1
Cgetarg...info...           integer :: npadr , integer(8) :: sarray%npadr 1
Cgetarg       call upack('array opos nclasp npadl npadr',sarray,opos,
Cgetarg      .  nclasp,npadl,npadr,0)
       
       opos=sarray%opos
       nclasp=sarray%nclasp
       npadl=sarray%npadl
       npadr=sarray%npadr

Cgetarg...info...           integer :: opgfsl , integer(8) :: sarray%opgfsl 1
Cgetarg       call upack('array opgfsl',sarray,opgfsl,0,0,0,0)
       
       opgfsl=sarray%opgfsl

      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas
      mordrn = 0
      ncplx = 1
C ... Clebsch-Gordan and related coefficents
Cgetarg...info...           integer :: ocg , integer(8) :: slat%ocg 1
Cgetarg...info...           integer :: ojcg , integer(8) :: slat%ojcg 1
Cgetarg...info...           integer :: oidxcg , integer(8) :: slat%oidxcg 1
Cgetarg...info...           integer :: ocy , integer(8) :: slat%ocy 1
Cgetarg       call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
       
       ocg=slat%ocg
       ojcg=slat%ojcg
       oidxcg=slat%oidxcg
       ocy=slat%ocy

C ... Global strux data
Cgetarg...info...           integer :: nkaps , integer(8) :: sstr%nkaps 1
Cgetarg...info...           real(8) :: kap2(n0) , real(8) :: sstr%kaps 6
Cgetarg...info...           real(8) :: rmaxs , real(8) :: sstr%rmax 1
Cgetarg...info...           real(8) :: rfit , real(8) :: sstr%rfit 1
Cgetarg...info...           integer :: amode , integer(8) :: sstr%amode 1
Cgetarg       call upack('str nkaps kaps rmax rfit amode',sstr,nkaps,kap2,rmaxs,
Cgetarg      .  rfit,amode)
       
       nkaps=sstr%nkaps
       i_copy_size=size(sstr%kaps) 
       call dcopy(i_copy_size,sstr%kaps,1,kap2,1) 
       rmaxs=sstr%rmax
       rfit=sstr%rfit
       amode=sstr%amode

Cgetarg...info...           integer :: lmaxw , integer(8) :: sstr%lmaxw 1
Cgetarg...info...           real(8) :: drwats , real(8) :: sstr%drwats 1
Cgetarg...info...           real(8) :: siinv(5) , real(8) :: sstr%iinv 5
Cgetarg...info...           integer :: loka , integer(8) :: sstr%loka 1
Cgetarg       call upack('str lmaxw drwats iinv loka',sstr,lmaxw,drwats,siinv,
Cgetarg      .  loka,0)
       
       lmaxw=sstr%lmaxw
       drwats=sstr%drwats
       i_copy_size=size(sstr%iinv) 
       call dcopy(i_copy_size,sstr%iinv,1,siinv,1) 
       loka=sstr%loka

Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: avw , real(8) :: slat%avw 1
Cgetarg       call upack('lat alat plat avw',slat,alat,plat,avw,0,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       avw=slat%avw


C --- Parameters determining kind and number of strux ---
      call isanrg(mod(amode,10),0,3,'ASASTR:','strux mode',.true.)
C     3rd generation uses Tank's tral conventions
      if (mod(amode,10) .eq. 2) then
        itral = 4
        lsdot = .true.
        nkapn = nkaps
        nkaps = 1
        nmto = 1
      elseif (mod(amode,10) .eq. 1 .or. mod(amode,10) .eq. 3) then
        itral = 0
        lsdot = .false.
        nkapn = 1
        nmto = 0
        if (mod(amode,10) .eq. 3) then
Cgetarg...info...           integer :: ivl , integer(8) :: sstr%ivl 1
Cgetarg...info...           real(8) :: rmxg0 , real(8) :: sstr%rmaxg 1
Cgetarg           call upack('str ivl rmaxg',sstr,ivl,rmxg0,0,0,0)
           
           ivl=sstr%ivl
           rmxg0=sstr%rmaxg

          call isanrg(ivl,0,2,'ASASTR:','vl mode',.true.)
          call fsanrg(rmxg0,0d0,1d0,0d0,'ASASTR:','rmaxg/rmax',.true.)
c            print *,' asastr: ivl, rmxg0 =',ivl, rmxg0  
        endif
C     2rd generation only for kappa=0
      else
        itral = 0
        nkaps = 1
        nkapn = 1
        kap2(1) = 0
        lsdot = .true.
        nmto = 0
C  ...  Find hcr in terms of alpha, for printout
        call hcr2a(sspec)
      endif

C --- Compare structure constants ---
      if (cmdopt('--chk',5,0,outs)) then
       
        print *, prgnam // 'comparing strx STR and STR1 ...'
        i = 16+8+4+1
        ltmp = strdif('STR','STR1')
        return
      endif

C --- Printout ---
      if (iprint() .ge. 10 .and. mod(amode,10) .eq. 2) then
        call awrit4('%N ASASTR: make S'//
     .    '%?#n#, Sdot## for'//
     .    '%?#n>1#%-1j %i energies%j# kap^2=%d#.'//
     .    '  Watson sphere%?#n<0#: none#%-1j lmax=%i#.',
     .    ' ',80,stdo,isw(lsdot),nkapn,kap2,lmaxw)
        i = awrite('%x%9fitral = %i:',outs,80,0,itral,xx,xx,xx,xx,
     .    xx,xx,xx)
        call mktrli(itral,outs(i+3:))
        call awrit0('%a',outs,80,-stdo)
      elseif (iprint() .ge. 10 .and. mod(amode,10) .le. 1) then
        call awrit4('%N ASASTR: make strux'//
     .    '%?#n#, Sdot## for'//
     .    '%?#n>1#%-1j %i energies%j# kap^2=%d#.'//
     .    '%?#n<0##%-1j%N%9fWatson sphere lmax=%i#',
     .    ' ',80,stdo,isw(lsdot),nkaps,kap2,lmaxw)
c to adjust later
      elseif (iprint() .ge. 10 .and. mod(amode,10) .eq. 3) then
        call awrit4('%N ASASTR: make strux'//
     .    '%?#n#, Sdot## for'//
     .    '%?#n>1#%-1j %i energies%j# kap^2=%d#.'//
     .    '%?#n<0##%-1j%N%9fWatson sphere lmax=%i#',
     .    ' ',80,stdo,isw(lsdot),nkaps,kap2,lmaxw)
      endif

      if (mod(amode,10) .gt. 0) then
        call awrit3('%9fEnergies: %n:1d a.u.  avw = %,6;6d',' ',80,stdo,
     .    nkaps*nkapn,kap2,avw)
      else
        call awrit1('         avw = %,6;6d',' ',80,stdo,avw)
      endif

C --- Default rmaxs, if input is not positive ---
      if (rmaxs .le. 0d0) then
        rmaxs = 2.7d0*avw
        call info5(10,0,0,'%9f'//
     .    'Use default rmaxs = %;3d a.u. = %;3d*avw = %;3d*alat',
     .    rmaxs,rmaxs/avw,rmaxs/alat,0,0)
      endif

      if (siinv(1) .ne. 0)
     .  call awrit3('%9fIterative inversion: ncut = %d  nit = %d'//
     .  '  tol = %g',' ',80,stdo,siinv(1),siinv(2),siinv(3))

      if (mod(amode,10) .eq. 2) then
        call isanrg(nkapn,2,5,' ','number of energies NEL',.true.)
      elseif (mod(amode,10) .eq. 1 .or. mod(amode,10) .eq. 3) then
        call isanrg(nkaps,1,2,' ','number of energies NEL',.true.)
      else
        call isanrg(nkaps,1,1,' ','number of energies NEL',.true.)
      endif

C --- Get neighbor table iax for each atom in the cluster ---
      if (lpbc .eq. 0) then
        i = 3
        j = -1
      elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
        i = 2
        j = 1
      else
        call rx('ASASTR: not implemented for lpbc>1')
      endif
C ... Make nttab,ontab,oiax
Cgetarg       mxcsiz = igets('str mxnbr',sstr)
       mxcsiz = int(sstr%mxnbr) 

      call pairs(nbas,nbasp,alat,plat,rmaxs/2,w(opos),
     .  -1,i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
C      call pairs(nbas,nbasp,1d0,plat,rmaxs*avw/alat/2,w(opos),
C     .  -1,i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
C ... Patch iax(6) for the padded basis
      if (nbasp.gt.nbas) call pairp6(nbas,npadl,npadr,w(oiax),w(ontab))
C ... Added pairs for the embedded cluster
      if (mordrn .eq. 1) then
        call rx('pairec commented out')
C        mxnbr = 2*rmaxs**3*nbasp
C        call redfi(oiax, niax*mxnbr)
C        call pairec(nbas,0,0,w(oipc),alat/alat,plate,w(obas),
C     .    w(orham),nttab,w(ontab),w(oiax),w(owk),mxcsiz)
C        call redfi(oiax, niax*nttab)
      endif
C ... Make basis-dimensioning-related parts of iax
      call sp2cls('spec lmxb',sspec,0,1,1,nspec,olmx)
      call mkiaxd(nttab,w(olmx),w(oips),w(oiax))
      call rlse(olmx)

C      print *, '!! for now, assign nds to nl**2'
       nds = nl**2

C --- Screening parameters for all kappa's ---
      i = nds*nbaspp*nkaps**2*nkapn
      call defdr(oalpha,-i)
      call defdr(oadot,-i)
      if (mod(amode,10) .eq. 2 .or. mod(amode,10) .eq. 0) then
        call defdr(otral,-4*nds*nbaspp*nkapn)
        call defdr(otrad,-4*nds*nbaspp*nkapn)
      elseif (mod(amode,10) .eq. 1) then
        call defdr(otral,1)
        call defdr(otrad,1)
      elseif (mod(amode,10) .eq. 3) then
        call defdr(otral,1)
        call defdr(otrad,1)
        call defi(ontabg, nbas)
c       call defi(okmx,   nspec)
c       call defdr(oehvl,  nspec*nl)
c       call defdr(orsm,  nspec*nl)
c       call defdr(ormt,  nspec*nl)
c ... read rsmh
c       do is = 1, nspec
c         call dpzero(orbp,2*n0*n0)
c         call upack('spec orbp lmxb',sspec,is,orbp,lmxb,0,0)
c         call dcopy(n0,orbp(1,1,1),1,rsmh(1,is),1)
c           print *, 'is, lmxb =',is,lmxb
c           print *, 'rsmh =',(rsmh(il,is),il=1,lmxb+1)
c       enddo
      else
        call isanrg(nkaps,1,1,' ','number of energies NEL',.true.)
      endif
      call defdr(orwats,nbaspp)

C ... Some temporary arrays needed to generate alpha
      call sp2cls('spec lmxb',sspec,0,1,1,nspec,olmx)
      call sp2cls('spec hcr',sspec,0,nl,1,nspec,ohcr)
C     Unpack RSMH
      if (mod(amode,10) .eq. 3) then
        call sp2cls('spec orbp',sspec,0,nl,1,nspec,orsmh)
        call sp2cls('spec kmxt',sspec,0,1,1,nspec,okmx)
        call sp2cls('spec rsma',sspec,0,1,1,nspec,orsma)
        if (ivl .ne. 0)  
     .    call sp2cls('spec ehvl',sspec,0,nl,1,nspec,oehvl)
      endif

C ... Second generation: screening parms defined from spec->alpha
      if (mod(amode,10) .eq. 0) then
C       Note: maadot and makalp should be merged
C       Need to extract alp from sspec, and pass to makalp
        call maadot(sspec,nl,nbaspp,amode,w(oips),avw,
     .    w(oalpha),w(oadot),w(otral),w(otrad))

      elseif (mod(amode,10) .eq. 1 .or. mod(amode,10) .eq. 3) then
        call makalp(nl,nds,nbas,nkaps,kap2,w(ohcr),100*loka+1,
     .    w(oips),w(olmx),w,w(oalpha))

        if (mod(amode,10) .eq. 3) then
c ... adjust ntab for Gaussians
        if (nbas .ne. nbasp) 
     .  call rx('ASASTR: padded basis is not implemented in mode 3')

        rmxg = rmxg0*rmaxs
        call pairg(nbas,0,-1,alat,plat,w(opos),rmxg/2,
     .  w(ontab),w(oiax),w(ontabg),nttabg)

        endif

C ... NMTO: screening parms defined from spec->hcr.  Uses OKA now
      elseif (mod(amode,10) .eq. 2) then
        call dscal(nl*nspec,1/avw,w(ohcr),1)
        call dscal(nkapn,avw**2,kap2,1)
        call mktra2(11,loka,nbaspp,w(oips),nl,w(olmx),avw,itral,kap2,
     .    nkapn,w(ohcr),1,w(otral),w(otrad),w(oalpha),w(oadot),w,w)
      endif

C --- Watson sphere radii ---
      if (lmaxw .ge. 0) then
        call mkrwat(alat/avw,w(opos),nbaspp,w(oiax),w(ontab),
     .    w(oips),nl,nttab,plat,w(ohcr),drwats,w(orwats))
C       call awrit2('%n:1d',' ',80,stdo,nbaspp,w(orwats))
      endif

C ... Count the number of inequivalent strx
      nitab = nttab
Cgetarg       if (lgors('str lequiv,1',sstr)) then
       if ( iand(1,int(sstr%lequiv)) .ne.0 ) then 

        nitab = 0
        do  i = 1, nbas
          if (nmto .eq. 0) then
            o1 = oalpha
            o2 = oadot
            ltmp = .true.
            nla = nl**2
          else
            o1 = ohcr
            o2 = 0
            ltmp = .false.
            nla = nl
          endif
          j = streqv(w(ontab),ltmp,nla,w(o1),w(o2),i,0,plat,
     .      w(opos),w(oiax),nitab)
        enddo
      endif

C --- Real-space screened structure constants ---
      i = nds**2*nkaps**2*nitab*nkapn
      call defdr(os,-i)
      if (lsdot) then 
        call defdr(osd,-i)
      else
        osd = 1
      endif

C ... Make the strux
      if (mod(amode,10) .eq. 1 .or. mod(amode,10) .eq. 3) then
        call strrs(nbas,npadl,npadr,alat,plat,w(opos),w(orwats),
     .  nl,
     .  kap2,
     .  nkaps,lmaxw,siinv,
     .  w(ocy),w(ocg),w(oidxcg),w(ojcg),lpgf,lsdot,nitab.ne.nttab,
     .  loka,
     .  w(ontab),w(oiax),w(oalpha),w(oadot),w(os),w(osd))
        if (mod(amode,10) .eq. 3) then
c ... also make the value-laplacian unit basis out of G0 and G1/Hs/Hs-dot  
          call defdr(osg,  4*nds**2*nttabg)
          call strg(ivl,nbas,nl,w(okmx),w(orsma),alat,plat,w(opos),
     .      w(oips),w(orsmh),w(ohcr),w(oehvl),
     .      w(ocy),w(ocg),w(oidxcg),w(ojcg),
     .      w(ontab),w(oiax),w(ontabg),w(osg))
        endif
      else
        call strscr(loka,nbas,npadl,npadr,alat/avw,plat,w(opos),
     .    w(orwats),nl,kap2,nkapn,lmaxw,siinv,w(ocy),w(ocg),w(oidxcg),
     .    w(ojcg),lpbc.eq.11,lsdot,nitab.ne.nttab,w(ontab),w(oiax),
     .    w(oalpha),w(oadot),w(otral),w(otrad),w(os),w(osd))
      endif
C ... Sort iax(7) by iax(1..5)
      call ppair5(1,nbasp,plat,w(opos),-1d0,w(ontab),w(oiax))

C --- Write strux to disk ---
      ckbas = cksumf(w(opos),3*nbasp)
C ... for now
      lio = 1 + 100*(ncplx-1)
      if (mod(amode,10) .eq. 1) lio = lio+20000
      if (mod(amode,10) .eq. 3) then 
        lio = lio+20000
        liog = lio
      endif 
      if (mod(amode,10) .eq. 2) lio = lio+10000
      call iostr(lio,'STR',nl,nbasp,max(nkaps,nkapn),kap2,itral,ckbas,
     .  lmaxw,nitab,oalpha,oiax,ontab,os)
      call fclose(fopno('STR'))
      if (lsdot) then
        call iostr(lio+800,'SDOT',nl,nbasp,max(nkaps,nkapn),kap2,itral,
     .    ckbas,lmaxw,nitab,oadot,oiax,ontab,osd)
      call fclose(fopno('SDOT'))
      else
        call dfclos(fopn('SDOT'))
      endif
      if (mod(amode,10) .eq. 3) then 
        call iosg(liog,'STRG',ivl,nl,nbas,ckbas,nttabg,
     .    oiax,ontab,ontabg,osg)
        call fclose(fopno('STRG'))
        call rlse(osg)
      endif
      call rlse(os)
      call rlse(olmx)

C --- Printout ---
Cgetarg       if (lgors('str lshow,1',sstr)) then
       if ( iand(1,int(sstr%lshow)) .ne.0 ) then 

        call query('V>=30 to display structure constants',-1,0)
        lio = lio - 1 + 8
        i = nkaps
        if (mod(amode,10) .eq. 2) i = nkapn
        if (mod(amode,10) .eq. 3) liog = lio 
        call iostr(lio,'STR',nl,nbasp,i,kap2,itral,ckbas,lmaxw,
     .             nitab,oalpha,oiax,ontab,os)
        call shostr(nds,nitab,nbasp,plat,w(opos),lio,w(oalpha),w(oiax),
     .             w(ontab),w(os),nkaps,nkapn,kap2,1,1d0)
        if (lsdot) then
        call query('V>=40 to display sdot',-1,0)
        call iostr(lio+800,'SDOT',nl,nbasp,i,kap2,itral,ckbas,lmaxw,
     .             nitab,oalpha,oiax,ontab,os)
        call shostr(nds,nitab,nbasp,plat,w(opos),lio+800,w(oalpha),
     .    w(oiax),w(ontab),w(os),nkaps,nkapn,kap2,1,1d0)
        endif
        if (mod(amode,10) .eq. 3) then 
          call iosg(liog,'STRG',ivl,nl,nbas,ckbas,nttabg,
     .      oiax1,ontab1,ontabg1,osg1)
          call fclose(fopno('STRG'))
          call shosg(ivl,nl*nl,nttabg,nbas,plat,w(opos),w(oiax1),
     .             w(ontab1),w(ontabg1),w(osg1),1,0.1d0)
        endif
      endif

C --- Plot screened orbital structure constants ---
      if (cmdopt('--plot',6,0,outs)) then
        call chkstr(sctrl,slat,sspec,sstr,sarray)
      elseif (cmdopt('--pltg',6,0,outs)) then
        call chksg(sctrl,slat,sspec,sstr,sarray)
      endif

      call rx0(prgnam)
Cgetarg       end
       end subroutine asastr 


c...deb
C       subroutine xxn(s,nds)
C       integer nds,s(nds),i
C       print *, (s(i),i=1,nds)
C       end
C       subroutine xxr(s,nds)
C       integer nds,i
C       double precision s(nds)
C       print *, (s(i),i=1,nds)
C       end
c...deb

