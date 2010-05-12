Cdelw1 var ['oag', '0', 'oag', '3*ngmx', 'real(8)', 'rv_w_', 'mksym']
Cdelw1 var ['osymgr', '0', 'osymgr', '9*ngmx', 'real(8)', 'rv_w_', 'mksym']
Cdelw1 var ['oipc', '0', 'oipc', 'nsite', 'integer', 'iv_w_', 'mksym']
Cdelw1 var ['oics', '0', 'oics', 'nsite', 'integer', 'iv_w_', 'mksym']
Cdelw1 var ['onrspc', '0', 'onrspc', 'nbas', 'integer', 'iv_w_', 'mksym']
Cdelw1 var ['oistab', '0', 'oistab', '(ngmx+1)*nbas', 'integer', 'iv_w_', 'mksym']
Cdelw1 var ['oips2', '0', 'oips2', 'ngmx*nbas', 'integer', 'iv_w_', 'mksym']
Cdelw1 var ['opos2', '0', 'opos2', '3*ngmx*nbas', 'real(8)', 'rv_w_', 'mksym']
Cdelw1 var ['oips', '-1', 'undef', '0', 'undef', 'x', 'mksym']
Cdelw1 var ['opos', '-1', 'undef', '0', 'undef', 'x', 'mksym']
Cdelw1 var ['onrc', '0', 'onrc', 'nspec', 'integer', 'iv_w_', 'mksym']
Cdelw1 var ['oclabl', '1', 'oclabl', 'nclspp', 'real(8)', 'rv_w_', 'mksym']
Cdelw1 var ['onrcp', '1', 'onrcp', 'nclspp', 'integer', 'iv_w_', 'mksym']
Cdelw1 var ['oipcp', '-1', 'undef', '0', 'undef', 'x', 'mksym']
Cdelw1 undel [['undef', 'wref'], 'oips', 'mksym']
Cdelw1 undel [['undef', 'wref'], 'opos', 'mksym']
Cdelw1 undel [['undef', 'wref'], 'oipcp', 'mksym']
Cdelw1 undel [['wref'], 'oipc', 'mksym']
Cdelw1 undel [['wref'], 'onrc', 'mksym']
Cdelw1 undel [['wref'], 'onrcp', 'mksym']
Cdelw1 undel [['wref'], 'oclabl', 'mksym']
Cdelw1 undel [['wref'], 'oics', 'mksym']
Cdelw1 undel [['wref'], 'oistab', 'mksym']
Cdelw1 undel [['wref'], 'oag', 'mksym']
Cdelw1 undel [['wref'], 'osymgr', 'mksym']
Cgetarg...info...           structure ['mksym', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['mksym', 'lat', 'slat']
Cgetarg...info...           structure ['mksym', 'array', 'sarray']
Cgetarg...info...           structure ['mksym', 'site', 'ssite']
Cgetarg...info...           structure ['pvsym2', 'site', 'ssite']
Cgetarg...info...           use_to_add ['mksym', 'm_struc_def']
Cgetarg...info...           use_to_add ['pvsym2', 'm_struc_def']
      subroutine mksym(mode,slabl,ssymgr,sctrl,slat,ssite,sarray)

      use m_struc_def  !Cgetarg

C- Setup for symmetry group
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : 1s digit
Ci           0  make space group only
Ci           1  add inversion to point group (see aginv below)
Ci           2  Same as 1, but make additionally ag,istab for extra
Ci              operations, using -g for rotation part; see Remarks
Ci           10s digit
Ci           0  do nothing about classes
Ci           1  split species into classes
Ci           2  Also assign class labels to each class
Ci           4  Assume class and species are the same.
Ci   slabl : species labels
Ci   ssymgr: string containing symmetry group generators.
Ci           if ssymgr contains 'find', mksym will add basis atoms as
Ci           needed to guarantee generators are valid, and generate
Ci           internally any additonal group operations needed to
Ci           complete the space group.
Ci
Cio Inputs/outputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nspec nsite lpgf
Ci     Stored:    nclass
Ci     Passed to: *
Cio  slat  :struct containing information about the lattice
Cio    Elts read: plat ldist dist
Cio    Stored:    npgrp nsgrp oistab oag osymgr
Cio          :nsgrp= number of space group operations.
Cio          :npgrp= number of point group operations.  It may differ
Cio                  from nsgrp if inversion is added.
Cio          :istab= table of site permutations for each group op
Cio          :ag   = translation part the group ops
Cio          :symgr= point group (rotation) part of each group op
Cio  ssite :struct containing site-specific information
Cio    Elts read: *
Ci     Stored:    class clabel
Cio  sarray:structure containing offsets to various arrays
Cio    Elts read: opos oips npadl npadr
Cio    Stored:    onrcp oipcp oics oipc oclabl nclasp onrc
Co           oipc,oipcp,onrc,oics and oclabl, and
Co           oipc:  pointer to class table, pad equiv to double pad
Co           oipcp: pointer to class table, pad nequiv to double pad
Co Outputs
Co  ssymgr:  may be enlarged if it contains 'find' initially
Co  slat     the following are set: npgrp nsgrp oistab oag osymgr
Cr Remarks
Cr   In certain cases the inversion operation may be added to the space
Cr   group, for purposes of k integration.  This is permissible when the
Cr   hamiltonian has the form h(-k) = h*(k).  In that case, the
Cr   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
Cr
Cr   Also, the Green's functions are related G(-k) = Gtranspose(k).
Cr   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
Cr   then G(-k) = Gtranspose(k), and the same (g,ag) information is
Cr   needed for either rotation.
Cu Updates
Cu   04 Jan 06 Enabled keyword 'point' in ssymgr
Cu             Generates extra group ops when generators create more
Cu             than the maximum expected
Cu   07 Oct 05 Put altered positions into site structure
Cu   03 Nov 01 Adapted to altered gensym
Cu   26 Jan 01 Added new '2' option to 1s digit of mode
Cu   31 May 00 revised meaning of mode
C ----------------------------------------------------------------------
C     implicit none
      character*8 slabl(1),ssymgr*(*)
      integer mode,nsgrp,npgrp
Cgetarg       double precision sctrl(1),slat(1),ssite(1),sarray(1)
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_site)::ssite(*)
      type(s_array)::sarray

C ... Heap
      integer w(1)
      common /w/ w
C Local variables
      logical T,F,cmdopt,a2bin,ltmp
      integer idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nbas,
     .  nbas0,nbasp,nsite,nclass,nclasp,nclspp,ngen,ngnmx,nspec,usegen,
     .  ldist,nggen,ngmx,stdo,nglob,incli !npadl,npadr,
Cdelw1 do not change  [['wref'], 'oclabl', 'mksym']
Cdelw1 do not change  [['wref'], 'oics', 'mksym']
Cdelw1 do not change  [['undef', 'wref'], 'oips', 'mksym']
Cdelw1 do not change  [['wref'], 'oistab', 'mksym']
Cdelw1 do not change  [['wref'], 'onrc', 'mksym']
Cdelw1 do not change  [['wref'], 'onrcp', 'mksym']
Cdelw1 do not change  [['undef', 'wref'], 'opos', 'mksym']
Cdelw1 do not change  [['wref'], 'oipc', 'mksym']
Cdelw1 do not change  [['undef', 'wref'], 'oipcp', 'mksym']
Cdelw1 do not change  [['wref'], 'oag', 'mksym']
Cdelw1 do not change  [['wref'], 'osymgr', 'mksym']
Cdelw1       integer oclabl,oics,oips,oistab,onrc,onrcp,
Cdelw1      .  onrspc,opos,oipc,oipcp,oiwk,oag,osymgr,opos2,oips2,aginv
      integer:: oclabl , oics , oips , oistab , onrc , onrcp , opos 
     . , oipc , oipcp , oiwk , oag , osymgr , aginv
      integer ,allocatable :: iv_w_onrspc(:)
      real(8) ,allocatable :: rv_w_opos2(:)
      integer ,allocatable :: iv_w_oips2(:)

      parameter (T=.true., F=.false., ngnmx=10)
      character*120 gens,strn*72
      double precision gen(9,ngnmx),plat(3,3),dist(3,3),qlat(3,3),xx,
     .  fptol

Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg...info...           integer :: nsite , integer(8) :: sctrl%nsite 1
Cgetarg       call upack('ctrl nbas nspec nsite',sctrl,nbas,nspec,nsite,0,0)

      nbas=sctrl%nbas
      nspec=sctrl%nspec
      nsite=sctrl%nsite

Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           integer :: ldist , integer(8) :: slat%ldist 1
Cgetarg...info...           real(8) :: dist(3,3) , real(8) :: slat%dist 9
Cgetarg       call upack('lat plat ldist dist',slat,plat,ldist,dist,0,0)

      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      ldist=slat%ldist
      i_copy_size=size(slat%dist)
      call dcopy(i_copy_size,slat%dist,1,dist,1)

      ngmx = 48
      stdo = nglob('stdo')
      if (iprint() .ge. 30) write(stdo,'(1x)')

C ... Re-entry when ngmx was increased
    5 continue

C ... For memory allocation of class data, assume initially nclass=nbas
Cdelw1 do not change oag because of ['wref']
      call defrr(oag,   3*ngmx)
Cdelw1 do not change osymgr because of ['wref']
      call defrr(osymgr,9*ngmx)
Cdelw1 do not change oipc because of ['wref']
      call defi(oipc,  nsite)
Cdelw1 do not change oics because of ['wref']
      call defi(oics,  nsite)
Cdelw1       call defi(onrspc,nbas)
      allocate(iv_w_onrspc(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_onrspc(:)=0

Cgetarg...info...           integer :: opos , integer(8) :: sarray%opos 1
Cgetarg...info...           integer :: oips , integer(8) :: sarray%oips 1
Cgetarg       call upack('array opos oips',sarray,opos,oips,0,0,0)

Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 138
Cdelw1 do not change opos because of ['undef', 'wref']
      opos=sarray%opos
Cdelw1 warning(1) , probably  oips  is not defined yet at linenumber= 139
Cdelw1 do not change oips because of ['undef', 'wref']
      oips=sarray%oips

C     call defi(oips,nbas)
C     call defrr(opos,3*nbas)
C     call spackv(10,'site spec',ssite,1,nbas,w(oips))
C     call spackv(10,'site pos',ssite,1,nbas,w(opos))

*     call prmx('basis',w(opos),3,3,nbas)
*     print *, (w(oips-1+i), i=1,nbas)

C ... Extract keywords from input string; copy rest to gens
C     keyword 'point' is appended to gens
      call words(ssymgr,ngen)
      j1 = 1
      idest = 1
      usegen = 2
      gens = ' '
      ltmp = .false.
      do  ig = 1, ngen
        call word(ssymgr,ig,j1,j2)
        if (ssymgr(j1:j2) .eq. 'find') then
          usegen = 0
        else if (ssymgr(j1:j2) .eq. 'point') then
          ltmp = .true.
        else
          call strncp(gens,ssymgr,idest,j1,j2-j1+2)
          idest = idest+j2-j1+2
        endif
      enddo
      if (ltmp) then
        call strncp(gens,'point',idest,1,5)
      endif

C --- Generate space group ---
      nbas0 = nbas
      if (cmdopt('--fixpos',8,0,strn)) then
        j1 = 8+1
        if (strn(9:13) .eq. ':tol=') then
          j1 = 13
        endif
        if (strn(9:9) .ne. ':' .or.
     .    .not. a2bin(strn,fptol,4,0,' ',j1,len(strn))) fptol = 1d-5
      else
        fptol = 0
      endif

C ... When generating the group the basis may become enlarged ...
C     copy larger files relevant arrays to larger workspace
Cdelw1 do not change oistab because of ['wref']
      call defi (oistab,(ngmx+1)*nbas)
C     call defi (oistb2,(ngmx+1)*nbas)
Cdelw1       call defi (oips2, ngmx*nbas)
      allocate(iv_w_oips2(ngmx*nbas))
Ctakao_ZeroClear_NotRequiered        if (ngmx*nbas<0) iv_w_oips2(:)=0

Cdelw1       call defdr(opos2, 3*ngmx*nbas)
      allocate(rv_w_opos2(3*ngmx*nbas))
Ctakao_ZeroClear_NotRequiered        if (3*ngmx*nbas<0) rv_w_opos2(:)=0.0d0

Cdelw1 warning(1) , probably  oips  is not defined yet at linenumber= 191
Cdelw1 do not change oips because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oips  is not defined yet at linenumber= 191
Cdelw1 do not change oips because of ['undef', 'wref']
Cdelw1       call icopy(nbas,w(oips),1,w(oips2),1)
      call icopy ( nbas , w ( oips ) , 1 , iv_w_oips2 , 1 )

Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 192
Cdelw1 do not change opos because of ['undef', 'wref']
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 192
Cdelw1 do not change opos because of ['undef', 'wref']
Cdelw1       call dcopy(3*nbas,w(opos),1,w(opos2),1)
      call dcopy ( 3 * nbas , w ( opos ) , 1 , rv_w_opos2 , 1 )

Cdelw1 do not change oag because of ['wref']
Cdelw1 do not change osymgr because of ['wref']
Cdelw1 do not change oistab because of ['wref']
Cdelw1 do not change oag because of ['wref']
Cdelw1 do not change osymgr because of ['wref']
Cdelw1 do not change oistab because of ['wref']
Cdelw1       call gensym(slabl,gens,usegen,T,F,fptol,F,nbas,nspec,ngmx,plat,
Cdelw1      .  plat,ldist,dist,w(opos2),w(oips2),w(onrspc),nsgrp,w(osymgr),
Cdelw1      .  w(oag),ngen,gen,ssymgr,nggen,isym,w(oistab))
      call gensym ( slabl , gens , usegen , t , f , fptol , f , nbas 
     . , nspec , ngmx , plat , plat , ldist , dist , rv_w_opos2 , iv_w_oips2 
     . , iv_w_onrspc , nsgrp , w ( osymgr ) , w ( oag ) , ngen , gen 
     . , ssymgr , nggen , isym , w ( oistab ) )

      if (nbas .gt. nbas0) call rxs('gensym: the basis was enlarged.',
     .  ' Check group operations.')
      if (nggen .gt. nsgrp) then
        call info2(1,0,0,' MKSYM (warning): generators create more than'
     .    //' ngmx=%i group ops ...',ngmx,0)
        ngmx = ngmx*16
Cdelw1 rlse name= oag old_list= oag osymgr oipc oics onrspc oistab oips2 opos2
Cdelw1 rlse name= oag new_list= (None)
Cdelw1 not deallocate oistab because of [wref]
Cdelw1 not deallocate oics because of [wref]
Cdelw1 not deallocate oipc because of [wref]
Cdelw1 not deallocate osymgr because of [wref]
Cdelw1 not deallocate oag because of [wref]
Cdelw1         call rlse(oag)
        if (allocated(rv_w_opos2)) deallocate(rv_w_opos2)
        if (allocated(iv_w_oips2)) deallocate(iv_w_oips2)
        if (allocated(iv_w_onrspc)) deallocate(iv_w_onrspc)
        call rlse(oag)

        goto 5
      endif
      if (fptol .ne. 0) then
Cgetarg         call spackv(11,'site pos',ssite,1,nbas,w(opos2))
        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbas
Cdelw1 warning(1) , probably  opos2  is not defined yet at linenumber= 209
Cdelw1          call spackv_array_copy_r8_r8('p',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos2))
          call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .    , i_copy_size , i_spackv + 1 - 1 , rv_w_opos2 )

        enddo

Cgetarg         call spackv(10,'site pos',ssite,1,nbas,w(opos))
        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbas
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 215
Cdelw1 do not change opos because of ['undef', 'wref']
          call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
        enddo

      endif
Cdelw1 rlse name= oips2 old_list= (None)

Cki#error ERROR, try to release name= oips2 ,but list=null at linenumber= 219 list= (None)
Cdelw1       call rlse(oips2)


C --- Add inversion to point group ---
      incli = -1
      npgrp = nsgrp
      if (mod(mode,10) .ne. 0) then
        ngen = ngen+1
        call dpzero(gen(1,ngen),9)
        gen(1,ngen) = -1
        gen(5,ngen) = -1
        gen(9,ngen) = -1
        call pshpr(iprint()-40)
Cdelw1 warning(1) , probably  osymgr  is not defined yet at linenumber= 231
Cdelw1 do not change osymgr because of ['wref']
        call grpgen(gen(1,ngen),1,w(osymgr),npgrp,ngmx)
        call poppr
        incli = npgrp-nsgrp
      endif
C     Not necessary, because it will be remade anyway
C     call redfi(oistab,npgrp*nbas)

C --- Printout of symmetry operations ---
      call info5(20,0,0,' MKSYM:  found %i space group operations'//
     .  '%?#n#; adding inversion generated %i ops#%j#'//
     .  '%?#(n==0)# ... includes inversion##',nsgrp,npgrp-nsgrp,npgrp,
     .  incli,0)
      if (incli .eq. -1) then
        call info0(30,0,0,'         no attempt to add inversion symmetry')
      endif
      if (mod(mode/10,10) .eq. 0) goto 100

C --- Split species into classes ---
Cdelw1 do not change onrc because of ['wref']
      call defi(onrc, nspec)
Cdelw1 warning(1) , probably  oipc  is not defined yet at linenumber= 250
Cdelw1 do not change oipc because of ['wref']
Cdelw1 warning(1) , probably  oips  is not defined yet at linenumber= 250
Cdelw1 do not change oips because of ['undef', 'wref']
      call icopy(nsite,w(oips),1,w(oipc),1)

Cdelw1 warning(1) , probably  oipc  is not defined yet at linenumber= 252
Cdelw1 do not change oipc because of ['wref']
Cdelw1 warning(1) , probably  oics  is not defined yet at linenumber= 252
Cdelw1 do not change oics because of ['wref']
Cdelw1 warning(1) , probably  oistab  is not defined yet at linenumber= 252
Cdelw1 do not change oistab because of ['wref']
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 252
Cdelw1 do not change opos because of ['undef', 'wref']
Cdelw1 do not change onrc because of ['wref']
      call splcls(mod(mode/10,10).eq.4,w(opos),nbas,nsgrp,w(oistab),
     .  nspec,slabl,nclass,w(oipc),w(oics),w(onrc))
ccccccccccc
      print *,'zzz nclass=',nclass
cccccccccc
Cgetarg...info...           integer :: npadl , integer(8) :: sarray%npadl 1
Cgetarg...info...           integer :: npadr , integer(8) :: sarray%npadr 1
Cgetarg       call upack('array npadl npadr',sarray,npadl,npadr,0,0,0)

c       npadl=sarray%npadl
c       npadr=sarray%npadr

      nclasp = nclass
      nclspp = nclass
      nbasp  = nbas !+ npadl + npadr

C ... (pgf) Extra classes for padding layers
Cgetarg       lpgf = igets('ctrl lpgf',sctrl)
C$$$       lpgf = int(sctrl%lpgf(1))
C$$$      if (lpgf .gt. 0 .and. lpgf .lt. 10) then
C$$$C   ... Remake ipc,ics for padding and double-padding layers
C$$$        call pgfpp(nbas,npadl,npadr,1,1,.true.,.false.,w(oipc),w(oipc))
C$$$        call redfi(oics,nsite)
C$$$        call defi(oiwk,nsite)
C$$$        call pvsym1(nbas,npadl+npadr,nclasp,w(oipc),w(oics))
C$$$        nclspp = nclasp
C$$$        call pvsym1(nbasp,npadl+npadr,nclspp,w(oipc),w(oics))
C$$$        if (iprint() .ge. 30) then
C$$$          call awrit2('         %i classes for embedding region + %i '//
C$$$     .      'for padding layers',' ',80,stdo,nclass,nclasp-nclass)
C$$$        endif
C$$$        call redfi(oics,   nclspp)
C$$$C   ... Separate ipc, ipcd (dble padded ic eqiv =sngl pad vs not equiv)
C$$$        call defi (oipcp,  nsite)
C$$$        call icopy(nsite,w(oipc),1,w(oipcp),1)
C$$$        call icopy(npadl+npadr,w(oipc+nbas),1,w(oipc+nbasp),1)
C$$$      else
C   ... Reallocate arrays as permanent arrays
Cdelw1 do not change oics because of ['wref']
      call redfi(oics,   nclspp)
Cdelw1 warning(1) , probably  oipc  is not defined yet at linenumber= 291
Cdelw1 do not change oipc because of ['wref']
Cdelw1 warning(1) , probably  oipcp  is not defined yet at linenumber= 291
Cdelw1 do not change oipcp because of ['undef', 'wref']
      oipcp = oipc
C$$$      endif

C ... Remake istab
      if (mod(mode,10) .eq. 2) then
Cdelw1 do not change oistab because of ['wref']
        call defi(oistab,  npgrp*nbas)
      else
Cdelw1 do not change oistab because of ['wref']
        call defi(oistab,  nsgrp*nbas)
      endif
      call dinv33(plat,1,qlat,xx)
C     call prmx('pos',w(opos),3,3,nbas)
C      call symtab(nbas,w,w(opos),w(osymgr),w(oag),-nsgrp,plat,qlat,w,
C     .  w(oistab))
C      call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
C      call symtbl(0,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
C     .  w(oistab))
C      call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 308
Cdelw1 do not change oag because of ['wref']
Cdelw1 warning(1) , probably  osymgr  is not defined yet at linenumber= 308
Cdelw1 do not change osymgr because of ['wref']
Cdelw1 do not change oistab because of ['wref']
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 308
Cdelw1 do not change opos because of ['undef', 'wref']
      call symtbl(1,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
     .  w(oistab))
C     Debugging check
C     call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
C     call symtbl(0,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
C     .   w(oistb2))
C     call shoist(w(oistb2),nbas,w(oag),w(osymgr),nsgrp)
C     call istbpm(w(oistab),nbas,nsgrp,w(oistb2))
C     call shoist(w(oistb2),nbas,w(oag),w(osymgr),nsgrp)
C     stop

C   ... For each i*grpop fix translation ag to be that of grpop
      if (mod(mode,10) .eq. 2) then
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 321
Cdelw1 do not change oag because of ['wref']
Cdelw1 warning(1) , probably  osymgr  is not defined yet at linenumber= 321
Cdelw1 do not change osymgr because of ['wref']
Cdelw1 do not change oistab because of ['wref']
        if (aginv(w(osymgr),nsgrp,npgrp,nbas,w(oag),w(oistab)) .ne. 0)
     .      call rx('bug in mksym')
      endif

C ... poke ipc into ssite
Cgetarg       call spackv(11,'site class',ssite,1,nbasp,w(oipc))
      i_copy_size=1;
      do i_spackv=1,nbasp
Cdelw1 warning(1) , probably  oipc  is not defined yet at linenumber= 329
Cdelw1 do not change oipc because of ['wref']
        call spackv_array_copy_i8_i('p',ssite(i_spackv)%class,i_copy_size,i_spackv+1-1,w(oipc))
      enddo

C --- Create nrc and class labels ---
      if (mod(mode/10,10) .eq. 2 .or. mod(mode/10,10) .eq. 4) then
Cdelw1 do not change oclabl because of ['wref']
        call defrr(oclabl, nclspp)
Cdelw1 do not change onrcp because of ['wref']
        call defi (onrcp,  nclspp)
Cdelw1 do not change oics because of ['wref']
Cdelw1 do not change oclabl because of ['wref']
Cdelw1 do not change onrcp because of ['wref']
Cdelw1 warning(1) , probably  oipcp  is not defined yet at linenumber= 336
Cdelw1 do not change oipcp because of ['undef', 'wref']
        call pvsym2(3,2*nbasp-nbas,nclspp,w(oics),w(oipcp),nspec,slabl,
     .    ssite,w(oclabl),w(onrcp))
C$$$     if (nclasp .eq. nclass) then !we only treat nclasp=nclass case. takao Apr2010
Cdelw1 do not change onrc because of ['wref']
Cdelw1 do not change onrcp because of ['wref']
        onrc = onrcp
C$$$     else
C$$$          call defi (onrc,   nclasp)
C$$$          call icopy(nclasp,w(onrcp),1,w(onrc),1)
C$$$          call iaxpy(nclasp-nclass,1,w(onrcp+nclasp),1,w(onrc+nclass),1)
C$$$C          call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
C$$$C          call awrit2('%n:1i',' ',100,6,nclspp,w(onrcp))
C$$$     endif
      else
Cdelw1 do not change oclabl because of ['wref']
        oclabl = 0
Cdelw1 do not change onrcp because of ['wref']
        onrcp = 0
      endif

C       call awrit2('%n:1i',' ',100,6,nsite,w(oipc))
C       call awrit2('%n:1i',' ',100,6,nsite,w(oipcp))
C       pause

C --- Cleanup: poke class and symmetry info into structures ---
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg       call pack1('ctrl nclass',sctrl,nclass)

      sctrl%nclass=nclass

Cgetarg...info...           integer :: onrcp , integer(8) :: sarray%onrcp 1
Cgetarg...info...           integer :: oipcp , integer(8) :: sarray%oipcp 1
Cgetarg       call pack5('array onrcp oipcp',sarray,onrcp,oipcp,0,0,0)

Cdelw1 do not change onrcp because of ['wref']
      sarray%onrcp=onrcp
Cdelw1 warning(1) , probably  oipcp  is not defined yet at linenumber= 367
Cdelw1 do not change oipcp because of ['undef', 'wref']
      sarray%oipcp=oipcp

Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg...info...           integer :: oipc , integer(8) :: sarray%oipc 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: onrc , integer(8) :: sarray%onrc 1
Cgetarg       call pack5('array oics oipc oclabl nclasp onrc',sarray,oics,oipc,
Cgetarg      .  oclabl,nclasp,onrc)

Cdelw1 do not change oics because of ['wref']
      sarray%oics=oics
Cdelw1 warning(1) , probably  oipc  is not defined yet at linenumber= 378
Cdelw1 do not change oipc because of ['wref']
      sarray%oipc=oipc
Cdelw1 do not change oclabl because of ['wref']
      sarray%oclabl=oclabl
      sarray%nclasp=nclasp
Cdelw1 do not change onrc because of ['wref']
      sarray%onrc=onrc

c     call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
  100 continue
Cgetarg...info...           integer :: npgrp , integer(8) :: slat%npgrp 1
Cgetarg...info...           integer :: nsgrp , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: oistab , integer(8) :: slat%oistab 1
Cgetarg...info...           integer :: oag , integer(8) :: slat%oag 1
Cgetarg...info...           integer :: osymgr , integer(8) :: slat%osymgr 1
Cgetarg       call pack5('lat npgrp nsgrp oistab oag osymgr',slat,npgrp,nsgrp,
Cgetarg      .  oistab,oag,osymgr)

      slat%npgrp=npgrp
      slat%nsgrp=nsgrp
Cdelw1 do not change oistab because of ['wref']
      slat%oistab=oistab
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 396
Cdelw1 do not change oag because of ['wref']
      slat%oag=oag
Cdelw1 warning(1) , probably  osymgr  is not defined yet at linenumber= 397
Cdelw1 do not change osymgr because of ['wref']
      slat%osymgr=osymgr


Cgetarg       end
Cdelw1 w_varlist remains: onrc oics oistab oclabl onrcp
Cdelw1 w_varlistundel: [oips opos oipcp oipc onrc onrcp oclabl oics oistab oag osymgr]
Cdelw1 w_varlist (undel), remains: [onrc oics oistab oclabl onrcp]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate onrcp because of [wref]
Cdelw1 not deallocate oclabl because of [wref]
Cdelw1 not deallocate oistab because of [wref]
Cdelw1 not deallocate oics because of [wref]
Cdelw1 not deallocate onrc because of [wref]
      end subroutine mksym

C$$$      subroutine pvsym1(nbas,npad,nclass,ipc,ics)

C$$$C- Fix classes for padding layers
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   nbas   number of basis atoms without padding layers
C$$$Ci   nclass number of classes without padding layers
C$$$Ci Outputs
C$$$Co   nclass updated
C$$$Co   ipc(nbas..nbas+npad) updated
C$$$Co   ics updated
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$      integer nbas,nclass,ipc(nbas),ics(nclass)
C$$$      integer i,npad,nc,is,ic,j
C$$$      if (npad .eq. 0) return

C$$$C ... ipc->-ipc to flag which classes have been renumbered
C$$$      do  5  i = nbas+1, nbas+npad
C$$$    5 ipc(i) = -ipc(i)

C$$$C ... For each site in padding layer, assign class
C$$$      nc = nclass
C$$$C ... Look for next occurrence of class not assigned
C$$$      do  10  i = nbas+1, nbas+npad

C$$$        if (ipc(i) .gt. 0) goto 10
C$$$C   ... Class ipc(i) need to be assigned:  do for each site in class
C$$$C       ic is (-)original class, is is species of original class
C$$$        ic = ipc(i)
C$$$        is = ics(-ic)

C$$$        nc = nc+1
C$$$        do  12  j = i, nbas+npad
C$$$   12   if (ipc(j) .eq. ic) ipc(j) = nc

C$$$C   ... Species index is same as that of original class
C$$$        ics(nc) = is

C$$$   10 continue

C$$$      nclass = nc
C$$$C      call awrit2('%n:1i',' ',100,6,nbas+npad,ipc)
C$$$C      call awrit2('%n:1i',' ',100,6,nclass,ics)
C$$$C      pause

C$$$Cgetarg       end
C$$$       end subroutine pvsym1

      subroutine pvsym2(mode,nbas,nclass,ics,ipc,nspec,slabl,ssite,
     .  dclabl,nrc)

      use m_struc_def  !Cgetarg

C- Create class labels from species labels (double precision format)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nosplt: T copy class and species
Ci     mode: 0 do nothing
Ci           1 create class labels clabl
Ci           2 create number of sites in each class nrc
Ci      ipc: for padding sites ib ipc(ib) = class
Ci   nclass: number of classes
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    clabel
Co Outputs:
Co   dclabl: class labels in double precision format
Co      nrc: number of sites in each class
Cu Updates
Cu   18 Dec 01 Packs class label into ssite->clabel
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nbas,nclass,nspec,ics(1),ipc(nbas),nrc(1)
Cgetarg       double precision dclabl(nclass),ssite(1)
      real(8):: dclabl(nclass)
      type(s_site)::ssite(*)

      character*8 slabl(nspec)
      integer ic,iclbsj,idx,is,ib
      character*8 clabl

C --- Make class labels from species labels ---
      if (mod(mode,2) .eq. 1) then
        do  10  is = 1, nspec
          do  12  idx = 1, nbas
            ic = iclbsj(is,ics,-nclass,idx)
            if (ic .lt. 0) goto 13
            call clabel(slabl,is,idx,clabl)
            call s8tor8(clabl,dclabl(ic))
   12     continue
   13     continue
   10   continue
      endif

      do  20  ib = 1, nbas
        ic = ipc(ib)
        call r8tos8(dclabl(ic),clabl)
Cgetarg         call spacks(1,'site clabel',ssite,clabl,ib,ib)
        do i_spacks=ib,ib
          call spacks_copy('p',ssite(i_spacks)%clabel,ib,ib,clabl,i_spacks)
        enddo

   20 continue

C      do  ib = 1, nbas
C        call spacks(0,'site clabel',ssite,clabl,ib,ib)
C        print *, ib, clabl
C      enddo
C      pause

C --- Create nrc ---
      if (mod(mode/2,2) .eq. 1) then
        call iinit(nrc,nclass)
        do  30  ib = 1, nbas
          ic = ipc(ib)
          nrc(ic) = nrc(ic)+1
   30   continue
      endif

C ... debugging : check class labels
C      call awrit2('%n:1i',' ',100,6,nclass,nrc)
C      do  50  ic = 1, nclass
C        call r8tos8(dclabl(ic),clabl)
C        print 333, ic, clabl
C  333   format(i4,1x,a)
C   50 continue
C      pause

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine pvsym2


      integer function aginv(g,npg,nsg,nbas,ag,istab)

C- Finds translation part for extra inversion-generated group ops
C ----------------------------------------------------------------------
Ci Inputs
Ci   g     :point group operations
Ci   npg   :number of point group operations
Ci   nsg   :number of space group operations
ci   nbas  :dimensions istab.  If zero, istab is not altered.
Co Outputs
Co   ag    :for each group op ig in range nsg+1..npg,
Co         :the corresponding group op jg is found such that
Co         :g(jg) = -g(ig), and ag(jg) is copied to ag(ig); see Remarks
Co   istab :istab(jg) is copied into istab(ig) just in the same way as
Co         :ag(jg) is copied into ag(ig).
Co         :The sign of istab(1) is set negative to indicate that this
Co         :operation is derived from a true space group operation.
Cr Remarks
Cr   In certain cases the inversion operation may be added to the space
Cr   group, for purposes of k integration.  This is permissible when the
Cr   hamiltonian has the form h(-k) = h*(k).  In that case, the
Cr   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
Cr
Cr   Also, the Green's functions are related G(-k) = Gtranspose(k).
Cr   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
Cr   then G(-k) = Gtranspose(k), and the same (g,ag) information is
Cr   needed for either rotation.  Let g~ = -g:
Cr     G0(g~^-1 -k) into G(k),
Cr
Cr
Cr   Routine aborts with aginv=ig if for group ig it cannot find
Cr   corresponding jg such that g(jg) = -g(ig)
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npg,nsg,nbas,istab(nbas,1)
      double precision g(3,3,npg),ag(3,npg)
C ... Local parameters
      integer ig,jg
      double precision ginv(3,3)
      logical grpeql

      do  ig = npg+1, nsg
        call dpcopy(g(1,1,ig),ginv,1,9,-1d0)
        do  jg = 1, nsg
          if (grpeql(g(1,1,jg),ginv)) then
            call dcopy(3,ag(1,jg),1,ag(1,ig),1)
            if (nbas .gt. 0) then
              call icopy(nbas,istab(1,jg),1,istab(1,ig),1)
              istab(1,ig) = -istab(1,jg)
            endif
C           call awrit2(' aginv mapping %,2i ->%,2i',' ',80,6,ig,jg)
            goto 10
          endif
        enddo
        aginv = ig
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
   10   continue
      enddo
      aginv = 0
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end function aginv


