Cchp1 cmdvar= integer,iv_p_,oipcp
Cchp1 cmdvar= integer,iv_p_,oipc
Cchp1 cmdvar= real(8),rv_p_,osymgr
Cchp1 cmdvar= real(8),rv_p_,oag
Cchp1 cmdvar= integer,iv_p_,onrcp
Cchp1 cmdvar= integer,iv_p_,onrc
Cchp1 cmdvar= integer,iv_p_,oips
Cchp1 cmdvar= real(8),iv_p_,oistab
Cchp1 cmdvar= integer,iv_p_,oics
Cchp1 cmdvar= real(8),rv_p_,oclabl
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
     .nbas0,nbasp,nsite,nclass,nclasp,nclspp,ngen,ngnmx,nspec,usegen,
     .ldist,nggen,ngmx,stdo,nglob,incli !npadl,npadr,
Cchp1       integer:: oclabl , oics , oips , oistab , onrc , onrcp , opos 
Cchp1      ., oipc , oipcp , oiwk , oag , osymgr , aginv
Cchp1        integer:: oics , oips , oistab , onrc , onrcp , opos , oipc , 
Cchp1      . oipcp , oiwk , oag , osymgr , aginv 
Cchp1         integer:: oips , oistab , onrc , onrcp , opos , oipc , oipcp 
Cchp1      .  , oiwk , oag , osymgr , aginv 
Cchp1          integer:: oips , onrc , onrcp , opos , oipc , oipcp , oiwk , 
Cchp1      .   oag , osymgr , aginv 
Cchp1           integer:: onrc , onrcp , opos , oipc , oipcp , oiwk , oag , osymgr 
Cchp1      .    , aginv 
Cchp1            integer:: onrcp , opos , oipc , oipcp , oiwk , oag , osymgr , 
Cchp1      .     aginv 
Cchp1             integer:: opos , oipc , oipcp , oiwk , oag , osymgr , aginv 
Cchp1              integer:: opos , oipc , oipcp , oiwk , osymgr , aginv 
Cchp1               integer:: opos , oipc , oipcp , oiwk , aginv 
Cchp1                integer:: opos , oipcp , oiwk , aginv 
                integer:: opos , oiwk , aginv 
                integer,pointer :: iv_p_oipcp(:)

               integer,pointer :: iv_p_oipc(:)

              real(8),pointer :: rv_p_osymgr(:)

             real(8),pointer :: rv_p_oag(:)

            integer,pointer :: iv_p_onrcp(:)

           integer,pointer :: iv_p_onrc(:)

          integer,pointer :: iv_p_oips(:)

         integer,pointer :: iv_p_oistab(:)

        integer,pointer :: iv_p_oics(:)

       real(8),pointer :: rv_p_oclabl(:)

      integer ,allocatable :: nrspc_iv(:)
      real(8) ,allocatable :: pos2_rv(:)
      integer ,allocatable :: ips2_iv(:)

      integer,allocatable:: iv_a_tmp(:)

      parameter (T=.true., F=.false., ngnmx=10)
      character*120 gens,strn*72
      double precision gen(9,ngnmx),plat(3,3),dist(3,3),qlat(3,3),xx,
     .fptol


      nbas=sctrl%nbas
      nspec=sctrl%nspec
      nsite=sctrl%nsite


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
Cchp1       call defrr(oag,   3*ngmx)
       allocate(rv_p_oag(abs(3*ngmx)))
       if (3*ngmx<0) rv_p_oag(:)=0.0d0

Cchp1       call defrr(osymgr,9*ngmx)
       allocate(rv_p_osymgr(abs(9*ngmx)))
       if (9*ngmx<0) rv_p_osymgr(:)=0.0d0

Cchp1       call defi(oipc,  nsite)
       allocate(iv_p_oipc(abs(nsite)))
       if (nsite<0) iv_p_oipc(:)=0

Cchp1       call defi(oics,  nsite)
       allocate(iv_p_oics(abs(nsite)))
       if (nsite<0) iv_p_oics(:)=0

      allocate(nrspc_iv(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) nrspc_iv(:)=0


      opos=sarray%opos
Cchp1       oips=sarray%oips
       iv_p_oips => sarray%iv_p_oips 


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
     .  .not. a2bin(strn,fptol,4,0,' ',j1,len(strn))) fptol = 1d-5
      else
        fptol = 0
      endif

C ... When generating the group the basis may become enlarged ...
C     copy larger files relevant arrays to larger workspace
Cchp1       call defi (oistab,(ngmx+1)*nbas)
       allocate(iv_p_oistab(abs((ngmx+1)*nbas)))
       if ((ngmx+1)*nbas<0) iv_p_oistab(:)=0.0d0

C     call defi (oistb2,(ngmx+1)*nbas)
      allocate(ips2_iv(ngmx*nbas))
Ctakao_ZeroClear_NotRequiered        if (ngmx*nbas<0) ips2_iv(:)=0

      allocate(pos2_rv(3*ngmx*nbas))
Ctakao_ZeroClear_NotRequiered        if (3*ngmx*nbas<0) pos2_rv(:)=0.0d0

Cchp1       call icopy ( nbas , w ( oips ) , 1 , ips2_iv , 1 )
       call icopy ( nbas , iv_p_oips , 1 , ips2_iv , 1 ) 


      call dcopy ( 3 * nbas , w ( opos ) , 1 , pos2_rv , 1 )

Cchp1       call gensym ( slabl , gens , usegen , t , f , fptol , f , nbas 
Cchp1      ., nspec , ngmx , plat , plat , ldist , dist , pos2_rv , ips2_iv 
Cchp1      ., nrspc_iv , nsgrp , w ( osymgr ) , w ( oag ) , ngen , gen 
Cchp1      ., ssymgr , nggen , isym , w ( oistab ) )
Cchp1        call gensym ( slabl , gens , usegen , t , f , fptol , f , nbas 
Cchp1      . , nspec , ngmx , plat , plat , ldist , dist , pos2_rv , ips2_iv 
Cchp1      . , nrspc_iv , nsgrp , w ( osymgr ) , w ( oag ) , ngen , gen , 
Cchp1      . ssymgr , nggen , isym , iv_p_oistab ) 
Cchp1         call gensym ( slabl , gens , usegen , t , f , fptol , f , nbas 
Cchp1      .  , nspec , ngmx , plat , plat , ldist , dist , pos2_rv , ips2_iv 
Cchp1      .  , nrspc_iv , nsgrp , w ( osymgr ) , rv_p_oag , ngen , gen , ssymgr 
Cchp1      .  , nggen , isym , iv_p_oistab ) 
         call gensym ( slabl , gens , usegen , t , f , fptol , f , nbas 
     .   , nspec , ngmx , plat , plat , ldist , dist , pos2_rv , ips2_iv 
     .   , nrspc_iv , nsgrp , rv_p_osymgr , rv_p_oag , ngen , gen , ssymgr 
     .   , nggen , isym , iv_p_oistab ) 




      if (nbas .gt. nbas0) call rxs('gensym: the basis was enlarged.',
     .' Check group operations.')
      if (nggen .gt. nsgrp) then
        call info2(1,0,0,' MKSYM (warning): generators create more than'
     .  //' ngmx=%i group ops ...',ngmx,0)
        ngmx = ngmx*16
        if (allocated(pos2_rv)) deallocate(pos2_rv)
        if (allocated(ips2_iv)) deallocate(ips2_iv)
        if (allocated(nrspc_iv)) deallocate(nrspc_iv)
Cchp1         call rlse(oag)
         if (associated(rv_p_oag)) deallocate(rv_p_oag)


        goto 5
      endif
      if (fptol .ne. 0) then
        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbas
          call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .    , i_copy_size , i_spackv + 1 - 1 , pos2_rv )

        enddo

        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbas
          call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
        enddo

      endif

Cki#error ERROR, try to release name= oips2 ,but list=null at linenumber= 219 list= (None)


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
Cchp1         call grpgen(gen(1,ngen),1,w(osymgr),npgrp,ngmx)
         call grpgen ( gen ( 1 , ngen ) , 1 , rv_p_osymgr , npgrp , ngmx 
     .   ) 

        call poppr
        incli = npgrp-nsgrp
      endif
C     Not necessary, because it will be remade anyway
C     call redfi(oistab,npgrp*nbas)

C --- Printout of symmetry operations ---
      call info5(20,0,0,' MKSYM:  found %i space group operations'//
     .'%?#n#; adding inversion generated %i ops#%j#'//
     .'%?#(n==0)# ... includes inversion##',nsgrp,npgrp-nsgrp,npgrp,
     .incli,0)
      if (incli .eq. -1) then
        call info0(30,0,0,'         no attempt to add inversion symmetry')
      endif
      if (mod(mode/10,10) .eq. 0) goto 100

C --- Split species into classes ---
Cchp1       call defi(onrc, nspec)
       allocate(iv_p_onrc(abs(nspec)))
       if (nspec<0) iv_p_onrc(:)=0

Cchp1       call icopy(nsite,w(oips),1,w(oipc),1)
Cchp1        call icopy ( nsite , iv_p_oips , 1 , w ( oipc ) , 1 ) 
        call icopy ( nsite , iv_p_oips , 1 , iv_p_oipc , 1 ) 



Cchp1       call splcls(mod(mode/10,10).eq.4,w(opos),nbas,nsgrp,w(oistab),
Cchp1      .nspec,slabl,nclass,w(oipc),w(oics),w(onrc))
Cchp1        call splcls ( mod ( mode / 10 , 10 ) .eq.4 , w ( opos ) , nbas 
Cchp1      . , nsgrp , w ( oistab ) , nspec , slabl , nclass , w ( oipc ) 
Cchp1      . , iv_p_oics , w ( onrc ) ) 
Cchp1         call splcls ( mod ( mode / 10 , 10 ) .eq.4 , w ( opos ) , nbas 
Cchp1      .  , nsgrp , iv_p_oistab , nspec , slabl , nclass , w ( oipc ) , 
Cchp1      .  iv_p_oics , w ( onrc ) ) 
Cchp1          call splcls ( mod ( mode / 10 , 10 ) .eq.4 , w ( opos ) , nbas 
Cchp1      .   , nsgrp , iv_p_oistab , nspec , slabl , nclass , w ( oipc ) , 
Cchp1      .   iv_p_oics , iv_p_onrc ) 
          call splcls ( mod ( mode / 10 , 10 ) .eq.4 , w ( opos ) , nbas 
     .    , nsgrp , iv_p_oistab , nspec , slabl , nclass , iv_p_oipc , 
     .    iv_p_oics , iv_p_onrc ) 




ccccccccccc
      print *,'zzz nclass=',nclass
cccccccccc

c       npadl=sarray%npadl
c       npadr=sarray%npadr

      nclasp = nclass
      nclspp = nclass
      nbasp  = nbas !+ npadl + npadr

C ... (pgf) Extra classes for padding layers
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
Cchp1       call redfi(oics,   nclspp)
       i_data_size=size(iv_p_oics); allocate(iv_a_tmp(i_data_size))
       iv_a_tmp=iv_p_oics; deallocate(iv_p_oics)
       i_data_size=min(i_data_size,nclspp); allocate(iv_p_oics(nclspp))
       iv_p_oics(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)

cki      oipcp = oipc
        iv_p_oipcp => iv_p_oipc
C$$$      endif

C ... Remake istab
      if (mod(mode,10) .eq. 2) then
Cchp1         call defi(oistab,  npgrp*nbas)
         allocate(iv_p_oistab(abs(npgrp*nbas)))
         if (npgrp*nbas<0) iv_p_oistab(:)=0.0d0

      else
Cchp1         call defi(oistab,  nsgrp*nbas)
         allocate(iv_p_oistab(abs(nsgrp*nbas)))
         if (nsgrp*nbas<0) iv_p_oistab(:)=0.0d0

      endif
      call dinv33(plat,1,qlat,xx)
C     call prmx('pos',w(opos),3,3,nbas)
C      call symtab(nbas,w,w(opos),w(osymgr),w(oag),-nsgrp,plat,qlat,w,
C     .  w(oistab))
C      call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
C      call symtbl(0,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
C     .  w(oistab))
C      call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
Cchp1       call symtbl(1,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
Cchp1      .w(oistab))
Cchp1        call symtbl ( 1 , fptol , nbas , w , w ( opos ) , w ( osymgr 
Cchp1      . ) , w ( oag ) , nsgrp , qlat , iv_p_oistab ) 
Cchp1         call symtbl ( 1 , fptol , nbas , w , w ( opos ) , w ( osymgr 
Cchp1      .  ) , rv_p_oag , nsgrp , qlat , iv_p_oistab ) 
         call symtbl ( 1 , fptol , nbas , w , w ( opos ) , rv_p_osymgr 
     .   , rv_p_oag , nsgrp , qlat , iv_p_oistab ) 



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
Cchp1         if (aginv(w(osymgr),nsgrp,npgrp,nbas,w(oag),w(oistab)) .ne. 0)
Cchp1      .  call rx('bug in mksym')
Cchp1          if ( aginv ( w ( osymgr ) , nsgrp , npgrp , nbas , w ( oag ) 
Cchp1      .   , iv_p_oistab ) .ne. 0 ) call rx ( 'bug in mksym' ) 
Cchp1           if ( aginv ( w ( osymgr ) , nsgrp , npgrp , nbas , rv_p_oag , 
Cchp1      .    iv_p_oistab ) .ne. 0 ) call rx ( 'bug in mksym' ) 
           if ( aginv ( rv_p_osymgr , nsgrp , npgrp , nbas , rv_p_oag , 
     .     iv_p_oistab ) .ne. 0 ) call rx ( 'bug in mksym' ) 



      endif

C ... poke ipc into ssite
      i_copy_size=1;
      do i_spackv=1,nbasp
Cchp1         call spackv_array_copy_i8_i('p',ssite(i_spackv)%class,i_copy_size,i_spackv+1-1,w(oipc))
         call spackv_array_copy_i8_i ( 'p' , ssite ( i_spackv ) %class 
     .   , i_copy_size , i_spackv + 1 - 1 , iv_p_oipc ) 

      enddo

C --- Create nrc and class labels ---
      if (mod(mode/10,10) .eq. 2 .or. mod(mode/10,10) .eq. 4) then
Cchp1         call defrr(oclabl, nclspp)
         allocate(rv_p_oclabl(abs(nclspp)))
         if (nclspp<0) rv_p_oclabl(:)=0.0d0

Cchp1         call defi (onrcp,  nclspp)
         allocate(iv_p_onrcp(abs(nclspp)))
         if (nclspp<0) iv_p_onrcp(:)=0

Cchp1         call pvsym2(3,2*nbasp-nbas,nclspp,w(oics),w(oipcp),nspec,slabl,
Cchp1      .  ssite,w(oclabl),w(onrcp))
Cchp1          call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , w ( oics ) , w 
Cchp1      .   ( oipcp ) , nspec , slabl , ssite , rv_p_oclabl , w ( onrcp ) 
Cchp1      .   ) 
Cchp1           call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , iv_p_oics , w ( 
Cchp1      .    oipcp ) , nspec , slabl , ssite , rv_p_oclabl , w ( onrcp ) ) 
Cchp1            call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , iv_p_oics , w ( 
Cchp1      .     oipcp ) , nspec , slabl , ssite , rv_p_oclabl , iv_p_onrcp ) 
            call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , iv_p_oics , iv_p_oipcp 
     .      , nspec , slabl , ssite , rv_p_oclabl , iv_p_onrcp ) 




C$$$     if (nclasp .eq. nclass) then !we only treat nclasp=nclass case. takao Apr2010
cki        onrc = onrcp
        iv_p_onrc => iv_p_onrcp
C$$$     else
C$$$          call defi (onrc,   nclasp)
C$$$          call icopy(nclasp,w(onrcp),1,w(onrc),1)
C$$$          call iaxpy(nclasp-nclass,1,w(onrcp+nclasp),1,w(onrc+nclass),1)
C$$$C          call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
C$$$C          call awrit2('%n:1i',' ',100,6,nclspp,w(onrcp))
C$$$     endif
      else
cki        oclabl = 0
        nullify(rv_p_oclabl)
cki        onrcp = 0
        nullify(iv_p_onrcp)
      endif

C       call awrit2('%n:1i',' ',100,6,nsite,w(oipc))
C       call awrit2('%n:1i',' ',100,6,nsite,w(oipcp))
C       pause

C --- Cleanup: poke class and symmetry info into structures ---

      sctrl%nclass=nclass


Cchp1       sarray%onrcp=onrcp
       sarray%iv_p_onrcp => iv_p_onrcp 

Cchp1       sarray%oipcp=oipcp
       sarray%iv_p_oipcp => iv_p_oipcp 



Cchp1       sarray%oics=oics
       sarray%iv_p_oics => iv_p_oics 

Cchp1       sarray%oipc=oipc
       sarray%iv_p_oipc => iv_p_oipc 

Cchp1       sarray%oclabl=oclabl
       sarray%rv_p_oclabl => rv_p_oclabl 

      sarray%nclasp=nclasp
Cchp1       sarray%onrc=onrc
       sarray%iv_p_onrc => iv_p_onrc 


c     call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
  100 continue

      slat%npgrp=npgrp
      slat%nsgrp=nsgrp
Cchp1       slat%oistab=oistab
       slat%iv_p_oistab => iv_p_oistab 

Cchp1       slat%oag=oag
       slat%rv_p_oag => rv_p_oag 

Cchp1       slat%osymgr=osymgr
       slat%rv_p_osymgr => rv_p_osymgr 



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
     .dclabl,nrc)

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
        return
   10   continue
      enddo
      aginv = 0
      end function aginv


