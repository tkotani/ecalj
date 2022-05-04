      module m_mksym
      use m_mksym_util,only: gensym,grpgen,symtbl
      public :: m_mksym_init,
! ibas ==> iclass=ipc(ibas),iv_a_oipc ==> ispec=ics(iclass),iv_a_oics
     &     rv_a_oag,   iv_a_oics, iv_a_oipc , lat_npgrp, lat_nsgrp,
     &     rv_a_osymgr,iv_a_oistab, ctrl_nclass, iclasstaf_, symops_af_,ag_af_,ngrpaf_
      
      integer , allocatable,protected ::  iv_a_oics (:)
      integer,  allocatable,protected ::  iv_a_oipc(:)
      real(8) , allocatable,protected ::  rv_a_oag (:)
      real(8) , allocatable,protected ::  rv_a_osymgr (:)
      integer , allocatable,protected ::  iv_a_oistab (:)
      integer,allocatable,protected:: iclasstaf_(:)
      real(8),allocatable,protected:: symops_af_(:,:,:), ag_af_(:,:)
      integer,protected:: lat_npgrp, lat_nsgrp, ctrl_nclass, ngrpaf_
      
      private 
      integer , allocatable,protected ::  iv_a_onrc (:)
      contains
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine m_mksym_init(prgnam) !latsym(prgnam)
      use m_lmfinit,only: v_ssite,nbas,sstrnsymg,ctrl_noinv,
     &     symgaf,iv_a_oips,slabl,mxspec,nsite,procid,master
! v_ssite%pos modified
      use m_lattic,only: rv_a_opos,m_lattic_init
!! Driver for calling mksymaf and mksym
!-------------      
!!note: Antiferro symmetric self-energy mode. Feb2021 recovered
!! To obtain self-consistency, it may be useful to keep AF condition during iteration.
!! We need to set SYMGRPAF (still developing...)
! sample for nio
!     SYMGRPAF i:(1,1,1) ! translation + inversion makes antiferro symmetry.
!     SYMGRP r3d   !this keeps spin axis. (probably larger symmetry allowed for SO=0)
!STRUC   ALAT={a} PLAT= 0.5 0.5 1.0  0.5 1.0 0.5  1.0 0.5 0.5
!        NBAS= 4  NSPEC=3
!SITE    ATOM=Niup POS=  .0   .0   .0   AF=1    <--- AF symmetric pair
!        ATOM=Nidn POS= 1.0  1.0  1.0   AF=-1   <---
!        ATOM=O POS=  .5   .5   .5
!        ATOM=O POS= 1.5  1.5  1.5
!------------------------------------------------------------
      implicit none
      character,intent(in)::  prgnam*(*)
      integer:: ibas,lc,j,iprint
      integer,parameter::recln=511
      logical ::cmdopt,ipr10=.false.
      character strn*(recln),strn2*(recln),outs(recln)
!! lattice set up      
      do ibas=1,nbas
        v_ssite(ibas)%pos = rv_a_opos(:,ibas)
      enddo
!! mksym mksymaf generate symmetry operations; split species into classes  ---
!! mksym may correct atomic position slightly (not tested well).
      strn = 'find'
      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
      if (cmdopt('--nosym',7,0,outs).or.cmdopt('--pdos',6,0,outs)) strn = ' '
!! when lmfgw, 1st digit of lc is zero--> no inversion added in mksym.
      lc = 10
      if (.not. prgnam == 'LMFGWD') lc = 12 ! add inversion for (not lmfgw)
      if (prgnam=='LMFA'.or.prgnam=='LMCHK') then
         lc=20
         if(ctrl_noinv /=1 ) lc = lc+1  !+1 means add inversion 
      endif
      if(len_trim(symgaf)>0) lc=12 ! inversion allowed for AF case. Good?
      if(.not.(prgnam=='LMFA'.or.prgnam=='LMCHK')) ipr10= iprint()>10 !this is only for master
      if(len_trim(symgaf)>0) then
        if(ipr10) then
           write(6,*)
           write(6,"(a)")       ' AF: ======================================== '
           write(6,"(a)")       ' AF: Antiferro mode: SPGGRAF='//trim(symgaf)
           write(6,"(a)")       ' AF:  (neglct waring in GENSYM) '
           do j=1,nbas
              write(6,"(a,2i3)") ' AF:  ibas,AF=',j,v_ssite(j)%iantiferro
           enddo
        endif
        strn2=trim(strn)//' '//trim(symgaf)
        call mksymaf(v_ssite,iv_a_oips,nbas,procid==master,strn2,lc,slabl,mxspec,nsite)
        if(ipr10) write(6,"(a)") ' AF: ===== end of AF section================= '
        if(ipr10) write(6,"(a)")       
      endif
      if(procid==master) call pshpr(60)
      call mksym(lc,slabl,strn,v_ssite,iv_a_oips) 
      if(procid==master) call poppr
      end subroutine m_mksym_init
      
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      subroutine mksymaf(v_ssite,iv_a_oips_in,nbas,imaster,strn2,lc,slabl,mxspec,nsite)!all input
!     mksymaf sets data above used in sugw.F. a little complicated...
      use m_struc_def,only: s_site
      intent(in)::       v_ssite,iv_a_oips_in,nbas,imaster,strn2,lc,slabl,mxspec,nsite
      logical::                         imaster
      character strn2*(*)
      character(8) :: slabl(mxspec) 
      type(s_site):: v_ssite(nbas)
      type(s_site),allocatable:: v_ssite2(:) !pos is rewritten in lmfp
      integer:: inumaf,iv_a_oips_in(nsite),j,k,nbas,  lc,       mxspec,nsite,ib
      integer,allocatable::iv_a_oips(:)
      allocate(v_ssite2(nbas),iv_a_oips(nsite))
      v_ssite2 = v_ssite
      iv_a_oips=iv_a_oips_in
      inumaf = 0
      do j=1,nbas
      do k=j,nbas
            if( v_ssite2(j)%iantiferro+v_ssite2(k)%iantiferro==0) then
               iv_a_oips(k) = iv_a_oips(j) !to drive mksymx
               inumaf=inumaf+1
               exit
            endif
      enddo
      enddo
      if(imaster) call pshpr(60)
      call mksym(lc,slabl,strn2,v_ssite2,iv_a_oips) !strn2 and v_ssite2 are used.
      if(imaster) call poppr()
      if(imaster) write(6,"(a)")' AF: mksym, generator= SYMGRP+SYMGRPAF= '//trim(strn2)
      ngrpaf_     = lat_nsgrp 
      allocate(iclasstaf_(nbas),symops_af_(3,3,ngrpaf_),ag_af_(3,ngrpaf_))
      call dcopy ( ngrpaf_ * 9 , rv_a_osymgr , 1 , symops_af_ , 1 )
      call dcopy ( ngrpaf_ * 3 , rv_a_oag ,    1 , ag_af_ , 1 )
      if(imaster) write(6,"(a,i3)") ' AF: ngrpaf=',ngrpaf_
      do ib=1,nbas
         iclasstaf_(ib)=v_ssite2(ib)%class
      enddo
      deallocate(v_ssite2)
      end subroutine
      
cSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs      
      subroutine mksym(mode,slabl,ssymgr,ssite,iv_a_oips)
      use m_ftox
      use m_struc_def
      use m_lattic,only: lat_plat,rv_a_opos!,lat_dist
      use m_lmfinit,only: ctrl_nbas,ctrl_nsite,ctrl_nspec,stdo!,lat_ldist
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
Cio    Stored:    onrcp oipcp oics oipc oclabl onrc
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
      implicit none
      type(s_site)::ssite(*)
      character*8 slabl(*),ssymgr*(*)
      integer mode,nsgrp,npgrp,ibas
      integer ::iwdummy ,iwdummy1(1)
      logical T,F,cmdopt,ltmp
      integer:: idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nbas,
     .nbas0,nsite,nclass,ngen,ngnmx,nspec,usegen,
     .nggen,ngmx,incli, oiwk , aginv  !ldist,
      integer ,allocatable :: nrspc_iv(:)
      real(8) ,allocatable :: pos2_rv(:,:)
      integer ,allocatable :: ips2_iv(:)
      integer,allocatable:: iv_a_tmp(:)
      parameter (T=.true., F=.false., ngnmx=10)
      character*120 gens,strn*72
      double precision gen(9,ngnmx),plat(3,3),qlat(3,3),xx !,dist(3,3)
      integer:: i_copy_size,i_data_size,i_spackv
      integer:: iv_a_oips(*)
      nbas =ctrl_nbas
      nspec=ctrl_nspec
      nsite=ctrl_nsite
      plat =lat_plat
c      ldist=lat_ldist           ! ldist mode need to be recovered if necessary
c      dist=lat_dist
c      i_copy_size=size(lat_dist)
c      call dcopy(i_copy_size,lat_dist,1,dist,1)
c      stdo = globalvariables%stdo
      ngmx = 48
C ... Re-entry when ngmx was increased
    5 continue
C ... For memory allocation of class data, assume initially nclass=nbas
      if(allocated(rv_a_oag)) deallocate(rv_a_oag,rv_a_osymgr,iv_a_oipc,iv_a_oics)
      allocate( rv_a_oag(3*ngmx)    )
      allocate( rv_a_osymgr(9*ngmx) )
      allocate( iv_a_oipc(nsite)  )
      allocate( iv_a_oics(nsite)  )
      allocate( nrspc_iv(nbas) )
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
c        else if (ssymgr(j1:j2) .eq. 'point') then
c          ltmp = .true.
        else
          call strncp(gens,ssymgr,idest,j1,j2-j1+2)
          idest = idest+j2-j1+2
        endif
      enddo
c      if (ltmp) then
c        call strncp(gens,'point',idest,1,5)
c      endif
C --- Generate space group ---
      nbas0 = nbas
      if (cmdopt('--fixpos',8,0,strn)) call Rx('fixpos is removed. outside of lmf')
C ... When generating the group the basis may become enlarged ...
      if(allocated(iv_a_oistab)) deallocate(iv_a_oistab) !lat_iv_a_oistab
      allocate(iv_a_oistab(abs((ngmx+1)*nbas)))
      allocate(ips2_iv(ngmx*nbas))
      allocate(pos2_rv(3,ngmx*nbas))
      call icopy ( nbas , iv_a_oips , 1 , ips2_iv , 1 ) 
      call dcopy ( 3 * nbas , rv_a_opos , 1 , pos2_rv , 1 )
      call gensym ( slabl , gens , usegen , t , f , f , nbas 
     .   , nspec , ngmx , plat , plat , pos2_rv , ips2_iv !, ldist , dist
     .   , nrspc_iv , nsgrp , rv_a_osymgr , rv_a_oag , ngen , gen , ssymgr 
     .   , nggen , isym , iv_a_oistab ) 
      if (nbas .gt. nbas0) call rxs('gensym: the basis was enlarged.',' Check group operations.')
      if (nggen .gt. nsgrp) then
        write(stdo,ftox)' MKSYM (warning): generators create more than ngmx=',ngmx,' group ops ...'
        ngmx = ngmx*16
        if (allocated(pos2_rv)) deallocate(pos2_rv)
        if (allocated(ips2_iv)) deallocate(ips2_iv)
        if (allocated(nrspc_iv)) deallocate(nrspc_iv)
        goto 5
      endif
C --- Add inversion to point group ---
      incli = -1
      npgrp = nsgrp
      if (mod(mode,10) .ne. 0) then
        ngen = ngen+1
        gen(:,ngen) = [-1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0]
        call pshpr(iprint()-40)
        call grpgen ( gen ( 1 , ngen ) , 1 , rv_a_osymgr , npgrp , ngmx  ) 
        call poppr
        incli = npgrp-nsgrp
      endif
C --- Printout of symmetry operations ---
      write(stdo,ftox)' MKSYM: found ',nsgrp,' space group operations'
      if(nsgrp/=npgrp) write(stdo,ftox)' MKSYM: adding inversion gives ',npgrp,' operations'
      if(incli .eq. -1) write(stdo,*)'  no attempt to add inversion symmetry'
      if(mod(mode/10,10) .eq. 0) goto 100
! Split species into classes : ibas ==> iclass=ipc(ibas) ==> ispec=ics(iclass)
      if(allocated(iv_a_onrc)) deallocate(iv_a_onrc)
       allocate(iv_a_onrc(abs(nspec)))
       call icopy ( nsite , iv_a_oips , 1 , iv_a_oipc , 1 )
       call splcls ( mod ( mode / 10 , 10 ) .eq.4 , rv_a_opos , nbas 
     .     , nsgrp , iv_a_oistab , nspec , slabl , nclass , iv_a_oipc , 
     .     iv_a_oics , iv_a_onrc ) 
C   ... Reallocate arrays as permanent arrays
       i_data_size=size(iv_a_oics); allocate(iv_a_tmp(i_data_size))
       iv_a_tmp=iv_a_oics; deallocate(iv_a_oics)
       i_data_size=min(i_data_size,nclass); allocate(iv_a_oics(nclass))
       iv_a_oics(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)
C ... Remake istab
       if (allocated(iv_a_oistab)) deallocate(iv_a_oistab)
       allocate(iv_a_oistab(abs(nsgrp*nbas)))
       call dinv33(plat,1,qlat,xx)
       call symtbl(1, nbas, iwdummy1, rv_a_opos , rv_a_osymgr, rv_a_oag, nsgrp, qlat, iv_a_oistab) 
C ... poke ipc into ssite
      do ibas=1,nbas
         ssite(ibas)%class = iv_a_oipc(ibas)
      enddo
C --- Cleanup: poke class and symmetry info into structures ---
      ctrl_nclass=nclass
  100 continue
      lat_npgrp=npgrp
      lat_nsgrp=nsgrp
      end subroutine mksym
      end module m_mksym
