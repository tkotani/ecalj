      subroutine mksym(mode,slabl,ssymgr,sctrl,slat,ssite,sarray)
      use m_struc_def  !Cgetarg
      use m_globalvariables
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
      integer ::iwdummy 
      logical T,F,cmdopt,a2bin,ltmp
      integer idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nbas,
     .nbas0,nbasp,nsite,nclass,nclasp,nclspp,ngen,ngnmx,nspec,usegen,
     .ldist,nggen,ngmx,stdo,nglob,incli !npadl,npadr,
                 integer:: oiwk , aginv 
      integer ,allocatable :: nrspc_iv(:)
      real(8) ,allocatable :: pos2_rv(:,:)
      integer ,allocatable :: ips2_iv(:)
      integer,allocatable:: iv_a_tmp(:)
      parameter (T=.true., F=.false., ngnmx=10)
      character*120 gens,strn*72
      double precision gen(9,ngnmx),plat(3,3),dist(3,3),qlat(3,3),xx,fptol

      nbas =sctrl%nbas
      nspec=sctrl%nspec
      nsite=sctrl%nsite
      plat =slat%plat
      ldist=slat%ldist
      i_copy_size=size(slat%dist)
      call dcopy(i_copy_size,slat%dist,1,dist,1)
      ngmx = 48
      stdo = globalvariables%stdo
      if (iprint() .ge. 30) write(stdo,'(1x)')
C ... Re-entry when ngmx was increased
    5 continue
C ... For memory allocation of class data, assume initially nclass=nbas
      if(allocated(slat%rv_a_oag)) then
        deallocate(slat%rv_a_oag,slat%rv_a_osymgr,sarray%iv_a_oipc,sarray%iv_a_oics)
      endif  
      allocate( slat%rv_a_oag(3*ngmx)    )
      allocate( slat%rv_a_osymgr(9*ngmx) )
      allocate( sarray%iv_a_oipc(nsite)  )
      allocate( sarray%iv_a_oics(nsite)  )
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
      if(allocated(slat%iv_a_oistab)) deallocate(slat%iv_a_oistab)
      allocate(slat%iv_a_oistab(abs((ngmx+1)*nbas)))
      allocate(ips2_iv(ngmx*nbas))
      allocate(pos2_rv(3,ngmx*nbas))
      call icopy ( nbas , sarray%iv_a_oips , 1 , ips2_iv , 1 ) 
      call dcopy ( 3 * nbas , slat%rv_a_opos , 1 , pos2_rv , 1 ) 
      call gensym ( slabl , gens , usegen , t , f , fptol , f , nbas 
     .   , nspec , ngmx , plat , plat , ldist , dist , pos2_rv , ips2_iv 
     .   , nrspc_iv , nsgrp , slat%rv_a_osymgr , slat%rv_a_oag , ngen , gen , ssymgr 
     .   , nggen , isym , slat%iv_a_oistab ) 

      if (nbas .gt. nbas0) call rxs('gensym: the basis was enlarged.',
     .' Check group operations.')
      if (nggen .gt. nsgrp) then
        call info2(1,0,0,' MKSYM (warning): generators create more than'
     .  //' ngmx=%i group ops ...',ngmx,0)
        ngmx = ngmx*16
        if (allocated(pos2_rv)) deallocate(pos2_rv)
        if (allocated(ips2_iv)) deallocate(ips2_iv)
        if (allocated(nrspc_iv)) deallocate(nrspc_iv)
        if (allocated(slat%rv_a_oag)) deallocate(slat%rv_a_oag)
        goto 5
      endif
      if (fptol .ne. 0) then
        do i_spackv=1,nbas
           ssite( i_spackv )%pos = pos2_rv(:,i_spackv)
          slat%rv_a_opos(3*(i_spackv-1)+1:3*(i_spackv-1)+3)=ssite(i_spackv)%pos
        enddo
      endif
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
        call grpgen ( gen ( 1 , ngen ) , 1 , slat%rv_a_osymgr , npgrp , ngmx  ) 
        call poppr
        incli = npgrp-nsgrp
      endif
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
      if(allocated(sarray%iv_a_onrc)) deallocate(sarray%iv_a_onrc)
       allocate(sarray%iv_a_onrc(abs(nspec)))
       call icopy ( nsite , sarray%iv_a_oips , 1 , sarray%iv_a_oipc , 1 ) 
       call splcls ( mod ( mode / 10 , 10 ) .eq.4 , slat%rv_a_opos , nbas 
     .     , nsgrp , slat%iv_a_oistab , nspec , slabl , nclass , sarray%iv_a_oipc , 
     .     sarray%iv_a_oics , sarray%iv_a_onrc ) 

c       npadl=sarray%npadl
c       npadr=sarray%npadr

      nclasp = nclass
      nclspp = nclass
      nbasp  = nbas !+ npadl + npadr

C   ... Reallocate arrays as permanent arrays
       i_data_size=size(sarray%iv_a_oics); allocate(iv_a_tmp(i_data_size))
       iv_a_tmp=sarray%iv_a_oics; deallocate(sarray%iv_a_oics)
       i_data_size=min(i_data_size,nclspp); allocate(sarray%iv_a_oics(nclspp))
       sarray%iv_a_oics(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)

C ... Remake istab
      if (mod(mode,10) .eq. 2) then
         if (allocated(slat%iv_a_oistab)) deallocate(slat%iv_a_oistab)
         allocate(slat%iv_a_oistab(abs(npgrp*nbas)))
         if (npgrp*nbas<0) slat%iv_a_oistab(:)=0.0d0
      else
         if (allocated(slat%iv_a_oistab)) deallocate(slat%iv_a_oistab)
         allocate(slat%iv_a_oistab(abs(nsgrp*nbas)))
         if (nsgrp*nbas<0) slat%iv_a_oistab(:)=0.0d0
      endif
      call dinv33(plat,1,qlat,xx)
           call symtbl ( 1 , fptol , nbas , iwdummy , slat%rv_a_opos , slat%rv_a_osymgr 
     .     , slat%rv_a_oag , nsgrp , qlat , slat%iv_a_oistab ) 
C   ... For each i*grpop fix translation ag to be that of grpop
      if (mod(mode,10) .eq. 2) then
           if ( aginv ( slat%rv_a_osymgr , nsgrp , npgrp , nbas , slat%rv_a_oag , 
     .     slat%iv_a_oistab ) .ne. 0 ) call rx ( 'bug in mksym' ) 
      endif
C ... poke ipc into ssite
      i_copy_size=1;
      do i_spackv=1,nbasp
         ssite ( i_spackv )%class = sarray%iv_a_oipc(i_spackv)
      enddo
C --- Create nrc and class labels ---
      if (mod(mode/10,10) .eq. 2 .or. mod(mode/10,10) .eq. 4) then
        allocate(sarray%rv_a_oclabl(nclspp))
        allocate(sarray%iv_a_onrcp(nclspp))
        sarray%rv_a_oclabl(:)='' !0.0d0
        sarray%iv_a_onrcp(:)=0
        call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , sarray%iv_a_oics , sarray%iv_a_oipc 
     .   , nspec , slabl , ssite , sarray%rv_a_oclabl , sarray%iv_a_onrcp ) 
        deallocate(sarray%iv_a_onrc)
        call move_alloc(from=sarray%iv_a_onrcp, to=sarray%iv_a_onrc)
      else
        if (allocated(sarray%rv_a_oclabl)) deallocate(sarray%rv_a_oclabl)
      endif
C --- Cleanup: poke class and symmetry info into structures ---
      sctrl%nclass=nclass
      sarray%nclasp=nclasp
  100 continue
      slat%npgrp=npgrp
      slat%nsgrp=nsgrp
      end subroutine mksym

!! -----------------------------------
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
      type(s_site)::ssite(*)
      character*8 slabl(nspec)
      integer ic,iclbsj,idx,is,ib
      character(8):: clabl,dclabl(nclass)
C --- Make class labels from species labels ---
      if (mod(mode,2) .eq. 1) then
        do  10  is = 1, nspec
          do  12  idx = 1, nbas
            ic = iclbsj(is,ics,-nclass,idx)
            if (ic .lt. 0) goto 13
            call clabel(slabl,is,idx,clabl)
c            call s8tor8(clabl,dclabl(ic))
            dclabl(ic)=clabl
   12     continue
   13     continue
   10   continue
      endif
      do  20  ib = 1, nbas
        ic = ipc(ib)
        ssite(ib)%clabel = dclabl(ic) !clabl
   20 continue
C --- Create nrc ---
      if (mod(mode/2,2) .eq. 1) then
        call iinit(nrc,nclass)
        do  30  ib = 1, nbas
          ic = ipc(ib)
          nrc(ic) = nrc(ic)+1
   30   continue
      endif
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
c      double precision g(3,3,npg),ag(3,npg)
      double precision g(3,3,*),ag(3,*)
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


