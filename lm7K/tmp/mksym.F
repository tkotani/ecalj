      subroutine mksym(mode,slabl,ssymgr,sctrl,slat,ssite,sarray)

      use m_struc_def           !Cgetarg
      use m_globalvariables
C     - Setup for symmetry group
C     ----------------------------------------------------------------------
C     i Inputs
C     i   mode  : 1s digit
C     i           0  make space group only
C     i           1  add inversion to point group (see aginv below)
C     i           2  Same as 1, but make additionally ag,istab for extra
C     i              operations, using -g for rotation part; see Remarks
C     i           10s digit
C     i           0  do nothing about classes
C     i           1  split species into classes
C     i           2  Also assign class labels to each class
C     i           4  Assume class and species are the same.
C     i   slabl : species labels
C     i   ssymgr: string containing symmetry group generators.
C     i           if ssymgr contains 'find', mksym will add basis atoms as
C     i           needed to guarantee generators are valid, and generate
C     i           internally any additonal group operations needed to
C     i           complete the space group.
C     i
C     io Inputs/outputs
C     i   sctrl :struct for program flow parameters; see routine uctrl
C     i     Elts read: nbas nspec nsite lpgf
C     i     Stored:    nclass
C     i     Passed to: *
C     io  slat  :struct containing information about the lattice
C     io    Elts read: plat ldist dist
C     io    Stored:    npgrp nsgrp oistab oag osymgr
C     io          :nsgrp= number of space group operations.
C     io          :npgrp= number of point group operations.  It may differ
C     io                  from nsgrp if inversion is added.
C     io          :istab= table of site permutations for each group op
C     io          :ag   = translation part the group ops
C     io          :symgr= point group (rotation) part of each group op
C     io  ssite :struct containing site-specific information
C     io    Elts read: *
C     i     Stored:    class clabel
C     io  sarray:structure containing offsets to various arrays
C     io    Elts read: opos oips npadl npadr
C     io    Stored:    onrcp oipcp oics oipc oclabl nclasp onrc
C     o           oipc,oipcp,onrc,oics and oclabl, and
C     o           oipc:  pointer to class table, pad equiv to double pad
C     o           oipcp: pointer to class table, pad nequiv to double pad
C     o Outputs
C     o  ssymgr:  may be enlarged if it contains 'find' initially
C     o  slat     the following are set: npgrp nsgrp oistab oag osymgr
C     r Remarks
C     r   In certain cases the inversion operation may be added to the space
C     r   group, for purposes of k integration.  This is permissible when the
C     r   hamiltonian has the form h(-k) = h*(k).  In that case, the
C     r   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
C     r
C     r   Also, the Green's functions are related G(-k) = Gtranspose(k).
C     r   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
C     r   then G(-k) = Gtranspose(k), and the same (g,ag) information is
C     r   needed for either rotation.
C     u Updates
C     u   04 Jan 06 Enabled keyword 'point' in ssymgr
C     u             Generates extra group ops when generators create more
C     u             than the maximum expected
C     u   07 Oct 05 Put altered positions into site structure
C     u   03 Nov 01 Adapted to altered gensym
C     u   26 Jan 01 Added new '2' option to 1s digit of mode
C     u   31 May 00 revised meaning of mode
C     ----------------------------------------------------------------------
C     implicit none
      character*8 slabl(1),ssymgr*(*)
      integer mode,nsgrp,npgrp
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_site)::ssite(*)
      type(s_array)::sarray

C     ... Heap
      integer ::iwdummy 

C     Local variables
      logical T,F,cmdopt,a2bin,ltmp
      integer idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nbas,
     .     nbas0,nbasp,nsite,nclass,nclasp,nclspp,ngen,ngnmx,nspec,usegen,
     .     ldist,nggen,ngmx,stdo,nglob,incli !npadl,npadr,
      integer:: oiwk , aginv 
c     kino Dec.8.2011:                   real(8),pointer :: rv_p_opos(:) =>NULL()

c     kino Dec.16.2011:                  integer,pointer :: iv_p_oipcp(:) =>NULL()

c     kino Dec.16.2011:                 integer,pointer :: iv_p_oipc(:) =>NULL()

c     kino Dec.14.2011:                real(8),pointer :: rv_p_osymgr(:) =>NULL()

c     kino Dec.14.2011:               real(8),pointer :: rv_p_oag(:) =>NULL()

c     kino Dec.16.2011:              integer,pointer :: iv_p_onrcp(:) =>NULL()

c     kino Dec.16.2011:             integer,pointer :: iv_p_onrc(:) =>NULL()

c     kino Dec.8.2011:            integer,pointer :: iv_p_oips(:) =>NULL()

c     kino Dec.14.2011:           integer,pointer :: iv_p_oistab(:) =>NULL()

c     kino Dec.14.2011:          integer,pointer :: iv_p_oics(:) =>NULL()

c     kino Dec.14.2011:         real(8),pointer :: rv_p_oclabl(:) =>NULL()

      integer ,allocatable :: nrspc_iv(:)
      real(8) ,allocatable :: pos2_rv(:)
      integer ,allocatable :: ips2_iv(:)

      integer,allocatable:: iv_a_tmp(:)

      parameter (T=.true., F=.false., ngnmx=10)
      character*120 gens,strn*72
      double precision gen(9,ngnmx),plat(3,3),dist(3,3),qlat(3,3),xx,
     .     fptol


      nbas=sctrl%nbas
      nspec=sctrl%nspec
      nsite=sctrl%nsite


      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      ldist=slat%ldist
      i_copy_size=size(slat%dist)
      call dcopy(i_copy_size,slat%dist,1,dist,1)

      ngmx = 48
C     hangenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      if (iprint() .ge. 30) write(stdo,'(1x)')

C     ... Re-entry when ngmx was increased
    5 continue

C     ... For memory allocation of class data, assume initially nclass=nbas
c     kino Dec.14.2011:         allocate(rv_p_oag(abs(3*ngmx)))
c     kino 12月.27.2011:         allocate(slat%rv_p_oag(abs(3*ngmx)))
      allocate(slat%rv_a_oag(abs(3*ngmx)))
c     kino Dec.14.2011:         if (3*ngmx<0) rv_p_oag(:)=0.0d0
c     kino 12月.27.2011:         if (3*ngmx<0) slat%rv_p_oag(:)=0.0d0
      if (3*ngmx<0) slat%rv_a_oag(:)=0.0d0

c     kino Dec.14.2011:         allocate(rv_p_osymgr(abs(9*ngmx)))
c     kino Dec.28.2011:         allocate(slat%rv_p_osymgr(abs(9*ngmx)))
      allocate(slat%rv_a_osymgr(abs(9*ngmx)))
c     kino Dec.14.2011:         if (9*ngmx<0) rv_p_osymgr(:)=0.0d0
c     kino Dec.28.2011:         if (9*ngmx<0) slat%rv_p_osymgr(:)=0.0d0
      if (9*ngmx<0) slat%rv_a_osymgr(:)=0.0d0

c     kino Dec.16.2011:         allocate(iv_p_oipc(abs(nsite)))
c     kino Dec.28.2011:         allocate(sarray%iv_p_oipc(abs(nsite)))
      allocate(sarray%iv_a_oipc(abs(nsite)))
c     kino Dec.16.2011:         if (nsite<0) iv_p_oipc(:)=0
c     kino Dec.28.2011:         if (nsite<0) sarray%iv_p_oipc(:)=0
      if (nsite<0) sarray%iv_a_oipc(:)=0

c     kino Dec.14.2011:         allocate(iv_p_oics(abs(nsite)))
c     kino Dec.27.2011:         allocate(sarray%iv_p_oics(abs(nsite)))
      allocate(sarray%iv_a_oics(abs(nsite)))
c     kino Dec.14.2011:         if (nsite<0) iv_p_oics(:)=0
c     kino Dec.27.2011:         if (nsite<0) sarray%iv_p_oics(:)=0
      if (nsite<0) sarray%iv_a_oics(:)=0

      allocate(nrspc_iv(nbas))


c     kino Dec.8.2011:         rv_p_opos => sarray%rv_p_opos 

c     kino Dec.8.2011:         iv_p_oips => sarray%iv_p_oips 


C     call defi(oips,nbas)
C     call defrr(opos,3*nbas)
C     call spackv(10,'site spec',ssite,1,nbas,w(oips))
C     call spackv(10,'site pos',ssite,1,nbas,w(opos))

*     call prmx('basis',w(opos),3,3,nbas)
*     print *, (w(oips-1+i), i=1,nbas)

C     ... Extract keywords from input string; copy rest to gens
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

C     --- Generate space group ---
      nbas0 = nbas
      if (cmdopt('--fixpos',8,0,strn)) then
         j1 = 8+1
         if (strn(9:13) .eq. ':tol=') then
            j1 = 13
         endif
         if (strn(9:9) .ne. ':' .or.
     .        .not. a2bin(strn,fptol,4,0,' ',j1,len(strn))) fptol = 1d-5
      else
         fptol = 0
      endif

C     ... When generating the group the basis may become enlarged ...
C     copy larger files relevant arrays to larger workspace
c     kino Dec.14.2011:         allocate(iv_p_oistab(abs((ngmx+1)*nbas)))
c     kino Dec.28.2011:         allocate(slat%iv_p_oistab(abs((ngmx+1)*nbas)))
      allocate(slat%iv_a_oistab(abs((ngmx+1)*nbas)))
c     kino Dec.14.2011:         if ((ngmx+1)*nbas<0) iv_p_oistab(:)=0.0d0
c     kino Dec.28.2011:         if ((ngmx+1)*nbas<0) slat%iv_p_oistab(:)=0.0d0
      if ((ngmx+1)*nbas<0) slat%iv_a_oistab(:)=0.0d0

C     call defi (oistb2,(ngmx+1)*nbas)
      allocate(ips2_iv(ngmx*nbas))

      allocate(pos2_rv(3*ngmx*nbas))

c     kino Dec.8.2011:         call icopy ( nbas , iv_p_oips , 1 , ips2_iv , 1 ) 
c     kino Dec.28.2011:         call icopy ( nbas , sarray%iv_p_oips , 1 , ips2_iv , 1 ) 
      call icopy ( nbas , sarray%iv_a_oips , 1 , ips2_iv , 1 ) 


c     kino Dec.8.2011:         call dcopy ( 3 * nbas , rv_p_opos , 1 , pos2_rv , 1 ) 
      call dcopy ( 3 * nbas , sarray%rv_p_opos , 1 , pos2_rv , 1 ) 


      call gensym ( slabl , gens , usegen , t , f , fptol , f , nbas 
     .     , nspec , ngmx , plat , plat , ldist , dist , pos2_rv , ips2_iv 
c     kino Dec.14.2011:       .   , nrspc_iv , nsgrp , rv_p_osymgr , rv_p_oag , ngen , gen , ssymgr 
c     kino Dec.14.2011:       .   , nrspc_iv , nsgrp , slat%rv_p_osymgr , rv_p_oag , ngen , gen , ssymgr 
c     kino 12月.27.2011:       .   , nrspc_iv , nsgrp , slat%rv_p_osymgr , slat%rv_p_oag , ngen , gen , ssymgr 
c     kino Dec.28.2011:       .   , nrspc_iv , nsgrp , slat%rv_p_osymgr , slat%rv_a_oag , ngen , gen , ssymgr 
     .     , nrspc_iv , nsgrp , slat%rv_a_osymgr , slat%rv_a_oag , ngen , gen , ssymgr 
c     kino Dec.14.2011:       .   , nggen , isym , iv_p_oistab ) 
c     kino Dec.28.2011:       .   , nggen , isym , slat%iv_p_oistab ) 
     .     , nggen , isym , slat%iv_a_oistab ) 




      if (nbas .gt. nbas0) call rxs('gensym: the basis was enlarged.',
     .     ' Check group operations.')
      if (nggen .gt. nsgrp) then
         call info2(1,0,0,' MKSYM (warning): generators create more than'
     .        //' ngmx=%i group ops ...',ngmx,0)
         ngmx = ngmx*16
         if (allocated(pos2_rv)) deallocate(pos2_rv)
         if (allocated(ips2_iv)) deallocate(ips2_iv)
         if (allocated(nrspc_iv)) deallocate(nrspc_iv)
c     kino Dec.14.2011:           if (associated(rv_p_oag)) deallocate(rv_p_oag)
c     kino 12月.27.2011:           if (associated(slat%rv_p_oag)) deallocate(slat%rv_p_oag)
         if (allocated(slat%rv_a_oag)) deallocate(slat%rv_a_oag)


         goto 5
      endif
      if (fptol .ne. 0) then
         i_copy_size=size(ssite(1)%pos)
         do i_spackv=1,nbas
            call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .           , i_copy_size , i_spackv + 1 - 1 , pos2_rv )

         enddo

         i_copy_size=size(ssite(1)%pos)
         do i_spackv=1,nbas
            call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
c     kino Dec.8.2011:       .     , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 
     .           , i_copy_size , i_spackv + 1 - 1 , sarray%rv_p_opos ) 

         enddo

      endif

C     ki#error ERROR, try to release name= oips2 ,but list=null at linenumber= 219 list= (None)


C     --- Add inversion to point group ---
      incli = -1
      npgrp = nsgrp
      if (mod(mode,10) .ne. 0) then
         ngen = ngen+1
         call dpzero(gen(1,ngen),9)
         gen(1,ngen) = -1
         gen(5,ngen) = -1
         gen(9,ngen) = -1
         call pshpr(iprint()-40)
c     kino Dec.14.2011:           call grpgen ( gen ( 1 , ngen ) , 1 , rv_p_osymgr , npgrp , ngmx 
c     kino Dec.28.2011:           call grpgen ( gen ( 1 , ngen ) , 1 , slat%rv_p_osymgr , npgrp , ngmx 
         call grpgen ( gen ( 1 , ngen ) , 1 , slat%rv_a_osymgr , npgrp , ngmx 
     .        ) 

         call poppr
         incli = npgrp-nsgrp
      endif
C     Not necessary, because it will be remade anyway
C     call redfi(oistab,npgrp*nbas)

C     --- Printout of symmetry operations ---
      call info5(20,0,0,' MKSYM:  found %i space group operations'//
     .     '%?#n#; adding inversion generated %i ops#%j#'//
     .     '%?#(n==0)# ... includes inversion##',nsgrp,npgrp-nsgrp,npgrp,
     .     incli,0)
      if (incli .eq. -1) then
         call info0(30,0,0,'         no attempt to add inversion symmetry')
      endif
      if (mod(mode/10,10) .eq. 0) goto 100

C     --- Split species into classes ---
c     kino Dec.16.2011:         allocate(iv_p_onrc(abs(nspec)))
      allocate(sarray%iv_p_onrc(abs(nspec)))
c     kino Dec.16.2011:         if (nspec<0) iv_p_onrc(:)=0
      if (nspec<0) sarray%iv_p_onrc(:)=0

c     kino Dec.8.2011:          call icopy ( nsite , iv_p_oips , 1 , iv_p_oipc , 1 ) 
c     kino Dec.16.2011:          call icopy ( nsite , sarray%iv_p_oips , 1 , iv_p_oipc , 1 ) 
c     kino Dec.28.2011:          call icopy ( nsite , sarray%iv_p_oips , 1 , sarray%iv_p_oipc , 1 ) 
c     kino Dec.28.2011:          call icopy ( nsite , sarray%iv_p_oips , 1 , sarray%iv_a_oipc , 1 ) 
      call icopy ( nsite , sarray%iv_a_oips , 1 , sarray%iv_a_oipc , 1 ) 



c     kino Dec.8.2011:             call splcls ( mod ( mode / 10 , 10 ) .eq.4 , rv_p_opos , nbas 
      call splcls ( mod ( mode / 10 , 10 ) .eq.4 , sarray%rv_p_opos , nbas 
c     kino Dec.14.2011:       .     , nsgrp , iv_p_oistab , nspec , slabl , nclass , iv_p_oipc , 
c     kino Dec.16.2011:       .     , nsgrp , slat%iv_p_oistab , nspec , slabl , nclass , iv_p_oipc , 
c     kino Dec.28.2011:       .     , nsgrp , slat%iv_p_oistab , nspec , slabl , nclass , sarray%iv_p_oipc , 
c     kino Dec.28.2011:       .     , nsgrp , slat%iv_a_oistab , nspec , slabl , nclass , sarray%iv_p_oipc , 
     .     , nsgrp , slat%iv_a_oistab , nspec , slabl , nclass , sarray%iv_a_oipc , 
c     kino Dec.14.2011:       .     iv_p_oics , iv_p_onrc ) 
c     kino Dec.16.2011:       .     sarray%iv_p_oics , iv_p_onrc ) 
c     kino Dec.27.2011:       .     sarray%iv_p_oics , sarray%iv_p_onrc ) 
     .     sarray%iv_a_oics , sarray%iv_p_onrc ) 





ccccccccccc
      print *,'zzz nclass=',nclass
cccccccccc

c     npadl=sarray%npadl
c     npadr=sarray%npadr

      nclasp = nclass
      nclspp = nclass
      nbasp  = nbas             !+ npadl + npadr

C     ... (pgf) Extra classes for padding layers
C$$$  lpgf = int(sctrl%lpgf(1))
C$$$  if (lpgf .gt. 0 .and. lpgf .lt. 10) then
C$$$  C   ... Remake ipc,ics for padding and double-padding layers
C$$$  call pgfpp(nbas,npadl,npadr,1,1,.true.,.false.,w(oipc),w(oipc))
C$$$  call redfi(oics,nsite)
C$$$  call defi(oiwk,nsite)
C$$$  call pvsym1(nbas,npadl+npadr,nclasp,w(oipc),w(oics))
C$$$  nclspp = nclasp
C$$$  call pvsym1(nbasp,npadl+npadr,nclspp,w(oipc),w(oics))
C$$$  if (iprint() .ge. 30) then
C$$$  call awrit2('         %i classes for embedding region + %i '//
C$$$  .      'for padding layers',' ',80,stdo,nclass,nclasp-nclass)
C$$$  endif
C$$$  call redfi(oics,   nclspp)
C$$$  C   ... Separate ipc, ipcd (dble padded ic eqiv =sngl pad vs not equiv)
C$$$  call defi (oipcp,  nsite)
C$$$  call icopy(nsite,w(oipc),1,w(oipcp),1)
C$$$  call icopy(npadl+npadr,w(oipc+nbas),1,w(oipc+nbasp),1)
C$$$  else
C     ... Reallocate arrays as permanent arrays
c     kino Dec.14.2011:         i_data_size=size(iv_p_oics); allocate(iv_a_tmp(i_data_size))
c     kino Dec.27.2011:         i_data_size=size(sarray%iv_p_oics); allocate(iv_a_tmp(i_data_size))
      i_data_size=size(sarray%iv_a_oics); allocate(iv_a_tmp(i_data_size))
c     kino Dec.14.2011:         iv_a_tmp=iv_p_oics; deallocate(iv_p_oics)
c     kino Dec.27.2011:         iv_a_tmp=sarray%iv_p_oics; deallocate(sarray%iv_p_oics)
      iv_a_tmp=sarray%iv_a_oics; deallocate(sarray%iv_a_oics)
c     kino Dec.14.2011:         i_data_size=min(i_data_size,nclspp); allocate(iv_p_oics(nclspp))
c     kino Dec.27.2011:         i_data_size=min(i_data_size,nclspp); allocate(sarray%iv_p_oics(nclspp))
      i_data_size=min(i_data_size,nclspp); allocate(sarray%iv_a_oics(nclspp))
c     kino Dec.14.2011:         iv_p_oics(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)
c     kino Dec.27.2011:         sarray%iv_p_oics(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)
      sarray%iv_a_oics(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)

c     ki      oipcp = oipc
c     kino Dec.16.2011:          iv_p_oipcp => iv_p_oipc
c     kino Dec.28.2011:          sarray%iv_p_oipcp => sarray%iv_p_oipc
C$$$  endif

C     ... Remake istab
      if (mod(mode,10) .eq. 2) then
c     kino Dec.14.2011:           allocate(iv_p_oistab(abs(npgrp*nbas)))
c     kino Dec.28.2011:           allocate(slat%iv_p_oistab(abs(npgrp*nbas)))
         if (allocated(slat%iv_a_oistab)) deallocate(slat%iv_a_oistab)
         allocate(slat%iv_a_oistab(abs(npgrp*nbas)))
c     kino Dec.14.2011:           if (npgrp*nbas<0) iv_p_oistab(:)=0.0d0
c     kino Dec.28.2011:           if (npgrp*nbas<0) slat%iv_p_oistab(:)=0.0d0
         if (npgrp*nbas<0) slat%iv_a_oistab(:)=0.0d0

      else
c     kino Dec.14.2011:           allocate(iv_p_oistab(abs(nsgrp*nbas)))
c     kino Dec.28.2011:           allocate(slat%iv_p_oistab(abs(nsgrp*nbas)))
         if (allocated(slat%iv_a_oistab)) deallocate(slat%iv_a_oistab)
         allocate(slat%iv_a_oistab(abs(nsgrp*nbas)))
c     kino Dec.14.2011:           if (nsgrp*nbas<0) iv_p_oistab(:)=0.0d0
c     kino Dec.28.2011:           if (nsgrp*nbas<0) slat%iv_p_oistab(:)=0.0d0
         if (nsgrp*nbas<0) slat%iv_a_oistab(:)=0.0d0

      endif
      call dinv33(plat,1,qlat,xx)
C     call prmx('pos',w(opos),3,3,nbas)
C     call symtab(nbas,w,w(opos),w(osymgr),w(oag),-nsgrp,plat,qlat,w,
C     .  w(oistab))
C     call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
C     call symtbl(0,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
C     .  w(oistab))
C     call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
c     kino Dec.8.2011:             call symtbl ( 1 , fptol , nbas , iwdummy , rv_p_opos , rv_p_osymgr 
c     kino Dec.14.2011:             call symtbl ( 1 , fptol , nbas , iwdummy , sarray%rv_p_opos , rv_p_osymgr 
c     kino Dec.28.2011:             call symtbl ( 1 , fptol , nbas , iwdummy , sarray%rv_p_opos , slat%rv_p_osymgr 
      call symtbl ( 1 , fptol , nbas , iwdummy , sarray%rv_p_opos , slat%rv_a_osymgr 
c     kino Dec.14.2011:       .     , rv_p_oag , nsgrp , qlat , iv_p_oistab ) 
c     kino Dec.14.2011:       .     , slat%rv_p_oag , nsgrp , qlat , iv_p_oistab ) 
c     kino 12月.27.2011:       .     , slat%rv_p_oag , nsgrp , qlat , slat%iv_p_oistab ) 
c     kino Dec.28.2011:       .     , slat%rv_a_oag , nsgrp , qlat , slat%iv_p_oistab ) 
     .     , slat%rv_a_oag , nsgrp , qlat , slat%iv_a_oistab ) 





C     Debugging check
C     call shoist(w(oistab),nbas,w(oag),w(osymgr),nsgrp)
C     call symtbl(0,fptol,nbas,w,w(opos),w(osymgr),w(oag),nsgrp,qlat,
C     .   w(oistb2))
C     call shoist(w(oistb2),nbas,w(oag),w(osymgr),nsgrp)
C     call istbpm(w(oistab),nbas,nsgrp,w(oistb2))
C     call shoist(w(oistb2),nbas,w(oag),w(osymgr),nsgrp)
C     stop

C     ... For each i*grpop fix translation ag to be that of grpop
      if (mod(mode,10) .eq. 2) then
c     kino Dec.14.2011:             if ( aginv ( rv_p_osymgr , nsgrp , npgrp , nbas , rv_p_oag , 
c     kino Dec.14.2011:             if ( aginv ( slat%rv_p_osymgr , nsgrp , npgrp , nbas , rv_p_oag , 
c     kino 12月.27.2011:             if ( aginv ( slat%rv_p_osymgr , nsgrp , npgrp , nbas , slat%rv_p_oag , 
c     kino Dec.28.2011:             if ( aginv ( slat%rv_p_osymgr , nsgrp , npgrp , nbas , slat%rv_a_oag , 
         if ( aginv ( slat%rv_a_osymgr , nsgrp , npgrp , nbas , slat%rv_a_oag , 
c     kino Dec.14.2011:       .     iv_p_oistab ) .ne. 0 ) call rx ( 'bug in mksym' ) 
c     kino Dec.28.2011:       .     slat%iv_p_oistab ) .ne. 0 ) call rx ( 'bug in mksym' ) 
     .        slat%iv_a_oistab ) .ne. 0 ) call rx ( 'bug in mksym' ) 



      endif

C     ... poke ipc into ssite
      i_copy_size=1;
      do i_spackv=1,nbasp
         call spackv_array_copy_i8_i ( 'p' , ssite ( i_spackv ) %class 
c     kino Dec.16.2011:       .   , i_copy_size , i_spackv + 1 - 1 , iv_p_oipc ) 
c     kino Dec.28.2011:       .   , i_copy_size , i_spackv + 1 - 1 , sarray%iv_p_oipc ) 
     .        , i_copy_size , i_spackv + 1 - 1 , sarray%iv_a_oipc ) 

      enddo

C     --- Create nrc and class labels ---
      if (mod(mode/10,10) .eq. 2 .or. mod(mode/10,10) .eq. 4) then
c     kino Dec.14.2011:           allocate(rv_p_oclabl(abs(nclspp)))
c     kino Dec.27.2011:           allocate(sarray%rv_p_oclabl(abs(nclspp)))
         allocate(sarray%rv_a_oclabl(abs(nclspp)))
c     kino Dec.14.2011:           if (nclspp<0) rv_p_oclabl(:)=0.0d0
c     kino Dec.27.2011:           if (nclspp<0) sarray%rv_p_oclabl(:)=0.0d0
         if (nclspp<0) sarray%rv_a_oclabl(:)=0.0d0

c     kino Dec.16.2011:           allocate(iv_p_onrcp(abs(nclspp)))
         allocate(sarray%iv_p_onrcp(abs(nclspp)))
c     kino Dec.16.2011:           if (nclspp<0) iv_p_onrcp(:)=0
         if (nclspp<0) sarray%iv_p_onrcp(:)=0

c     kino Dec.14.2011:              call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , iv_p_oics , iv_p_oipcp 
c     kino Dec.16.2011:              call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , sarray%iv_p_oics , iv_p_oipcp 
c     kino Dec.27.2011:              call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , sarray%iv_p_oics , sarray%iv_p_oipcp 
c     kino Dec.28.2011:              call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , sarray%iv_a_oics , sarray%iv_p_oipcp 
c     kino Dec.28.2011:              call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , sarray%iv_a_oics , sarray%iv_p_oipc 
         call pvsym2 ( 3 , 2 * nbasp - nbas , nclspp , sarray%iv_a_oics , sarray%iv_a_oipc 
c     kino Dec.14.2011:       .      , nspec , slabl , ssite , rv_p_oclabl , iv_p_onrcp ) 
c     kino Dec.16.2011:       .      , nspec , slabl , ssite , sarray%rv_p_oclabl , iv_p_onrcp ) 
c     kino Dec.27.2011:       .      , nspec , slabl , ssite , sarray%rv_p_oclabl , sarray%iv_p_onrcp ) 
     .        , nspec , slabl , ssite , sarray%rv_a_oclabl , sarray%iv_p_onrcp ) 




C$$$  if (nclasp .eq. nclass) then !we only treat nclasp=nclass case. takao Apr2010
c     ki        onrc = onrcp
c     kino Dec.16.2011:          iv_p_onrc => iv_p_onrcp
         sarray%iv_p_onrc => sarray%iv_p_onrcp
C$$$  else
C$$$  call defi (onrc,   nclasp)
C$$$  call icopy(nclasp,w(onrcp),1,w(onrc),1)
C$$$  call iaxpy(nclasp-nclass,1,w(onrcp+nclasp),1,w(onrc+nclass),1)
C$$$  C          call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
C$$$  C          call awrit2('%n:1i',' ',100,6,nclspp,w(onrcp))
C$$$  endif
      else
c     ki        oclabl = 0
c     kino Dec.14.2011:          nullify(rv_p_oclabl)
c     kino Dec.27.2011:          nullify(sarray%rv_p_oclabl)
         if (allocated(sarray%rv_a_oclabl)) deallocate(sarray%rv_a_oclabl)
c     ki        onrcp = 0
c     kino Dec.16.2011:          nullify(iv_p_onrcp)
         nullify(sarray%iv_p_onrcp)
      endif

C     call awrit2('%n:1i',' ',100,6,nsite,w(oipc))
C     call awrit2('%n:1i',' ',100,6,nsite,w(oipcp))
C     pause

C     --- Cleanup: poke class and symmetry info into structures ---

      sctrl%nclass=nclass


c     kino Dec.16.2011:         sarray%iv_p_onrcp => iv_p_onrcp 

c     kino Dec.16.2011:         sarray%iv_p_oipcp => iv_p_oipcp 



c     kino Dec.14.2011:         sarray%iv_p_oics => iv_p_oics 

c     kino Dec.16.2011:         sarray%iv_p_oipc => iv_p_oipc 

c     kino Dec.14.2011:         sarray%rv_p_oclabl => rv_p_oclabl 

      sarray%nclasp=nclasp
c     kino Dec.16.2011:         sarray%iv_p_onrc => iv_p_onrc 


c     call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
 100  continue

      slat%npgrp=npgrp
      slat%nsgrp=nsgrp
c     kino Dec.14.2011:         slat%iv_p_oistab => iv_p_oistab 

c     kino Dec.14.2011:         slat%rv_p_oag => rv_p_oag 

c     kino Dec.14.2011:         slat%rv_p_osymgr => rv_p_osymgr 



      end subroutine mksym

C$$$  subroutine pvsym1(nbas,npad,nclass,ipc,ics)

C$$$  C- Fix classes for padding layers
C$$$  C ----------------------------------------------------------------------
C$$$  Ci Inputs
C$$$  Ci   nbas   number of basis atoms without padding layers
C$$$  Ci   nclass number of classes without padding layers
C$$$  Ci Outputs
C$$$  Co   nclass updated
C$$$  Co   ipc(nbas..nbas+npad) updated
C$$$  Co   ics updated
C$$$  C ----------------------------------------------------------------------
C$$$  C     implicit none
C$$$  integer nbas,nclass,ipc(nbas),ics(nclass)
C$$$  integer i,npad,nc,is,ic,j
C$$$  if (npad .eq. 0) return

C$$$  C ... ipc->-ipc to flag which classes have been renumbered
C$$$  do  5  i = nbas+1, nbas+npad
C$$$  5 ipc(i) = -ipc(i)

C$$$  C ... For each site in padding layer, assign class
C$$$  nc = nclass
C$$$  C ... Look for next occurrence of class not assigned
C$$$  do  10  i = nbas+1, nbas+npad

C$$$  if (ipc(i) .gt. 0) goto 10
C$$$  C   ... Class ipc(i) need to be assigned:  do for each site in class
C$$$  C       ic is (-)original class, is is species of original class
C$$$  ic = ipc(i)
C$$$  is = ics(-ic)

C$$$  nc = nc+1
C$$$  do  12  j = i, nbas+npad
C$$$  12   if (ipc(j) .eq. ic) ipc(j) = nc

C$$$  C   ... Species index is same as that of original class
C$$$  ics(nc) = is

C$$$  10 continue

C$$$  nclass = nc
C$$$  C      call awrit2('%n:1i',' ',100,6,nbas+npad,ipc)
C$$$  C      call awrit2('%n:1i',' ',100,6,nclass,ics)
C$$$  C      pause

C$$$  Cgetarg       end
C$$$  end subroutine pvsym1

      subroutine pvsym2(mode,nbas,nclass,ics,ipc,nspec,slabl,ssite,
     .     dclabl,nrc)

      use m_struc_def           !Cgetarg

C     - Create class labels from species labels (double precision format)
C     ----------------------------------------------------------------------
C     i Inputs:
C     i   nosplt: T copy class and species
C     i     mode: 0 do nothing
C     i           1 create class labels clabl
C     i           2 create number of sites in each class nrc
C     i      ipc: for padding sites ib ipc(ib) = class
C     i   nclass: number of classes
C     i   ssite :struct for site-specific information; see routine usite
C     i     Elts read: *
C     i     Stored:    clabel
C     o Outputs:
C     o   dclabl: class labels in double precision format
C     o      nrc: number of sites in each class
C     u Updates
C     u   18 Dec 01 Packs class label into ssite->clabel
C     ----------------------------------------------------------------------
C     implicit none
      integer mode,nbas,nclass,nspec,ics(1),ipc(nbas),nrc(1)
      real(8):: dclabl(nclass)
      type(s_site)::ssite(*)

      character*8 slabl(nspec)
      integer ic,iclbsj,idx,is,ib
      character*8 clabl

C     --- Make class labels from species labels ---
      if (mod(mode,2) .eq. 1) then
         do  10  is = 1, nspec
            do  12  idx = 1, nbas
               ic = iclbsj(is,ics,-nclass,idx)
               if (ic .lt. 0) goto 13
               call clabel(slabl,is,idx,clabl)
               call s8tor8(clabl,dclabl(ic))
 12         continue
 13         continue
 10      continue
      endif

      do  20  ib = 1, nbas
         ic = ipc(ib)
         call r8tos8(dclabl(ic),clabl)
         do i_spacks=ib,ib
            call spacks_copy('p',ssite(i_spacks)%clabel,ib,ib,clabl,i_spacks)
         enddo

 20   continue

C     do  ib = 1, nbas
C     call spacks(0,'site clabel',ssite,clabl,ib,ib)
C     print *, ib, clabl
C     enddo
C     pause

C     --- Create nrc ---
      if (mod(mode/2,2) .eq. 1) then
         call iinit(nrc,nclass)
         do  30  ib = 1, nbas
            ic = ipc(ib)
            nrc(ic) = nrc(ic)+1
 30      continue
      endif

C     ... debugging : check class labels
C     call awrit2('%n:1i',' ',100,6,nclass,nrc)
C     do  50  ic = 1, nclass
C     call r8tos8(dclabl(ic),clabl)
C     print 333, ic, clabl
C     333   format(i4,1x,a)
C     50 continue
C     pause

      end subroutine pvsym2


      integer function aginv(g,npg,nsg,nbas,ag,istab)

C     - Finds translation part for extra inversion-generated group ops
C     ----------------------------------------------------------------------
C     i Inputs
C     i   g     :point group operations
C     i   npg   :number of point group operations
C     i   nsg   :number of space group operations
c     i   nbas  :dimensions istab.  If zero, istab is not altered.
C     o Outputs
C     o   ag    :for each group op ig in range nsg+1..npg,
C     o         :the corresponding group op jg is found such that
C     o         :g(jg) = -g(ig), and ag(jg) is copied to ag(ig); see Remarks
C     o   istab :istab(jg) is copied into istab(ig) just in the same way as
C     o         :ag(jg) is copied into ag(ig).
C     o         :The sign of istab(1) is set negative to indicate that this
C     o         :operation is derived from a true space group operation.
C     r Remarks
C     r   In certain cases the inversion operation may be added to the space
C     r   group, for purposes of k integration.  This is permissible when the
C     r   hamiltonian has the form h(-k) = h*(k).  In that case, the
C     r   eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
C     r
C     r   Also, the Green's functions are related G(-k) = Gtranspose(k).
C     r   Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
C     r   then G(-k) = Gtranspose(k), and the same (g,ag) information is
C     r   needed for either rotation.  Let g~ = -g:
C     r     G0(g~^-1 -k) into G(k),
C     r
C     r
C     r   Routine aborts with aginv=ig if for group ig it cannot find
C     r   corresponding jg such that g(jg) = -g(ig)
C     ----------------------------------------------------------------
C     implicit none
C     ... Passed parameters
      integer npg,nsg,nbas,istab(nbas,1)
c     double precision g(3,3,npg),ag(3,npg)
      double precision g(3,3,*),ag(3,*)
C     ... Local parameters
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
C     call awrit2(' aginv mapping %,2i ->%,2i',' ',80,6,ig,jg)
               goto 10
            endif
         enddo
         aginv = ig
         return
 10      continue
      enddo
      aginv = 0
      end function aginv


