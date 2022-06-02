module m_mksym_util
  use m_lgunit,only:stdo
  public gensym,grpgen,symtbl
  private
  real(8):: fptol=0d0
contains
  !$$$      subroutine pvsym2(mode,nbas,nclass,ics,ipc,nspec,slabl,ssite,
  !$$$     .dclabl,nrc)
  !$$$      use m_struc_def  !Cgetarg
  !$$$C- Create class labels from species labels (double precision format)
  !$$$C ----------------------------------------------------------------------
  !$$$Ci Inputs:
  !$$$Ci   nosplt: T copy class and species
  !$$$Ci     mode: 0 do nothing
  !$$$Ci           1 create class labels clabl
  !$$$Ci           2 create number of sites in each class nrc
  !$$$Ci      ipc: for padding sites ib ipc(ib) = class
  !$$$Ci   nclass: number of classes
  !$$$Ci   ssite :struct for site-specific information; see routine usite
  !$$$Ci     Elts read: *
  !$$$Ci     Stored:    clabel
  !$$$Co Outputs:
  !$$$Co   dclabl: class labels in double precision format
  !$$$Co      nrc: number of sites in each class
  !$$$Cu Updates
  !$$$Cu   18 Dec 01 Packs class label into ssite->clabel
  !$$$C ----------------------------------------------------------------------
  !$$$      implicit none
  !$$$      integer mode,nbas,nclass,nspec,ics(1),ipc(nbas),nrc(1)
  !$$$      type(s_site)::ssite(*)
  !$$$      character*8 slabl(nspec)
  !$$$      integer ic,iclbsj,idx,is,ib
  !$$$      character(8):: clabl,dclabl(nclass)
  !$$$c$$$C --- Make class labels from species labels ---
  !$$$c$$$      if (mod(mode,2) .eq. 1) then
  !$$$c$$$        do  10  is = 1, nspec
  !$$$c$$$          do  12  idx = 1, nbas
  !$$$c$$$            ic = iclbsj(is,ics,-nclass,idx)
  !$$$c$$$            if (ic .lt. 0) goto 13
  !$$$c$$$            call clabel(slabl,is,idx,clabl)
  !$$$c$$$c            call s8tor8(clabl,dclabl(ic))
  !$$$c$$$            dclabl(ic)=clabl
  !$$$c$$$   12     continue
  !$$$c$$$   13     continue
  !$$$c$$$   10   continue
  !$$$c$$$      endif
  !$$$c$$$      do  20  ib = 1, nbas
  !$$$c$$$        ic = ipc(ib)
  !$$$c$$$        ssite(ib)%clabel = dclabl(ic) !clabl
  !$$$c$$$   20 continue
  !$$$C --- Create nrc ---
  !$$$      if (mod(mode/2,2) .eq. 1) then
  !$$$        call iinit(nrc,nclass)
  !$$$        do  30  ib = 1, nbas
  !$$$          ic = ipc(ib)
  !$$$          nrc(ic) = nrc(ic)+1
  !$$$   30   continue
  !$$$      endif
  !$$$      end subroutine pvsym2

  subroutine gensym(slabl,gens,usegen,lcar,lfix,lsmall,nbas, &
       nspec,ngmx,plat,platcv,bas,ips,nrspec,ng,g, &
       ag,ngen,gen,nwgens,nggen,isym,istab)
    !- Generate the space group !,ldist,dist
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   slabl: name of the different species.
    !i   gens:  a list of generators, in symbolic representation
    !i          NB: this list is not required; see Remarks.
    !i   usegen:0 Find any additional group operations for this basis.
    !i          1 Also, extra basis atoms are added as needed to guarantee
    !i            the group operations created from gens are valid.
    !i          2 Do neither 0 nor 1.
    !i   lcar:  (not used)
    !i          T express ag,positions in cartesian coordinates
    !i          F express in units of conventional unit cell
    !i   lfix:  T: do not rotate or shift lattice
    !i   fptol: >0:Adjust positions slightly, rendering them as exactly
    !i          possible consistent with the symmetry group.  Any sites
    !i          within a lattice vector of tol are considered to be
    !i          at the same point.
    !i   nspec: number of classes, atoms in same class are symmetry-related
    !i   plat:  primitive lattice vectors (scaled by alat)
    !i   platcv:Used to scale translation part of generators,
    !i         :when translation part specified as a multiple of
    !i         :lattice vectors.  Can be same as plat but
    !i         :primitive lattice vectors of "conventional unit cell"
    !i         :are sometimes used to specify these translations, e.g.
    !i         :when generated from spacegroup data in some books.
    !i   ldist: lattice deformation matrix key; see lattdf
    !i   dist:  lattice deformation matrix; see lattdf
    ! o Inputs/Outputs (altered only if usegen=F)
    ! o  nbas:  On input, number of atoms in the basis
    ! o         On output nbas may be enlarged, depending symops and usegen
    ! o  bas:   basis vectors
    ! o         On output bas may be enlarged, depending symops and usegen
    ! o  ips:   the jth atom belongs to spec ips(j)
    ! o         On output ips may be enlarged, depending symops and usegen
    !o Outputs:
    !o   istab: site ib is transformed into istab(ib,ig) by operation ig
    !o   g:     symmetry operation matrix (assumed dimensioned >=ngmx)
    !o   ag:    symmetry operation vector (assumed dimensioned >=ngmx)
    !o   ... The following are generated if usegen=F
    !o   isym:  numbers characterizing the symmetry of lattice and crystal
    !o          isym(1) produces index for underlying lattice (see symlat)
    !o   lsmall:if T: a smaller unit cell can be found
    !o   nrspec:number of atoms in the ith class
    !o   ng:    number of group operations
    !o   ngen:  number of symmetry generators
    !o   gen:   generators in matrix form
    !o   nwgens:generators in ascii form
    !o   nggen :number of group ops generated by generators.
    !o         :Usually nggen=ng; however nggen can exceed ng if
    !o         :supercell is artificial -> extra translations; see groupg
    !l Local variables
    !l   modes: 0 -> sgroup compares point and space groups
    !l          1 -> sgroup compares point groups only
    !l          (set with 'points' keyword in gens)
    !r Remarks:
    !r   gensym generates the space group, using the following prescription:
    !r     1.  Any generators supplied from input gens
    !r         are checked for consistency with the underlying lattice.
    !r     2.  The space group is made from these generators.
    !r     3.  if usegen<2, missing basis atoms are added to make
    !r         the basis consistent with the supplied symmetry.
    !r     4.  nrspec is created
    !r     ... Unless usegen is 0, nothing more is done
    !r     5.  The point group of the underlying lattice without the
    !r         basis is generated.
    !r     6.  The full space group is generated from the point group
    !r     7.  A set of generators for this group is created
    !r   This program was adapted from the Stuttgart ASA version lmto-46.
    !b Bugs:
    !b   auto symmetry finder can fail with supercells, where extra
    !b   group operations include the same point group but inequivalent
    !b   translations.  Solution: have symcry call sgroup to see
    !b   if the space group is enlarged.  If so, space group should be
    !b   enlarged.
    !u Updates
    !u   04 Jan 06 Enabled keyword 'point' in ssymgr, returns if ng>ngmx
    !u   13 Dec 03 Uses platcv when scaling translation part of symgrp
    !u   05 Apr 03 Call sgroup looking only for point group ops;see bugs
    !u   03 Nov 01 Shortened argument list, eliminating duplicate bas,ips
    ! ----------------------------------------------------------------------
    implicit none
    ! Passed parameters:
    integer :: nbas,isym(*),istab(nbas,*),nspec,ngen,ngmx, &
         ng,nrspec(nspec),usegen,ips(nbas),nggen !,ldist
    double precision :: plat(9),platcv(9),g(9,*),ag(3,*),bas(3,nbas)!,fptol=0d0
    character(8) ::  slabl(*), gens*(*), nwgens*(*)
    logical :: lcar,lfix
    integer:: i , j , ibas , ic , iprint , ngnmx , igen , mxint , modes,ig
    real(8) ,allocatable :: wk_rv(:)
    double precision :: qlat(3,3),vol,platt(9)
    logical :: lsmall!,latvec
    parameter(ngnmx=10)
    double precision :: gen(9,ngnmx),agen(3,ngnmx)
    integer ::iwdummy,iwdummy1(1)
    !      stdo = lgunit(1)
    call rxx(.not. lcar, 'gensym not implemented lcar')
    call rxx(lfix,  'gensym not implemented lfix')
    call rxx(lsmall,'gensym not implemented lsmall')
    call dinv33(plat,1,qlat,vol) !Reciprocal lattice vectors ---
    ! Symmetry group as given by input generators ---
    nwgens = gens
    modes = 0
    call words(gens,ngen)
    if (ngen > 0) then
       call word(gens,ngen,i,j)
       if (gens(i:j) == 'point') then
          j = i-1
          modes = 1
       endif
       call psymop(gens(1:j),platcv,gen,agen,ngen)
       nwgens = gens(1:j)
    endif
    ! ... Rotate the generators
    !      call pshpr(iprint()-11)
    !      call lattdf ( - ldist , dist , plat , 0 , iwdummy , ngen , gen )
    !      call poppr
    do  10  igen = 1, ngen
       call grpprd(gen(1,igen),plat,platt)
       !       call dmpy(gen(1,igen),3,1,plat,3,1,platt,3,1,3,3,3)
       if ( .NOT. latvec(3,1d-5,qlat,platt)) &
            call fexit(-1,111,' Exit -1 GENSYM: '// &
            'generator %i imcompatible with underlying lattice',igen)
10  enddo

    ! ... Set up space group (g,ag,ng) given point group generators gen
    call sgroup(10+modes,gen,agen,ngen,g,ag,nggen,ngmx,qlat)
    ng = min(nggen,ngmx)
    if (nggen > ngmx) return

    !$$$C --- Add new atoms to the basis according to symmetry ---
    !$$$      if (usegen .lt. 2) then
    !$$$        call addbas(fptol,bas,slabl,ips,nbas,ng,qlat,g,ag)
    !$$$      endif
    ! ... Make nrspec ... i should be nspec
    i = mxint(nbas,ips)
    if (i /= nspec .AND. iprint() > 0) &
         call awrit2(' GENSYM (warning) %i species supplied but only '// &
         '%i spec used ...%N%8fpossible errors in class data',' ',120,6, &
         nspec,i)
    nspec = i
    call iinit(nrspec,nspec)
    do  22  ibas = 1, nbas
       ic = ips(ibas)
       nrspec(ic) = nrspec(ic)+1
22  enddo
    ! --- check if unit cell is the smallest possible one (not implemented)
    !      call chkcel(alat,bas,csym,ips,isym,lsmall,nbas,nspec,
    !     .            nrspec,plat,qlat)

    ! --- Complete the space group ---
    if (usegen == 0) then
       !       call rotlat(alat,bas,csym,isym,lfix,nbas,plat,qlat)
       !   ... Symmetry of lattice without basis
       call symlat(plat,ng,g,isym(1))
       !   ... Symmetry of lattice with basis
       allocate(wk_rv(3*nbas))
       call symcry (bas , wk_rv , ips , nbas , nspec , nrspec &
            , ng , plat , qlat , g , ag , istab )
       if (allocated(wk_rv)) deallocate(wk_rv)
       if (ng > ngmx) return
       nwgens = ' '
       call groupg(0,ng,g,ag,plat,ngen,gen,agen,nwgens,nggen)
    else
       nggen = 0
       call symtbl ( 0 , nbas , iwdummy1 , bas , g , ag , ng , qlat , istab )
    endif
    ! --- Make istab2 ---
    !      do  i = 1, ng
    !        do  ibas = 1, nbas
    !          ic = istab(ibas,i)
    !          istab2(ic,i) = ibas
    !        enddo
    !      enddo
    if (iprint() >= 80 .AND. ng > 1) then
       write(stdo,*)'  ib  istab ...'
       do  i = 1, nbas
          write(stdo,"(i4,':',48i3)") i, (istab(i,ig), ig=1,ng)
       enddo
       write(stdo,*)' GENSYM: site permutation table for group operations ...'
    endif
    ! --- Adjust basis to conform with symops to numerical precision ---
    !      if (fptol .gt. 0) call fixpos(bas,nbas,fptol,ng,plat,g,ag,istab)
    !      endif
  end subroutine gensym

  !$$$      subroutine addbas(tol,bas,clabl,ips,nbas,ngrp,qlat,g,ag)
  !$$$C- Adds the missing basis atoms to get the right symmetry
  !$$$C ----------------------------------------------------------------------
  !$$$Ci Inputs:
  !$$$Ci   clabl :name of the different inequivalent atom
  !$$$Ci   ips:the jth atom belongs to spec ips(j)
  !$$$Ci   ngrp  :number of group operations
  !$$$Ci   qlat  :primitive translation vectors in reciprocal space
  !$$$Ci   g     :symmetry operation matrix
  !$$$Ci   symops:symmetry operation symbol
  !$$$Ci   ag    :symmetry operation vector (dimensionless)
  !$$$Cio Inputs/Output:
  !$$$Cio  bas   :basis vectors (dimensionless)
  !$$$Cio         on output list has been completed by the new positions
  !$$$Cio  nbas  :number of atoms in the basis
  !$$$Cr Remarks:
  !$$$Cr For each atom, the symmetry-related positions are generated.
  !$$$Cr If this position is empty a new atom is added.
  !$$$Cr At the end, a check is made to ensure that
  !$$$Cr atoms of different speces do not occupy the same positions.
  !$$$C ----------------------------------------------------------------------
  !$$$C     implicit none
  !$$$C Passed parameters:
  !$$$      integer ips(*),nbas,ngrp
  !$$$      double precision qlat(3,3),bas(3,*),g(9,*),ag(3,*),tol
  !$$$      character*8 clabl(*)
  !$$$C Local parameters:
  !$$$      integer i,isop,ipr,lgunit,nbasnw,ibas,jbas,ic,jc,m,novlp,stdo
  !$$$      double precision bast(3),dbas(3),tol1
  !$$$      logical latvec
  !$$$      character*50 sg
  !$$$C External calls:
  !$$$      external  daxpy,dcopy,dmpy,errmsg,iprint,lgunit
  !$$$C Intrinsic functions:
  !$$$      intrinsic  idnint
  !$$$      character(200)::aaa
  !$$$
  !$$$      tol1 = tol
  !$$$      if (tol .eq. 0) tol1 = 1d-5
  !$$$      stdo = lgunit(1)
  !$$$
  !$$$C --- If no atom of same spec at the rotated site, add it ---
  !$$$      call getpr(ipr)
  !$$$      nbasnw = nbas
  !$$$      do  10  ibas = 1, nbas
  !$$$        do  20  isop = 2, ngrp
  !$$$          ic = ips(ibas)
  !$$$          call dmpy(g(1,isop),3,1,bas(1,ibas),3,1,bast,3,1,3,1,3)
  !$$$          call daxpy(3,1.d0,ag(1,isop),1,bast,1)
  !$$$          do  30  jbas = 1, nbasnw
  !$$$            jc = ips(jbas)
  !$$$            if (ic .eq. jc) then
  !$$$              do  32  m = 1, 3
  !$$$                dbas(m) = bast(m)-bas(m,jbas)
  !$$$   32         continue
  !$$$              if (latvec(1,tol1,qlat,dbas)) goto 22
  !$$$            endif
  !$$$   30     continue
  !$$$          nbasnw = nbasnw+1
  !$$$          if (ipr .ge. 50) then
  !$$$            if (nbasnw .eq. nbas+1) write(stdo,304)
  !$$$            call asymop(g(1,isop),ag(1,isop),' ',sg)
  !$$$            call skpblb(sg,len(sg),i)
  !$$$            write(stdo,303) clabl(ic),ibas,nbasnw,sg(1:i+1)
  !$$$          endif
  !$$$          ips(nbasnw) = ic
  !$$$          call dcopy(3,bast,1,bas(1,nbasnw),1)
  !$$$   22     continue
  !$$$   20   continue
  !$$$   10 continue
  !$$$
  !$$$C --- Printout ---
  !$$$      if (nbasnw > nbas) then
  !$$$        if (ipr >= 10) then
  !$$$c          call awrit2('%N ADDBAS: The basis was enlarged from %i'//
  !$$$c     .    ' to %i sites%N         The additional sites are: %N',
  !$$$          write(stdo,"(/,' ADDBAS: The basis was enlarged from ',i0,
  !$$$     &     ' to ',i0,' sites',/,'         The additional sites are: ')") nbas,nbasnw
  !$$$          write(stdo,301) (clabl(ips(ibas)),(bas(i,ibas),i=1,3), ibas=nbas+1,nbasnw)
  !$$$          write(stdo,'(a)') ' '
  !$$$        endif
  !$$$        nbas = nbasnw
  !$$$      else
  !$$$        if (ipr>40) write(stdo,*)' ADDBAS: basis is already complete --- no sites added'
  !$$$      endif
  !$$$
  !$$$C --- Error whether atoms are sitting on same position ---
  !$$$      novlp = 0
  !$$$      do  40  ibas = 1, nbas
  !$$$        ic = ips(ibas)
  !$$$        do  50  jbas = 1, ibas-1
  !$$$          jc = ips(jbas)
  !$$$          do  52  m = 1, 3
  !$$$            dbas(m) = bas(m,ibas)-bas(m,jbas)
  !$$$   52     continue
  !$$$          if (latvec(1,tol1,qlat,dbas)) then
  !$$$            write(stdo,400) ibas,clabl(ic),(bas(m,ibas),m=1,3),
  !$$$     .      jbas,clabl(jc),(bas(m,jbas),m=1,3)
  !$$$            novlp = novlp+1
  !$$$          endif
  !$$$   50   continue
  !$$$   40 continue
  !$$$      if (novlp>0) then
  !$$$         write(aaa,"(' ADDBAS: basis has ',i0,' overlapping site(s)')")novlp
  !$$$         call rx(aaa)
  !$$$      endif
  !$$$c     call fexit(-1,111,' Exit -1 ADDBAS: basis has %i overlapping site(s)',novlp)
  !$$$  301 format(8x,'ATOM=',a4,1x,'POS=',3f12.7)
  !$$$  303 format(10x,a4,2x,i3,2x,'-> ',i3,3x,5x,a)
  !$$$  304 format(/' ADDBAS: Spec   Atom  New_atom     Operation'/9x,35('-'))
  !$$$  400 format(/' ADDBAS: atom ',i3,', SPEC ',a4,' POS=',3f9.5,
  !$$$     ./'     and atom ',i3,', SPEC ',a4,' POS=',3f9.5,
  !$$$     ./'     are at the same positions. '/)
  !$$$      end
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine sgroup(mode,gen,agen,ngen,g,ag,ng,ngmx,qb)
    !      use m_lmfinit,only: stdo
    !- Sets up space group given generators (gen,agen).
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :a compound set of switches
    !i         :1s digit
    !i         :0 two groups compare to equal when both their point
    !i         :  and space parts compare equal
    !i         :1 two groups compare to equal when their point
    !i         :  group compares equal.  This eliminates
    !i         :  space groups that that have the same point group
    !i         :  but differing translational symmetry, which can
    !i         :  occur for artifically large supercells
    !i         :10s digit
    !i         :0 if ng>ngmx, abort with error message
    !i         :1 if ng>ngmx, return with ng=ngmx+1
    !i   gen   :rotation part of generators of the group
    !i   agen  :translation part of space group generator
    !i   ngen  :number of generators
    !i   ngmx  :maximum allowed number of group operations
    !i   qb    :vectors of a microcell in the Brillouin zone
    !o Outputs
    !o   g     :point group operations
    !o   ag    :translation part of space group
    !o   ng    :number of group operations
    !r Remarks
    !r   Operations are defined as (g,a)(v):=g*v+a
    !r   where g is a (3x3) matrix, a is a vector.
    !r   Always returns the identity operation as one group operation
    !u Updates
    !u   04 Jan 06 Added 10s digit mode
    !u   14 Mar 03 Added mode
    ! ----------------------------------------------------------------------
    implicit none
    integer:: mode,ngen,ng,ngmx
    integer:: ipr,igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,is,nnew,n,m,mode0,mode1
    real(8):: gen(9,ngen),g(9,ngmx),qb(3,3),agen(3,ngen),ag(3,ngmx), &
         h(9),hh(9),e(9),sig(9),asig(3),ah(3),ahh(3),ae(3)
    !      logical :: spgeql
    character:: sout*80,sg*35
    data e/1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/, ae/0d0,0d0,0d0/
    call getpr(ipr)
    sout = ' '
    call spgcop(e,ae,g,ag)
    ng = 1
    mode0 = mod(mode,10)
    mode1 = mod(mode/10,10)
    ! --- For each generator, do ---
    do  80  igen = 1, ngen
       call spgcop(gen(1,igen),agen(1,igen),sig,asig)
       ! --- Extend the group by all products with sig ----
       do  9  ig = 1, ng
          if (spgeql(mode0,g(1,ig),ag(1,ig),sig,asig,qb)) then
             if (ipr > 30) call awrit2(' Generator %i already in group '// &
                  'as element %i',' ',80,stdo,igen,ig)
             !        write(stdo,650) igen,ig
             !  650   format(' generator',i3,'  is already in group as element',i3)
             goto 80
          endif
9      enddo
       ! ... Determine order (= power of sig that restores unit operation)
       call spgcop(sig,asig,h,ah)
       do  1  itry = 1, 100
          iord = itry
          if (spgeql(mode0,h,ah,e,ae,qb)) goto 2
          call spgprd(sig,asig,h,ah,h,ah)
1      enddo
       ! ... Products of type  g1 sig**p g2
2      nnow = ng
       if (ipr >= 40) call awrit2('%a  %i is %i,',sout,80,0,igen,iord)
       do  8  j = 1, ng
          call spgcop(g(1,j),ag(1,j),h,ah)
          do  10  ip = 1, iord-1
             call spgprd(sig,asig,h,ah,h,ah)
             do  11  i = 1, ng
                call spgprd(g(1,i),ag(1,i),h,ah,hh,ahh)
                do  12  k = 1, nnow
                   if ( spgeql(mode0,g(1,k),ag(1,k),hh,ahh,qb) ) goto 11
12              enddo
                !         call asymop(hh,ahh,' ',sg)
                !         write(stdo,'('' sgroup adding'',i3,2x,a)') nnow+1,sg
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                call spgcop(hh,ahh,g(1,nnow),ag(1,nnow))
11           enddo
10        enddo
          if (j == 1) n2 = nnow
8      enddo
       ! ... Products with more than one sandwiched sigma-factor
       m1 = ng+1
       m2 = nnow
       do  20  is = 2, 50
          nnew = 0
          do 211 n = ng+1,n2
             do 21  m = m1, m2
                call spgprd(g(1,n),ag(1,n),g(1,m),ag(1,m),h,ah)
                do  22  k = 1, nnow
                   if (spgeql(mode0,g(1,k),ag(1,k),h,ah,qb)) goto 21
22              enddo
                nnew = nnew+1
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                call spgcop(h,ah,g(1,nnow),ag(1,nnow))
21           enddo
211       enddo
          m1 = m2+1
          m2 = nnow
          if (nnew == 0) goto 25
20     enddo
25     continue
       ng = nnow
80  enddo
    ! --- Printout ---
    if (ipr >= 30) then
       if (sout /= ' ' .AND. ipr >= 60) call awrit0 &
            (' Order of generator'//sout//'%a%b',' ',80,stdo)
       call awrit2(' SGROUP: %i symmetry operations from %i '// &
            'generators',' ',80,stdo,ng,ngen)
       if (ipr >= 60 .AND. ng > 1) then
          write(stdo,'('' ig  group op'')')
          do  60  ig = 1, ng
             call asymop(g(1,ig),ag(1,ig),' ',sg)
             write(stdo,'(i4,2x,a)') ig,sg
60        enddo
       endif
    endif
    return
99  continue
    if (mode1 == 0) call rx1( &
         'SGROUP: ng greater than ngmx=%i: probably bad translation',ngmx)
    !      call info2(1,0,0,
    !     .  ' SGROUP (warning) ng greater than ngmx=%i ... exiting',ngmx,0)
    ng = ngmx+1
  end subroutine sgroup

  subroutine spgprd(g1,a1,g2,a2,g,a)
    !     implicit none
    double precision :: &
         g1(3,3),g2(3,3),g(3,3),sum,a1(3),a2(3),a(3),h(3,3),ah(3)
    integer :: i,j,k
    do 101 i=1,3
       do 10 j=1,3
          sum=0d0
          do 11 k=1,3
             sum=sum+g1(i,k)*g2(k,j)
11        enddo
          h(i,j)=sum
10     enddo
101 enddo
    do j=1,3
       do i=1,3
          g(i,j)=h(i,j)
       enddo
    enddo
    do  i=1,3
       ah(i)=a1(i)
       do  j=1,3
          ah(i)=ah(i)+g1(i,j)*a2(j)
       enddo
    enddo
    do 14 i=1,3
       a(i)=ah(i)
14  enddo
    return
  end subroutine spgprd

  subroutine spgcop(g,ag,h,ah)
    integer :: i
    double precision :: h(9),g(9),ag(3),ah(3)
    do 10 i=1,9
       h(i)=g(i)
       if (dabs(h(i)) < 1.d-10) h(i)=0d0
10  enddo
    do 11 i=1,3
       ah(i)=ag(i)
       if (dabs(ah(i)) < 1.d-10) ah(i)=0d0
11  enddo
  end subroutine spgcop

  subroutine gpfndx(g,ag,ia,ja,pos,nrc,rb,qb)
    !- Finds atom ja which is transformed into ia by group operation g,ag.
    !     implicit none
    integer :: ia,ja
    double precision :: g(3,3),ag(3),pos(3,1),d(3),rb(3,3),qb(3,3)
    integer :: ka,nrc,m,k
    !     integer mode(3)
    !      mode(1) = 2
    !      mode(2) = 2
    !      mode(3) = 2
    ja = 0
    do  11  ka = 1, nrc
       do  m = 1, 3
          d(m) = ag(m) - pos(m,ia)
          do  k = 1, 3
             d(m) = d(m) + g(m,k)*pos(k,ka)
          enddo
       enddo
       call shorbz(d,d,rb,qb)
       if (abs(d(1))+abs(d(2))+abs(d(3)) < 1d-4) then
          ja = ka
          return
       endif
11  enddo
  end subroutine gpfndx

  subroutine groupg(mode,ng,g,ag,plat,ngen,gen,agen,gens,ngout)
    !- Finds a set of generators for the symmetry group
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   mode  :0 two groups compare to equal when both their point
    !i         :  and space parts compare equal
    !i         :1 two groups compare to equal when their point
    !i         :  group compares equal.  This eliminates
    !i         :  space groups that that have the same point group
    !i         :  but differing translational symmetry, which can
    !i         :  occur for artifically large supercells.
    !i   plat  :primitive translation vectors in real space
    !i   qlat  :primitive translation vectors in reciprocal space
    !i   g:symmetry operation symbol
    !i   ng:number of symmetry operations as supplied by the generators
    !o Outputs:
    !o   gen,ngen:generators, and number needed to produce g
    !o   ngout :number of group ops generated by (gen,ngen)
    !o         :usually ngout=ng unless artificial translations
    !o   gens  :ascii representation of generators
    !r Remarks:
    !r   The smallest set of generators is sought.
    !r   This subroutine performs the inverse function of sgroup.
    !u Updates
    !u   09 Jul 08 Extra check to find new generators beyond
    !u             the given ones.
    !u   12 May 07 Always returns gens, independent of verbosity
    !u   04 Jan 06 Returns ngout
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters:
    integer :: mode,ngen,ng,ngout
    double precision :: plat(9)
    double precision :: gen(9,*),agen(3,*),g(9,*),ag(3,*)
    character*(*) gens
    ! Local parameters:
    integer :: imax,isop,ngloc,ngmax,iprint,ngen0,ngmx
    integer :: i1,i2,j1,j2
    parameter (ngmx=48*64)
    character(100) :: sg,sg1,sout,sout2
    double precision :: gloc(3,3,ngmx),agloc(3,ngmx),qlat(9),xx,vec(3),vol

    ! --- Starting number of group ops ---
    call dinv33(plat,1,qlat,vol)
    call pshpr(1)
    ngen0 = ngen
    call sgroup(0,gen,agen,ngen,gloc,agloc,ngout,ngmx,qlat)

10  continue
    ! --- Do until enough generators added to make whole group ---
    if (ngout < ng) then
       !   ... Run through all symops, choosing whichever adds the most ops
       imax = 0
       ngmax = 0
       do  isop = 1, ng
          call dcopy(9,g(1,isop),1,gen(1,ngen+1),1)
          call dcopy(3,ag(1,isop),1,agen(1,ngen+1),1)
          !         call pshpr(61)
          call sgroup(mode,gen,agen,ngen+1,gloc,agloc,ngloc,ngmx,qlat)
          !         call poppr
          if (ngloc > ngmax) then
             imax = isop
             ngmax = ngloc
             ngout = ngloc
          endif
       enddo
       ngen = ngen+1
       call dcopy(9,g(1,imax),1,gen(1,ngen),1)
       call dcopy(3,ag(1,imax),1,agen(1,ngen),1)
       goto 10
    endif

    !     One last pass in case extra generators
    if ( .TRUE. ) then
       !   ... Run through all symops, choosing whichever adds the most ops
       imax = 0
       ngmax = ngout
       do  isop = 1, ng
          call dcopy(9,g(1,isop),1,gen(1,ngen+1),1)
          call dcopy(3,ag(1,isop),1,agen(1,ngen+1),1)
          !         call pshpr(61)
          call sgroup(mode,gen,agen,ngen+1,gloc,agloc,ngloc,ngmx,qlat)
          !         call poppr
          if (ngloc > ngmax) then
             imax = isop
             ngmax = ngloc
             ngout = ngloc
          endif
       enddo
       if (ngout > ngmax) then
          ngen = ngen+1
          call dcopy(9,g(1,imax),1,gen(1,ngen),1)
          call dcopy(3,ag(1,imax),1,agen(1,ngen),1)
       endif
    endif

    call poppr

    ! --- Create gens, optionally printout ---
    !     if (iprint() .ge. 20)  then
    if (ngen0 == 0) then
       call info0(20,0,0,' GROUPG: the following '// &
            'are sufficient to generate the space group:')
    else
       call info2(20,0,0,' GROUPG: %i generator(s) were added to '// &
            'complete the group%?#n#:',ngen-ngen0,ngen-ngen0)
    endif
    sout = ' '
    sout2 = ' '
    do  20  isop = 1, ngen
       call asymop(gen(1,isop),agen(1,isop),':',sg)
       call awrit0('%a '//sg,sout(9:),len(sout)-9,0)
       call dcopy(3,agen(1,isop),1,vec,1)
       call dgemm('N','N',1,3,3,1d0,agen(1,isop),1,qlat,3,0d0,vec,1)
       call asymop(gen(1,isop),vec,'::',sg1)
       call word(sg1,1,i1,i2)
       call shorbz(vec,vec,plat,qlat)
       call asymop(gen(1,isop),vec,'::',sg)
       call word(sg,1,j1,j2)
       if (i2-i1 < j2-j1) sg = sg1
       call awrit0('%a '//sg,sout2(9:),len(sout2)-9,0)
       ! rint *,'len(sout2)=',isop,len(sout2)
20  enddo
    if (ngen > ngen0 .AND. iprint() >= 20) then
       write(stdo,"(' Generator(cart): ', a)") trim(adjustl(sout))  ! call awrit0('%a',sout,len(sout),-stdo)
       write(stdo,"(' Generator(frac): ', a)") trim(adjustl(sout2)) ! call awrit0('%a',sout2,len(sout2),-stdo)
    endif
    gens = sout2
    !     endif
    if (ngout > ng) then
       call info2(20,0,0, &
            '%9f(warning) %i group ops supplied but generators create'// &
            ' %i ops',ng,ngout)
    endif
  end subroutine groupg
  !$$$!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  !$$$      subroutine fixpos(pos,nbas,tol,ng,plat,g,ag,istab)
  !$$$C- Adjusts site positions to agree with given symmmetry
  !$$$C  to machine precision
  !$$$C ----------------------------------------------------------------------
  !$$$Ci Inputs:
  !$$$Ci   pos:   basis vectors (scaled by alat)
  !$$$Ci   nbas:  number of atoms in the basis
  !$$$Ci   tol:   largest separation at which atoms are considered coincident
  !$$$Ci   ng:    number of symmetry operations
  !$$$Ci   plat:  primitive lattice vectors (scaled by alat)
  !$$$Ci   g,ag:  point and translation group operators
  !$$$Ci   istab: atom transformation table; see symtab
  !$$$Co Outputs:
  !$$$Co   pos:   basis vectors are adjusted.
  !$$$Cr Remarks:
  !$$$Cr   Generally atoms of the same class do not sit exactly on
  !$$$Cr   symmetry-related positions. In this subroutine each atomic
  !$$$Cr   position is replaced by the average of the position itself and
  !$$$Cr   the generated positions of the atoms of the same class.
  !$$$C ----------------------------------------------------------------------
  !$$$      implicit none
  !$$$C Passed parameters:
  !$$$      integer nbas,ng,istab(nbas,ng)
  !$$$      double precision pos(3,*),plat(*),g(9,*),ag(3,*),tol
  !$$$      integer ibas,jbas,m,ig,lgunit
  !$$$      double precision dbas(3),bast(3),sum,tol2,qlat(3,3),vol,ddot
  !$$$      double precision sdpos(3,nbas)
  !$$$      tol2 = 2*tol
  !$$$      call dpzero(sdpos,3*nbas)
  !$$$      sum = 0
  !$$$      call dinv33(plat,1,qlat,vol)
  !$$$      do  10  ibas = 1, nbas
  !$$$        do  20  ig = 1, ng
  !$$$          jbas = istab(ibas,ig)
  !$$$          call dmpy(g(1,ig),3,1,pos(1,ibas),3,1,bast,3,1,3,1,3)
  !$$$          do  30  m = 1, 3
  !$$$            dbas(m) = bast(m) + ag(m,ig) - pos(m,jbas)
  !$$$   30     continue
  !$$$  333     format(a,3f12.6)
  !$$$          call shorbz(dbas,dbas,plat,qlat)
  !$$$c         print 334, 'output dbas', dbas
  !$$$C     ... Debugging check
  !$$$          sum = sum + abs(dbas(1))+abs(dbas(2))+abs(dbas(3))
  !$$$          if (abs(dbas(1)) .gt. tol2 .or. abs(dbas(2)) .gt. tol2.or.
  !$$$     .    abs(dbas(3)) .gt. tol2) call fexit(-1,111,
  !$$$     .    'Exit -1 FIXPOS: positions incompatible with symgrp:'//
  !$$$     .    '  dpos=%d',max(dbas(1),dbas(2),dbas(3)))
  !$$$          if (abs(dbas(1)) .gt. tol .or. abs(dbas(2)) .gt. tol .or.
  !$$$     .    abs(dbas(3)) .gt. tol) call awrit4(
  !$$$     .    ' FIXPOS (warning): sites %i,%i incompatible '//
  !$$$     .    'with grp op %i:  dpos=%d',' ',80,lgunit(1),
  !$$$     .    ibas,jbas,ig,max(dbas(1),dbas(2),dbas(3)))
  !$$$  334     format(a,3f18.12)
  !$$$          call daxpy(3,1d0,dbas,1,sdpos(1,jbas),1)
  !$$$   20   continue
  !$$$   10 continue
  !$$$      sum = dsqrt(ddot(3*nbas,sdpos,1,sdpos,1)/3/nbas)
  !$$$      call daxpy(3*nbas,1d0/ng,sdpos,1,pos,1)
  !$$$      call awrit1(' FIXPOS: shifted site positions by average %;3g',' ',
  !$$$     .80,lgunit(1),sum/ng)
  !$$$      end
  !!
  subroutine grpgen(gen,ngen,symops,ng,ngmx)
    !- Generate all point symmetry operations from the generation group
    ! ----------------------------------------------------------------
    !i Inputs
    !i   gen,ngen,ngmx
    !i   ng  (12 Sep 96): if>0 , add symops to the ng already in list.
    !o Outputs
    !o   symops,ng
    !r Remarks
    !r   This works for point groups only and is set up for integer
    !r   generators.
    ! ----------------------------------------------------------------
    !     implicit none
    integer :: ngen,ng,ngmx
    double precision :: gen(9,ngen),symops(9,ngmx)
    double precision :: h(9),hh(9),e(9),sig(9),ae(3)
    integer :: igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,n,m
    integer :: ipr
    !      logical grpeql
    character(80) :: sout
    data e /1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/, ae/0d0,0d0,0d0/

    call getpr(ipr)
    sout = ' '
    call grpcop(e,symops)
    if (ng < 1) ng = 1
    do  80  igen = 1, ngen
       call grpcop(gen(1,igen),sig)
       ! ---   Extend the group by all products with sig ---
       do  9  ig = 1, ng
          if (grpeql(symops(1,ig),sig) .AND. ipr > 30) &
               call awrit2(' Generator %i already in group as element %i', &
               ' ',80,stdo,igen,ig)
          if (grpeql(symops(1,ig),sig)) goto 80
9      enddo

       ! ---   Determine order ---
       call grpcop(sig,h)
       do  1  itry = 1, 100
          iord = itry
          if (grpeql(h,e)) goto 2
          call grpprd(sig,h,h)
1      enddo
       ! --- Products of type  g1 sig**p g2 ---
2      nnow = ng
       if(ipr >= 40) call awrit2('%a  %i is %i,',sout,80,0,igen,iord)
       do  8  j = 1, ng
          call grpcop(symops(1,j),h)
          do  10  ip = 1, iord-1
             ! ... h = sig**ip
             call grpprd(sig,h,h)
             do  11  i = 1, ng
                ! ... hh = symops_i sig**ip
                call grpprd(symops(1,i),h,hh)
                do  12  k = 1, nnow
                   if ( grpeql(symops(1,k),hh) ) goto 11
12              enddo
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                call grpcop(hh,symops(1,nnow))
                !              print 333, (symops(k,nnow), k=1,9), nnow
                !  333         format(9f12.6,i3)
11           enddo
10        enddo
          if (j == 1) n2 = nnow
8      enddo

       ! --- Products with more than one sandwiched sigma-factor ---
       m1 = ng+1
       m2 = nnow
       do  20  i = 2, 50
          do  121  n = ng+1, n2
             do  21  m = m1, m2
                call grpprd(symops(1,n),symops(1,m),h)
                do  22  k = 1, nnow
                   if (grpeql(symops(1,k),h)) goto 21
22              enddo
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                call grpcop(h,symops(1,nnow))
21           enddo
121       enddo
          if (m2 == nnow) goto 25
          m1 = m2 + 1
          m2 = nnow
20     enddo
25     continue
       ng = nnow
80  enddo

    ! --- Printout ---
    if (ipr >= 30) then
       if (sout /= ' ' .AND. ipr >= 60) call awrit0 &
            (' Order of generator'//sout//'%a%b',' ',80,stdo)
       call awrit2(' GRPGEN: %i symmetry operations from %i '// &
            'generator(s)',' ',80,stdo,ng,ngen)
    endif
    if (ipr >= 80 .AND. ng > 1) then
       write(stdo,'('' ig  group op'')')
       do  60  ig = 1, ng
          call asymop(symops(1,ig),ae,' ',sout)
          write(stdo,'(i4,2x,a)') ig,sout(1:35)
60     enddo
    endif

    !      if (ipr .ge. 110) then
    !        print *, 'group operations:'
    !        call ywrm(0,' ',1,i1mach(2),'(5f12.6)',symops,1,9,9,ng)
    !      endif
    return
99  call rx('GRPGEN: too many elements')
  end subroutine grpgen
  subroutine grpcop(g,h)
    !- Copy matrix
    !     implicit none
    double precision :: h(9),g(9)
    integer :: i
    do  10  i = 1, 9
       h(i) = g(i)
10  enddo
  end subroutine grpcop

  subroutine grpprd(g1,g2,g1xg2)
    !- Returns the product of two point group operations
    !     implicit none
    double precision :: g1(3,3),g2(3,3),g1xg2(3,3),h(3,3),sum
    integer :: i,j,k
    do   i = 1, 3
       do    j = 1, 3
          sum = 0d0
          do  11  k = 1, 3
             sum = sum + g1(i,k)*g2(k,j)
11        enddo
          h(i,j) = sum
       enddo
    enddo
    do   j = 1, 3
       do   i = 1, 3
          g1xg2(i,j) = h(i,j)
       enddo
    enddo
  end subroutine grpprd
  !! -------------------------------------
  subroutine psymop(t,plat,g,ag,ng)
    !- Parse symbolic representation of symmetry group operations
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   t,nt  string of symmetry operations, separated by spaces
    !i   plat  lattice vectors that scale translation part ag
    !i         (if needed, i.e. if translation specified by '::')
    !o Outputs:
    !o   g,ng  group op (3x3 matrix) for each input, and number
    !r Remarks:
    !r   Symbols have two parts, first the point group part, followed
    !r   By an optional translation.  The point group part has the form
    !r   O(nx,ny,nz) where O is one of M, I or Rj for mirror, inversion
    !r   and j-fold rotations, respectively, and nx,ny,nz are a triplet
    !r   of indices specifying the axis of operation.
    !r   (nx,ny,nz) is one of (1,0,0), (0,1,0), (0,0,1) and (1,1,1),
    !r   it can be abbreviated as x,X, y,Y, z,Z and d,D, respectively.
    !r   Also permissible are products, eg I*R4X.
    !r   The translation is also of the form (n1,n2,n3)
    !r   Example: the following input
    !r     R3D(0,0,0) Mx R2(1/2,sqrt(3)/2,0)(pi,0,0) my*i'
    !r   is nonsensical, but permissible and generates four group ops.
    !r   10 Jan 1997 now generates g=transpose of prior versions.
    ! ----------------------------------------------------------------------
    !     implicit none
    character*(*) t
    double precision :: plat(3,3),g(9,1),h(9),hh(9),ag(3,1),vec(3)
    integer :: nt,ng,i
    logical :: flgp
    character*1:: leftp='('
    ! --- Do until no more symbolic representation, do ---
    nt = len(t)
    ng = 0
    i = 0
90  call skipbl(t,nt,i)
    if (i >= nt) return
    ng = ng+1
    call parsop(t,i,g(1,ng))
    if (t(i+1:i+1) == '*') then
       i = i+1
       call parsop(t,i,h)
       call grpprd(g(1,ng),h,hh)
       !       call dmpy(g(1,ng),3,1,h,3,1,hh,3,1,3,3,3)
       !       call dvcpy(hh,1,g(1,ng),1,9)
       call dcopy(9,hh,1,g(1,ng),1)
    endif
    call dpzero(ag(1,ng),3)
    ! ... Compatibility with old :T(x,y,z)
    if (t(i+1:i+2) == ':T' .OR. t(i+1:i+2) == ':t') i=i+2
    ! ... Compatibility with ::(x,y,z)
    flgp = .false.
    if (t(i+1:i+2) == '::') then
       flgp = .true.
       i=i+2
    elseif (t(i+1:i+1) == ':') then
       i=i+1
    endif
    if (t(i+1:i+1) == leftp) then
       if ( .NOT. parsvc(-1,t,i,ag(1,ng))) &
            call fexit(-1,111,' Exit -1 PSYMOP: '// &
            'failed to parse translation, ig=%i',ng)
       if (flgp) then
          call dcopy(3,ag(1,ng),1,vec,1)
          call grpop(vec,ag(1,ng),plat,1)
          !         call dgemm('N','N',3,1,3,1d0,plat,3,vec,3,0d0,ag(1,ng),3)
       endif
    endif
    goto 90
  end subroutine psymop
  subroutine parsop(t,i,a)
    !- Parse string for a point group operator
    double precision :: v(3),sp,c,s,pi2,a(3,3),ddot
    character(1) :: t(0:*)
    !      logical parsvc
    integer :: i,j,k,nrot,iii
    pi2 = 8*datan(1d0)
    if (t(i) == 'r' .OR. t(i) == 'R') then
       i = i+1
       read(t(i),'(i1)',err=99) nrot
       i = i+1
       if ( .NOT. parsvc(-1,t,i,v)) goto 99
       sp = ddot(3,v,1,v,1)
       sp = 1d0/dsqrt(sp)
       do  14  k = 1, 3
          v(k) = v(k)*sp
14     enddo
       c = dcos(pi2/nrot)
       s = dsin(pi2/nrot)
       do  16  k = 1, 3
          do  15  j = 1, 3
             a(k,j) = (1-c)*v(j)*v(k)
15        enddo
          a(k,k) = a(k,k) + c
16     enddo
       a(2,1) = a(2,1) + s*v(3)
       a(1,3) = a(1,3) + s*v(2)
       a(3,2) = a(3,2) + s*v(1)
       a(1,2) = a(1,2) - s*v(3)
       a(3,1) = a(3,1) - s*v(2)
       a(2,3) = a(2,3) - s*v(1)
    else if (t(i) == 'm' .OR. t(i) == 'M') then
       i = i+1
       if ( .NOT. parsvc(-1,t,i,v)) goto 99
       sp = ddot(3,v,1,v,1)
       do  11  j = 1, 3
          do  12  k = 1, 3
             a(j,k) = -2.d0*v(k)*v(j)/sp
12        enddo
          a(j,j) = a(j,j) + 1d0
11     enddo
    else if (t(i) == 'i' .OR. t(i) == 'I') then
       i = i+1
       !       call dvcpy(0d0,0,a,1,9)
       !       call dvcpy(-1d0,0,a,4,3)
       call dpzero(a,9)
       a(1,1) = -1
       a(2,2) = -1
       a(3,3) = -1
    else if (t(i) == 'e' .OR. t(i) == 'E') then
       i = i+1
       !       call dvcpy(0d0,0,a,1,9)
       !       call dvcpy(-1d0,0,a,4,3)
       call dpzero(a,9)
       a(1,1) = 1
       a(2,2) = 1
       a(3,3) = 1
    else
       goto 99
    endif
    return
99  print *, 'PARSOP: parse error at ',(t(iii),iii = 0,i),'  ...'
    call fexit(-1,119,' ',0)
  end subroutine parsop
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine symlat(platcp,ngrp,grp,isym)
    !- Generates the (point) symmetry operations of the lattice
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   platcp:lattice vectors of most compact primitive unit cell
    !o Outputs:
    !o   ngrp  :number of allowed symmetry operations
    !o   grp   :symmetry operation matrix
    !o   isym  :index to lattice type, calculated from ngrp:
    !o          ngrp   isym    name
    !o                  0     shouldn't happen
    !o            2     1     triclinic
    !o            4     2     monoclinic
    !o            8     3     orthorhombic
    !o           16     4     tetragonal
    !o           12     5     rhombohedral
    !o           24     6     hexagonal
    !o           48     7     cubic
    !r Remarks:
    !r   symlat analyzes the primitive translations of the bravais
    !r   lattice in order to supply the symmetry operations of the lattice.
    !r   It gives the number ngrp of allowed operations as well as
    !r   these operations themselves.
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ngrp,isym
    double precision :: platcp(3,3),grp(9,*)
    integer :: i,iprint,ltmax,ll1,m,m1,m2,m3,mm,nrot(4)
    parameter(ltmax=3,ll1=ltmax*2+1)
    double precision :: platt(9),qlatcp(3,3),mat(9),vecg(3),vol
    logical :: lirr
    character(12) :: csym1(0:7)
    data nrot /2,3,4,6/
    data csym1 /'indefinite','triclinic','monoclinic','orthorhombic', &
         'tetragonal','rhombohedral','hexagonal','cubic'/
    mm(i,m) = ltmax-(mod(i,ll1**m)-mod(i,ll1**(m-1)))/ll1**(m-1)
    call dinv33(platcp,1,qlatcp,vol)
    ! --- Start out with E and I ---
    ngrp = 2
    call csymop(-1,grp(1,1),.false.,1,[0d0,0d0,0d0])
    call csymop(-1,grp(1,2), .true.,1,[0d0,0d0,0d0])
    ! --- Find all possible rotation axes ---
    do  10  i = 0, (ll1**3-1)/2-1
       m1 = mm(i,1)
       m2 = mm(i,2)
       m3 = mm(i,3)
       lirr = .true.
       do  12  m = 2, ll1
          lirr = lirr.and.(mod(m1,m).ne.0.or.mod(m2,m).ne.0.or. &
               mod(m3,m).ne.0)
12     enddo
       if (lirr) then
          do  14  m = 1, 3
             vecg(m) = m1*platcp(m,1) + m2*platcp(m,2) + m3*platcp(m,3)
14        enddo

          do  16  m = 1, 4
             !       ... Matrix for this symmetry operation
             call csymop(-1,mat,.false.,nrot(m),vecg)
             call grpprd(mat,platcp,platt)
             !           call dmpy(mat,3,1,platcp,3,1,platt,3,1,3,3,3)
             !       ... Add it and i*symop, if allowed
             if (latvec(3,1d-5,qlatcp,platt)) then
                call csymop(-1,grp(1,ngrp+1),.false.,nrot(m),vecg)
                call csymop(-1,grp(1,ngrp+2),.true. ,nrot(m),vecg)
                ngrp = ngrp+2
                if (m /= 1) then
                   call csymop(-1,grp(1,ngrp+1),.false.,-nrot(m),vecg)
                   call csymop(-1,grp(1,ngrp+2),.true. ,-nrot(m),vecg)
                   ngrp = ngrp+2
                endif
             endif
16        enddo
       endif
10  enddo
    isym = 0
    if (ngrp == 2) isym=1
    if (ngrp == 4) isym=2
    if (ngrp == 8) isym=3
    if (ngrp == 16) isym=4
    if (ngrp == 12) isym=5
    if (ngrp == 24) isym=6
    if (ngrp == 48) isym=7
    if (iprint() >= 30) call awrit1(' SYMLAT: Bravais system is ' &
         //csym1(isym)//'%a with %i symmetry operations.',' ',80,stdo,ngrp)
  end subroutine symlat
  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
  subroutine symcry(bas,bast,ipc,nbas,nclass,nrclas, &
       ng,plat,qlat,g,ag,istab)
    !- Generates the symmetry ops of the crystal from those of the lattice
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   tol:   tol for which atoms are considered to be at the same site
    !i          use 0 for symtbl to pick internal default
    !i   bas   :basis vectors (scaled by alat)
    !i   ipc   :the jth atom belongs to class ipc(j)
    !i   nbas  :number of atoms in the basis
    !i   nclass:number of classes, atoms in same class are symmetry-related
    !i   nrclas:number of atoms in the ith class
    !i   plat  :primitive lattice vectors (scaled by alat)
    !i   qlat  :primitive translation vectors in reciprocal space
    !i   bast  :work array of same dimension as bas
    ! o Inputs/Outputs:
    ! o  ng    :number of allowed symmetry operations (see Remarks)
    ! o         on input  number of symmetry operations of the lattice
    ! o         on output number of symmetry operations of the crystal
    ! o  g     :symmetry operation matrices
    !o Outputs:
    !o   ag    :symmetry operation vector
    !o   istab :site ib is transformed into istab(ib,ig) by operation ig
    !r Remarks:
    !r   symcry finds the subset of the allowed ng point operations of the
    !r   lattice without a basis (see symlat.f) that are valid for
    !r   the crystal.
    !r
    !r   This routine is based on ATFTMT written by Worlton and Warren,
    !r   CPC 3, 88 (1972).
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters:
    integer :: nbas,ng,ipc(nbas),nclass,istab(nbas,ng),nrclas(nclass)
    double precision :: plat(9,*),qlat(*),bas(3,*),bast(3,*), &
         g(9,*),ag(3,*)
    !      real(8):: tol=0d0
    ! Local parameters:
    integer :: ibas,ic,iclbsj,icmin,ig,ipr,jbas,kbas,kc, &
         m,mbas,nj,nm,ng0
    !     integer mode(3)
    double precision :: dbas(3),tol0,tol1
    parameter (tol0=1d-5)
    !      logical latvec
    character sg*35
    call getpr(ipr)
    tol1 = fptol
    if(fptol==0d0) tol1 = tol0
    ! --- Find the class with minimum number of atoms ---
    icmin = 1
    do  5  ic = 1, nclass
       if (nrclas(ic) < nrclas(icmin) .AND. nrclas(ic) > 0) icmin = ic
5   enddo
    ibas = iclbsj(icmin,ipc,nbas,1)

    ! --- For each group op, see whether it only shifts basis by some T ---
    ng0 = ng
    ng = 0
    do  30  ig = 1, ng0
       !   ... Rotate the basis by g
       call dmpy(g(1,ig),3,1,bas,3,1,bast,3,1,3,nbas,3)
       do  20  nj = 1, nrclas(icmin)
          jbas = iclbsj(icmin,ipc,nbas,nj)
          !     ... This is a candidate for translation ag
          do  22  m = 1, 3
             ag(m,ng+1) = bas(m,jbas)-bast(m,ibas)
22        enddo
          call shorbz(ag(1,ng+1),ag(1,ng+1),plat,qlat)
          !          mode(1) = 2
          !          mode(2) = 2
          !          mode(3) = 2
          !          call shorps(1,plat,mode,ag(1,ng+1),ag(1,ng+1))
          !     ... See whether candidate works for all sites; also make istab
          do  10  kbas = 1, nbas
             kc = ipc(kbas)
             do  12  nm = 1, nrclas(kc)
                mbas = iclbsj(kc,ipc,nbas,nm)
                do  14  m = 1,3
                   dbas(m) = bas(m,mbas)-bast(m,kbas)-ag(m,ng+1)
14              enddo
                if (latvec(1,tol1,qlat,dbas)) then
                   istab(kbas,ng+1) = mbas
                   goto 10
                endif
12           enddo
             !       ... Candidate not valid
             if (ipr >= 90) then
                call asymop(g(1,ig),ag(1,ng+1),' ',sg)
                call awrit1(' symcry: excluded candidate ig=%,2i  '//sg &
                     //'%a',' ',80,stdo,ig)
             endif
             goto 20
10        enddo

          !     --- Valid ag found; add g to list ---
          ng = ng+1
          if (ig > ng) call dcopy(9,g(1,ig),1,g(1,ng),1)
          if (ipr >= 70) then
             call asymop(g(1,ng),ag(1,ng),' ',sg)
             call awrit1(' symcry: accepted candidate ig=%,2i  '//sg &
                  //'%a',' ',80,stdo,ig)
          endif
          goto 30
20     enddo
30  enddo
    if (ipr >= 30) call awrit2(' SYMCRY: crystal invariant under '// &
         '%i symmetry operations for tol=%;3g',' ',80,stdo,ng,tol1)
    if (ipr >= 60 .AND. ng > 1) then
       write(stdo,'('' ig  group op'')')
       do  60  ig = 1, ng
          call asymop(g(1,ig),ag(1,ig),' ',sg)
          write(stdo,'(i4,2x,a)') ig,sg
60     enddo
    endif
  end subroutine symcry
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine asymop(grp,ag,asep,sg)
    !- Generate the symbolic representation of a group operation
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i  grp,ag :  space group rotation + translation matrix
    !i  asep:     If first character is 'f', rot. and trans. vectors
    !i            are not converted into small algebraic expressions.
    !i            If first char. is 'f', The second-* character are
    !i            used for the separators.
    !o Outputs:
    !o  sg  :  symbolic representation of group op
    !b Bugs
    !b  No check is made on the length of sg
    ! ----------------------------------------------------------------------
    !     implicit none
    double precision :: grp(3,3),ag(3)
    character*(*) sg,asep
    ! Local variables
    double precision :: vecg(3),dasum,tiny
    integer :: nrot,ip,isw,awrite,i1,i2,fmtv
    logical :: li!,parsvc
    parameter(tiny=1d-4)
    ! --- Get consitutents of grp ---
    call csymop(1,grp,li,nrot,vecg)
    ! --- Rotational part ---
    i1 = 1
    fmtv = 0
    if (asep(1:1) == 'f') then
       fmtv = 4
       i1 = 2
    endif
    sg = ' '
    if (nrot == 1) then
       sg = 'i*i'
       ip = 3
       if (li) sg = 'i'
       if (li) ip = 1
    else
       if (li .AND. nrot == 2) then
          sg = 'm'
          ip = 1
       else
          !          ip = awrite('%?#n#i*##r%i',sg,len(sg),0,isw(li),nrot,0,0,0,0,0,0)
          if(  li   ) sg='i*r'//char(48+nrot)
          if( .NOT. li) sg='r'//char(48+nrot)
          ip =len(trim(sg))
       endif
       call rxx(.not. parsvc(2+fmtv,sg,ip,vecg),'bug in asymop')
    endif
    ! --- Translational part ---
    if (dasum(3,ag,1) > tiny) then
       if (asep(i1:i1) /= ' ') then
          call nword(asep,1,i1,i2)
          sg(ip+1:) = asep(i1:i2)
          ip = ip+i2-i1+1
       endif
       call rxx(.not. parsvc(1+fmtv,sg,ip,ag),'bug in asymop')
    endif
  end subroutine asymop
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine csymop(iopt,grp,li,nrot,vecg)
    !- Decomposes a group operation into its consitutents, or vice-versa
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   iopt  := -1 to convert (nrot,vecg,li) to grp
    !i          =  1 to convert grp to to (nrot,vecg,li)
    !o Inputs/Outputs:
    ! o grp   :group operation matrix
    ! o li    :if T: inversion or rotoinversion
    ! o nrot  :rotation angle = 2*pi/nrot
    ! o vecg  :rotation axis
    !r Remarks
    !r   for nrot > 2 the matrix is non-symmetric and the rotation
    !r   axis can be calculated from the antisymmetric part.
    !r   For nrot = 2 this not possible.  However, the squared vector
    !r   components are given by:  mat(i,i) = 2 v_i * v_i - 1.
    !r   This is used for the largest component. The others are taken
    !r   from: mat(i,j) = 2 v_i * v_j for i ne j.  This way we also
    !r   get the right phases between the components.
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters:
    integer :: nrot,iopt
    double precision :: vecg(3),grp(3,3)
    logical :: li
    ! Local parameters:
    integer :: i,idamax,j,in
    double precision :: costbn,detop,ddet33,dnrm2,sinpb3,tiny,twopi,vfac, &
         wk(9),sintbn,omcos,ddot
    parameter(tiny=1d-5)
    ! External calls:
    external daxpy,dcopy,ddet33,dpzero,dnrm2,dscal,idamax,ddot

    twopi = 8*datan(1d0)
    ! --- Make grp from (nrot,vecg,li) ---
    if (iopt == -1) then
       call dpzero(grp,9)
       in = iabs(nrot)
       if (in <= 0 .OR. in == 5 .OR. in > 6) &
            call fexit(-1,111,'%N Exit -1 CSYMOP: '// &
            'abs(nrot) must 1,2,3,4 or 6, but is %i',in)
       if (in == 1) then
          call dcopy(3,1d0,0,grp,4)
       else
          sintbn = dsin(twopi/nrot)
          costbn = dcos(twopi/nrot)
          omcos  = 1d0-costbn
          call rxx(dnrm2(3,vecg,1).lt.tiny, &
               'CSYMOP: zero rotation vector')
          call dscal(3,1/sqrt(ddot(3,vecg,1,vecg,1)),vecg,1)
          grp(1,1) = omcos*vecg(1)*vecg(1) + costbn
          grp(1,2) = omcos*vecg(1)*vecg(2) - sintbn*vecg(3)
          grp(1,3) = omcos*vecg(1)*vecg(3) + sintbn*vecg(2)
          grp(2,1) = omcos*vecg(2)*vecg(1) + sintbn*vecg(3)
          grp(2,2) = omcos*vecg(2)*vecg(2) + costbn
          grp(2,3) = omcos*vecg(2)*vecg(3) - sintbn*vecg(1)
          grp(3,1) = omcos*vecg(3)*vecg(1) - sintbn*vecg(2)
          grp(3,2) = omcos*vecg(3)*vecg(2) + sintbn*vecg(1)
          grp(3,3) = omcos*vecg(3)*vecg(3) + costbn
       endif
       if (li) call dscal(9,-1d0,grp(1,1),1)
       ! --- Make (nrot,vecg,li) from grp ---
    else if (iopt == 1) then
       ! ... Require |determinant=1|
       call dinv33(grp,0,wk,detop)
       if (dabs(dabs(detop)-1.0d0) > tiny) &
            call fexit(-1,111,'%N Exit -1 ASYMOP: '// &
            'determinant of group op must be +/- 1, but is %d',detop)
       detop = dsign(1.d0,detop)
       !   ... li is T if to multiply by inversion
       li = detop .lt. 0d0
       !   ... Multiply operation grp with detop to guarantee pure rotation
       call dscal(9,detop,grp(1,1),1)
       !   --- Calculate rotation angle from the normalization of v ---
       !       sum_i grp(i,i) = sum_i (1-cos) v_i*v_i + 3*cos = 1 + 2 * cos
       !       costbn = -0.5d0
       !       call daxpy(3,0.5d0,grp(1,1),4,costbn,0)
       costbn = 0.5d0*(-1 + grp(1,1) + grp(2,2) + grp(3,3))

       if (dabs(costbn-1d0) < tiny) then
          nrot = 1
          call dpzero(vecg,3)
       else
          !     ... See Remarks
          nrot = idnint(twopi/dacos(dmax1(-1d0,costbn)))
          if (nrot == 2) then
             do  10  i = 1, 3
                vecg(i) = 0.5d0*(grp(i,i)+1.0d0)
10           enddo
             j = idamax(3,vecg,1)
             if (vecg(j) < 0d0) &
                  call fexit2(-1,111,' Exit -1 ASYMOP:  bad component %i'// &
                  ' of operation.  Diagonal element = %d',j,grp(j,j))
             vecg(j) = dsqrt(vecg(j))
             vfac = 0.5d0/vecg(j)
             do  12  i = 1, 3
                if (i /= j) vecg(i) = vfac*grp(i,j)
12           enddo
          else
             vecg(1) = grp(3,2)-grp(2,3)
             vecg(2) = grp(1,3)-grp(3,1)
             vecg(3) = grp(2,1)-grp(1,2)
          endif

          !     --- Renormalize at least one component to 1 ---
          !         to allow for abbreviations as 'D', 'X', 'Y' or 'Z'
          sinpb3 = dsqrt(.75d0)
          if (dabs((sinpb3-dabs(vecg(1)))*(sinpb3-dabs(vecg(2)))* &
               (sinpb3-dabs(vecg(3)))) > tiny) then
             do  20  j = 3, 1,-1
                vfac = dabs(vecg(j))
                if(vfac > tiny) call dscal(3,1.d0/vfac,vecg,1)
20           enddo
          endif
       endif
       call dscal(9,detop,grp(1,1),1)
    endif
  end subroutine csymop
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine symtbl(mode,nbas,ipc,pos,g,ag,ng,qlat,istab)
    !- Make symmetry transformation table for posis atoms; check classes
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :1st digit
    !i         :0  site ib is transformed into istab(ib,ig) by grp op ig
    !i         :1  site istab(i,ig) is transformed into site i by grp op ig
    !i         :   (NB: old routine gpfndx used this convention)
    !i         :10s digit
    !i         :1  check atom classes
    !i   tol   :tol for which atoms are considered to be at the same site
    !i         :use 0 for symtbl to pick internal default
    !i   nbas  :size of basis
    !i   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
    !i   pos   :pos(i,j) are Cartesian coordinates of jth atom in basis
    !i   g     :point group operations
    !i   ag    :translation part of space group
    !i   ng    :number of group operations
    !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
    !o Outputs
    !o   istab :table of site permutations for each group op; see mode
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer :: nbas,ng,mode
    integer :: ipc(1),istab(nbas,1)
    double precision :: pos(3,1),g(9,1),ag(3,1),qlat(9)!,tol
    ! Local variables
    integer :: ib,ic,ig,jb,jc,mode1,mode10 !,oiwk
    double precision :: tol0,tol1
    character(200)::aaa
    parameter (tol0=1d-5)
    !      real(8):: tol=0d0
    integer,allocatable:: w_oiwk(:)
    if (ng == 0) return
    mode1 = mod(mode,10)
    mode10 = mod(mode/10,10)
    tol1 = fptol !tol
    if(fptol==0d0) tol1 = tol0
    !     --- Make atom transformation table ---
    do  20  ig = 1, ng
       do  10  ib = 1, nbas
          call grpfnd(tol1,g,ag,ig,pos,nbas,qlat,ib,jb)
          if (jb == 0) then
             !     .    call fexit2(-1,111,' Exit -1 SYMTBL: no map for atom '//
             !     .    'ib=%i, ig=%i',ib,ig)
             write(aaa,"('SYMTBL: no map for atom ib=',i0,' ig=',i0)") ib,ig
             call rx(aaa)
          endif
          if (mode10 /= 0) then
             ic = ipc(ib)
             jc = ipc(jb)
             if (ic /= jc) then
                write(aaa,"('SYMTBL: site ',i0,' not in same class as mapped site ',i0,', ig=',i0)") &
                     ib,jb,ig
                call rx(aaa)
                !               call fexit3(-1,111,' Exit -1 SYMTBL: '//
                !     .      'site %i not in same class as mapped site %i, ig=%i',
                !     .      ib,jb,ig)
             endif
          endif
          if (mode1 == 0) then
             istab(ib,ig) = jb
          else
             istab(jb,ig) = ib
          endif
10     enddo
20  enddo
    if (mode10 == 0) return
    ! --- Check atom classes ---
    allocate(w_oiwk(nbas))
    do  50  ib = 1, nbas
       ic = ipc(ib)
       w_oiwk=0
       do  30  ig = 1, ng
          w_oiwk(istab(ib,ig)) = 1
30     enddo
       do  40  jb = 1, nbas
          if (w_oiwk(jb) == 1) goto 40
          jc = ipc(jb)
          if (ic == jc) call fexit2(-1,111,' Exit -1 SYMTBL:  '// &
               'site ib=%i in same class as inequivalent site jb=%i',ib,jb)
40     enddo
50  enddo
    deallocate(w_oiwk)
  end subroutine symtbl

  subroutine istbpm(istab,nbas,ng,istab2)
    !- Makes inverse of istab
    integer :: nbas,ng
    integer :: istab(nbas,ng),istab2(nbas,ng)
    integer :: ib,ig,ibp
    do  ig = 1, ng
       do  ib = 1, nbas
          ibp = istab(ib,ig)
          istab2(ibp,ig) = ib
       enddo
    enddo
  end subroutine istbpm
  logical function parsvc(iopt,t,ip,v)
    !- Parses string, converting to a vector, or vice-versa
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   iopt   -1 to convert string t to vector,
    !i           1 to convert vector v to string t
    !i           2 same as 1, but use shorthand 'D', 'X', 'Y', 'Z'.
    !i   ip      position in t first char (0 for first char)
    !i Inputs/Outputs:
    ! o  t       string representation of vector (see Remarks)
    ! o  v       vector
    !o Outputs:
    !o   ip      position after last char: t(ip-1) is last char
    !o   parsvc  T if parse was successful; F if not.
    !r Remarks
    !r   The general string represntation has the form '(a,b,c)'
    !r   where a, b and c are expressions.
    !r   The following shorthand notations are allowed:
    !r   'D' for  (1,1,1)
    !r   'X', 'Y', 'Z' for (1,0,0), (0,1,0), (0,0,1).
    !b Bugs
    !b   no check is made on the length of t
    ! ----------------------------------------------------------------------
    !     implicit none
    integer :: ip
    double precision :: v(3)
    character*(1) t(0:*)
    ! Local variables
    double precision :: tiny,x,y,z,d
    character rchr*9, sout*50
    integer :: itrm,a2vec,awrite,ix(3),ich,iopt,m,i,iz,id
    logical :: lveq0(3),lveq1(3),a2bin
    parameter (tiny=1d-4)

    data rchr /'(XxYyZzDd'/

    ! --- Convert t to vec ---
    if (iopt == -1) then
       call dpzero(v,3)
       ich = 0
       call chrps2(t(ip),rchr,len(rchr),ich,ich,itrm)
       !   ... First char not in rchr: not a recognizable vector
       parsvc = .false.
       if (itrm == 0) return
       !   ... '('
       if (itrm == 1) then
          ip = ip+1
          ! this doesn't handle complicated cases
          !         parsvc = a2vec(t,ip+100,ip,4,',)',2,3,3,ix,v) .eq. 3
          if ( .NOT. a2bin(t,v,4,0,',',ip,-1)) return
          if ( .NOT. a2bin(t,v,4,1,',',ip,-1)) return
          if ( .NOT. a2bin(t,v,4,2,')',ip,-1)) return
          parsvc = .true.
          return
       else
          !         ... 'd'
          if (itrm >= 8) then
             v(1) = 1
             v(2) = 1
             v(3) = 1
             !         ... 'x' 'y' or 'z'
          else
             v(itrm/2) = 1
          endif
          parsvc = .true.
          ip = ip+1
       endif
       return

       ! --- Convert vec to t ---
    elseif (mod(iopt,4) == 1 .OR. mod(iopt,4) == 2) then
       parsvc = .true.
       do    i = 1, 3
          lveq0(i) = dabs(v(i)) .lt. tiny
          lveq1(i) = dabs(v(i)-1) .lt. tiny
       enddo
       !   ... Exclude shorthand
       if (mod(iopt,4) == 1) then
          lveq0(3) = .false.
          lveq1(3) = .false.
       endif
       if     (lveq1(1) .AND. lveq1(2) .AND. lveq1(3)) then
          t(ip) = 'd'
       elseif (lveq1(1) .AND. lveq0(2) .AND. lveq0(3)) then
          t(ip) = 'x'
       elseif (lveq0(1) .AND. lveq1(2) .AND. lveq0(3)) then
          t(ip) = 'y'
       elseif (lveq0(1) .AND. lveq0(2) .AND. lveq1(3)) then
          t(ip) = 'z'
       else
          t(ip) = '('
          do  12  i = 1, 3
             m = awrite('%x%;7d%?#n<>3#,#)#', &
                  sout,len(sout),0,v(i),i,i,i,i,i,i,i)
             x = v(i)
             y = v(i)*dsqrt(3d0)*4
             if (abs(x) > tiny .AND. dabs(dabs(x)-1) > tiny &
                  .AND. iopt <= 4) then
                if (abs(1/x-nint(1/x)) < tiny) then
                   if (x > 0) m = awrite('%x1/%d%?#n<>3#,#)#', &
                        sout,len(sout),0,1/x,i,i,i,i,i,i,i)
                   if (x < 0) m = awrite('%x-1/%d%?#n<>3#,#)#', &
                        sout,len(sout),0,-1/x,i,i,i,i,i,i,i)
                elseif (abs(y-nint(y)) < tiny .AND. &
                     abs(y) > .5d0) then
                   d = 12
                   iz = nint(abs(y))
                   id = 12
                   if (iz == 2) id = 6
                   if (iz == 4) id = 3
                   if (iz == 6) id = 2
                   if (iz == 12) id = 1
                   y = y/(12/id)
                   m = awrite('%x%?#n==1#%j#%d*#sqrt(3)/%i%?#n<>3#,#)#', &
                        sout,len(sout),0,nint(y),y,id,i,i,i,i,i)
                endif
             endif
             call strncp(t(ip+1),sout,1,1,m)
             ip = ip+m
12        enddo
       endif
       ip = ip+1
    endif
  end function parsvc
  !      subroutine fmain
  !      implicit none
  !      double precision v(3)
  !      integer ip
  !      logical parsvc,lsw
  !      character *50 t


  !      call dpzero(v,3)
  !      v(1) = 1
  !      ip = 1
  !      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
  !      lsw = parsvc(1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      call dpzero(v,3)
  !      v(2) = 1
  !      ip = 1
  !      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
  !      lsw = parsvc(1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      call dpzero(v,3)
  !      v(3) = 1
  !      ip = 1
  !      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
  !      lsw = parsvc(1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      call dpzero(v,3)
  !      v(1) = 1
  !      v(2) = 1
  !      v(3) = 1
  !      ip = 1
  !      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
  !      lsw = parsvc(1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      call dpzero(v,3)
  !      v(1) = 1d0/3
  !      v(2) = .25d0
  !      v(3) = .3d0
  !      ip = 1
  !      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
  !      lsw = parsvc(1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      call dpzero(v,3)
  !      v(1) = 1d0/3
  !      v(2) = .25d0
  !      v(3) = .3d0
  !      ip = 1
  !      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
  !      lsw = parsvc(5,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1;7d '//t//'%a '//t(ip:ip+1),
  !     .  ' ',80,6,lsw,ip,v)

  !      call dpzero(v,3)
  !      v(1) = 5*sqrt(3d0)/12
  !      v(2) = -sqrt(3d0)/2
  !      v(3) = sqrt(3d0)/6
  !      ip = 1
  !      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
  !      lsw = parsvc(1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1;7d '//t//'%a '//t(ip:ip+1),
  !     .  ' ',95,6,lsw,ip,v)

  !      t = ' xa'
  !      ip = 1
  !      lsw = parsvc(-1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      t = ' yb'
  !      ip = 1
  !      lsw = parsvc(-1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      t = ' zc'
  !      ip = 1
  !      lsw = parsvc(-1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      t = ' Dd'
  !      ip = 1
  !      lsw = parsvc(-1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      t = ' (.1,2-1,pi/4)e'
  !      ip = 1
  !      lsw = parsvc(-1,t,ip,v)
  !      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
  !     .  80,6,lsw,ip,v)

  !      end


  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine grpfnd(tol,g,ag,ig,pos,nbas,qlat,ia,ja)
    !- Find index to site ja site into which g,ag transforms site ia
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   tol   :tolerance in site positions
    !i   g     :rotation part of space group
    !i   ag    :translation part of space group
    !i   ig    :which group operation
    !i   pos   :basis vectors
    !i   nbas  :size of basis
    !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
    !i   ia    :site for which to find equivalent by group op (g,ag)
    !i         :If ia<0, absolute value of ia is used and additionally
    !i         :point group operation -g is used.
    !o Outputs
    !o   ja    :site that ia is transformed into by (g,ag)
    !o         :i.e. R(ja) = g(ig) R(ia) + ag(ig)
    !o         :if zero, no equivalent site was found.
    !u Updates
    !u   26 Jan 01  Add ability to operate with -g (ia<0)
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer :: ia,ja,ig
    double precision :: g(3,3,ig),ag(3,ig),pos(3,1),qlat(3,3),tol
    ! Local parameters
    double precision :: d(3),d2(3)
    !      logical latvec
    integer :: ka,nbas,m,k
    ka = iabs(ia)
    do    m = 1, 3
       d(m) = ag(m,ig)
       do    k = 1, 3
          d(m) = d(m) + g(m,k,ig)*pos(k,ka)
       enddo
    enddo
    if (ia < 0) call dscal(3,-1d0,d,1)

    ja = 0
    do  10  ka = 1, nbas
       d2(1) = d(1) - pos(1,ka)
       d2(2) = d(2) - pos(2,ka)
       d2(3) = d(3) - pos(3,ka)
       if (latvec(1,tol,qlat,d2)) then
          ja = ka
          return
       endif
10  enddo
  end subroutine grpfnd

  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  logical function spgeql(mode,g1,a1,g2,a2,qb)
    !- Determines whether space group op g1 is equal to g2
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 -> space group is compared
    !i         :1 -> only point group is compared
    !i   g1,a1 :first space group
    !i   g2,a2 :second space group
    !i   qb    :reciprocal lattice vectors
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    integer :: mode
    double precision :: g1(9),g2(9),a1(3),a2(3),qb(3,3)
    integer :: m,iq,iac
    double precision :: c,ca,dc
    spgeql=.true.
    do 10 m=1,9
       if (dabs(g1(m)-g2(m)) > 1.d-5) then
          spgeql=.false.
          return
       endif
10  enddo
    if (mode == 1) return
    do 20 iq=1,3
       c=(a1(1)-a2(1))*qb(1,iq)+(a1(2)-a2(2))*qb(2,iq) +(a1(3)-a2(3))*qb(3,iq)
       ca=dabs(c)
       iac=ca+0.5d0
       dc=ca-iac
       if (dabs(dc) > 1.d-5) then
          spgeql=.false.
          return
       endif
20  enddo
    return
  end function spgeql
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  logical function grpeql(g1,g2)
    !- Checks if G1 is equal to G2
    !     implicit none
    double precision :: g1(9),g2(9),dabs,tol,x1,x2
    parameter (tol = 1d-8)
    logical :: ddif
    integer :: i
    ddif(x1,x2) = dabs(x1-x2) .gt. tol
    grpeql = .false.
    do  10  i = 1, 9
       if (ddif(g1(i),g2(i))) return
10  enddo
    grpeql = .true.
  end function grpeql
  ! ssssssssssssssssssssssssssssssssssssssss
  SUBROUTINE GRPOP(V,V1,G,I)
    !     implicit none
    double precision :: G(3,3,*),V(3),V1(3)
    integer :: i
    V1(1) = G(1,1,I)*V(1) + G(1,2,I)*V(2) + G(1,3,I)*V(3)
    V1(2) = G(2,1,I)*V(1) + G(2,2,I)*V(2) + G(2,3,I)*V(3)
    V1(3) = G(3,1,I)*V(1) + G(3,2,I)*V(2) + G(3,3,I)*V(3)
    RETURN
  END SUBROUTINE GRPOP
  ! sssssssssssssssssssssssssssssssssssssssssssssss
  logical function latvec(n,tol,qlat,vec)
    !- Checks whether a set of vectors are lattice vectors
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   n     :number of vectors
    !i   tol   :tolerance
    !i   qlat  :primitive translation vectors in reciprocal space
    !i   vec   :double-precision vector
    !o Outputs:
    !o   latvec:T if all vectors are lattice vectors within spec'd tol
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters:
    integer :: n
    double precision :: qlat(3,3),vec(3,n),tol
    ! Local parameters:
    integer :: i,m
    double precision :: vdiff

    latvec = .false.
    do  10  i = 1, n
       do  20  m = 1, 3
          vdiff = vec(1,i)*qlat(1,m) + &
               vec(2,i)*qlat(2,m) + &
               vec(3,i)*qlat(3,m)
          if (dabs(vdiff-dnint(vdiff)) > tol) return
20     enddo
10  enddo
    latvec = .true.
  end function latvec

end module m_mksym_util

