!>Subroutines used in m_mksym
module m_mksym_util 
  use m_lgunit,only:stdo
  use m_ftox
  public mksym,mptauof,rotdlmm
  private
  real(8),parameter:: toll=1d-4,tiny=1d-4,epsr=1d-12
  integer,parameter:: ngnmx=10
  integer,parameter:: ngmx = 48
contains
  subroutine mksym(modeAddinversion,slabl,ssymgr,iv_a_oips, iclass,nclass,npgrp,nsgrp,rv_a_oag,rv_a_osymgr,iv_a_oics,iv_a_oistab)! Setup symmetry group. Split species into classes, Also assign class labels to each class
    use m_lmfinit,only: nbas,stdo,nspec           
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat,rv_a_opos
    use m_mpitk,only: master_mpi
    implicit none
    intent(in)::   modeAddinversion,slabl,ssymgr,iv_a_oips
    intent(out)::                                            iclass,nclass,npgrp,nsgrp,rv_a_oag,rv_a_osymgr,iv_a_oics,iv_a_oistab
    !i modeAddinversion  : 
    !i           =0  Not add inversion
    !i           =1  Add inversion to point group. Make additionally ag,istab for extra operations, using -g for rotation part; see Remarks
    !i slabl : species labels
    !i ssymgr: string containing symmetry group generators.
    !i           if ssymgr contains 'find', mksym will add basis atoms as
    !i           needed to guarantee generators are valid, and generate
    !i           internally any additonal group operations needed to
    !i           complete the space group.
    !r Remarks
    !r   In certain cases the inversion operation may be added to the space group, for purposes of k integration.  This is permissible when the
    !r   hamiltonian has the form h(-k) = h*(k).  In that case, the eigenvectors z(k) of h(k) are related to z(-k) as z(-k) = z*(k).
    !r
    !r   Also, the Green's functions are related G(-k) = Gtranspose(k). Thus if g is a space group operation rotating G0(g^-1 k) into G(k),
    !r   then G(-k) = Gtranspose(k), and the same (g,ag) information is needed for either rotation.
    integer :: modeAddinversion,nsgrp,npgrp,ibas,iwdummy1(1),idest,ig,iprint,igets,isym(10),j1,j2,lpgf,nclass,ngen, nggen,incli
    integer:: iv_a_oips(nbas),iclass(nbas),ifind, iv_a_oics(nbas),iv_a_oistab(ngmx*nbas)
    character(8) :: slabl(*),ssymgr*(*)
    character(120) :: gens
    real(8) :: gen(3,3,ngnmx), rv_a_oag(3,ngmx),rv_a_osymgr(3,3,ngmx)
    integer,allocatable ::  iv_a_onrc (:), iv_a_oipc(:) 
    logical:: symfind
    ifind = index(ssymgr,'find')
    gens = merge(ssymgr(1:ifind-1)//' '//ssymgr(ifind+4:),ssymgr, ifind>0)
    if(master_mpi) write(stdo,*)' Generators except find: ',trim(gens)
    symfind = ifind>0
    call gensym(slabl,gens,symfind,nbas,nspec,ngmx,plat,plat,rv_a_opos(:,1:nbas),iv_a_oips, & !Generate space group ops
         nsgrp, rv_a_osymgr,rv_a_oag, ngen,gen,ssymgr, nggen,isym,iv_a_oistab)
    if(nggen>ngmx) call rx('mksym: nggen>ngmx')
    incli = -1
    npgrp = nsgrp
    if(modeAddinversion /= 0) then !Add inversion to point group
       ngen = ngen+1
       gen(:,:,ngen) = reshape([-1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0],[3,3])
       call pshpr(iprint()-40)
       call grpgen(gen(1,1,ngen),1, rv_a_osymgr,npgrp, ngmx)
       call poppr
       incli = npgrp-nsgrp
    endif  ! Printout of symmetry operations !    if(master_mpi) write(stdo,ftox)'  mksym: found ',nsgrp,' space group operations'
    if(master_mpi.and.nsgrp/=npgrp) write(stdo,ftox) &
         '    adding inversion gives',npgrp,' operations for generating k points; enforce real for dmatu for LDA+U'
    if(master_mpi.and.incli == -1) write(stdo,*)'  no attempt to add inversion symmetry'
    allocate(iv_a_onrc(nspec))
    allocate(iv_a_oipc,source=iv_a_oips(1:nbas))
    call splcls(rv_a_opos,nbas,nsgrp,iv_a_oistab,nspec,slabl,nclass,iv_a_oipc,iv_a_oics,iv_a_onrc) !Split species into classes
    !                                                   ibas ==> iclass=ipc(ibas) ==> ispec=ics(iclass)
    !  allocate(iv_a_oistab(nsgrp*nbas))
    call symtbl(1, nbas, rv_a_opos , rv_a_osymgr, rv_a_oag, nsgrp, qlat, iv_a_oistab)
    iclass(1:nbas)=iv_a_oipc(1:nbas) 
  end subroutine mksym

  subroutine gensym(slabl,gens,symfind,nbas,nspec,ngmx,plat,platcv,bas,ips, ng,g,ag, ngen,gen,nwgens,nggen,isym,istab) !Generate the space group ops
    !i Inputs:
    !i   slabl: name of the different species.
    !i   gens:  a list of generators, in symbolic representation  NB: this list is not required; see Remarks.
    !i   symfind:=T Find any additional group operations for this basis.
    !i   nspec: number of classes, atoms in same class are symmetry-related
    !i   plat:  primitive lattice vectors (scaled by alat)
    !i   platcv:(=plat usually). Used to scale translation part of generators, when translation part specified as a multiple of
    !i         :lattice vectors.  Can be same as plat but primitive lattice vectors of "conventional unit cell"
    !i         :are sometimes used to specify these translations, e.g. when generated from spacegroup data in some books.
    !i  nbas:  number of atoms in the basis
    !i  bas:   basis vectors
    !i  ips:   the jth atom belongs to spec ips(j)
    !o Outputs:
    !o   istab: site ib is transformed into istab(ib,ig) by operation ig
    !o   g:     symmetry operation matrix (assumed dimensioned >=ngmx)
    !o   ag:    symmetry operation vector (assumed dimensioned >=ngmx)
    !o   ... The following are generated if symfind=T
    !o   isym:  numbers characterizing the symmetry of lattice and crystal. isym is index for underlying lattice (see symlat)
    !o   ng:    number of group operations
    !o   ngen:  number of symmetry generators
    !o   gen:   generators in matrix form
    !o   nwgens:generators in ascii form
    !o   nggen :number of group ops generated by generators.
    !o         :Usually nggen=ng; however nggen can exceed ng if supercell is artificial -> extra translations; see groupg
    !r Remarks:
    !r   gensym generates the space group, using the following prescription:
    !r     1.  Any generators supplied from input gens
    !r         are checked for consistency with the underlying lattice.
    !r     2.  The space group is made from these generators.
    !r     For symfind=T, we do follwoings
    !r       * The point group of the underlying lattice without the basis is generated.
    !r       * The full space group is generated from the point group
    !r       * A set of generators for this group is created
    !r   This program was adapted from the Stuttgart ASA version lmto-46.
    implicit none
    logical :: symfind
    integer :: nbas,isym(*),istab(nbas,*),nspec,ngen,ngmx, ng,nrspec(nspec),ips(nbas),nggen 
    integer:: i , j , ibas , ic , iprint, igen , mxint , ig,iwdummy1(1)
    real(8) :: plat(3,3),platcv(3,3),g(3,3,*),ag(3,*),bas(3,nbas), qlat(3,3),vol,platt(3,3),gen(3,3,ngmx),agen(3,ngmx)
    character(8) ::  slabl(*), gens*(*), nwgens*(*),xn
    character(100) :: sg,sg1,sout,sout2
    call dinv33(plat,1,qlat,vol) !Reciprocal lattice vectors --- 
    nwgens = gens
    ngen=0 !initialization buf obata added at 2023-5-12
    if(len_trim(gens) > 0) then
       call psymop(trim(gens),platcv,gen,agen,ngen)     ! Symmetry group as given by input generators ---
       nwgens = trim(gens)
    endif
    do igen = 1, ngen
       if(.NOT.latvec(3,toll,qlat,matmul(gen(:,:,igen),plat))) &
            call rx('GENSYM: imcompatible with lattice generators. igen='//trim(xn(igen)))
    enddo !write(6,*)'goto sgroup',ngen  !  do i=1,ngen     !       write(stdo,ftox)i,' gen agen=',ftof(gen(:,i),2),' ',ftof(agen(:,i),3)!  enddo
    call sgroup(gen,agen,ngen, g,ag,ng, ngmx,qlat) ! Set up space group (g,ag,ng) for  generators (gen,agen,ngen)
    !                  write(stdo,ftox)'end of sgroup ngen nggen ngmx',ngen,nggen,ngmx,ftof(ag(:,1),3)
    if (maxval(ips)/= nspec.AND.iprint() > 0)write(stdo,ftox)' GENSYM (warning)',nspec,'species supplied but only',i,'spec used ...'
    nspec = maxval(ips) 
    nrspec= [(count(ips(1:nbas)==ic),ic=1,nspec)] ! nrspec:number of atoms in the ith class
    if (symfind) then !Complete the space group. SYMGRP find. symfind=T
       call symlat(plat,ng,g,isym(1)) !Symmetry of lattice without bas
       call symcry(bas, ips,nbas,nspec,nrspec,ng,plat,qlat,g,ag,istab ) ! Symmetry of lattice with bas
       if (ng > ngmx) return
       nwgens = ' '
       groupgblock:block 
         !o Outputs: 
         !o   gen,ngen:generators, and number needed to produce g
         !o   nggen :number of group ops generated by (gen,ngen)
         !o         :usually nggen=ng unless artificial translations
         !o   nwgens  :ascii representation of generators
         integer :: imax,isop,ngloc,ngmax,iprint,ngen0,i1,i2,j1,j2,icount,ngmax2
         real(8)::gloc(3,3,ngmx),agloc(3,ngmx),xx,vec(3)
         call pshpr(1)
         ngen0 = ngen
         call sgroup(gen,agen,ngen, gloc,agloc,nggen, ngmx,qlat)
         icount=0
         ngmax2=0
         do ! a set of generators (gen,agen) to maximize number of spc operations
            icount=icount+1
            imax = 0
            ngmax = 0
            do  isop = 1, ng!   ... Run through all symops, choosing whichever adds the most ops
               gen(:,:,ngen+1)= g(:,:,isop) 
               agen(:,ngen+1) = ag(:,isop)  
               call sgroup(gen,agen,ngen+1, gloc,agloc,ngloc,ngmx,qlat)
               if (ngloc > ngmax) then
                  imax = isop
                  ngmax = ngloc
                  nggen = ngloc
               endif
            enddo
            if(ngmax>ngmax2) then
               ngmax2=ngmax
            else
               exit
            endif
            ngen = ngen+1
            gen(:,:,ngen)= g(:,:,imax) 
            agen(: ,ngen)= ag(: ,imax) 
            if(iprint()>0) write(stdo,ftox)'   Enlarging ngen=',ngen,' ng nggen=',ng,nggen
         enddo
         call poppr
         if(iprint()>0.and. ngen0 == 0) then
            write(stdo,ftox)' groupg: the following are sufficient to generate the space group:'
         elseif(iprint()>0) then
            write(stdo,ftox)' groupg:',ngen-ngen0,'generator(s) were added to complete the group:'
         endif
         sout = ''
         sout2 = ''
         do isop = 1, ngen
            call asymop(gen(1,1,isop),agen(1,isop),':',sg) !cartesian
            sout=trim(sout)//' '//trim(sg)
            call asymop(gen(1,1,isop),matmul(agen(1:3,isop),qlat(:,:)),'::',sg1) !fractional
            sout2=trim(sout2)//' '//trim(sg1) 
         enddo
         if (ngen > ngen0 .AND. iprint() >= 20) then
            write(stdo,"('  Generators:  trans(cart)= ', a)")trim(adjustl(sout)) 
            write(stdo,"('  Generators:: trans(frac)= ', a)")trim(adjustl(sout2))
         endif
         nwgens = sout2
         if(iprint()>0.and.nggen>ng) write(stdo,ftox)' (warning) Enlarged ng=',nggen,' > lattice ng=',ng
         ng =nggen
         g(:,:,1:ng)= gloc(:,:,1:ng)
         ag(:,1:ng) = agloc(:,1:ng)
       endblock groupgblock
       call symtbl(0,nbas,bas,g,ag,ng,qlat,istab )
    else
       call symtbl(0,nbas,bas,g,ag,ng,qlat,istab )
       nggen = ng
    endif
    if(iprint()>0) then
       write(stdo,"(' gensym: ig group ops (:vector means translation in cartesian)')")
       do  ig = 1, ng
          call asymop(g(:,:,ig),ag(1,ig),':',sg)
          write(stdo,'(i5,2x,a)') ig,trim(sg)
       enddo
       write(stdo,"(a)")' gensym: site permutation table for group operations ...'
       write(stdo,"('  ib/ig:',48i3)")  [(ig,ig=1,ng)]
       do i = 1, nbas
          write(stdo,"(i7,':',48i3)") i,(istab(i,ig), ig=1,ng)
       enddo
    endif
  end subroutine gensym
  subroutine sgroup(gen,agen,ngen, g,ag,ng, ngmx,qb)  !Get space group ops for given generators (gen,agen,ngen)
    !i Inputs
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
    !r   Operations are defined as (g,a)(v):=g*v+a where g is a (3x3) matrix, a is a vector.
    !r   Always returns the identity operation as one group operation
    implicit none
    integer:: ngen,ng,ngmx
    integer:: ipr,igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,is,nnew,n,m
    real(8):: gen(3,3,ngen),g(3,3,ngmx),qb(3,3),agen(3,ngen),ag(3,ngmx),h(3,3),hh(3,3),sig(3,3),asig(3),ah(3),ahh(3)
    character:: sout*80,sg*35
    real(8),parameter:: e(9)=[1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0],ae(3)=0d0
    character(8):: xt,xn
    call getpr(ipr)
    sout = ' '
    ag=0d0
    g(:,:,1)=reshape(e,[3,3])
    ng = 1
    igenloop: do  80  igen = 1, ngen !For each generator, do ---  write(stdo,ftox)'do80',igen,ftof(gen(1:9,igen),3),'  ',ftof(agen(:,igen),3)
       call spgcop(gen(:,:,igen),agen(:,igen),sig,asig)
       do ig = 1, ng ! --- Extend the group by all products with sig ----
          if (spgeql(g(:,:,ig),ag(:,ig),sig,asig,qb)) then
             if (ipr > 30) write(stdo,ftox)' Generator ',igen,' already in group as element ',ig
             cycle
          endif
       enddo
       call spgcop(sig,asig,h,ah) !Determine order (= power of sig that restores unit operation)
       do  itry = 1, 100
          iord = itry
          if (spgeql(h,ah,e,ae,qb)) exit
          call spgprd(sig,asig,h,ah,h,ah)
       enddo
       nnow = ng
       do     j = 1, ng
          call spgcop(g(:,:,j),ag(:,j),h,ah)
          do ip = 1, iord-1
             call spgprd(sig,asig,h,ah,h,ah)
             do  i = 1, ng
                call spgprd(g(:,:,i),ag(:,i),h,ah,hh,ahh) ! ... Products of type  g1 sig**p g2
                if(any([(spgeql(g(:,:,k),ag(:,k),hh,ahh,qb),k=1,nnow)])) cycle
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                call spgcop(hh,ahh,g(:,:,nnow),ag(:,nnow))
             enddo
          enddo
          if (j == 1) n2 = nnow
       enddo
       m1 = ng+1
       m2 = nnow
       do      is = 2, 50
          nnew = 0
          do    n = ng+1,n2
             do m = m1, m2
                call spgprd(g(:,:,n),ag(:,n),g(:,:,m),ag(:,m),h,ah) ! ... Products with more than one sandwiched sigma-factor
                if(any([(spgeql(g(:,:,k),ag(:,k),h,ah,qb),k=1,nnow)])) cycle
                nnew = nnew+1
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                call spgcop(h,ah,g(:,:,nnow),ag(:,nnow))
             enddo
          enddo
          m1 = m2+1
          m2 = nnow
          if (nnew == 0) exit
       enddo
       ng = nnow
80  enddo igenloop
    if(ipr >= 30) then
       if (sout /= ' ' .AND. ipr >= 60) write(stdo,ftox)' Order of generator'//trim(sout)
       write(stdo,ftox)' sgroup: ',ng,'symmetry operations from',ngen,'generators'
       if (ipr >= 60 .AND. ng > 1) then
          write(stdo,'('' sgroup: ig op'')')
          do  ig = 1, ng
             call asymop(g(1,1,ig),ag(1,ig),':',sg)
!             write(6,*)' g=',g(:,1,ig)
!             write(6,*)'   ',g(:,2,ig)
!             write(6,*)'   ',g(:,3,ig)
!             write(6,*)'ag=',ag(:,ig)
             write(stdo,'(5x,i4,2x,a)') ig,sg
          enddo
       endif
!       stop 'qqqqqqqqqqqqqqq'
    endif
    return
99  continue
    call rx('SGROUP: ng='//trim(xn(ng))//' > '//trim(xn(ngmx))//' probably bad translation')
  end subroutine sgroup
  subroutine grpgen(gen,ngen,symops,ng,ngmx) !Generate all point symmetry operations from the generation group
    use m_ftox
    !i   gen,ngen,ngmx
    !i   if ng>0, add symops to the ng already in list.
    !o Outputs: symops,ng
    !r Remarks   This works for point groups only and is set up for integer  generators.
    implicit none
    integer :: ngen,ng,ngmx
    real(8) :: gen(3,3,ngen),symops(3,3,ngmx), h(3,3),hh(3,3),sig(3,3)
    integer :: igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,n,m, ipr
    character(80) :: sout
    character(8):: xn
    real(8),parameter:: ee(9)=[1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],e(3,3)= reshape(ee,[3,3]),ae(3)=0d0
    call getpr(ipr)
    sout = ' '
    symops(:,:,1)= e
    if (ng < 1) ng = 1
    igenloop: do  80  igen = 1, ngen
       sig = gen(:,:,igen) !  Extend the group by all products with sig ---
       do ig = 1, ng 
          if (grpeql(symops(:,:,ig),sig) .AND. ipr > 30)  write(stdo,ftox)' Generator ',igen,' already in group as element',ig
          if (grpeql(symops(:,:,ig),sig)) goto 80
       enddo
       h=sig
       do  itry = 1, 100
          iord = itry
          if (grpeql(h,e)) exit
          h=matmul(sig,h) 
       enddo
       nnow = ng
       if(ipr >= 40) write(stdo,ftox) trim(sout),' ',igen,' is',iord
       do j = 1, ng !Products of type  g1 sig**p g2 ---
          h = symops(:,:,j) 
          do   ip = 1, iord-1  
             h = matmul(sig,h) ! h = sig**ip
             do i = 1, ng    
                hh = matmul(symops(:,:,i),h) ! hh = symops_i sig**ip
                if(any([(grpeql(symops(:,:,k),hh),k=1,nnow)])) cycle
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                symops(:,:,nnow)=hh 
             enddo
          enddo
          if (j == 1) n2 = nnow
       enddo
       m1 = ng+1
       m2 = nnow
       do       i = 2, 50 ! --- Products with more than one sandwiched sigma-factor ---
          do    n = ng+1, n2
             do m = m1, m2
                h= matmul(symops(:,:,n),symops(:,:,m)) 
                if(any([(grpeql(symops(:,:,k),h),k=1,nnow)])) cycle
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                symops(:,:,nnow)=h 
             enddo
          enddo
          if (m2 == nnow) exit
          m1 = m2 + 1
          m2 = nnow
       enddo
       ng = nnow
80  enddo igenloop
    if( ipr >= 30) then
       if(sout /= ' ' .AND. ipr >= 60) write(stdo,ftox)' Order of generator '//trim(sout)
       write(stdo,ftox)'GRPGEN:',ng,'symmetry operations from',ngen,'generator(s)'
       if(ipr >= 80 .AND. ng > 1) then
          write(stdo,'('' ig  group op'')')
          do  ig = 1, ng
             call asymop(symops(1,1,ig),ae,' ',sout)
             write(stdo,'(i4,2x,a)') ig,trim(sout)
          enddo
       endif
    endif
    return
99  continue
    call rx('GRPGEN: too many elements nnow ngmx='//trim(xn(nnow))//' '//trim(xn(ngmx)))
  end subroutine grpgen
  subroutine psymop(t,plat, g,ag,ng) !- Get symmetry group operations(g,ag,ng) from symbolic representation t
    !i Inputs:
    !i   t     string of symmetry operations, separated by spaces
    !i   plat  lattice vectors that scale translation part ag  (if needed, i.e. if translation specified by '::')
    !o Outputs:
    !o   g,ng  group op (3x3 matrix) for each input, and number
    !r Remarks:
    !r   Symbols have two parts, first the point group part, followed by an optional translation.  The point group part has the form
    !r   O(nx,ny,nz) where O is one of M, I or Rj for mirror, inversion and j-fold rotations, respectively, and nx,ny,nz are a triplet
    !r   of indices specifying the axis of operation.   (nx,ny,nz) is one of (1,0,0), (0,1,0), (0,0,1) and (1,1,1),
    !r   it can be abbreviated as x,X, y,Y, z,Z and d,D, respectively.  Also permissible are products, eg I*R4X.
    !r   The translation is also of the form (n1,n2,n3) :
    !r  CAUTION!: 2023: we do not allowe math operations in parenthesis. Give numerical number 6digits (if you like to recover, touch ctrl2ctrlp.py)
    implicit none
    character(*):: t
    real(8) :: plat(3,3),g(3,3,*),h(3,3),hh(3,3),ag(3,1),vec(3)
    integer :: nt,ng,i
    logical :: flgp
    character*1:: leftp='(' 
    nt = len(t)
    ng = 0
    i = 0
    do ! Do until no more symbolic representation, do ---
       call skipbl(t,nt,i)
       i=i+1 !  write(6,*)'  psymop:start###'//trim(t(i:)),' i=',i
       if (i >= nt) return
       ng = ng+1
       call parsop(t,i,g(1,1,ng)) !  write(6,*)'  Endof parsop1 i=',i,trim(t(i:))
       if (t(i:i) == '*') then
          i = i+1
          call parsop(t,i,h)
          g(:,:,ng)=matmul(g(1:3,1:3,ng),h) 
       endif
       ag(:,ng)=0d0 
       flgp = .false. 
       if (t(i:i+1) == '::') then !fractional
          flgp = .true.
          i=i+2
       elseif (t(i:i) == ':') then !cartesian 
          i=i+1
       endif
       if (t(i:i) == leftp) then
          if( .NOT. parsvc(t,i,ag(1:3,ng))) call rxi('psymop: failed to parse translation ig=',ng)
          if(flgp) ag(:,ng) = matmul(plat(:,:),ag(1:3,ng)) 
       endif
    enddo
  end subroutine psymop
  subroutine parsop(t,i,a)    !- Parse string for a point group operator
    real(8) :: v(3),sp,c,s,pi2,a(3,3),ddot
    character(*) :: t
    character(8):: xn
    real(8),parameter:: e(9)=[1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0]
    integer :: i,j,k,nrot,iii !    write(*,*)"parsopinput=@@@"//trim(t)//"@@@",i
    pi2 = 8*datan(1d0)
    if (t(i:i) == 'r' .OR. t(i:i) == 'R') then !rotation
       i = i+1
       read(t(i:i),'(i1)',err=99) nrot
       i = i+1
       if( .NOT. parsvc(t,i,v)) goto 99
       v(1:3) = v(1:3)/sum(v**2)**.5
       c = dcos(pi2/nrot)
       s = dsin(pi2/nrot)
       do  k = 1, 3
          a(k,1:3) = (1-c)*v(1:3)*v(k)
          a(k,k) = a(k,k) + c
       enddo
       a(2,1) = a(2,1) + s*v(3)
       a(1,3) = a(1,3) + s*v(2)
       a(3,2) = a(3,2) + s*v(1)
       a(1,2) = a(1,2) - s*v(3)
       a(3,1) = a(3,1) - s*v(2)
       a(2,3) = a(2,3) - s*v(1)
    elseif (t(i:i) == 'm' .OR. t(i:i) == 'M') then !mirror
       i = i+1
       if( .NOT. parsvc(t,i,v)) goto 99
       sp = sum(v**2)
       do  j = 1, 3
          a(j,:) = -2d0*v(:)*v(j)/sp
          a(j,j) = a(j,j) + 1d0
       enddo
    elseif (t(i:i) == 'i' .OR. t(i:i) == 'I') then !inversion
       i = i+1
       a=-reshape(e,[3,3])
    elseif (t(i:i) == 'e' .OR. t(i:i) == 'E') then !identity
       i = i+1
       a= reshape(e,[3,3])
    else
       goto 99
    endif
    return
99  continue
    call rx('PARSOP: parse error at '//trim(xn(i))//'t(i:)='//trim(t(i:)))
  end subroutine parsop
  subroutine symlat(platcp,ngrp,grp,isym)  !- Generates the point symmetry operations of the lattice
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
    !r   symlat analyzes the primitive translations of the bravais lattice in order to supply the symmetry operations of the lattice.
    !r   It gives the number ngrp of allowed operations as well as these operations themselves.
    implicit none
    integer :: ngrp,isym
    real(8) :: platcp(3,3),grp(9,*)
    integer :: i,iprint,ltmax,ll1,m,m1,m2,m3,mm
    parameter(ltmax=3,ll1=ltmax*2+1)
    real(8) :: platt(3,3),qlatcp(3,3),mat(3,3),vecg(3),vol
    logical :: lirr
    character(12),parameter:: csym1(0:7)=[character(12):: &
         'indefinite','triclinic','monoclinic','orthorhombic','tetragonal','rhombohedral','hexagonal','cubic']
    integer,parameter:: nrot(4)=[2,3,4,6],ngtab(7)=[2,4,8,16,12,24,48]
    mm(i,m) = ltmax-(mod(i,ll1**m)-mod(i,ll1**(m-1)))/ll1**(m-1)
    call dinv33(platcp,1,qlatcp,vol)
    ngrp = 2 ! --- Start out with E and I ---
    call csymop(grp(1,1),.false.,1,[0d0,0d0,0d0])
    call csymop(grp(1,2), .true.,1,[0d0,0d0,0d0])
    do  10  i = 0, (ll1**3-1)/2-1 !Find all possible rotation axes ---
       m1 = mm(i,1); m2 = mm(i,2); m3 = mm(i,3)
       if( all([((mod(m1,m)/=0.or.mod(m2,m)/=0.or.mod(m3,m)/=0),m=2,ll1)]) ) then
          vecg(:) = matmul(platcp,[m1,m2,m3]) 
          do  16  m = 1, 4 
             call csymop(mat,.false.,nrot(m),vecg) ! Matrix for this symmetry operation
             platt=matmul(mat,platcp)
             if (latvec(3,toll,qlatcp,platt)) then !       ... Add it and i*symop, if allowed
                call csymop(grp(1,ngrp+1),.false.,nrot(m),vecg)
                call csymop(grp(1,ngrp+2),.true. ,nrot(m),vecg)
                ngrp = ngrp+2
                if (m /= 1) then
                   call csymop(grp(1,ngrp+1),.false.,-nrot(m),vecg)
                   call csymop(grp(1,ngrp+2),.true. ,-nrot(m),vecg)
                   ngrp = ngrp+2
                endif
             endif
16        enddo
       endif
10  enddo
    isym = findloc(ngtab,value=ngrp,dim=1)
    if(iprint()>=30) write(stdo,ftox)' symlat: Bravais system is '//csym1(isym)//' with',ngrp,'symmetry operations.'
  end subroutine symlat
  subroutine symcry(bas,ipc,nbas,nclass,nrclas, ng,plat,qlat,g,ag,istab) ! Generates the symmetry ops of the crystal from those of the lattice
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
    ! io Inputs/Outputs:
    ! io  ng    :number of allowed symmetry operations (see Remarks)
    ! io         on input  number of symmetry operations of the lattice
    ! io         on output number of symmetry operations of the crystal
    ! io  g     :symmetry operation matrices
    !o Outputs:
    !o   ag    :symmetry operation vector
    !o   istab :site ib is transformed into istab(ib,ig) by operation ig
    !r Remarks:
    !r   symcry finds the subset of the allowed ng point operations of the
    !r   lattice without a basis (see symlat.f) that are valid for
    !r   the crystal.
    !r
    !r   This routine is based on ATFTMT written by Worlton and Warren,  CPC 3, 88 (1972).
    ! ----------------------------------------------------------------------
    implicit none
    integer :: nbas,ng,ipc(nbas),nclass,istab(nbas,ng),nrclas(nclass)
    real(8) :: plat(3,3),qlat(3,3),bas(3,nbas),bast(3,nbas), g(3,3,*),ag(3,*)
    integer :: ibas,ic,icmin,ig,ipr,jbas,kbas,kc, m,mbas,nj,nm,ng0
    real(8) :: dbas(3)
    character(135):: sg
    real(8):: rfrac(3)
    call getpr(ipr)
    icmin = 1
    do ic = 1, nclass !Find the class with minimum number of atoms ---
       if(nrclas(ic) < nrclas(icmin) .AND. nrclas(ic) > 0) icmin = ic
    enddo
    ibas = iclbsjx(ipc,nbas, icmin,1)
    ! --- For each group op, see whether it only shifts basis by some T ---
    ng0 = ng
    ng = 0
    do  30  ig = 1, ng0 !  write(stdo,ftox)'check initial generator ig=',ig
       bast= matmul(g(:,:,ig),bas)!   ... Rotate the basis by g
       do  20  nj = 1, nrclas(icmin)
          jbas = iclbsjx(ipc,nbas, icmin,nj)
          ag(:,ng+1) = bas(:,jbas)-bast(:,ibas) ! This is a candidate for translation ag
          rfrac = matmul(ag(:,ng+1)-epsr,qlat)
          ag(:,ng+1) = matmul(plat,rfrac -nint(rfrac)+epsr)
          do  10  kbas = 1, nbas
             kc      = ipc(kbas)
             do  nm  = 1, nrclas(kc)
                mbas = iclbsjx(ipc,nbas, kc,nm)
                dbas = bas(:,mbas)-bast(:,kbas)-ag(:,ng+1)
                if (latvec(1,toll,qlat,dbas)) then
                   istab(kbas,ng+1) = mbas
                   goto 10
                endif
             enddo
             call asymop(g(:,:,ig),ag(1,ng+1),' ',sg) !write(stdo,ftox)' symcry: excluded candidate ig=',ig,' op=',trim(sg) !Candidate not valid
             goto 20
10        enddo
          ng = ng+1 !  --- Valid ag found; add g to list ---
          if (ig > ng) g(:,:,ng)=g(:,:,ig) 
          exit
20     enddo
30  enddo
    if(ipr>=30)write(stdo,ftox)' symcry: crystal invariant under',ng,'following symmetry operations for tol=',ftof(toll)
!    if (ipr >= 60 .AND. ng > 1) then
!       write(stdo,'('' -- ig group op: symcry'')')
!       do  ig = 1, ng
!          call asymop(g(:,:,ig),ag(1,ig),' ',sg)
!          write(stdo,'(i5,2x,a)') ig,trim(sg)
!       enddo
!    endif
  end subroutine symcry
  subroutine asymop(grpin,ag,asep,sg)  ! Generate the symbolic representation sg of a group operation 
    !i  grpin,ag :  space group rotation + translation matrix
    !i  asep: 
    !o  sg  :  symbolic representation of group op
    implicit none
    real(8) :: grp(3,3),ag(3),vecg(3),grpin(3,3),costbn,detop,ddet33,dnrm2,sinpb3,vfac,wk(9)
    character(*):: sg,asep
    integer :: nrot,ip,isw,i1,i2,fmtv,llen,i,idamax,j,in
    logical :: li
    real(8),parameter:: twopi = 8*datan(1d0)
    character(8):: xn
    grp=grpin
    call dinv33(grp,0,wk,detop)
    if(dabs(dabs(detop)-1d0)>tiny) call rx('Exit -1 ASYMOP: determinant of group op must be +/- 1, but is '//trim(ftof(detop)))
    detop = dsign(1d0,detop) !sign of determinant
    li = detop<0d0 !   ... li is T if to multiply by inversion
    grp= detop*grp !Multiply operation grp with detop to guarantee pure rotation 
    costbn = 0.5d0*(-1 + grp(1,1) + grp(2,2) + grp(3,3))
    if (dabs(costbn-1d0) < tiny) then
       nrot = 1
       vecg=0d0 
    else
       nrot = idnint(twopi/dacos(dmax1(-1d0,costbn)))
       if (nrot == 2) then
          vecg = 0.5d0*[((grp(i,i)+1.0d0),i=1,3)]
          j = idamax(3,vecg,1)
          if(vecg(j) < 0d0)call rx('ASYMOP: bad operation j='//trim(xn(j))//'. Diagonal element is '//ftof(grp(j,j)))
          vecg(j) = dsqrt(vecg(j))
          vfac = 0.5d0/vecg(j)
          do i = 1, 3
             if (i /= j) vecg(i) = vfac*grp(i,j)
          enddo
       else
          vecg=[grp(3,2)-grp(2,3), grp(1,3)-grp(3,1), grp(2,1)-grp(1,2)]
       endif
       sinpb3 = dsqrt(.75d0) 
       if (dabs((sinpb3-dabs(vecg(1)))*(sinpb3-dabs(vecg(2)))*(sinpb3-dabs(vecg(3)))) > tiny) then
          do  j = 3, 1,-1 !Renormalize at least one component to 1 to allow for abbreviations as 'D', 'X', 'Y' or 'Z'
             vfac = dabs(vecg(j))
             if(vfac > tiny) vecg=1d0/vfac*vecg
          enddo
       endif
    endif
    sg=''
!    write(stdo,ftox)'nrotnnnnn',nrot,li,'vecg',vecg,'ag=',ag
    if(nrot == 1) then ! Rotational part
       sg = merge('i','e',li) 
       ip=len(trim(sg))+1
    else
       if(li.and.nrot==2) then
          sg='m'
       else   
          sg=merge('i*','  ',li)//'r'//char(48+nrot)
       endif
       ip=len(trim(sg))+1
       call rxx(.not. parsvc2(.true.,sg,ip,vecg),'bug in asymop 2')!rotation axis
    endif
    sg=adjustl(sg)
    if(sum(abs(ag))>tiny) then !Translational part added
!       print *,'sg=',sg,'asep=',trim(asep)
       if(asep(1:1)/=' ') sg=trim(sg)//trim(asep)
       ip=len(trim(sg))+1
       call rxx(.not. parsvc2(.false.,sg,ip,ag),'bug in asymop 1')
     endif  
  end subroutine asymop
  subroutine csymop(grp,li,nrot,vecg) !Convert (nrot,vecg,li) to grp
    use m_ftox
    ! i li    :if T: inversion or rotoinversion
    ! i nrot  :rotation angle = 2*pi/nrot
    ! i vecg  :rotation axis
    ! o grp   :group operation matrix
    !r Remarks
    !r   for nrot > 2 the matrix is non-symmetric and the rotation axis can be calculated from the antisymmetric part.
    !r   For nrot = 2 this not possible.  However, the squared vector components are given by:  mat(i,i) = 2 v_i * v_i - 1.
    !r   This is used for the largest component. The others are taken from: mat(i,j) = 2 v_i * v_j for i ne j.  This way we also
    !r   get the right phases between the components.
    implicit none
    integer :: nrot,iopt, i,idamax,j,in
    real(8) :: vecg(3),grp(3,3), costbn,detop,ddet33,dnrm2,sinpb3,vfac, wk(9),sintbn,omcos,ddot
    logical :: li
    character(8):: xt
    real(8),parameter:: twopi = 8*datan(1d0)
    in = iabs(nrot)
    if(in <= 0 .OR. in == 5 .OR. in > 6)call rx('CSYMOP: abs(nrot) must 1,2,3,4 or 6, but is '//trim(xt(in)))
    if(in == 1) then
       grp=reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[3,3])
    else
       call rxx(dnrm2(3,vecg,1).lt.tiny, 'CSYMOP: zero rotation vector')
       sintbn = sin(twopi/nrot)
       costbn = cos(twopi/nrot)
       omcos  = 1d0-costbn
       vecg= 1d0/sum(vecg**2)**.5 *vecg 
       grp(1,1:3) = omcos*vecg(1)*vecg(:) + [costbn,         -sintbn*vecg(3),  sintbn*vecg(2)]
       grp(2,1:3) = omcos*vecg(2)*vecg(:) + [sintbn*vecg(3),          costbn, -sintbn*vecg(1)]
       grp(3,1:3) = omcos*vecg(3)*vecg(:) + [-sintbn*vecg(2), sintbn*vecg(1),  costbn]
    endif
    if(li) grp=-grp 
  end subroutine csymop
  subroutine symtbl(mode,nbas,pos,g,ag,ng,qlat,istab) ! Make symmetry transformation table for posis atoms; check classes
    !i Inputs
    !i   mode  :1st digit
    !i         :0  site ib is transformed into istab(ib,ig) by grp op ig
    !i         :1  site istab(i,ig) is transformed into site i by grp op ig
    !i   nbas  :size of basis
    !i   pos   :pos(i,j) are Cartesian coordinates of jth atom in basis
    !i   g     :point group operations
    !i   ag    :translation part of space group
    !i   ng    :number of group operations
    !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
    !o Outputs  istab :table of site permutations for each group op; see mode
    implicit none
    integer :: nbas,ng,mode, istab(nbas,1),ib,ic,ig,jb,jc,ka
    real(8) :: pos(3,1),g(3,3,1),ag(3,1),qlat(9)
    character(8):: xt,xn
    do    ig = 1, ng !Make atom transformation table ---
       do ib = 1, nbas          
          jb=findloc( [(latvec(1,toll,qlat, matmul(g(:,:,ig),pos(:,ib))+ag(:,ig)-pos(:,ka)), ka=1,nbas)],dim=1,value=.true.)!ib is mapped to jb by g,ag 
          if(jb == 0) call rx("SYMTBL: no map for atom ib="//trim(xn(ib))//" ig="//trim(xn(ig)))
          if (mode == 0) then
             istab(ib,ig) = jb
          else
             istab(jb,ig) = ib
          endif
       enddo
    enddo
  end subroutine symtbl
  logical function parsvc2(modex,t,ip,v)  
    implicit none
    logical:: modex
    integer :: ip
    real(8) :: v(3)
    character(*) :: t
    real(8) :: x,y,z,d
    character sout*50, add*1,soutx*50
    integer :: itrm,ix(3),ich,iopt,m,i,iz,id,mx !,awrite !,a2vec
    character(9),parameter:: rchr='(XxYyZzDd'
    parsvc2 = .true.
    t(ip:ip)=' '
    if(modex) then
       if( all(abs(v(:)-[1d0,1d0,1d0])<tiny) ) t(ip:ip)='d'
       if( all(abs(v(:)-[1d0,0d0,0d0])<tiny) ) t(ip:ip)='x'
       if( all(abs(v(:)-[0d0,1d0,0d0])<tiny) ) t(ip:ip)='y'
       if( all(abs(v(:)-[0d0,0d0,1d0])<tiny) ) t(ip:ip)='z'
    endif
    if(t(ip:ip)==' ') then
       t(ip:ip)='('
       do i = 1, 3
          write(sout,ftox) ftof(v(i))
          if(abs(nint(v(i))-v(i))<1d-6) write(sout,ftox) nint(v(i))
          sout = trim(adjustl(sout))//merge(')',',',i==3)
          m= len_trim(sout)
          t(ip+1:ip+m)=trim(sout) 
          ip = ip+m
       enddo
    endif
    ip = ip+1
  end function parsvc2
  logical function parsvc(t,ip,v)  !Converting ascii t to vector v.  We handle only v=(num,num,num),d,x,y,z, where num is wo math operations.
    !io  ip     position in t first char.  out: position ofr last char
    !i  t       string representation of vector (see Remarks)
    !o  v       vector
    !o  parsvc  T if parse was successful; F if not.
    !r  The following shorthand notations are allowed: 'D'(1,1,1), 'X', 'Y', 'Z' for (1,0,0), (0,1,0), (0,0,1).
    !note: 2023 we removed math operation. If you need to recover, touch ctrlp2ctrl.py.
    implicit none
    integer :: itrm,ix(3),ich,iopt,m,i,iz,id,mx,ip,mmm
    real(8) :: v(3)
    character(*) :: t
    character(1):: rkey
    v=0d0
    rkey=t(ip:ip)
    parsvc = .false.
    if(scan(t(ip:ip),'(XxYyZzDd')==0) return !if no translation
    if (rkey=='(') then ! we expect t(ip:)='(..., ..., ...)
       ip = ip+1
       mmm=index(t(ip:),')')
       if(mmm<=0) call rx('cannot find right parensis for vector')
       if(verify(t(ip:ip+mmm-2),'0123456789.Dde-, ')>0) call rx('trans.part of symops should be numerical w/o math operations'&
            //trim(t(ip:ip+mmm-2)))
       read(t(ip:ip+mmm-2),*) v 
       ip=ip+mmm-2+2
    else
       if(rkey=='d' .or. rkey=='D') v = [1d0,1d0,1d0]
       if(rkey=='x' .or. rkey=='X') v = [1d0,0d0,0d0]
       if(rkey=='y' .or. rkey=='Y') v = [0d0,1d0,0d0]
       if(rkey=='z' .or. rkey=='Z') v = [0d0,0d0,1d0]
       ip = ip+1 
    endif
    parsvc = .true.
  end function parsvc
  subroutine splcls(bas,nbas,ng,istab,nspec,slabl,nclass,ipc, ics,nrclas) !- Splits species into classes
    use m_lgunit,only:stdo
    !i   bas,nbas: dimensionless basis vectors, and number
    !i   nspec:    number of species
    !io   ipc:      on input, site j belongs to species ipc(j)
    !i   slabl:    on input, slabl is species label
    !i   ng:       number of group operations
    !i   istab:    site ib is transformed into istab(ib,ig) by grp op ig
    !o  Outputs:
    !o   ipc:      site j belongs to class ipc(j)
    !o   ics:      class i belongs to species ics(i)
    !o   nclass:   number of classes
    !o   nrclas:   number of classes per each species
    implicit none
    logical :: nosplt
    integer :: nbas,nspec,nclass,ng,istab(nbas,ng),ipc(nbas), ics(nbas),nrclas(nspec)
    real(8) :: bas(3,*)
    character(8) :: slabl(*)
    integer :: ib,ic,icn,ig,jb,m,i,is,ipr,idx,ispec,j
    logical :: lyetno
    character(80) :: outs,clabl=''
    call getpr(ipr)
    nclass = nspec
    ics = [(i,i=1,nspec)]
    ic = 1
    do while(ic <= nclass) 
       is = ics(ic)
       ib = iclbsjx(ipc,nbas, ic,1)
       if (ib == 0) goto 11 !   ... No sites of this class ... skip
       lyetno = .true.
       do 20  jb = 1, nbas !For each basis atom in this class, do
          if (ipc(jb) == ic) then !class of jb
             if(  any(istab(ib,1:ng) == jb).or.&                    !If there is a g mapping ib->jb, sites are equivalent
                  any(istab(jb,1:ng) == ib).or.&                    !If there is a g mapping jb->ib, sites are equivalent
                  any([(istab(istab(ib,ig),ig)== jb,ig=1,ng)]).or.&   !If there is a g mapping ib->kb,jb, sites are equivalent
                  any([(istab(istab(jb,ig),ig)== ib,ig=1,ng)])) cycle !If there is a g mapping jb->kb,ib, sites are equivalent
             if (lyetno) then !If the classes haven't been split yet, do so
                nclass = nclass+1
                icn  =  nclass
                ics(icn) = is
                nrclas(is) = nrclas(is)+1
                lyetno = .false.
             endif
             if(nclass > nbas) call rx('splcls:  problem with istab')
             icn  =  nclass
             ipc(jb)=  icn !class index
          endif
20     enddo
11     continue
       ic = ic + 1
    enddo
    if(ipr>=30) then
       write(stdo,"(a)")' splcls:  ibas iclass ispec label(ispec)'
       do j=1,nbas
          ic   = ipc(j) !class
          ispec= ics(ic)!spec
          write(stdo,"(a,3i6,a)")"       ",j,ic,ispec,'     '//trim(slabl(ispec))
       enddo
    endif
  end subroutine splcls
  subroutine spgcop(g,ag,h,ah)
    real(8):: h(9),g(9),ag(3),ah(3)
    h = merge(0d0, g, dabs(g) <1d-8)
    ah= merge(0d0,ag, dabs(ag)<1d-8)
  end subroutine spgcop
  subroutine spgprd(g1,a1,g2,a2,g,a)
    implicit none
    real(8) :: g1(3,3),g2(3,3),g(3,3),sum,a1(3),a2(3),a(3),h(3,3),ah(3)
    integer :: i,j,k
    h=matmul(g1,g2)
    g=h !tk does not know why g=matmul(g1,g2) fails for gfortran gcc9.4.0 2023march
    ah=a1+matmul(g1,a2)
    a=ah
  end subroutine spgprd
  logical function spgeql(g1,a1,g2,a2,qb) ! Determines whether space group op g1 is equal to g2
    implicit none !i      g1,a1 :first space group,  g2,a2 :second space group, qb:reciprocal lattice vectors
    integer :: m,iq,iac
    real(8) :: g1(9),g2(9),a1(3),a2(3),qb(3,3),adiff(3)
    adiff = matmul(a1-a2,qb)
    spgeql= all([dabs(g1-g2),abs(adiff-nint(adiff))]<toll)
  end function spgeql
  logical function grpeql(g1,g2)    !- Checks if G1 is equal to G2
    implicit none
    real(8):: g1(9),g2(9)
    grpeql = all(dabs(g1-g2)<toll)
  end function grpeql
  logical function latvec(n,tol,qlat,vec) ! Checks whether a set of vec(1:3,n) are lattice vectors
    implicit none
    integer:: n
    real(8):: qlat(3,3),vec(3,n),tol, vdiff(n,3)
    vdiff  = matmul(transpose(vec(:,:)),qlat(:,:))
    latvec = all(reshape(abs(vdiff-nint(vdiff)),[n*3]) < tol)
  end function latvec
  integer function iclbsjx(ipc,nbas, ic,nrbas) !the nrbas-th atom belonging to class ic (ipc(ibas)==ic)
    implicit none
    integer :: ic,nbas,ipc(nbas),nrbas,ib,ibas
    iclbsjx = findloc( [(count([(ipc(ib)==ic,ib=1,ibas)]),ibas=1,nbas)], value=nrbas,dim=1)
  end function iclbsjx
  subroutine mptauof(symops,ng,plat,nbas,bas, iclass,miat,tiat,invg,delta,afmode) !- Mapping of atomic sites by points group operations.
    use  m_lmfinit,only: iantiferro
    !i  Input
    !i     symops(1,ng),ng,plat,nbas,bas(3,nbas)
    !i     iclass(nbas); denote class for each atom
    !o  Output
    !o    miat(ibas  ,ig); ibas-th atom is mapped to miat-th atom, by the ig-th
    !o    points group operation.  Origin is (0,0,0).
    !o    tiat(k,ibas,ig);
    !o    delta : shifting vector for non-symmorphic group.
    !o            r' = matmul (am, r) + delta
    !r  Remarks
    !r
    !r (1) The ibas-th atom (position at bas(k,ibas) ) is mapped to
    !r
    !r    bas( k,miat(ibas,ig) )+ tiat(k,ibas,ig), k=1~3.
    !r
    !r (2) tiat= unit translation
    implicit none
    integer :: ng,nbas, miat(nbas,ng),iclass(nbas),invg(ng), &
         nbmx, nsymx, ig,igd,i,j,ibas,mi,i1,i2,i3
    double precision :: SYMOPS(9,ng),plat(3,3), &
         tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas), &
         tr1,tr2,tr3,ep, dd1,dd2,dd3,t1,t2,t3
    integer::  iprintx=0
    integer :: ires(3, nbas, ng)
    integer:: ib1,ib2
    real(8) ::tran(3),delta(3,ng)
    logical:: cmdopt0
    logical,optional:: afmode
    data ep/1.0d-3/
    if(iprintx>=46) write(6,*)'MPTAUOf: search miat tiat for wave function rotation'
    do 10 ig=1,ng
       do igd=1,ng
          ! seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ) <= ep .AND. &
               abs( symops(2,ig)-symops(4,igd) ) <= ep .AND. &
               abs( symops(3,ig)-symops(7,igd) ) <= ep .AND. &
               abs( symops(4,ig)-symops(2,igd) ) <= ep .AND. &
               abs( symops(5,ig)-symops(5,igd) ) <= ep .AND. &
               abs( symops(6,ig)-symops(8,igd) ) <= ep .AND. &
               abs( symops(7,ig)-symops(3,igd) ) <= ep .AND. &
               abs( symops(8,ig)-symops(6,igd) ) <= ep .AND. &
               abs( symops(9,ig)-symops(9,igd) ) <= ep  ) then
             invg(ig)=igd
             goto 16
          endif
       enddo
16     continue
       do i=1,3
          do j=1,3
             am(i,j)=symops(i+3*(j-1),ig)
          enddo
       enddo
       do 120 ib1=1,nbas ! trial shift vector tran
          do 121 ib2=1,nbas
             tran =  bas(:,ib2)  - matmul(am,bas(:,ib1))
             if(present(afmode)) then
                if(iantiferro(ib1)==0) cycle
                if(iantiferro(ib2)==0) cycle
                if(iantiferro(ib1)+iantiferro(ib2)/=0) cycle
             endif
             do 30 ibas=1,nbas
                !bb1=matmul(am,bas(:,ibas))+trans
                b1=am(1,1)*bas(1,ibas)+am(1,2)*bas(2,ibas)+am(1,3)*bas(3,ibas) +tran(1)
                b2=am(2,1)*bas(1,ibas)+am(2,2)*bas(2,ibas)+am(2,3)*bas(3,ibas) +tran(2)
                b3=am(3,1)*bas(1,ibas)+am(3,2)*bas(2,ibas)+am(3,3)*bas(3,ibas) +tran(3)
                do 40 mi=1,nbas
                   if( iclass(mi) /= iclass(ibas) ) cycle
                   do  i1=-3,3
                      do  i2=-3,3
                         do  i3=-3,3
                            dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                            dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                            dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )
                            t1 = b1 - (bas(1,mi)+dd1)
                            t2 = b2 - (bas(2,mi)+dd2)
                            t3 = b3 - (bas(3,mi)+dd3)
                            if(abs(t1) <= ep .AND. abs(t2) <= ep .AND. abs(t3) <= ep) go to 60
                         enddo
                      enddo
                   enddo
40              enddo
                goto 121 ! seach failed, Not found mi and dd1. Try next (tr).
60              continue
                miat(ibas,ig)  = mi
                tiat(1,ibas,ig)= dd1
                tiat(2,ibas,ig)= dd2
                tiat(3,ibas,ig)= dd3
                ires(1,ibas,ig)= i1
                ires(2,ibas,ig)= i2
                ires(3,ibas,ig)= i3
30           enddo
             goto 21 ! When the do-30 loop has been completed, we get out of do-20 loop
121       enddo
120    enddo
       call rx('mptauof: Can not find miat and tiat')
21     continue
       delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000
       !- have gotten the translation-> check write --------------------
       if(iprintx >= 46) then
          write(6,4658)tran
4658      format('  Obtained translation operation=',3d12.4)
          do 123  ibas=1,nbas
             write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig), &
                  tiat(2,ibas,ig), tiat(3,ibas,ig), &
                  ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
150          format(' iiiiibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
123       enddo
       endif
10  enddo
  end subroutine mptauof

  subroutine rotdlmm(symops,ng,nl ,dlmm) ! Generate rotation matrix D^l_{m,m'} for L-representaiton,
    !  corresponding to points group operations.
    !i symops(9,ng),ng; point ops.
    !i nl; num.of l =lmax+1
    !o dlmm(2*nl-1,2*nl-1,0:nl-1,ng,2); D^l_{m,m'}. Indexes are for Real harmonics.
    !r  dlmmc is used as work area about 200kbyte used for  s,p,d,f -> nl=4
    !-----------------------------------------------------------------
    implicit double precision (a-h,o-z)
    integer:: is,i,ig,ikap,j,l,m,m1,m2,m3,md,mx,ix, ng,nl
    double precision :: SYMOPS(9,ng), am(3,3) ,fac1,fac2
    double precision :: dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
    double precision :: det,osq2
    complex(8):: msc(0:1,2,2), mcs(0:1,2,2),dum(2),&
         dlmmc(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
    complex(8),parameter:: Img=(0d0,1d0)
    integer:: debugmode
    real(8):: ep=1d-3 !ep was 1d-8 before feb2013
    do 10 ig =1,ng
       do  i=1,3
          do  j=1,3
             am(i,j) = symops(i+3*(j-1),ig)
          enddo
       enddo
       ! calculate determinant(signature)
       det= am(1,1)*am(2,2)*am(3,3) &
            -am(1,1)*am(3,2)*am(2,3) &
            -am(2,1)*am(1,2)*am(3,3) &
            +am(2,1)*am(3,2)*am(1,3) &
            +am(3,1)*am(1,2)*am(2,3) &
            -am(3,1)*am(2,2)*am(1,3)
       if(abs(abs(det)-1d0) >= 1d-10) then
          print *,' rotdlmm: det/=1 ig and det=',ig,det
          stop
       endif
       ! seek Euler angle   print *,' goto cbeta',ig,det
       cbeta = am(3,3)/det
       ! added region correction so as to go beyond domain error for functions, dsqrt and acos.
       if(abs(cbeta-1d0) <= 1d-6) cbeta= 1d0
       if(abs(cbeta+1d0) <= 1d-6) cbeta=-1d0
       beta = dacos(cbeta) ! beta= 0~pi
       sbeta= sin(beta)
       if(sbeta <= 1.0d-6) then
          calpha= 1d0
          salpha= 0d0
          alpha = 0d0
          cgamma= am(2,2)/det
          sgamma= am(2,1)/det
       else
          salpha =  am(2,3)/sbeta/det
          calpha =  am(1,3)/sbeta/det
          sgamma =  am(3,2)/sbeta/det
          cgamma = -am(3,1)/sbeta/det
       endif
       co2 = dcos(beta/2d0)
       so2 = dsin(beta/2d0)
       if(abs(calpha-1.0d0) <= 1.0d-6) calpha= 1.0d0
       if(abs(calpha+1.0d0) <= 1.0d-6) calpha=-1.0d0
       if(abs(cgamma-1.0d0) <= 1.0d-6) cgamma= 1.0d0
       if(abs(cgamma+1.0d0) <= 1.0d-6) cgamma=-1.0d0
       alpha=dacos(calpha)
       if(salpha < 0d0) alpha=-alpha
       gamma=dacos(cgamma)
       if(sgamma < 0d0) gamma=-gamma  !print *,'alpha beta gamma det=',alpha,beta,gamma,det
       do l =  0, nl-1
          do md= -l, l
             do m = -l, l
                !  from 'Ele theo. ang. mom. by M. E. Rose 5th 1967 Wisley and Sons.  p.52 (4.13)
                fac1 = dsqrt( igann(l+m)*igann(l-m)*igann(l+md)*igann(l-md) )
                fac2 = 0d0
                do ikap=0,2*l
                   if(l-md-ikap >= 0 .AND. l+m-ikap >= 0 &
                        .AND. ikap+md-m >= 0) then
                      add= dble((-1)**ikap)/( igann(l-md-ikap)*igann(l+m-ikap) &
                           *igann(ikap+md-m)*igann(ikap) )
                      if(2*l+m-md-2*ikap /= 0) add=add*co2**(2*l+m-md-2*ikap)
                      if(md-m+2*ikap /= 0)     add=add*(-so2)**(md-m+2*ikap)
                      fac2 = fac2+add
                   endif
                enddo
                ! l-th rep. is odd or even according to (det)**l
                dlmmc(md,m,l,ig) = fac1*fac2*det**l* cdexp( -Img*(alpha*md+gamma*m) )
             enddo
          enddo
       enddo
       am(1,1)= cos(beta)*cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)
       am(1,2)=-cos(beta)*cos(alpha)*sin(gamma)-sin(alpha)*cos(gamma)
       am(1,3)= sin(beta)*cos(alpha)
       am(2,1)= cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
       am(2,2)=-cos(beta)*sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)
       am(2,3)= sin(beta)*sin(alpha)
       am(3,1)=-sin(beta)*cos(gamma)
       am(3,2)= sin(beta)*sin(gamma)
       am(3,3)= cos(beta)
       if(abs(am(1,1)*det-symops(1,ig))>ep .OR. &
            abs(am(2,1)*det-symops(2,ig))>ep .OR. &
            abs(am(3,1)*det-symops(3,ig))>ep .OR. &
            abs(am(1,2)*det-symops(4,ig))>ep .OR. &
            abs(am(2,2)*det-symops(5,ig))>ep .OR. &
            abs(am(3,2)*det-symops(6,ig))>ep .OR. &
            abs(am(1,3)*det-symops(7,ig))>ep .OR. &
            abs(am(2,3)*det-symops(8,ig))>ep .OR. &
            abs(am(3,3)*det-symops(9,ig))>ep) then
          print *,' rotdlmm: not agree. symgrp and one by eular angle'
          stop
       endif
       if(debugmode()>9) then
          print *;print *;print *,' **** group ops no. ig=', ig
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
          print *,' by Eular angle '
          write(6,1731)am(1,1)*det,am(1,2)*det,am(1,3)*det
          write(6,1731)am(2,1)*det,am(2,2)*det,am(2,3)*det
          write(6,1731)am(3,1)*det,am(3,2)*det,am(3,3)*det
       endif
1731   format (' ',3f9.4)
10  enddo
    ! conversion to real rep. Belows are from csconvs
    !  msc mcs conversion matrix generation 2->m 1->-m for m>0
    osq2 = 1d0/sqrt(2d0)
    do m = 0,1
       Msc(m,1,:)= osq2*[complex(8):: (-1d0)**m, -Img*(-1d0)**m] !spherical to real(cubic)
       Msc(m,2,:)= osq2*[complex(8)::       1d0,            Img]
       Mcs(m,1,:)= osq2*[complex(8):: (-1d0)**m,      1d0]     !inverse
       Mcs(m,2,:)= osq2*[complex(8):: Img*(-1d0)**m, -Img]
    enddo
    converttoreal:do 123 is=1,ng ! convert to real rep.
       llooop:do 23   l =0,nl-1
          do  m2=-l,l
             do  m1= 1,l
                mx    = mod(m1,2)
                dum= [dlmmc(m2, m1,l,is), dlmmc(m2,-m1,l,is)]
                dlmmc(m2,  m1,l,is)= sum(dum(:)*msc(mx,:,1))
                dlmmc(m2, -m1,l,is)= sum(dum(:)*msc(mx,:,2))
             enddo
          enddo
          do m2=  1,l
             do m1= -l,l
                mx=mod(m2,2)
                dum= [dlmmc( m2, m1,l,is),dlmmc(-m2, m1,l,is)]
                dlmmc( m2, m1,l,is)= sum(mcs(mx,1,:)*dum(:))
                dlmmc(-m2, m1,l,is)= sum(mcs(mx,2,:)*dum(:))
             enddo
          enddo
          do m2=-l,l
             do m1=-l,l
                dlmm(m2,m1,l,is)=dreal( dlmmc(m2,m1,l,is) )
                if( abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12 ) &
                     call rx(' rotdlmm: abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12')
             enddo
          enddo
          if( .FALSE. ) then
             print *; print *,'  points ops  ig, l=', is,l,' cubic   '
             do m2=-l,l
                write(6,"(28f10.5)")( dreal(dlmmc (m2, m1,l,is) ), m1=-l,l)
                !    &    , ( dimag(dlmmc (m2, m1,l,is) ), m1=-l,l),( dlmm(m2, m1,l,is), m1=-l,l)
             enddo
          endif
23     enddo llooop
123 enddo converttoreal
    if(debugmode()>1) print *,' end of rotdlmm'
  end subroutine rotdlmm
  !--------------------------------------------
  double precision function igann(i)
    integer:: i,ix
    igann  = 1d0
    do ix =1,i
       igann=igann*dble(ix)
    enddo
  end function igann
end module m_mksym_util
