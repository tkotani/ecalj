c$$$      subroutine clsset(mode,nbas,pos,ncl,clssl,cllst,clord,clp) !,slat
c$$$      use m_lattic,only: lat_plat
c$$$      use m_lmfinit,only: lat_plate


c$$$      use m_struc_def  !Cgetarg
c$$$
c$$$C- Groups all sites in a crystal into clusters
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   mode  :1s digit determines how clssl is made.
c$$$Ci           0  clssl is input
c$$$Ci           1  divide lattice into microcells of size plate
c$$$Ci          10s digit set rationalizes cluster index numbering
c$$$Ci   slat  :struct containing information about the lattice
c$$$Ci   nbas  :size of basis (input)
c$$$Ci   pos   :basis vectors (input)
c$$$Co Outputs
c$$$Co   ncl   :number of clusters
c$$$Co   clssl :cluster index: site ib belongs to cluster clssl(ib)
c$$$Co   cllst :list of sites belonging to each cluster.
c$$$Co          The range of sites is clp(1,ic)+1 ... clp(1,ic+1)
c$$$Co   clord :permutation table of indices that sorts basis by cluster
c$$$Co   clp   :CL-dependent indices.  clsset returns for each cluster:
c$$$Co          1: cumulative number of basis atoms in CL=1 ... this CL
c$$$Co          (the following are not made here)
c$$$Co          2: not used now
c$$$Co          3: lower matrix dimension for diagonal part of this CL
c$$$Co          4: downfolded matrix dimension for diagonal part of this CL
c$$$Co          5: matrix dimension for CL, incl. off-diagonal
c$$$Co          6: not used now
c$$$Cr Remarks
c$$$Cu Updates
c$$$C ----------------------------------------------------------------------
c$$$      implicit none
c$$$      integer nclp
c$$$      parameter (nclp=9)
c$$$      integer mode,nbas,clssl(nbas),cllst(nbas),clord(nbas),clp(nclp,1)
c$$$      real(8):: pos(3,nbas)
c$$$c      type(s_lat)::slat
c$$$
c$$$C Local variables
c$$$      logical lrat,pvclst
c$$$      integer:: lgunit , m , k , i , ib , ipr , ipose(3) , i1mach , 
c$$$     .ndig , ndigi , ifac , ncl , idmx , ngmx , ngmn , ncol , nskip 
c$$$     ., isw , ix(3) , j , i1 , i2
c$$$      real(8) ,allocatable :: pos_rv(:)
c$$$      real(8) ,allocatable :: pose_rv(:)
c$$$      integer ,allocatable :: w_iv(:)
c$$$
c$$$      double precision plat(3,3),plate(3,3),qlat(3,3),qlate(3,3)
c$$$      double precision vol,volm
c$$$      character*80 outs,fmt
c$$$      integer:: ocllst
c$$$      integer,allocatable :: ncltb_rv(:)
c$$$      real(8) ,allocatable :: clctr_rv(:)
c$$$
c$$$C heap:
c$$$      integer ::iwdummy
c$$$
c$$$
c$$$C ... Setup
c$$$
c$$$c      i_copy_size=size(lat_plat)
c$$$c      call dcopy(i_copy_size,lat_plat,1,plat,1)
c$$$      plat=lat_plat
c$$$      lrat = mod(mode/10,10) .ne. 0
c$$$      call getpr(ipr)
c$$$
c$$$C --- Divide total cell into microcells ---
c$$$      if (mod(mode,10) .eq. 0) then
c$$$      elseif (mod(mode,10) .eq. 1) then
c$$$
c$$$c        i_copy_size=size(lat_plate)
c$$$c        call dcopy(i_copy_size,lat_plate,1,plate,1)
c$$$         plate=lat_plate
c$$$
c$$$C   ... Determine supercell size
c$$$        call dinv33(plat,1,qlat,vol)
c$$$        call dinv33(plate,1,qlate,volm)
c$$$
c$$$        print 352
c$$$  352   format(/15x,'Plat',31x,'Mlat')
c$$$        print 350, ((plat(m,k),m=1,3),(plate(m,k),m=1,3),k=1,3)
c$$$  350   format(3f10.5,5x,3f10.5)
c$$$        call awrit3('   Cell vol: %,3d   Cluster vol: %;3d'//
c$$$     .  '  ratio: %;3g',
c$$$     .  ' ',80,lgunit(1),abs(vol),abs(volm),abs(vol/volm))
c$$$
c$$$        ndig = int(dlog(dble(i1mach(9)))/dlog(2d0)/3)
c$$$        ifac = 2**ndig
c$$$        if (ifac .eq. 1024) ifac = 1000
c$$$        allocate(pose_rv(3*nbas))
c$$$
c$$$        ix(1) = 1
c$$$        ix(2) = 1
c$$$        ix(3) = 1
c$$$        call shorps(nbas,plat,ix,pos,pos)
c$$$c       call prmx('shortened pos',pos,3,3,nbas)
c$$$        allocate(pos_rv(3*nbas))
c$$$
c$$$        call dcopy ( 3 * nbas , pos , 1 , pos_rv , 1 )
c$$$
c$$$C   ... Shift pos off axis by a tiny amount, avoiding boundary problems
c$$$*       call daxpy(3*nbas,1d0,1d-5,0,w(opos),1)
c$$$c       call prmx('plat',plat,3,3,3)
c$$$c       call prmx('plate',plate,3,3,3)
c$$$C   ... w(opos) = positions as multiples of plat
c$$$        call prjpos ( 2 , ix , iwdummy , plat , nbas , pos , pos_rv )
c$$$
c$$$C       call prmx('pos, units of plat',w(opos),3,3,nbas)
c$$$C   ... Positions as multiples of plate
c$$$c        call prjpos ( 63 , ix , plat , plate , nbas , pos_rv , pos_rve )
c$$$
c$$$C       call prmx('pos, units of plate',w(opose),3,3,nbas)
c$$$C   ... index as ifac**2*pose(1) + ifac*pose(2) + pose(3), integer parts
c$$$        do  20  ib = 1, nbas
c$$$          call discop ( pose_rv , ipose , 3 , 3 * ib - 2 , 1 , 0 )
c$$$
c$$$          clssl(ib) = ipose(1)*ifac*ifac + ipose(2)*ifac + ipose(3)
c$$$   20   continue
c$$$
c$$$C   ... Attempt to combine clusters
c$$$   30   continue
c$$$        call ngroup(nbas,clssl,clord,ncl,ngmn,ngmx)
c$$$C        if (ncl .lt. nbas) then
c$$$C          if (pvclst(nbas,w(opose),clssl,clord)) goto 30
c$$$C        endif
c$$$        if (allocated(pos_rv)) deallocate(pos_rv)
c$$$        if (allocated(pose_rv)) deallocate(pose_rv)
c$$$
c$$$      else
c$$$        call rx1('clsset: unrecognized mode %i',mode)
c$$$      endif
c$$$
c$$$C --- Count and sort cluster list ---
c$$$      call ngroup(nbas,clssl,clord,ncl,ngmn,ngmx)
c$$$
c$$$C --- Rationalize clusters and printout
c$$$C     call awrit2(' clsset: %n:1i',' ',255,lgunit(1),nbas,clssl)
c$$$      idmx = clssl(clord(nbas))
c$$$      ndigi = max(int(dlog10(dble(max(idmx,1))))+2,5)
c$$$      ndig = int(dlog10(dble(nbas))) + 2
c$$$      if (ipr .ge. 20)
c$$$     .call awrit7('%N CLSSET: %i clusters.'//
c$$$     .'%?#n#%b (id rationalized).##'//
c$$$     .'  Min, max, average sizes= %i, %i, %;1d'//
c$$$     .'%?#n>=30#%N Cluster%(n-4)fid   sites ...',
c$$$     .' ',120,lgunit(1),ncl,isw(lrat),ngmn,ngmx,dble(nbas)/ncl,
c$$$     .ipr,ndigi)
c$$$      i2 = 0
c$$$      outs = ' '
c$$$      clp(1,1) = 0
c$$$      do  40  i = 1, ncl
c$$$        i1 = i2+1
c$$$        call nxtgrp(nbas,clssl(clord(i1)),clssl,clord,i1,i2)
c$$$        clp(1,i+1) = clp(1,i) + i2-i1+1
c$$$        if (lrat) clssl(clord(i1)) = i
c$$$        if (ipr .ge. 30) call awrit4('  %,ni %,ni:',outs,len(outs),
c$$$     .  0,ndig,i,ndigi,clssl(clord(i1)))
c$$$        k = 1
c$$$        do  42  j = i1, i2
c$$$          if (lrat) clssl(clord(j)) = i
c$$$          if (ipr .ge. 30)
c$$$     .    call awrit2('%a%,ni',outs,len(outs),0,ndig+k,clord(j))
c$$$          k=0
c$$$   42   continue
c$$$        call awrit0(outs,' ',-80,lgunit(1))
c$$$   40 continue
c$$$
c$$$C ... Make cllst
c$$$      allocate(ncltb_rv(ncl+1))
c$$$      allocate(clctr_rv(3*nbas))
c$$$
c$$$      call clsctr ( nbas , pos , clssl , clord , ncltb_rv , cllst 
c$$$     ., clctr_rv )
c$$$
c$$$      if (allocated(clctr_rv)) deallocate(clctr_rv)
c$$$      if (allocated(ncltb_rv)) deallocate(ncltb_rv)
c$$$
c$$$
c$$$C ... Print re-ordering of sites
c$$$C     call awrit2(' clord: %n:1i',' ',255,lgunit(1),nbas,clord)
c$$$C     call awrit2(' clssl: %n:1i',' ',255,lgunit(1),nbas,clssl)
c$$$      allocate(w_iv(nbas))
c$$$
c$$$      call ivprm ( 1 , nbas , clssl , w_iv , clord , .false. )
c$$$
c$$$      call awrit1(' ... sites reordered by cluster%?#n>30#:',' ',80,
c$$$     .lgunit(1),ipr)
c$$$      if (ipr .gt. 30) then
c$$$        idmx = clssl(clord(nbas))
c$$$        ndigi = max(int(dlog10(dble(idmx)))+2,5)
c$$$        call awrit3(' %%,%ii%%,%ii%%,%ii',fmt,len(fmt),0,ndig+1,ndig+1,
c$$$     .  ndigi+1)
c$$$        call word(fmt,1,i,k)
c$$$        ncol = 84/(2*ndig+ndigi+7)
c$$$        nskip = nbas/ncol
c$$$        if (ncol*(nbas/ncol) .ne. nbas) nskip = nskip+1
c$$$        outs = ' Old New   id'
c$$$        call arrprt ( outs , fmt ( i:k ) , 'iIi' , nbas , nskip , ncol 
c$$$     .   , nskip , '  | ' , clord , iwdummy , w_iv , iwdummy , iwdummy 
c$$$     .   , iwdummy , iwdummy , iwdummy )
c$$$
c$$$
c$$$        if (allocated(w_iv)) deallocate(w_iv)
c$$$
c$$$
c$$$      endif
c$$$
c$$$C ... Re-order clssl
c$$$      call icopy ( nbas , w_iv , 1 , clssl , 1 )
c$$$
c$$$C     call awrit2(' clssl: %n:1i',' ',255,lgunit(1),nbas,clssl)
c$$$
c$$$      end subroutine clsset
c$$$
c$$$      logical function pvclst(nbas,pose,clssl,clord)
c$$$
c$$$C- Attempt to combine clusters
c$$$C     implicit none
c$$$      integer nbas,clssl(nbas),clord(nbas)
c$$$      double precision pose(3,nbas)
c$$$      integer i1,i2,i,il,ib,ic,lgunit
c$$$      integer j1,j2,j,jl,jb,jc
c$$$      double precision tol
c$$$      parameter (tol=1d-5)
c$$$
c$$$      pvclst = .false.
c$$$    2 i2 = 0
c$$$      do  10  il = 1, nbas
c$$$        i1 = i2+1
c$$$        call nxtgrp(nbas,clssl(clord(i1)),clssl,clord,i1,i2)
c$$$        if (i1 .gt. i2) return
c$$$
c$$$        j2 = i2
c$$$        do  20  jl = il, nbas
c$$$          j1 = j2+1
c$$$          call nxtgrp(nbas,clssl(clord(j1)),clssl,clord,j1,j2)
c$$$          if (j1 .gt. j2) goto 20
c$$$C     ... Check whether every member in i close to every member in j
c$$$          do  12  i = i1, i2
c$$$            ib = clord(i)
c$$$            do  22  j = j1, j2
c$$$              jb = clord(j)
c$$$              if (abs(pose(1,ib)-pose(1,jb)) .gt. 1+tol .or.
c$$$     .        abs(pose(2,ib)-pose(2,jb)) .gt. 1+tol .or.
c$$$     .        abs(pose(3,ib)-pose(3,jb)) .gt. 1+tol) goto 20
c$$$   22       continue
c$$$   12     continue
c$$$          ic = clssl(clord(i1))
c$$$          jc = clssl(clord(j1))
c$$$          pvclst = .true.
c$$$C     ... Merge groups i and j
c$$$          call awrit2(' merge id=%i with id=%i',' ',80,lgunit(1),ic,jc)
c$$$          do  14  i = i1, i2
c$$$            clssl(clord(i)) = jc
c$$$   14     continue
c$$$C     ... resort the group
c$$$          call ivshel(1,nbas,clssl,clord,.true.)
c$$$          do  15  i = 1, nbas
c$$$            clord(i) = clord(i)+1
c$$$   15     continue
c$$$          goto 2
c$$$   20   continue
c$$$
c$$$   10 continue
c$$$
c$$$      call rx('bug in pvclst')
c$$$
c$$$      end function pvclst
c$$$
c$$$      subroutine clsctr(nbas,bas,clssl,clord,ncltab,cllst,ctr)
c$$$
c$$$C- Makes cluster list and center of gravity
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   nbas  :size of basis (input)
c$$$Ci   bas   :basis vectors, in units of alat (input)
c$$$Ci   clssl :cluster index: site ib belongs to cluster clssl(ib) (clsset.f)
c$$$Co Outputs
c$$$Co   ncltab:ncltab(ic)=# sites in clsst preceding cluster ic
c$$$Co   cllst :cllst(ncltab(ic+1)..ncltab(ic+1)) are the list of
c$$$Co          sites in cluster ic.
c$$$Co   ctr   :positions of cluster centers, for each site
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer nbas,clssl(nbas),clord(nbas),ncltab(1),cllst(nbas)
c$$$      double precision bas(3,nbas),ctr(3,nbas)
c$$$C ... Local parameters
c$$$      integer ncl,ib,kb,icl,kcl,ncli,ix,lgunit,ipr
c$$$      double precision cg
c$$$
c$$$C --- Make cluster list cllst ---
c$$$      icl = clssl(clord(1))
c$$$      ncl = 0
c$$$      ncltab(1) = 0
c$$$      do  10  ib = 1, nbas
c$$$        kb = clord(ib)
c$$$        kcl = clssl(kb)
c$$$        cllst(ib) = kb
c$$$        if (kcl .ne. icl) then
c$$$          ncl = ncl+1
c$$$          ncltab(ncl+1) = ib-1
c$$$          icl = kcl
c$$$        endif
c$$$   10 continue
c$$$      ncl = ncl+1
c$$$      ncltab(ncl+1) = nbas
c$$$
c$$$C --- Make cluster center of gravity ---
c$$$      do  20  ix = 1, 3
c$$$      do  20  icl = 1, ncl
c$$$        cg = 0
c$$$        ncli = ncltab(icl+1) - ncltab(icl)
c$$$        cg = 0
c$$$        do  24  ib = ncltab(icl)+1, ncltab(icl+1)
c$$$          kb = cllst(ib)
c$$$          cg = cg + bas(ix,kb)
c$$$   24   continue
c$$$        cg = cg/ncli
c$$$        do  26  ib = ncltab(icl)+1, ncltab(icl+1)
c$$$          kb = cllst(ib)
c$$$          ctr(ix,kb) = cg
c$$$   26   continue
c$$$   20 continue
c$$$
c$$$C --- Printout ---
c$$$      call getpr(ipr)
c$$$      if (ipr .lt. 40) return
c$$$      if (ipr .gt. 40) then
c$$$        print '(1x)'
c$$$        print 332
c$$$  332   format(/'  site',9x,'------- pos --------',14x,
c$$$     .  '------- pos-cg --------')
c$$$      endif
c$$$      do  30  icl = 1, ncl
c$$$        ncli = ncltab(icl+1) - ncltab(icl)
c$$$c       print 332, icl, ncli, (ctr(ix,kb), ix=1,3)
c$$$        kb = cllst(ncltab(icl)+1)
c$$$        call awrit3(' cluster%,3i  has %i elements.  %33pcenter ='//
c$$$     .  ' %3;11,6D',' ',80,lgunit(1),icl,ncli,ctr(1,kb))
c$$$        if (ipr .gt. 40) then
c$$$          do  34  ib = ncltab(icl)+1, ncltab(icl+1)
c$$$            kb = cllst(ib)
c$$$            print 333, kb, (bas(ix,kb), ix=1,3),
c$$$     .      (bas(ix,kb)-ctr(ix,kb), ix=1,3)
c$$$  333       format(i4,2x,3f11.6,3x,3f11.6)
c$$$   34     continue
c$$$        endif
c$$$   30 continue
c$$$
c$$$      end subroutine clsctr
c$$$
c$$$
