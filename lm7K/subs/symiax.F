Cdelw1 var ['oiaxib', '0', 'oiaxib', '-ldiaxb*nbas*niax', 'integer', 'iv_w_', 'symiax']
Cdelw1 var ['ontabb', '0', 'ontabb', 'nbas', 'integer', 'iv_w_', 'symiax']
      subroutine symiax(mode,plat,nbas,pos,g,ag,ng,ntab,iax,nttab,mxcsz)
C- Render a neighbor list compatible with symmetry operations
C ----------------------------------------------------------------------
Ci Inputs:
Ci   mode  :1s digit
Ci         :0 enlarge table to make it compatible with symops
Ci         :1 reduce  table to make it compatible with symops
Ci         :10s digit
Ci         :0 Add or subtract pairs found from the group operations
Ci         :  2..ng but do not require that each ri-rj has a
Ci         :  matching pair rj-ri.
Ci         :1 also guarantee that any ri-rj in the table has a
Ci         :  matching pair rj-ri
Ci         :100s digit
Ci         :0 do not sort table
Ci         :1 sort table by increasing length
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nbas  :size of basis
Ci   pos   :basis vectors
Ci         :size
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   ng    :number of group operations.
Ci         :Setting ng=0 and mode=1 has the effect of not checking
Ci         :for symmetry but just purging table of unused entries
Ci         :In this case, pos,g,ag are not used
Cio Inputs/Outputs
Cio  ntab  :number of pairs associated with each site.
Cio        :This array may be altered by the program
Cio  iax   : neighbor table containing pair information,
Cio        :in fixed-basis format (see iax2fd for conversion from
Cio        :standard form to fixed-basis format)
Cio        :This array may be altered by the program
Cio  nttab :On input, maximum allowed number of pairs in iax table
Cio        :(which depends on how it was dimensioned).
Cio        :On output, true number of pairs in iax table
Co Outputs
Cl   mxcsz :size of largest cluster
Cl   grow  :change in number of entries in iax table
Cr Remarks
Cr   To make table compatible with symmetry operations, extra
Cr   elements may be added to the existing table (mode 0), or elements
Cr   which don't have corresponding pairs in the table for each
Cr   rotation in the group are eliminated (mode 1)
Cu Updates
Cu   26 Sep 03 Bug fix: use grpfnd instead of gpfndx
Cu   29 Aug 02 First created
C ---------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer niax,nbas,ntab(nbas+1),ng,nttab,mode,mxcsz
      parameter (niax=10)
      integer iax(niax,nttab)
      double precision pos(3,*),g(9,*),ag(3,*),plat(3,3)
C ... Local parameters
      integer ldiaxb,ib,grow
Cdelw1       integer oiaxib,ontabb
      integer ,allocatable :: iv_w_oiaxib(:)
      integer ,allocatable :: iv_w_ontabb(:)

C ... Heap
      integer w(1)
      common /w/ w

C ... Copy iax to fixed-basis format
C     Assume no cluster in new iax table will not exceed 8x old one
      mxcsz = 0
      do  ib = 1, nbas
        mxcsz = max(mxcsz,ntab(ib+1)-ntab(ib))
      enddo
      ldiaxb = 8*mxcsz
Cdelw1       call defi(oiaxib,-ldiaxb*nbas*niax)
      allocate(iv_w_oiaxib(abs(-ldiaxb*nbas*niax)))
      if (-ldiaxb*nbas*niax<0) iv_w_oiaxib(:)=0

Cdelw1       call defi(ontabb,nbas)
      allocate(iv_w_ontabb(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_ontabb(:)=0

C     call yprm('iax before copy',0,iax,0,10,10,nttab)
Cdelw1       call iax2fd(0,nbas,ntab,iax,w(ontabb),w(oiaxib),ldiaxb)
      call iax2fd ( 0 , nbas , ntab , iax , iv_w_ontabb , iv_w_oiaxib 
     . , ldiaxb )


C ... Update the table
      if (ng .gt. 0) then
Cdelw1       call symia0(mode,plat,nbas,pos,g,ag,ng,w(ontabb),w(oiaxib),ldiaxb,
Cdelw1      .  grow)
        call symia0 ( mode , plat , nbas , pos , g , ag , ng , iv_w_ontabb 
     . , iv_w_oiaxib , ldiaxb , grow )

        ib = ntab(nbas+1)+grow
        if (ib .gt. nttab) call rxi(
     .  'symiax: not enough memory for new iax table: need dim',ib)
      endif

C ... Copy iax from fixed-basis format
Cdelw1       call iax2fd(1,nbas,ntab,iax,w(ontabb),w(oiaxib),nttab)
      call iax2fd ( 1 , nbas , ntab , iax , iv_w_ontabb , iv_w_oiaxib 
     . , nttab )

      mxcsz = 0
      do  ib = 1, nbas
        mxcsz = max(mxcsz,ntab(ib+1)-ntab(ib))
      enddo

C     debugging: second pass shouldn't change table
C      call iax2fd(0,nbas,ntab,iax,w(ontabb),w(oiaxib),ldiaxb)
C      call symia0(mode,plat,nbas,pos,g,ag,ng,w(ontabb),w(oiaxib),ldiaxb,
C     .  grow)
C      call iax2fd(1,nbas,ntab,iax,w(ontabb),w(oiaxib),nttab)

Cdelw1 rlse name= oiaxib old_list= oiaxib ontabb
Cdelw1 rlse name= oiaxib new_list= (None)
Cdelw1       call rlse(oiaxib)
      if (allocated(iv_w_ontabb)) deallocate(iv_w_ontabb)
      if (allocated(iv_w_oiaxib)) deallocate(iv_w_oiaxib)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

      subroutine symia0(mode,plat,nbas,pos,g,ag,ng,ntabib,iaxib,ldiaxb,
     .  grow)
C- Kernel for symiax. (Same as symiax but uses fixed-basis form for iax)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   mode  :1s digit
Ci         :0 enlarge table to make it compatible with symops
Ci         :1 reduce  table to make it compatible with symops
Ci         :10s digit
Ci         :0 Add or subtract pairs found from the group operations
Ci         :  2..ng but do not require that each ri-rj has a
Ci         :  matching pair rj-ri.
Ci         :1 also guarantee that any ri-rj in the table has a
Ci         :  matching pair rj-ri
Ci         :100s digit
Ci         :0 do not sort table
Ci         :1 sort table by increasing length
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nbas  :size of basis
Ci   pos   :basis vectors
Ci   ldiaxb:third dimension in iaxib, and maximum allowed cluster size
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   ng    :number of group operations
Cio Inputs/Outputs
Cio  ntabib:(Input) number of pairs associated with each site.
Cio        :This array may be altered by the program
Cio  iaxib :(Input) neighbor table containing pair information,
Cio        :in fixed-basis format (see iax2fd for conversion from
Cio        :standard form to fixed-basis format)
Cio        :This array may be altered by the program
Cl Local variables
Cl         :
Cr Remarks
Cr   To make table compatible with symmetry operations, extra
Cr   elements may be added to the existing table (mode 0), or elements
Cr   which don't have corresponding pairs in the table for each
Cr   rotation in the group are eliminated (mode 1)
Cu Updates
Cu   29 Aug 02  First created
C ---------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer niax,nbas,ntabib(nbas),ng,ldiaxb,mode,grow
      parameter (niax=10)
      integer iaxib(niax,nbas,ldiaxb)
      double precision pos(3,*),g(9,*),ag(3,*),plat(3,3)
C Local variables:
      integer ib,jb,ibas,isite,ig,ibp,jbp,ksite,ipr,ix,i,n,
     .  change,nttab,nttabn,isum,mode0,mode1,mode2
      double precision qlat(3,3),v(3),rv(3),vk(3),dx,fuzz
C ... ix-th component of connecting vector of pair i
      dx(ix,ibas,i) = pos(ix,iaxib(2,ibas,i)) - pos(ix,iaxib(1,ibas,i))
     .  + plat(ix,1)*iaxib(3,ibas,i)
     .  + plat(ix,2)*iaxib(4,ibas,i)
     .  + plat(ix,3)*iaxib(5,ibas,i)

C --- Setup ---
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (ng .le. 0) return
      fuzz = 1d-5
      call getpr(ipr)
      call dinv33(plat,1,qlat,v)
C     stdo = lgunit(1)
      nttab = isum(nbas,ntabib,1)
      grow = 0
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
C     call info(30,0,0,' symiax:  adjusting neighbor table for '//
C    .  'compatibility with symmmetry group',0,0)

C     call tcn('symiax')

C --- For each symmetry operation, do ---
      do  10  ig = 1, ng

        if (ig .eq. 1 .and. mode1 .eq. 0) goto 10

C       Re-entry point for restart in case some change made to table
   12   continue
        change = 0

C   --- For each pair in the neighbor list, do ---
        do  20  ib = 1, nbas
          do  22  isite = 1, ntabib(ib)

C         ib = iaxib(1,ib,isite)
            jb = iaxib(2,ib,isite)
            if (iaxib(1,ib,isite) .eq. 0 .or. jb .eq. 0) goto 22

            if (mode1 .ne. 1 .or. ig .ne. 1) then
C       ... Find ibp and jbp = sites that ib,jb are rotated into
C           call gpfndx(g(1,ig),ag(1,ig),ib,ibp,pos,nbas,plat,qlat)
C           call gpfndx(g(1,ig),ag(1,ig),jb,jbp,pos,nbas,plat,qlat)
              call grpfnd(fuzz,g,ag,ig,pos,nbas,qlat,ib,ibp)
              call grpfnd(fuzz,g,ag,ig,pos,nbas,qlat,jb,jbp)
              if (jbp .eq. 0) call rx('bug in symiax')
C       ... Original connecting vector v, and the rotated one rv
              v(1) = dx(1,ib,isite)
              v(2) = dx(2,ib,isite)
              v(3) = dx(3,ib,isite)
              call rotpnt(v,rv,g(1,ig))
C            if (ib .eq. 1 .and. isite .eq. 3) then
C              print *, ig, sngl(v)
C            endif
            else
C       ... Special case (ig=1) look for ri-rj inverse of rj-ri
              ibp = jb
              jbp = ib
              v(1) = dx(1,ib,isite)
              v(2) = dx(2,ib,isite)
              v(3) = dx(3,ib,isite)
              rv(1) = -v(1)
              rv(2) = -v(2)
              rv(3) = -v(3)
            endif

C     ... Find pair corresponding to ibp,jbp and rv
            do  40  ksite = 1, ntabib(ibp)

C       ... This is the corresponding pair only if jb of this pair = jbp
              if (iaxib(1,ibp,ksite) .ne. ibp) goto 40
              if (iaxib(2,ibp,ksite) .ne. jbp) goto 40
C       ... Also the connecting vector of ksite must be rv
              vk(1) = dx(1,ibp,ksite)
              vk(2) = dx(2,ibp,ksite)
              vk(3) = dx(3,ibp,ksite)
              if (abs(rv(1)-vk(1)) + abs(rv(2)-vk(2)) + abs(rv(3)-vk(3))
     .          .lt. fuzz) goto 42
   40       continue
C         No vector found.  Add a new element to the iax table
            if (mode0 .eq. 0) then
              do  ix = 1, 3
                rv(ix) = rv(ix) - (pos(ix,jbp) - pos(ix,ibp))
              enddo
              call dgemm('T','N',3,1,3,1d0,qlat,3,rv,3,0d0,vk,3)
              n = ntabib(ibp)+1
              grow = grow+1
              do  ix = 1, 3
                i = idnint(vk(ix))
C             if (abs(i-vk(ix)) .gt. fuzz) then
C               call rotpnt(v,rv,g(1,ig))
C               print *, 'oops'
C               rv(1) = rv(1) - (pos(1,jbp) - pos(1,ibp))
C               rv(2) = rv(2) - (pos(2,jbp) - pos(2,ibp))
C               rv(3) = rv(3) - (pos(3,jbp) - pos(3,ibp))
C             endif
                if (abs(i-vk(ix)) .gt. fuzz) call rx('bug in symiax')
                iaxib(2+ix,ibp,n) = i
C             call ivset(iaxib(1,n),6,10,0)
              enddo
              iaxib(1,ibp,n) = ibp
              iaxib(2,ibp,n) = jbp
              ntabib(ibp) = n
C           debugging check
              do  ix = 1, 3
                rv(ix) = rv(ix) + (pos(ix,jbp) - pos(ix,ibp))
                vk(ix) = dx(ix,ibp,n)
                if (abs(rv(ix)-vk(ix)) .gt. fuzz) call rx('bug in symiax')
              enddo
              change = 1
              call info8(70,0,0,' ig=%i isite=%i, ib,jb=%i %i add map to '
     .       //'site=%i, ibp,jbp=%i %i',ig,isite,ib,jb,ksite,ibp,jbp,0)
C         No vector found.  Remove element from the iax table
            else
              iaxib(1,ib,isite) = 0
              grow = grow-1
              change = -1
              call info8(70,0,0,' ig=%i isite=%i, ib,jb=%i %i remove: no '
     .        //'site=%i, ibp,jbp=%i %i',ig,isite,ib,jb,ksite,ibp,jbp,0)
            endif
            goto 22

C     ... Connecting vector found for rotated point
   42       continue
            if (ipr .ge. 110)
     .      call info8(90,0,0,' ig=%i isite=%i, ib,jb=%i %i  maps to '//
     .      'site=%i, ibp,jbp=%i %i',ig,isite,ib,jb,ksite,ibp,jbp,0)
   22     continue
   20   continue

C       Restart this loop if new pair added
        if (change .ne. 0) goto 12
   10 continue

      nttabn = isum(nbas,ntabib,1)
      if (grow .lt. 0) nttabn = nttabn + grow
      if (nttabn .gt. nttab) then
        call info5(30,0,0,' symiax: enlarged neighbor table from %i'//
     .    ' to %i pairs (%i symops)',nttab,nttabn,ng,0,0)
      elseif (nttabn .lt. nttab) then
        call info5(30,0,0,' symiax: reduced neighbor table from %i'//
     .    ' to %i pairs (%i symops)',nttab,nttabn,ng,0,0)
      else
        call info(30,0,0,' symiax: '
     .  //'neighbor table is compatible with symops',0,0)
      endif

C ... Sort enlarged table by increasing length
      if (grow .gt. 0 .and. mode2 .ne. 0)
     .  call symia1(plat,nbas,pos,ntabib,iaxib,ldiaxb)
C     call tcx('symiax')

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

      subroutine symia1(plat,nbas,pos,ntabib,iaxib,ldiaxb)
C- Kernel for symiax.  Sort iaxib table by increasing length
C ----------------------------------------------------------------------
Ci Inputs
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nbas  :size of basis
Ci   pos   :basis vectors
Ci   ntabib:number of pairs associated with each site.
Cio Inputs/Outputs
Cio  iaxib :(Input) neighbor table containing pair information,
Cio        :in fixed-basis format (see iax2fd for conversion from
Cio        :standard form to fixed-basis format)
Cio        :Array is reordered by increasing length
Co Outputs
Cr Remarks
Cr
Cu Updates
Cu   02 Sep 02 First created
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer niax,nbas,ntabib(nbas),ldiaxb
      parameter (niax=10)
      integer iaxib(niax,nbas,ldiaxb)
      double precision pos(3,*),plat(3,3)
C Local variables:
      integer ib,ibas,ix,i,is,iprm(ldiaxb),iwk(niax,ldiaxb)
      double precision dx,rtab(3,ldiaxb)

C ... ix-th component of connecting vector of pair i
      dx(ix,ibas,i) = pos(ix,iaxib(2,ibas,i)) - pos(ix,iaxib(1,ibas,i))
     .  + plat(ix,1)*iaxib(3,ibas,i)
     .  + plat(ix,2)*iaxib(4,ibas,i)
     .  + plat(ix,3)*iaxib(5,ibas,i)

      do  ib = 1, nbas
        do  is = 1, ntabib(ib)
          rtab(1,is) = dx(1,ib,is)
          rtab(2,is) = dx(2,ib,is)
          rtab(3,is) = dx(3,ib,is)
C         print 333, is,iaxib(1,ib,is),iaxib(2,ib,is),
C    .      iaxib(3,ib,is),iaxib(4,ib,is),iaxib(5,ib,is),
C    .      dx(1,ib,is),dx(2,ib,is),dx(3,ib,is),
C    .      sqrt(dx(1,ib,is)**2+dx(2,ib,is)**2+dx(3,ib,is)**2)
        enddo
C 333   format(i4,2x,5i4,2x,3f12.6,2x,f12.6)
        call dvheap(3,ntabib(ib),rtab,iprm,1d-4,111)
        do  is = 1, ntabib(ib)
          i = iprm(is)
          call icopy(niax,iaxib(1,ib,i),1,iwk(1,is),1)
        enddo
        do  is = 1, ntabib(ib)
          call icopy(niax,iwk(1,is),1,iaxib(1,ib,is),1)
C         print 333, is,iaxib(1,ib,is),iaxib(2,ib,is),
C    .      iaxib(3,ib,is),iaxib(4,ib,is),iaxib(5,ib,is),
C    .      dx(1,ib,is),dx(2,ib,is),dx(3,ib,is),
C    .      sqrt(dx(1,ib,is)**2+dx(2,ib,is)**2+dx(3,ib,is)**2)
        enddo
      enddo

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

