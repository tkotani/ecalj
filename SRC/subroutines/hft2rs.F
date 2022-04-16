      module m_hft2rs
      public:: hft2rs,hft2rsinit
      private
      contains
      subroutine hft2rsinit(n1,n2,n3,k1,k2,k3,qoff,isp,nsp,nbas,g,ag,
     .     nsgrp,range,offH,plat,pos,iv_a_ontab,iv_a_oiax,ndimh) !,ndhrs, hrs)
      use m_lmfinit,only: stdo
      use m_pairs,only: pairs,ppair1
      implicit none
      integer mode,idim,jdim,nbas,n1,n2,n3,k1,k2,k3,isp,nsp,nsgrp,i_data_size
      integer ia1,ia2,ndimh
      integer nkap0,n0H,niax
      parameter (nkap0=3,n0H=5,niax=10)
      integer offH(nbas+1)
      integer,allocatable :: iv_a_oiax(:)
      integer,allocatable :: iv_a_ontab(:)
      double precision plat(3,3),pos(3,*),range,g(9,*),ag(3,*),qoff(3)
      double complex hq(1)
      real(8):: hrs(1,1,1,1)    !ndhrs,ndhrs,nsp,1)
      integer:: ndhrs=1
C ... Local parameters
      integer i,ipr,iset,mode0,mode1,mode2,mode3,mxcsiz,nttab,nviax,nvtot
      integer ,allocatable :: iwk_iv(:),iv_a_tmp(:)
      double precision dmx(2)
      logical:: isanrg, l_dummy_isanrg
      integer ::iwdummy
      call tcn('hft2rsinit')
      call getpr(ipr)
C --- Mode 0 (setup) make and reduce iax table ---c     if (mode0 .eq. 0) then
      write(6,*) ' hft2rs: make neighbor table for r.s. hamiltonian using range = alat*',range
!  ... pair table for rs hamiltonian
      call pshpr(min(ipr-10,30))
      mxcsiz = 0
      call pairs ( nbas , nbas , 1d0 , plat , ( / range / 2 / ) , pos
     .     , ( / - 1 / ) , 3 , - 1 , (/iwdummy/) , nttab , iv_a_ontab , iv_a_oiax 
     .     , mxcsiz )
      call poppr
C     ... Standard Bloch transform doesn't require mapping to mesh
C     ... Mark all members in iax table with vector in sfz (for FFT)
      nviax = 0
      nvtot = 0
      mode = 0 !1000
      call hft2r1 ( mode , n1 , n2 , n3 , k1 , k2 , k3 , isp , nsp 
     .     , 1 , nbas , 1 , nbas , offh , ndimh , ndimh , hq , plat , pos
     .     , iv_a_ontab , iv_a_oiax , qoff , ndhrs , hrs , nviax , nvtot  , dmx )
      write(stdo,"(' hft2rs: found ',i8,' connecting vectors out of ',i8,' possible for FFT')")nviax,nvtot
!  ... Reduce iax table to elements on mesh in common with pairc
      call hft2r3 ( nbas  , iv_a_ontab , iv_a_oiax , nttab , mxcsiz )
!   ... Add pairs so that any ri-rj has corresponding rj-ri
      i = 110
      nttab = nttab*3
      i_data_size=size(iv_a_oiax)
      allocate(iv_a_tmp(i_data_size))
      iv_a_tmp=iv_a_oiax
      deallocate(iv_a_oiax)
      i_data_size=min(i_data_size,niax*nttab)
      allocate(iv_a_oiax(niax*nttab))
      iv_a_oiax(:i_data_size)=iv_a_tmp(:i_data_size)
      deallocate(iv_a_tmp)
      call pshpr(ipr)
      call symiax ( i , plat , nbas , pos , g , ag , nsgrp , iv_a_ontab, iv_a_oiax , nttab , mxcsiz )
      call poppr
C     ... Fill out iax(6,7,10)
      call pshpr(0)
      call ppair1 ( 0 , 1 , nbas , nbas , [- 1] , 1d0 , plat , pos ,
     .     [range / 2] , nttab , iv_a_ontab , iv_a_oiax , mxcsiz )
      call poppr
C     ... iax is now complete ... reallocate array
      i_data_size=size(iv_a_oiax); allocate(iv_a_tmp(i_data_size))
      iv_a_tmp=iv_a_oiax; deallocate(iv_a_oiax)
      i_data_size=min(i_data_size,niax*nttab); allocate(iv_a_oiax(niax*nttab))
      iv_a_oiax(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)
c$$$  if (ipr .gt. 40) then
c$$$  call pshpr(ipr-10)
c$$$  call ppair1 ( 20 , 1 , nbas , nbas , [- 1] , 1d0 , plat , pos ,
c$$$  .      [range / 2] , nttab , iv_a_ontab , iv_a_oiax , mxcsiz )
c$$$  call poppr
c$$$  endif
      call tcx('hft2rsinit')
      end subroutine hft2rsinit
      
      subroutine hft2rs(n1,n2,n3,k1,k2,k3,qoff,isp,nsp,nbas,g,ag,
     .nsgrp,range,offH,ia1,ia2,idim,jdim,hq,plat,pos,iv_a_ontab,iv_a_oiax,ndhrs, hrs)
C- Fourier transform of a matrix to real-space representation
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode    :See also Remarks
Ci           :1s digit specifies what is created.
Ci           :  0 (setup) create tables iax and ntab
Ci           :  1 (FT)    convert hq to h(T) by FFT
Ci           :  These modes are distinct so the caller can
Ci           :  allocate the appropriate memory for hrs.
Ci           :  Note: for converting hq to h(T) by standard FT,
Ci              see Remarks
Ci           :10s digit specifies whether hrs is real
Ci           :  0 No assumption made about hrs being real
Ci           :  1 Assume hrs is real (correct if hq is hermitian)
Ci           :100s digit specifies method for FT:
Ci           :  0 convert hq to h(T) by FFT.
Ci           :    hq is supplied on a uniform mesh of q-points
Ci           :    It is converted to R.S. by fft and copied to hrs.
Ci           :  1 convert hq to h(T) standard inverse Bloch transform
Ci           :    Note: hft2rs only does the setup for the
Ci           :    standard inverse Bloch transform; see Remarks
Ci           :1000s digit concerns generation of neighbor table:
Ci           :      (applicable only when 1s digit mode is zero)
Ci           :      See Note 6 in Remarks.
Ci           :  0 table is not altered by symmetry considerations
Ci           :  1 table is enlarged so that any pair connecting ri-rj
Ci           :    has a corresponding entry rj-ri.
Ci           :  2 table is reduced so that any pair connecting ri-rj
Ci           :    has a corresponding entry rj-ri.
Ci           :10000s digit specifies indexing of hrs_RL,R'L'
Ci           :  0 poke h_RL,R'L' (T) into hrs_RL,R'L'
Ci           :  1 poke h_RL,R'L' (T) into hrs_R'L',RL
Ci           :    (used when h(T) constructed for rj-ri rather than ri-rj)
Ci   n1,..n3 :number divisions in QP mesh
Ci   k1,..k3 :dimensions hq
Ci   qoff    :offset of gamma-point = q(1,1,1)
Ci   isp     :current spin channel (1 or 2): poke hq into hrs(isp)
Ci           :NB: hq has no spin index
Ci   nsp     :2 for spin-polarized case, otherwise 1
Ci   nbas    :size of basis
Ci   g       :point group operations
Ci   ag      :translation part of space group
Ci   nsgrp   :number of group operations
Ci   range   :
Ci   offH :Offsets to hamiltonian matrix (makidx.f)
Ci   idim    :hamiltonian dimension of ia1..ia2
Ci   jdim    :hamiltonian dimension of ib1..ib2
Ci   hq      :hamiltonian on uniform mesh of q-points
Ci           :or translation vectors (related to q-points by FT)
Ci   plat    :primitive lattice vectors, in units of alat
Ci   pos     :basis vectors
Ci   ndhrs   :leading dimensions of hrs: must be at least as large as
Ci           :the total number of orbitals of any atom.
Cio Inputs/Outputs
Cio  w(ontab):ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Cio           If 1s digit mode=0, memory is allocated for ntab and array
Cio           is created.  If 1s digit mode>1, ntab is an input.
Cio  w(oiax) :neighbor table containing pair information (pairc.f)
Cio           If 1s digit mode=0, memory is allocated for iax and array
Cio           is created.  If 1s digit mode>1, iax is an input.
Co   hrs     :real-space h in neighbor-list format (see bloch.f)
Co           :hrs is created when 1s digit mode>0 only.
Cr Remarks
Cr   hft2rs is designed to generate a real-space h_ij(T) from h_ij(q).
Cr   h_ij(T) is kept in standard tight-binding format, which entails
Cr   the following:
Cr     1. a neighbor table iax that contains a list of all connecting
Cr        vectors connecting sites i and j for the r.s. hq.
Cr     2. h_ij(T) is stored as a list of matrices; one matrix
Cr        connects all orbitals at site i to all orbitals at site j+T.
Cr     3. the forward Bloch transform of h_ij(T) should generate
Cr        h_ij (q); this routine generates h_ij(T) in the form that
Cr        routine Bloch will make h_ij(q) from h_ij(T).
Cr     4. The inverse Bloch transform may be done by FFT or by
Cr        a standard transform.  In the former case, the forward
Cr        transform from the hq(T) computed here should be exact on
Cr        the mesh of q-points from which hq(T) was made.
Cr        For the standard inverse transform it need not.
Cr     5. This routine is called once in a setup mode to create the
Cr        pair table (1s digit mode=0).  For the inverse transform by
Cr        FFT it is called a second time to do the mapping.  hft2rs will
Cr        not do the inverse transform using the standard approach,
Cr        because of complications connected with symmetrization.
Cr
Cr        For FFT, you will want to use the following procedure:
Cr          a) call hft2rs in setup mode (1s digit mode = 0)
Cr             Note: after this call, you may want to adjust the
Cr             neigbor table to accomodate symmetry operations;
Cr             see Note 7 below.
Cr          b) Create hq(q) in the full BZ.  If you have only hq(q)
Cr             for irreducible points, see hamfb3.f
Cr          c) call hft2rs again, with 1s digit mode = 1
Cr             At this point the FT is complete.  However it may
Cr             have someone undesirable properties as regards
Cr             symmetry. You may want to symmetrize the result;
Cr             see Note 6 below.
Cr        The standard inverse transform differs in the following ways:
Cr          * the inverse transform need to be exact
Cr          * the inverse transform uses hq(q) at the irreducible
Cr            qpoints
Cr        For this method use the following procedure:
Cr          a) call hft2rs in setup mode (1s digit mode = 0)
Cr             Note: after this call, you will want to adjust the
Cr             neigbor table to accomodate symmetry operations;
Cr             see Note 6 below.
Cr          b) make an unsymmetrized inverse Bloch transforming
Cr             by calling ibloch for each irreducible qp.
Cr          c) Symmetrize the inverse transform by calling rsmsym.
Cr
Cr   6.  Symmetry.  The pair table may be at variance with the symmetry
Cr       of the system.  This can happen for two reasons: points that
Cr       are symmetry-related may not all find their way into the pair
Cr       table because their connecting vector is within a tight
Cr       tolerance of the r.s. cutoff radius, or in the case of the FFT,
Cr       the list of points is reduce to make a list whose size
Cr       corresponds to the mesh of q-points.  You can cause hft2rs to
Cr       make the table compatible with symmetry operations by setting
Cr       the 1000s digit of mode in the setup call.  The symmetry
Cr       operations are passed as the triplet (g,ag,nsgrp)
Cr       hft2rs actually changes the table in two ways:
Cr
Cr       * by adding to the table all missing vectors that are
Cr         symmetry-related to vectors already in the table
Cr
Cr       * by adding to the table all vectors not in the table that
Cr         correspond to the inverse of vectors in the table.
Cr
Cr       Thus, even if you tell hft2rs to change the table with only
Cr       on group operation, it may affect the table.
Cr
Cr       Note, however, the FFT inverse will only poke hq(q) into
Cr       hq(T) for mesh of T it chooses.  If you want the transform
Cr       to be symmetrize, call rsmsym to do it.
Cr
Cr   Notes on the meaning of row and column dimensions:
Cr     ib=field(hq=gf)   or augmentation(hq=ham)
Cr     jb=source(hif=gf) or basis(hq=ham)
Cu Updates
Cu   25 Sep 04 Handles offset BZ mesh
Cu   30 Mar 03 New switch for copying h(T) to transpose hrs
Cu    3 Oct 02 redesigned
Cu   23 Jun 02 First written
C ----------------------------------------------------------------------
      implicit none
      integer mode,idim,jdim,nbas,n1,n2,n3,k1,k2,k3,isp,nsp,ndhrs,nsgrp,i_data_size
      integer ia1,ia2
      integer nkap0,n0H,niax
      parameter (nkap0=3,n0H=5,niax=10)
      integer offH(nbas+1)
      integer,allocatable :: iv_a_oiax(:)
      integer,allocatable :: iv_a_ontab(:)
      double precision plat(3,3),pos(3,*),range,g(9,*),ag(3,*),qoff(3)
      double complex hq(k1,k2,k3,idim,jdim,nsp)
      real(8):: hrs(ndhrs,ndhrs,nsp,1)
      integer i,ipr,iset,mode0,mode1,mode2,mode3,mxcsiz,nttab,nviax,nvtot
      integer ,allocatable :: iwk_iv(:),iv_a_tmp(:)
      double precision dmx(2)
      logical:: isanrg, l_dummy_isanrg
      integer ::iwdummy
      call tcn('hft2rs')
      mode= 1 
      call getpr(ipr)
      iset = 0
      call fftz3(hq,n1,n2,n3,k1,k2,k3,idim*jdim,iset,-1)
      call hft2r1 ( mode , n1 , n2 , n3 , k1 , k2 , k3 , isp , nsp 
     .    , ia1 , ia2 , 1 , nbas , offh , idim , jdim , hq , plat , pos
     .    , iv_a_ontab , iv_a_oiax , qoff , ndhrs , hrs , nviax , nvtot 
     .    , dmx )
      call info5(30,0,0,' hft2rs created hrs:  ndhrs=%i  '//
     .  'max Re(hrs) = %;3g  max Im(hrs) = %;3g',ndhrs,dmx,dmx(2),0,0)
      call tcx('hft2rs')
      end subroutine hft2rs
      
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
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
      integer ,allocatable :: iaxib_iv(:)
      integer ,allocatable :: ntabb_iv(:)
C ... Copy iax to fixed-basis format
C     Assume no cluster in new iax table will not exceed 8x old one
      mxcsz = 0
      do  ib = 1, nbas
        mxcsz = max(mxcsz,ntab(ib+1)-ntab(ib))
      enddo
      ldiaxb = 8*mxcsz
      allocate(iaxib_iv(ldiaxb*nbas*niax))
      iaxib_iv=0
      allocate(ntabb_iv(nbas))
      call iax2fd ( 0 , nbas , ntab , iax , ntabb_iv , iaxib_iv , ldiaxb )
C ... Update the table
      if (ng .gt. 0) then
        call symia0 ( mode , plat , nbas , pos , g , ag , ng , ntabb_iv 
     .  , iaxib_iv , ldiaxb , grow )
        ib = ntab(nbas+1)+grow
        if (ib .gt. nttab) call rxi('symiax: not enough memory for new iax table: need dim',ib)
      endif
C ... Copy iax from fixed-basis format
      call iax2fd ( 1 , nbas , ntab , iax , ntabb_iv , iaxib_iv , nttab )
      mxcsz = 0
      do  ib = 1, nbas
        mxcsz = max(mxcsz,ntab(ib+1)-ntab(ib))
      enddo
      if (allocated(ntabb_iv)) deallocate(ntabb_iv)
      if (allocated(iaxib_iv)) deallocate(iaxib_iv)
      end

      subroutine symia0(mode,plat,nbas,pos,g,ag,ng,ntabib,iaxib,ldiaxb,
     .grow)
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
      double precision pos(3,*),g(3,3,*),ag(3,*),plat(3,3)
C Local variables:
      integer ib,jb,ibas,isite,ig,ibp,jbp,ksite,ipr,ix,i,n,
     .change,nttab,nttabn,isum,mode0,mode1,mode2
      double precision qlat(3,3),v(3),rv(3),vk(3),dx,fuzz
C ... ix-th component of connecting vector of pair i
      dx(ix,ibas,i) = pos(ix,iaxib(2,ibas,i)) - pos(ix,iaxib(1,ibas,i))
     &     + sum(plat(ix,1:3)*iaxib(3:5,ibas,i))
c     .+ plat(ix,1)*iaxib(3,ibas,i)+ plat(ix,2)*iaxib(4,ibas,i)+ plat(ix,3)*iaxib(5,ibas,i)
C --- Setup ---
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
              call grpfnd(fuzz,g,ag,ig,pos,nbas,qlat,ib,ibp)
              call grpfnd(fuzz,g,ag,ig,pos,nbas,qlat,jb,jbp)
              if (jbp .eq. 0) call rx('bug in symiax')
C       ... Original connecting vector v, and the rotated one rv
              v = [(dx(i,ib,isite),i=1,3)]
              rv= matmul(g(:,:,ig),v) !call rotpnt(v,rv,g(1,ig))
            else
C       ... Special case (ig=1) look for ri-rj inverse of rj-ri
              ibp = jb
              jbp = ib
              v = [(dx(i,ib,isite),i=1,3)]
              rv(:) = -v(:)
            endif
C     ... Find pair corresponding to ibp,jbp and rv
            do  40  ksite = 1, ntabib(ibp)
C       ... This is the corresponding pair only if jb of this pair = jbp
              if (iaxib(1,ibp,ksite) .ne. ibp) goto 40
              if (iaxib(2,ibp,ksite) .ne. jbp) goto 40
C       ... Also the connecting vector of ksite must be rv
              vk(:) = [(dx(i,ibp,ksite),i=1,3)]
              if (sum(abs(rv(:)-vk(:)))<fuzz) goto 42 ! + abs(rv(2)-vk(2)) + abs(rv(3)-vk(3))
   40       continue
C         No vector found.  Add a new element to the iax table
            if (mode0 .eq. 0) then
              rv(:) = rv(:) - (pos(:,jbp) - pos(:,ibp))
              call dgemm('T','N',3,1,3,1d0,qlat,3,rv,3,0d0,vk,3)
              n = ntabib(ibp)+1
              grow = grow+1
              do  ix = 1, 3
                i = idnint(vk(ix))
                if (abs(i-vk(ix)) .gt. fuzz) call rx('bug in symiax')
                iaxib(2+ix,ibp,n) = i
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
     .        //'site=%i, ibp,jbp=%i %i',ig,isite,ib,jb,ksite,ibp,jbp,0)
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
     .  ' to %i pairs (%i symops)',nttab,nttabn,ng,0,0)
      elseif (nttabn .lt. nttab) then
        call info5(30,0,0,' symiax: reduced neighbor table from %i'//
     .  ' to %i pairs (%i symops)',nttab,nttabn,ng,0,0)
      else
        call info(30,0,0,' symiax: '
     .  //'neighbor table is compatible with symops',0,0)
      endif
C ... Sort enlarged table by increasing length
      if (grow .gt. 0 .and. mode2 .ne. 0)
     .call symia1(plat,nbas,pos,ntabib,iaxib,ldiaxb)
C     call tcx('symiax')
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
     .+ plat(ix,1)*iaxib(3,ibas,i)
     .+ plat(ix,2)*iaxib(4,ibas,i)
     .+ plat(ix,3)*iaxib(5,ibas,i)

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
      end subroutine
      end module m_hft2rs

      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine hft2r1(mode,n1,n2,n3,k1,k2,k3,isp,nsp,ia1,ia2,ib1,ib2,
     .offH,idim,jdim,hij,plat,pos,ntab,iax,qoff,ndhrs,hrs,nviax,nvtot,
     .dmx)
C- Copy a R.S. hamiltonian subblock on a uniform mesh of translation
C- vectors to a representation with a neighbor list,
C  or just mark entries in the table that correspond to mesh points
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode    :1s digit specifies what is created.
Ci           :  0 marks entries in row 8 of iax table; see iax below
Ci           :    In this mode, hij and hrs are not used
Ci           : >0 copy hij to hrs
Cixxx           :10s digit specifies whether hrs is real
Ci           :  0 No assumption made about hrs being real
Cixxx           :  1 Assume hrs is real (correct if hq is hermitian)
Ci           :10000s digit specifies indexing of hrs_RL,R'L'
Ci           :  0 poke h(T) into hrs_RL,R'L'
Ci           :  1 poke h(T) into hrs_R'L',RL
Ci           :    (used when h(T) constructed for rj-ri rather than ri-rj)
Ci   n1,..n3 :number divisions in QP mesh
Ci   isp     :current spin channel (1 or 2): poke hij into hrs(isp)
Ci           :NB: hij has no spin index
Ci   nsp     :2 for spin-polarized case, otherwise 1
Ci   ia1,ia2 :range of field point sites for which hij is tabulated
Ci   ib1,ib2 :range of source point sites for which hij is tabulated
Ci   offH    :Offsets to hamiltonian matrix (makidx.f)
Ci   idim    :hamiltonian dimension of ia1..ia2
Ci   jdim    :hamiltonian dimension of ib1..ib2
Ci   hij     :hamiltonian on uniform mesh of q-points
Ci   plat    :primitive lattice vectors, in units of alat
Ci   pos     :basis vectors
Ci   ntab    :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax     :neighbor table containing pair information (pairc.f)
Ci   ndhrs   :leading dimensions of hrs: must be at least as large as
Ci           :the total number of orbitals of any atom.
Co Outputs
Co   hrs     :diagonal GF for this principal layer
Co   nviax   :number of connecting vectors found in iax table added to nviax
Co   nvtot   :total number of connecting vectors added to nvtot
Co   dmx     :dmx(1) = largest real element found in h
Co           :dmx(2) = largest imaginary element found in h
Cr Remarks
Cr   ib=field(hij=gf)  or augmentation(hij=ham)
Cr   jb=source(hif=gf) or basis(hij=ham)
Cu Updates
Cu   01 May 06 dmx is properly global maximum
Cu   30 Mar 03 New switch for copying h(T) to transpose hrs
Cu   23 Jun 02 Adapted from gfq2rs.f
C ----------------------------------------------------------------------
      implicit none
      integer,parameter :: nkap0=3,n0H=5,niax=10
      integer:: mode,idim,jdim,ia1,ia2,ib1,ib2,n1,n2,n3,k1,k2,k3,isp,nsp,ndhrs,nviax,nvtot,
     & id,jd,i,ia,ib,ngmx,ng,is1,is2,nlma,nlmb,offi,offj,nhit,nmiss,
     & offH(*),iax(niax,1),ntab(ia2),wdummy
      real(8):: plat(3,3),pos(3,*),dmx(2),qoff(3),gmax,qlat(3,3),tau(3),dmx1(2)
      complex(8):: hij(k1,k2,k3,idim,jdim),hrs(ndhrs,ndhrs,nsp,1)
      real(8) ,allocatable :: rtab_rv(:), gv_rv(:)
      integer ,allocatable :: kv_iv(:)
      logical:: ltrans, isanrg, l_dummy_isanrg,debug=.false.
      id = offH(ia2+1) - offH(ia1)
      jd = offH(ib2+1) - offH(ib1)
      l_dummy_isanrg=isanrg(idim,id,id,'hft2rs:1 ','idim',.true.)
      l_dummy_isanrg=isanrg(jdim,jd,jd,'hft2rs:2 ','jdim',.true.)
      ltrans = .true. !mod(mode/10000,10) .ne. 0
      dmx(1) = 0
      dmx(2) = 0
      offi = 0
      do  10  ia = ia1, ia2
        offj = 0
        do  12  ib = ib1, ib2
          is1 = ntab(ib)+1
          is2 = ntab(ib+1)
          nlma = offH(ia+1) - offH(ia)
          nlmb = offH(ib+1) - offH(ib)
C --- Make list for a dimensionless lattice (alat=1) ---
          call dinv33(plat,1,qlat,gmax)
          call dpscop(pos,tau,3,3*ib-2,1,1d0)
          call dpsadd(tau,pos,3,1,3*ia-2,-1d0)
          ngmx = n1*n2*n3
          allocate(gv_rv(ngmx*3),kv_iv(ngmx*3))
          call pshpr(80*0)
          i = 208 + 400*0
          call gvlst2(1d0, qlat, -tau , n1,n2,n3, 0d0,1d5,0,i, ngmx,ng, kv_iv,gv_rv,wdummy,wdummy)
          call poppr
C --- For each connecting vector, mark iax table or poke into hrs ---
            call hft2r2 ( mode , k1 , k2 , k3 , ia , ib , jdim , idim , iax 
     .      , qoff , is1 , is2 , ndhrs , isp , nsp , nlmb , nlma , hij ( 
     .      1 , 1 , 1 , 1 + offj , 1 + offi ) , gv_rv , kv_iv , ng 
     .      , tau , plat , pos ,  hrs , hrs , nhit , nmiss ,   dmx1 )
          dmx(1) = max(dmx(1),dmx1(1))
          dmx(2) = max(dmx(2),dmx1(2))
          nviax = nviax + nhit
          nvtot = nvtot + nhit + nmiss
          if (allocated(kv_iv)) deallocate(kv_iv)
          if (allocated(gv_rv)) deallocate(gv_rv)
          offj = offj + offH(ib+1) - offH(ib)
   12   continue
        offi = offi + offH(ia+1) - offH(ia)
   10 continue
      end
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine hft2r2(mode,n1,n2,n3,ia,ib,idim,jdim,iax,qoff,is1,is2,
     .ndhrs,isp,nsp,nlma,nlmb,hij,gv,kv,ng,tau,plat,pos,hrs,hrsr,
     .nhit,nmiss,dmx)
C- For each connecting vector in a list taken from a uniform mesh of
C- translation vectors find corresponding place r.s neighbor table, and:
C  either mark iax table, or copy into hrs h generated from uniform mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode    :1s digit specifies what is created.
Ci           :  0 marks entries in row 8 of iax table; see iax below
Ci           :    In this mode, hij and hrs are not used
Ci           : >0 copy hij to hrs
Cixxx           :10s digit specifies whether hrs is real
Ci           :  0 No assumption made about hrs being real
Cixxx           :  1 Assume hrs is real
Ci           :Other digits are not used.
Ci   n1,..n3:number divisions in QP mesh
Cr   ia    :site index for field point (or row or augmentation dim)
Ci   ia    :site index for source point (or col or basis dimension)
Ci   idim  :augmentation dimension of hij
Ci   jdim  :basis dimension of hij
Ci   is1,is2:Look in iax table for pairs is1..is2
Ci   ndhrs :leading dimensions of hrs: must be at least as large as
Ci         :the total number of orbitals of any atom.
Ci   isp   :current spin channel (1 or 2): poke hij into hrs(isp)
Ci         :NB: hij has no spin index
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlma  :number of orbitals in field  (row) dimension
Ci   nlmb  :number of orbitals in source (col) dimension
Ci   hij   :R.S hamiltonian on uniform mesh of translation vectors
Ci   gv    :list of r.s. lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f) mapping gv
Ci         :to the uniform mesh.  Not used if 1s digit mode is zero
Ci   ng    :number of G-vectors
Ci   tau   :offset pos(ib)-pos(ia) added to lattice vectors gv to make
Ci         :true connecting vectors.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Cio Inputs/Outputs
Cio   iax  :neighbor table containing pair information (pairc.f)
Cio        :Output if 1s digit mode is zero:
Cio        :iax(8,k) set to -1 for pairs k that correspond
Cio        :         to a connecting vecor
Cio        :if 1s digit mode is >0, iax is not changed.
Co   hrs   :see hrsr
Co   hrsr  :copy hij to hrs or hrsr for pairs found in aix table.
Co         :hij copied EITHER into hrs (complex) OR hrsr (real)
Co         :depending on 10s digit mode
Co         :Also copy only takes place if 1s digit mode>0
Co   nhit  :number of connecting vectors found in aix table
Co   nmiss :number of connecting vectors not found in aix table
Co   dmx   :dmx(1) = largest real element found in h
Co         :dmx(2) = largest imaginary element found in h
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Mar 03 New switch for copying h(T) to transpose hrs
Cu   23 Jun 02  Adapted from pgfq2r
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax,idim,jdim,ia,ib,n1,n2,n3,ndhrs,isp,nsp,nlma,nlmb,is1,is2,mode,nhit,nmiss
      parameter (niax=10)
      integer ng,kv(ng,3),iax(niax,1)
      double precision tau(3),gv(ng,3),pos(3,*),plat(3,3)
      double precision hrsr(ndhrs,ndhrs,nsp,is2),dmx(2),qoff(3)
      double complex hrs(ndhrs,ndhrs,nsp,is2),hij(n1,n2,n3,idim,jdim)
      logical ltrans
      integer ig,i1,i2,i3,is,ix,offx,irow,icol,mode0,mode1,j,k
      double precision v(3),tol,dimx,dimy,TdotK,twopi
      double complex phase
      parameter (tol=1d-5)
      logical:: debug=.false.
      real(8):: rtab(3,is2-is1+1)
      mode0 = mod(mode,10)
c      mode1 = 0 !mod(mode/10,10)
      ltrans = .true. !mod(mode/10000,10) .ne. 0
c      if(debug) print *,'hft2r2: mode0 mode1 ltrans=',mode0,mode1,ltrans
      dimx = 0
      dimy = 0
      offx = is1-1
      twopi = 8*datan(1d0)
C --- Make -(connecting vector) = rtab = (pos(ia=field)-pos(ib=src)) ---
C     NB : Connecting vector is pos(src)-pos(field) = -rtab
      do   is = is1, is2
      do   ix = 1, 3
        rtab(ix,is-offx) = pos(ix,ia) - pos(ix,ib)
     .    + plat(ix,1)*iax(3,is) + plat(ix,2)*iax(4,is) + plat(ix,3)*iax(5,is)
      enddo
      enddo
C --- Loop through list of R vectors ---
      nhit = 0
      nmiss = 0
      do  10  ig = 1, ng, 1
        if (mode0 .ne. 0) then
          i1 = kv(ig,1)
          i2 = kv(ig,2)
          i3 = kv(ig,3)
        endif
        v(:) = -(gv(ig,:) - tau(:))
C   ... If a match in the iax table is found, copy hij there
        do  20  is = is1, is2
          if (ia .ne. iax(2,is)) cycle 
          if (sum(abs(v(:)+rtab(:,is-offx))) .lt. tol) then
c            if (abs(v(2)+rtab(2,is-offx)) .lt. tol) then
c              if (abs(v(3)+rtab(3,is-offx)) .lt. tol) then
                nhit = nhit+1
                if (mode0 .eq. 0) then
                  iax(8,is) = -1
                else
                      do  icol = 1, nlmb
                        do  irow = 1, nlma
                          hrs(icol,irow,isp,is) = hij(i1,i2,i3,irow,icol) 
                          dimx = max(dimx,abs(dble(hij(i1,i2,i3,irow,icol))))
                          dimy = max(dimy,abs(dimag(hij(i1,i2,i3,irow,icol))))
                        enddo
                      enddo
                endif
                goto 30
c              endif
c            endif
          endif
   20   continue
        nmiss = nmiss+1
   30   continue
   10 continue
      dmx(1) = dimx
      dmx(2) = dimy
      end

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine hft2r3(ib2,ntab,iax,nttab,mxcsiz)
C- Reduce iax table to those marked
C ----------------------------------------------------------------------
Ci Inputs
Cixx   ib1   :compact iax starting at site ib1
Ci   ib2   :compact iax ending at site ib2
Cio Inputs/Outputs
Cio  ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Cio  iax   :neighbor table containing pair information (pairc.f)
Cio        :Pairs for which iax(8,is) is set to -1 are retained;
Cio        :the remaining pairs are purged from the list.
Cio        :On output iax(8,*) is reset to zero.
Co   nttab :number of pairs in reduced table
Co   mxcsiz:maximum cluster size (number of pairs connected to any site)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   23 Jun 02  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer,parameter:: niax=10
      integer::ib1=1,ib2,ntab(*),iwk((niax+1)*mxcsiz),nttab,mxcsiz
      integer iax(niax,1)
      integer is1,is2,ib,nstart,i,nib
C --- Reorder iax table site-by-site ---
      do  ib = ib1, ib2
        is1 = ntab(ib)+1
        is2 = ntab(ib+1)
        nstart = is2-is1+1
C       Copy iax(8) to iwk(1+nstart)
c        call icopy(nstart,iax(8,is1),niax,iwk(1+nstart),1)
        iwk(1+nstart:2*nstart)=iax(8,is1:is1+nstart-1)
C       Sort, using iwk(1) as permutation array
        call ivheap(1,nstart,iwk(1+nstart),iwk,101)
C       Rearrange iax, putting marked ones at top
        call ivprm(niax,nstart,iax(1,is1),iwk(1+nstart),iwk,.true.)
      enddo
C --- Condense iax table and ntab ---
      mxcsiz = 0
      nttab = ntab(ib1)
      do  ib = ib1, ib2
C       Count nib = number of marked pairs for this site
        is1 = ntab(ib)+1
        is2 = ntab(ib+1)
        nib = 0
        do  i = is1, is2
          nib = i-is1
          if (iax(8,i) .eq. 0) goto 20
        enddo
        nib = is2-is1+1
   20   continue
        mxcsiz = max(mxcsiz,nib)
C       Append iax for this ib to the end of now shortened table
        if (is1 .ne. nttab+1) then
          call icopy(niax*nib,iax(1,is1),1,iax(1,nttab+1),1)
        endif
        ntab(ib) = nttab
        nttab = nttab+nib
      enddo
      ntab(ib2+1) = nttab
C --- Zero out iax(6..10) ---
      do  i = ntab(ib1)+1, nttab
        iax(6,i) = 0
        iax(7,i) = 0
        iax(8,i) = 0
        iax(9,i) = 0
        iax(10,i) = 0
      enddo
      end subroutine hft2r3
c$$$cccccccccccccccccccccccccccccccccccccc
c$$$      subroutine rotpnt(p,q,g)
c$$$C     implicit none
c$$$      double precision p(3),q(3),g(3,3),h(3)
c$$$      integer i,j
c$$$      do  i=1,3
c$$$        h(i)=0.d0
c$$$      do  j=1,3
c$$$         h(i)=h(i)+g(i,j)*p(j)
c$$$      enddo
c$$$      enddo
c$$$      do 2 i=1,3
c$$$        q(i)=h(i)
c$$$    2 continue
c$$$      end  subroutine
c$$$ 
