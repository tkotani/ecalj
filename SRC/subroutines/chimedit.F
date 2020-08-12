c$$$ commented out on 22mar2015

c$$$      subroutine chimedit(sopts,mode,sctrl,ssite,sspec,slat,spot,sbz,
c$$$     .nbas,nat,nspec)
c$$$
c$$$      use m_struc_def  !Cgetarg
c$$$      use m_globalvariables
c$$$C- Magnetic Linear response editor
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs/Outputs
c$$$Ci   sopts :command options performed automatically, before reading
c$$$Ci         :from standard input
c$$$Ci   mode  :0 longitudinal linear response, lmf
c$$$Ci         :1 transverse linear response, lmf
c$$$Ci   sctrl :struct for program flow parameters; see routine uctrl
c$$$Ci     Elts read: zbak lrs
c$$$Ci     Stored:    lrs
c$$$Ci     Passed to:
c$$$Ci   ssite :struct for site-specific information; see routine usite
c$$$Ci     Elts read:
c$$$Ci     Stored:
c$$$Ci     Passed to: dcopy
c$$$Ci   sspec :struct for species-specific information; see routine uspec
c$$$Ci     Elts read:
c$$$Ci     Stored:
c$$$Ci     Passed to: dcopy
c$$$Ci   slat  :struct for lattice information; see routine ulat
c$$$Ci     Elts read: nabc
c$$$Ci     Stored:
c$$$Ci     Passed to: dcopy rdovfa
c$$$Ci   spot  :struct for information about the potential; see routine upot
c$$$Ci     Elts read: nlml nlma
c$$$Ci     Stored:
c$$$Ci     Passed to: dcopy rdovfa
c$$$Ci   sbz   :struct for the Brillouin Zone; see routine ubz
c$$$Ci     Elts read:
c$$$Ci     Stored:
c$$$Ci     Passed to:
c$$$Ci   nbas  :size of basis
c$$$Ci   nat   :number atoms in basis with augmentation sites
c$$$Ci         :Note: if nat<nbas, there is a requirement that
c$$$Ci         :lmxa>-1 for nat sites, and
c$$$Ci         :and lmxa=-1 for nbas-nat sites
c$$$Co Outputs
c$$$Co   chiedit never returns.
c$$$Co   rst file can be written.
c$$$Cr Remarks
c$$$Cr   The density consists of a smooth part (smrho) plus
c$$$Cr   nbas atom-centered densities inside the MT spheres.
c$$$Cr   Their sum is the full charge density.
c$$$Cr   The local density is represented as the difference of the
c$$$Cr   two valence components in orhoat, plus the core density.
c$$$Cr   Density in the MT spheres:
c$$$Cr      mesh parameters rmt,nr,a;
c$$$Cr      total density rho (times r**2) to lmxl;
c$$$Cr      a spherical potential v0 defining the wave functions within rmt
c$$$Cr      pnu and idmod to lmxa
c$$$Cr   Smooth density
c$$$Cr      real part of complex*16 array smrho contains the density
c$$$Cr      k1,k2,k3 are the physical dimensions of the array
c$$$Cr      n1,n2,n3 are the dimensions of the mesh.
c$$$Cl Local variables
c$$$Cu Updates
c$$$Cu   06 Dec 08 First created
c$$$C  ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   nx       : cutoff number of w points in calculation, such  X0(q,w)
c$$$Ci   cutw     : cutoff value(in meV)  of w in calculation of X(q,w)
c$$$Ci   nbloch   : dimension of product basis
c$$$Ci   ngbin    : dimension of mix basis
c$$$Ci   nq       : total number of q points in calculation of X(q,w)
c$$$Ci   qpi      : vector of current q,   X(q,w)
c$$$Ci   qgbin    : vector of all q points,   qgbin(3,nq0i)
c$$$Ci   imbas    : magnetic atom index
c$$$Ci   momsite  : magnetic moment m
c$$$Ci   svec     : <B_I|m_a>  svec(nbloch,nmbas)
c$$$Ci   mmnorm   : |m|=sqrt(<m|B> <B|m> )=norm(svec) ; this is not m.
c$$$Ci            : ~e_a(r)=m_a(r) / sqrt(\int m_a(r)**2 dr)  then <~e_a|~e_a>=1
c$$$Ci            : <B_i|~e_a>=svec/mmnorm
c$$$Ci   biubj    : <B|U|B>  stoner matrix U(nbloch,nbloch)
c$$$Ci   zxq      : <B|X0|B>   chi0 matrix zxq(nbloch,nbloch,nx) denpend on (q,w)
c$$$Ci   chi0     : <~e|X0|~e> = <m|B><B|X0|B><B|m> / norm(<m_i|m_i>) /norm(<m_j|m_j>)
c$$$Ci            :            = <m|x0|m'>/ norm(<m_i|m_i>) /norm(<B|m_j>)
c$$$Ci            :            = <m_a|x0|m_a'>/(sqrt(<m_i|m_i> sqrt(<m_j|m_j>)
c$$$Ci            :     |~e_a> = |m_a>/<m_a|m_a>
c$$$Ci   gbvec    : <e^iqr|B+>  projection of |e^iqr> on mix basis set,  gbvec(ngbin,nq0i)
c$$$Ci   gbvec0   : <e^iqr|B>   subset(prod. bas. part) of gbvec,  gbvec0(nbloch,nq0i)
c$$$Ci   eiqrm    : eiqrm_i = <~e_i|e^{iqr}> =  <M_i|eiqr>/sqrt(<M_i|M_i>)
c$$$Ci            : sum( dconjg(gbvec(1:nbloch))*svec(1:nbloch,imb) )/mmnorm(imb)
c$$$Ci   freq     : frequency w in unit of ryberg
c$$$Ci   freq_meV : frequency w in unit of meV
c$$$Ci   mx2m     : <m|x|m>  from rigid spin approximation
c$$$Ci   mx0m     : <m|x0|m> from rigid spin approximation
c$$$Co Outputs:
c$$$Co   mim      : <m|i|m>  stoner parameter from rigid spin approximation
c$$$Co   x0meanx  : <m|B><B|X0|B><B|m>
c$$$Co   mUm      : <m|B><B|U|B><B|m>
c$$$Co   chike    : <m|B><B|X|B><B|m>
c$$$Co   iqrXiqr  : <e^{iqr}|B><B|X|B><B|e^{iqr}>
c$$$
c$$$
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character sopts*(*)
c$$$      integer mode,nbas,nat,nspec,ifi,n0
c$$$      parameter (n0=10)
c$$$      type(s_ctrl)::sctrl
c$$$      type(s_site)::ssite(*)
c$$$      type(s_spec)::sspec(*)
c$$$      type(s_lat)::slat
c$$$      type(s_pot)::spot
c$$$      type(s_bz)::sbz
c$$$
c$$$C ... Local parameters
c$$$      integer nglob,fopna,fopng,a2vec,havechi0,iq,nq,nsgrp
c$$$      integer i,j,l,im,m,n,j1,j2,js1,js2,nw,ne,iinear,lqp
c$$$
c$$$      logical lnsave,lsopts,isanrg,llshft(3)
c$$$      integer nkxyz(3),mxkp,nk1,nk2,nk3,k1,k2,k3
c$$$      equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
c$$$      integer nmag,nmagf,magat(100),magatf(100),ix(10),nmbas
c$$$ckino Dec.8.2011:            real(8),pointer :: rv_p_opos(:) =>NULL()
c$$$
c$$$ckino Dec.8.2011:           real(8),pointer :: rv_p_og(:) =>NULL()
c$$$
c$$$ckino Dec.8.2011:          real(8),pointer :: rv_p_oag(:) =>NULL()
c$$$
c$$$ckino Dec.8.2011:         integer,pointer :: iv_p_oistab(:) =>NULL()
c$$$
c$$$      integer ,allocatable :: gstar_iv(:)
c$$$      real(8) ,allocatable :: wgt_rv(:)
c$$$
c$$$      double precision xx,qpi(3),ddot,plat(3,3),rb(3,3),qb(3,3),
c$$$     .pos(3,nbas)
c$$$      character dc*1, fn*120, outs*150, strn*120
c$$$      integer,allocatable:: ipq(:,:,:)
c$$$      real(8),allocatable:: qp(:,:),qfbz(:,:),emesh(:),wk(:)
c$$$     .,momsite(:), mmnorm(:)
c$$$     .,chi0_r(:,:,:,:),chi0_i(:,:,:,:), eiqrm_r(:,:), eiqrm_i(:,:)
c$$$      complex(8),allocatable:: chi0(:,:,:,:), eiqrm(:,:)
c$$$      integer iy,iz, i1,i2,i3,i4
c$$$C     Local parameters-product basis
c$$$      integer nbloch, natom,nmbas2, iqxini,iqxend,nw_i,nw2
c$$$     .,ifchipm,nx, iw
c$$$      integer,allocatable :: imbas(:),nblocha(:)
c$$$      real(8),allocatable:: svec(:,:)
c$$$      complex(8),allocatable:: zzr(:,:),zxq(:,:,:,:),x0meanx(:,:)
c$$$
c$$$C     For offset q mesh
c$$$      integer is(3),lshft(3),ifac(3)
c$$$C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
c$$$c      double precision qk
c$$$c      integer jj1,jj2,jj3,k
c$$$c      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
c$$$c     .                    (jj2*ifac(2)-1)*qb(k,2) +
c$$$c     .                    (jj3*ifac(3)-1)*qb(k,3)
c$$$C ... Heap
c$$$      integer ::iwdummy
c$$$
c$$$
c$$$Changenglob      if (nglob('nsp') .ne. 2)
c$$$      if (globalvariables%nsp .ne. 2)
c$$$     .call rx('chim editor for only for magnetic systems')
c$$$
c$$$C --- Setup ---
c$$$
c$$$      i_copy_size=size(slat%plat)
c$$$      call dcopy(i_copy_size,slat%plat,1,plat,1)
c$$$ckino Dec.8.2011:         rv_p_opos => slat%rv_p_opos
c$$$
c$$$
c$$$
c$$$      nsgrp=slat%nsgrp
c$$$ckino Dec.8.2011:         iv_p_oistab => slat%iv_p_oistab
c$$$
c$$$ckino Dec.8.2011:         rv_p_og => slat%rv_p_osymgr
c$$$
c$$$ckino Dec.8.2011:         rv_p_oag => slat%rv_p_oag
c$$$
c$$$
c$$$
c$$$C     should be from sgw?
c$$$
c$$$      i_copy_size=size(sbz%nkabc)
c$$$      call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
c$$$
c$$$ckino Dec.8.2011:         call dcopy ( 3 * nbas , rv_p_opos , 1 , pos , 1 )
c$$$ckino Jan.04.2012:         call dcopy ( 3 * nbas , slat%rv_p_opos , 1 , pos , 1 )
c$$$      call dcopy ( 3 * nbas , slat%rv_a_opos , 1 , pos , 1 )
c$$$
c$$$      call pshpr(0)
c$$$      llshft = .false.
c$$$      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
c$$$      call poppr
c$$$      lshft = 0
c$$$
c$$$C ... Defaults
c$$$      fn = 'rst1'
c$$$C     havechi0 = 1 when chi0(PP) read in, 2 when chi0(PB) read in
c$$$      havechi0 = 0
c$$$C     nmag = # magnetic sites
c$$$      nmag = 0
c$$$C     number of k-points
c$$$      nq = 0
c$$$C     true if chi saved on disk
c$$$      lnsave = .false.
c$$$C     lqp = has to do with specification qp.
c$$$C     0=>nothing specified, 1=>nq specifed, 2=>qp specified, 4=>qp=irr
c$$$      lqp = 0
c$$$
c$$$      dc = sopts(1:1)
c$$$      if (dc .ne. ' ') then
c$$$        print 301
c$$$  301   format(//' Entering the magnetic response function editor. ',
c$$$     .  'Parsing command-line options ...')
c$$$        lsopts = .true.
c$$$        js2 = 0
c$$$      else
c$$$        print 302
c$$$  302   format(//' Welcome to the magnetic response function editor.  ',
c$$$     .  'Enter ''?'' to see options.')
c$$$        lsopts = .false.
c$$$      endif
c$$$      if (mode .ne. 1) call rx('chimedit not ready for mode ne 1')
c$$$
c$$$C ... Return here to resume parsing for arguments
c$$$   10 continue
c$$$      if (lsopts) then
c$$$        js2 = js2+1
c$$$        if (js2 .gt. len(sopts)) then
c$$$          lsopts = .false.
c$$$          goto 10
c$$$        endif
c$$$        if (sopts(js2:js2) .eq. dc) goto 10
c$$$        js1 = min(len(sopts),js2)
c$$$        call nwordg(sopts,0,dc,1,js1,js2)
c$$$        if (js2 .lt. js1) lsopts = .false.
c$$$      endif
c$$$
c$$$C 306 format(' Failed to parse string ',a,' ... try again.')
c$$$  100 continue
c$$$#if PRTNOCR
c$$$      print '(/'' Option : '',$)'
c$$$#else
c$$$      print '(/'' Option : '')'
c$$$#endif
c$$$      outs = ' '
c$$$      if (lsopts) then
c$$$        print '(a)', trim(sopts(js1:js2))
c$$$        outs = sopts(js1:js2)
c$$$      else
c$$$        read(*,'(a150)') outs
c$$$      endif
c$$$      call locase(outs)
c$$$
c$$$C ... Parse and execute the next command
c$$$      if (.false.) then
c$$$
c$$$      elseif (outs .eq. ' ') then
c$$$        print 304
c$$$  304   format(' Enter ''q'' to exit, ''a'' to abort',
c$$$     .  ' ''?'' to see menu')
c$$$        goto 10
c$$$
c$$$C ... Number of k-points
c$$$      elseif (outs(1:4) .eq. 'new ') then
c$$$
c$$$        call words(outs,nw)
c$$$        if (nw .ne. 2) goto 98
c$$$        call word(outs,2,j1,j2)
c$$$        if (allocated(qp)) deallocate(qp)
c$$$        if (outs(j1:j2) .eq. 'irr') then
c$$$          mxkp = nk1*nk2*nk3
c$$$          if (allocated(qfbz)) deallocate(qfbz)
c$$$          allocate(qfbz(3,mxkp))
c$$$          allocate(gstar_iv(abs(-mxkp-1)))
c$$$          if (-mxkp-1<0) gstar_iv(:)=0
c$$$
c$$$          gstar_iv = - 2
c$$$
c$$$          if (allocated(ipq)) deallocate(ipq)
c$$$          allocate(ipq(nk1,nk2,nk3))
c$$$          allocate(wgt_rv(abs(-mxkp)))
c$$$          if (-mxkp<0) wgt_rv(:)=0.0d0
c$$$
c$$$          call info0(20,1,0,' q-points in full BZ ...')
c$$$          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
c$$$     .     , 0 , ipq , qfbz , wgt_rv , nq , mxkp , 0 , 0 )
c$$$
c$$$
c$$$          call dpzero ( wgt_rv , mxkp )
c$$$
c$$$          call info2(20,0,0,' Irr. qp ...',0,0)
c$$$          allocate(qp(3,nq))
c$$$ckino Dec.8.2011:             call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og
c$$$ckino Dec.28.2011:             call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_p_osymgr
c$$$          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_a_osymgr 
c$$$     .     , nsgrp , ipq , qp , wgt_rv , nq , mxkp , gstar_iv , 0 )
c$$$
c$$$
c$$$          lqp = 4
c$$$          call rx('irr still in progress')
c$$$        else
c$$$          j = 0
c$$$          j = a2vec(outs(j1:),len(outs(j1:)),j,2,', ',2,-3,1,ix,nq)
c$$$          allocate(qp(3,nq))
c$$$          lqp = 1
c$$$        endif
c$$$        call info2(0,0,0,'%4p... new chi0:  %i k-points',nq,0)
c$$$        nmag = 0
c$$$        havechi0 = 0
c$$$        lnsave = .true.
c$$$
c$$$
c$$$C --- Read chi0 ---
c$$$      elseif (outs(1:5) .eq. 'read ') then
c$$$
c$$$        if (nq .le. 0) then
c$$$          call info0(0,0,0,'%6p... "new" required before read')
c$$$          goto 98
c$$$        endif
c$$$
c$$$        call words(outs,nw)
c$$$        if (nw .lt. 2) goto 98
c$$$        call word(outs,2,j1,j2)
c$$$
c$$$C   --- Kotani style, full matrix ---
c$$$c run  fe --pr29 '--chimedit~new 1~read tkpb'
c$$$        if (outs(j1:j2) .eq. 'tkpb') then
c$$$c          call info0(0,0,0,'%6p... not ready for tkpb')
c$$$          do  iq = 1, nq
c$$$            fn = 'ChiPM0000.fmat'
c$$$            write(fn(6:9),'(i4)') iq       ! idummy
c$$$            do  j = 6, 9
c$$$              if (fn(j:j) .eq. ' ') fn(j:j) = '0'
c$$$            enddo
c$$$            ifi = fopng(trim(fn),-1,5)
c$$$            rewind ifi
c$$$c             open(333, file='ChiPM0001.fmat',form='unformatted' )
c$$$            read(ifi) nbloch,natom,nmbas, iqxini,iqxend,nw_i,nw2
c$$$            write(6,308) nbloch,natom,nmbas, iqxini,iqxend, nw_i,nw2
c$$$ 308        format (i4,6i8)
c$$$            if (natom.ne.nbas) call rx0('natom ne nbas')
c$$$
c$$$C           First qp: allocate and assign
c$$$
c$$$            nx=200
c$$$            if (iq .eq. 1) then
c$$$              if (allocated(momsite)) deallocate(momsite)
c$$$              if (allocated(mmnorm )) deallocate(mmnorm)
c$$$              allocate(imbas(nmbas), momsite(nmbas), mmnorm(nmbas)
c$$$     &        ,nblocha(1:natom),svec(1:nbloch,1:nmbas)
c$$$     &        ,zzr(nbloch,1) )
c$$$              allocate( zxq(nbloch,nbloch,nq,nx) )
c$$$              allocate( emesh(nx))
c$$$            endif
c$$$
c$$$            read(ifi) imbas(1:nmbas),momsite(1:nmbas),mmnorm(1:nmbas)
c$$$            read(ifi) nblocha(1:natom),svec(1:nbloch,1:nmbas)
c$$$            read(ifi) zzr(1:nbloch,1)
c$$$
c$$$            do iw=1,nx
c$$$              read(ifi)  qpi,emesh(iw), zxq(1:nbloch,1:nbloch,iq,iw)
c$$$            enddo
c$$$            qp(:,iq) = qpi
c$$$
c$$$C===read in matrix X0
c$$$          enddo
c$$$
c$$$          call stonerpb(nq,nx,nmbas,nbloch,qp,momsite,mmnorm,emesh,zxq)
c$$$
c$$$
c$$$          goto 98
c$$$
c$$$
c$$$C   --- Kotani style, rigid spin approximation ---
c$$$        elseif (outs(j1:j2) .eq. 'tkrs') then
c$$$
c$$$C     ... For each k-point, do
c$$$          do  iq = 1, nq
c$$$C           Make Kotani-style file name
c$$$            fn = 'ChiPM0000.nlfc.mat'
c$$$            write(fn(6:9),'(i4)') iq       ! idummy
c$$$            do  j = 6, 9
c$$$              if (fn(j:j) .eq. ' ') fn(j:j) = '0'
c$$$            enddo
c$$$
c$$$C           Open file, read number and list of magnetic sites
c$$$            ifi = fopng(trim(fn),-1,1)
c$$$            rewind ifi
c$$$            read(ifi,*) nmagf       ! number of magnetic atoms, 'nmbas
c$$$            call info2(0,0,0,'%4p... reading file '//trim(fn)//
c$$$     .      ', %i magnetic sites',nmagf,0)
c$$$            if (nmag .gt. 100) call rx('increase size of magat')
c$$$            read(ifi,*) magatf(1:nmagf)  ! magnetic atom index  'imbas
c$$$
c$$$C           Sanity check
c$$$            if (nmagf .lt. nmag) then
c$$$              call info2(0,0,0,
c$$$     .        '%8pabort: %i magnetic sites sought but file '//
c$$$     .        'contains only %i',nmag,nmagf)
c$$$              call fclr(' ',ifi)
c$$$              goto 10
c$$$            endif
c$$$
c$$$C           If magnetic sites not specified, take from file
c$$$            if (iq .eq. 1) then
c$$$              nmag = nmagf
c$$$              magat(1:nmagf) = magatf(1:nmagf)      ! magat- first q point
c$$$            endif
c$$$C           Sanity check for subsequent qp
c$$$            if (isanrg(nmag,nmagf,nmagf,'        abort: ',
c$$$     .      'nmag',.false.)) goto 10
c$$$            do  j = 1, max(nmag,nmagf)
c$$$              if (magat(j) .ne. magatf(j)) then
c$$$                call info0(0,0,0,
c$$$     .          '%8pabort: magnetic site list does not match file:')
c$$$                print 345, ' sought:',magat(1:nmag)
c$$$                print 345, ' file:',  magat(1:nmagf)
c$$$  345           format(a10,100i4)
c$$$                call fclr(' ',ifi)
c$$$                goto 10
c$$$              endif
c$$$            enddo
c$$$
c$$$C           Get energy mesh ; put into wk
c$$$            allocate(wk(100000))
c$$$C           Skip next 3 lines
c$$$
c$$$            read(ifi,*) xx; read(ifi,*) xx; read(ifi,*) xx   !neglect the first 3 lines
c$$$            j = 0
c$$$            do while (.true.)
c$$$              read(ifi,*,end=30,err=30) qpi,wk(j+1)   ! readin q and omega
c$$$              j = j + 1
c$$$            enddo
c$$$   30       continue
c$$$            if (iq .eq. 1) then
c$$$              if(sum(abs(qpi(:))).ne.0d0) call rx('1st /= Gamma pnt')
c$$$            endif
c$$$
c$$$C           First qp: allocate and assign
c$$$            if (iq .eq. 1) then
c$$$              nmbas=nmag
c$$$              ne = j
c$$$              allocate(emesh(ne))
c$$$              emesh(1:ne) = wk(1:ne)
c$$$              allocate( momsite(nmagf), mmnorm(nmagf) )
c$$$              allocate( eiqrm(nmagf,nq),
c$$$     .        eiqrm_r(nmagf,nq),eiqrm_i(nmagf,nq) )
c$$$              allocate(chi0(nmagf,nmagf,nq,ne)
c$$$     .        ,chi0_r(nmagf,nmagf,nq,ne), chi0_i(nmagf,nmagf,nq,ne) )
c$$$
c$$$c              allocate( ,e1(nmagf),e2(nmagf))
c$$$            endif
c$$$C           Subsequent qp: assignments and sanity checks
c$$$            if (ne .ne. j) then
c$$$              call info5(0,0,0,'%8pabort, qp %i:  expected %i '//
c$$$     .        'energy points from but read %i',iq,ne,j,0,0)
c$$$              goto 10
c$$$            endif
c$$$            call daxpy(ne,-1d0,emesh,1,wk,1)
c$$$            if (ddot(ne,wk,1,wk,1) .gt. 1d-10) then
c$$$              call info2(0,0,0,'%8pabort, qp %i:  energy '//
c$$$     .        'mesh does not match first qp',iq,0)
c$$$              goto 10
c$$$            endif
c$$$            deallocate(wk)
c$$$            qp(:,iq) = qpi
c$$$            call info2(0,0,0,'%8pread %i energy points, qp=%3;11,6D',
c$$$     .      ne,qpi)
c$$$
c$$$C           Read chi0
c$$$            rewind ifi
c$$$            read(ifi,*) j; read(ifi,*) j
c$$$            read(ifi,*) momsite(1:nmagf)
c$$$            read(ifi,*) mmnorm(1:nmagf)
c$$$            read(ifi,*) (eiqrm_r(iy,iq),eiqrm_i(iy,iq),iy=1,nmbas)
c$$$            eiqrm(:,iq) = dcmplx(eiqrm_r(:,iq),eiqrm_i(:,iq) ) ! <e(iqr)|m>.
c$$$            do j = 1, ne
c$$$              read(ifi,*,end=40,err=40) qpi,xx , ( (chi0_r(iy,iz,iq,j)
c$$$     &        ,chi0_i(iy,iz,iq,j),iy=1,nmbas),iz=1,nmbas)
c$$$              chi0(:,:,iq,j)= dcmplx(chi0_r(:,:,iq,j),chi0_i(:,:,iq,j))
c$$$            enddo
c$$$
c$$$ 40         continue
c$$$
c$$$C           Cleanup for this qp
c$$$            call fclr(' ',ifi)
c$$$          enddo   ! end of iq loop
c$$$
c$$$c     check output
c$$$c          write(6,"(255i5)") nmagf
c$$$c          write(6,"(255i5)") magat(1:nmagf)
c$$$c          write(6,"(255d23.15)") momsite(1:nmagf)
c$$$c          write(6,"(255d23.15)")  mmnorm(1:nmagf)
c$$$c          write(6,*)
c$$$c          do  iq = 1, nq
c$$$c             write(6,"(255d23.15)")  eiqrm(1:nmagf,iq)
c$$$c          enddo
c$$$c          write(6,*)
c$$$c          do  iq = 1, 1
c$$$c             do j =1,ne
c$$$c                write(6,"(255d23.15)") emesh(j), chi0(:,:,iq,j)
c$$$c             enddo
c$$$c          enddo
c$$$
c$$$          do i =1,nmagf
c$$$            if (abs(momsite(i)).lt.1d-3) call rx('non-magnetic site ?')
c$$$          enddo
c$$$
c$$$          call stonerrsa(nq,ne,nmagf,qp,momsite,mmnorm,eiqrm,emesh,chi0)
c$$$
c$$$          do 45 iq = 1, nq
c$$$          do 45 j =1,ne
c$$$          do 45 i1 = 1,nmbas
c$$$          do 45 i2 = 1,nmbas
c$$$            chi0(i1,i2,iq,j)=              !D_{i,j}=<ebar_i|X0|ebar_j>
c$$$     .            chi0(i1,i2,iq,j)*momsite(i1)*momsite(i2)
c$$$     .            /mmnorm(i1) /mmnorm(i2)
c$$$ 45       continue
c$$$
c$$$
c$$$c omg
c$$$
c$$$
c$$$
c$$$c     ./lmf lsmo56 --chimedit
c$$$c
c$$$C         chi0 has been read: cleanup
c$$$          havechi0 = 1
c$$$
c$$$        else
c$$$          goto 98
c$$$
c$$$        endif
c$$$
c$$$C ... Specify list of magnetic sites
c$$$      elseif (outs(1:9) .eq. 'magsites ') then
c$$$
c$$$        if (havechi0 .eq. 0) then
c$$$          call info0(0,0,0,'%6p... "read" required before magsites')
c$$$          goto 98
c$$$        endif
c$$$
c$$$        call words(outs,nw)
c$$$        if (nw .ne. 2) goto 98
c$$$        call word(outs,2,j1,j2)
c$$$        call mkils0(outs(j1:j2),nmag,j)
c$$$        if (nmag .le. 0) then
c$$$          call info0(0,0,0,'%6p... Bad or null list : '//outs(j1:j2))
c$$$          nmag = 0
c$$$          goto 98
c$$$        endif
c$$$        if (nmag .gt. 100) call rx('increase size of magat')
c$$$        call mkilst(outs(j1:j2),nmag,magat)
c$$$        call ilst2a(magat,nmag,strn)
c$$$        call info2(0,0,0,'%3p... %i magnetic site%-1j%?#n==1##s#:  '//
c$$$     .  trim(strn),nmag,0)
c$$$C       Sanity check
c$$$        do  i = 1, nmag
c$$$          j = iinear(nmagf,magat(i),magatf,1)
c$$$          if (magat(i) .ne. magatf(j)) then
c$$$            call ilst2a(magatf,nmagf,strn)
c$$$            call info(0,0,0,'%7pabort, site %i is not among chi0 '//
c$$$     .      'list: '//trim(strn)//' .. restore chi0 list',magat(i),0)
c$$$            nmag = nmagf
c$$$            magat(1:nmagf) = magatf(1:nmagf)
c$$$            goto 10
c$$$          endif
c$$$        enddo
c$$$        call rx('magsites still in progress')
c$$$
c$$$C ... show
c$$$      elseif (outs(1:5) .eq. 'show ') then
c$$$        if (havechi0 .eq. 0) then
c$$$          call info2(0,0,0,' ... no chi0 read, '//
c$$$     .    '%?#n==0#no k-points specified#'//
c$$$     .    '%-1jwaiting to read %i k-points',
c$$$     .    nq,0)
c$$$        elseif (havechi0 .eq. 1) then
c$$$          call info5(0,0,0,' ... chi0 read, '//
c$$$     .    '%i site%-1j%?#n==1##s#:  nq=%i  ne=%i  emax=%;4d Ry',
c$$$     .    nmag,nq,ne,emesh(ne),0)
c$$$        else
c$$$          call rx('not ready for show')
c$$$        endif
c$$$
c$$$C ... Save
c$$$      elseif (outs(1:5) .eq. 'save ' .or. outs(1:6) .eq. 'savea ') then
c$$$        call rx('not ready for save')
c$$$C        lbin = outs(1:5) .eq. 'save '
c$$$C        lnsave = .false.
c$$$
c$$$C ... abort
c$$$      elseif (outs(1:2) .eq. 'a ') then
c$$$        call rx0('aborting chi editor ... no file written')
c$$$
c$$$C ... quit
c$$$      elseif (outs(1:2) .eq. 'q '. or. outs(1:5) .eq. 'quit ') then
c$$$        if (lnsave .and. havechi0 .gt. 0) then
c$$$          print '('' chipm file not saved ... really quit?'')'
c$$$          read(*,'(a150)') outs
c$$$          call locase(outs)
c$$$          if (.not. (outs(1:1) .eq. 'y' .or. outs(1:1) .eq. 'q'))
c$$$     .    goto 10
c$$$        endif
c$$$        call rx0('exit chi editor')
c$$$
c$$$C ... help
c$$$      elseif (outs .eq. '?') then
c$$$        print 310
c$$$        print 311
c$$$  310   format(
c$$$     .  ' Select one of these options:'/
c$$$     .  t4,'new nk|irr',t21,
c$$$     .  'New chi0: specify number of k-points to read.'/t21,
c$$$     .  'Optional irr => irreducible points'/
c$$$     .  t4,'magsites list',t21,'specify magnetic sites in basis'//
c$$$     .  t4,'read tk|tkrs',t21,
c$$$     .  'read chi0.  new must be input first.'/t21,
c$$$     .  'tk   => Kotani style, full matrix'/t21,
c$$$     .  'tkrs => Kotani style, rigid spin approximation')
c$$$
c$$$  311   format(/
c$$$     .  t4,'show',t21, 'Show summary information')
c$$$
c$$$      else
c$$$        print '(1x,''unrecognized option: '',a)', trim(outs)
c$$$
c$$$      endif
c$$$      goto 10
c$$$
c$$$   98 call info0(0,0,0,' chimedit:  improper usage of '//trim(outs)//
c$$$     .' ... nothing done')
c$$$      goto 10
c$$$
c$$$      if (allocated(wgt_rv)) deallocate(wgt_rv)
c$$$      if (allocated(gstar_iv)) deallocate(gstar_iv)
c$$$
c$$$      end subroutine chimedit
c$$$
c$$$
