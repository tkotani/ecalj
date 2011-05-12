      subroutine chimedit(sopts,mode,sctrl,ssite,sspec,slat,spot,sbz,
     .nbas,nat,nspec)

      use m_struc_def  !Cgetarg
      use m_globalvariables
C- Magnetic Linear response editor
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Ci   sopts :command options performed automatically, before reading
Ci         :from standard input
Ci   mode  :0 longitudinal linear response, lmf
Ci         :1 transverse linear response, lmf
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: zbak lrs
Ci     Stored:    lrs
Ci     Passed to:
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read:
Ci     Stored:
Ci     Passed to: dcopy
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:
Ci     Stored:
Ci     Passed to: dcopy
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc
Ci     Stored:
Ci     Passed to: dcopy rdovfa
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: nlml nlma
Ci     Stored:
Ci     Passed to: dcopy rdovfa
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read:
Ci     Stored:
Ci     Passed to:
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Co Outputs
Co   chiedit never returns.
Co   rst file can be written.
Cr Remarks
Cr   The density consists of a smooth part (smrho) plus
Cr   nbas atom-centered densities inside the MT spheres.
Cr   Their sum is the full charge density.
Cr   The local density is represented as the difference of the
Cr   two valence components in orhoat, plus the core density.
Cr   Density in the MT spheres:
Cr      mesh parameters rmt,nr,a;
Cr      total density rho (times r**2) to lmxl;
Cr      a spherical potential v0 defining the wave functions within rmt
Cr      pnu and idmod to lmxa
Cr   Smooth density
Cr      real part of complex*16 array smrho contains the density
Cr      k1,k2,k3 are the physical dimensions of the array
Cr      n1,n2,n3 are the dimensions of the mesh.
Cl Local variables
Cu Updates
Cu   06 Dec 08 First created
C  ----------------------------------------------------------------------
Ci Inputs
Ci   nx       : cutoff number of w points in calculation, such  X0(q,w)
Ci   cutw     : cutoff value(in meV)  of w in calculation of X(q,w)
Ci   nbloch   : dimension of product basis
Ci   ngbin    : dimension of mix basis
Ci   nq       : total number of q points in calculation of X(q,w)
Ci   qpi      : vector of current q,   X(q,w)
Ci   qgbin    : vector of all q points,   qgbin(3,nq0i)
Ci   imbas    : magnetic atom index
Ci   momsite  : magnetic moment m
Ci   svec     : <B_I|m_a>  svec(nbloch,nmbas)
Ci   mmnorm   : |m|=sqrt(<m|B> <B|m> )=norm(svec) ; this is not m.
Ci            : ~e_a(r)=m_a(r) / sqrt(\int m_a(r)**2 dr)  then <~e_a|~e_a>=1
Ci            : <B_i|~e_a>=svec/mmnorm
Ci   biubj    : <B|U|B>  stoner matrix U(nbloch,nbloch)
Ci   zxq      : <B|X0|B>   chi0 matrix zxq(nbloch,nbloch,nx) denpend on (q,w)
Ci   chi0     : <~e|X0|~e> = <m|B><B|X0|B><B|m> / norm(<m_i|m_i>) /norm(<m_j|m_j>)
Ci            :            = <m|x0|m'>/ norm(<m_i|m_i>) /norm(<B|m_j>)
Ci            :            = <m_a|x0|m_a'>/(sqrt(<m_i|m_i> sqrt(<m_j|m_j>)
Ci            :     |~e_a> = |m_a>/<m_a|m_a>
Ci   gbvec    : <e^iqr|B+>  projection of |e^iqr> on mix basis set,  gbvec(ngbin,nq0i)
Ci   gbvec0   : <e^iqr|B>   subset(prod. bas. part) of gbvec,  gbvec0(nbloch,nq0i)
Ci   eiqrm    : eiqrm_i = <~e_i|e^{iqr}> =  <M_i|eiqr>/sqrt(<M_i|M_i>)
Ci            : sum( dconjg(gbvec(1:nbloch))*svec(1:nbloch,imb) )/mmnorm(imb)
Ci   freq     : frequency w in unit of ryberg
Ci   freq_meV : frequency w in unit of meV
Ci   mx2m     : <m|x|m>  from rigid spin approximation
Ci   mx0m     : <m|x0|m> from rigid spin approximation
Co Outputs:
Co   mim      : <m|i|m>  stoner parameter from rigid spin approximation
Co   x0meanx  : <m|B><B|X0|B><B|m>
Co   mUm      : <m|B><B|U|B><B|m>
Co   chike    : <m|B><B|X|B><B|m>
Co   iqrXiqr  : <e^{iqr}|B><B|X|B><B|e^{iqr}>


C     implicit none
C ... Passed parameters
      character sopts*(*)
      integer mode,nbas,nat,nspec,ifi,n0
      parameter (n0=10)
      type(s_ctrl)::sctrl
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_pot)::spot
      type(s_bz)::sbz

C ... Local parameters
      integer nglob,fopna,fopng,a2vec,havechi0,iq,nq,nsgrp
      integer i,j,l,im,m,n,j1,j2,js1,js2,nw,ne,iinear,lqp

      logical lnsave,lsopts,isanrg,llshft(3)
      integer nkxyz(3),mxkp,nk1,nk2,nk3,k1,k2,k3
      equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
      integer nmag,nmagf,magat(100),magatf(100),ix(10),nmbas
          real(8),pointer :: rv_p_opos(:) =>NULL()

         real(8),pointer :: rv_p_og(:) =>NULL()

        real(8),pointer :: rv_p_oag(:) =>NULL()

       integer,pointer :: iv_p_oistab(:) =>NULL()

      integer ,allocatable :: gstar_iv(:)
      real(8) ,allocatable :: wgt_rv(:)

      double precision xx,qpi(3),ddot,plat(3,3),rb(3,3),qb(3,3),
     .pos(3,nbas)
      character dc*1, fn*120, outs*150, strn*120
      integer,allocatable:: ipq(:,:,:)
      real(8),allocatable:: qp(:,:),qfbz(:,:),emesh(:),wk(:)
     .,momsite(:), mmnorm(:)
     .,chi0_r(:,:,:,:),chi0_i(:,:,:,:), eiqrm_r(:,:), eiqrm_i(:,:)
      complex(8),allocatable:: chi0(:,:,:,:), eiqrm(:,:)
      integer iy,iz, i1,i2,i3,i4
C     Local parameters-product basis
      integer nbloch, natom,nmbas2, iqxini,iqxend,nw_i,nw2
     .,ifchipm,nx, iw
      integer,allocatable :: imbas(:),nblocha(:)
      real(8),allocatable:: svec(:,:)
      complex(8),allocatable:: zzr(:,:),zxq(:,:,:,:),x0meanx(:,:)

C     For offset q mesh
      integer is(3),lshft(3),ifac(3)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
c      double precision qk
c      integer jj1,jj2,jj3,k
c      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
c     .                    (jj2*ifac(2)-1)*qb(k,2) +
c     .                    (jj3*ifac(3)-1)*qb(k,3)
C ... Heap
       integer ::iwdummy 


Changenglob      if (nglob('nsp') .ne. 2)
      if (globalvariables%nsp .ne. 2)
     .call rx('chim editor for only for magnetic systems')

C --- Setup ---

      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
       rv_p_opos => slat%rv_p_opos 



      nsgrp=slat%nsgrp
       iv_p_oistab => slat%iv_p_oistab 

       rv_p_og => slat%rv_p_osymgr 

       rv_p_oag => slat%rv_p_oag 



C     should be from sgw?

      i_copy_size=size(sbz%nkabc)
      call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)

       call dcopy ( 3 * nbas , rv_p_opos , 1 , pos , 1 ) 

      call pshpr(0)
      llshft = .false.
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
      call poppr
      lshft = 0

C ... Defaults
      fn = 'rst1'
C     havechi0 = 1 when chi0(PP) read in, 2 when chi0(PB) read in
      havechi0 = 0
C     nmag = # magnetic sites
      nmag = 0
C     number of k-points
      nq = 0
C     true if chi saved on disk
      lnsave = .false.
C     lqp = has to do with specification qp.
C     0=>nothing specified, 1=>nq specifed, 2=>qp specified, 4=>qp=irr
      lqp = 0

      dc = sopts(1:1)
      if (dc .ne. ' ') then
        print 301
  301   format(//' Entering the magnetic response function editor. ',
     .  'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the magnetic response function editor.  ',
     .  'Enter ''?'' to see options.')
        lsopts = .false.
      endif
      if (mode .ne. 1) call rx('chimedit not ready for mode ne 1')

C ... Return here to resume parsing for arguments
   10 continue
      if (lsopts) then
        js2 = js2+1
        if (js2 .gt. len(sopts)) then
          lsopts = .false.
          goto 10
        endif
        if (sopts(js2:js2) .eq. dc) goto 10
        js1 = min(len(sopts),js2)
        call nwordg(sopts,0,dc,1,js1,js2)
        if (js2 .lt. js1) lsopts = .false.
      endif

C 306 format(' Failed to parse string ',a,' ... try again.')
  100 continue
#if PRTNOCR
      print '(/'' Option : '',$)'
#else
      print '(/'' Option : '')'
#endif
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif
      call locase(outs)

C ... Parse and execute the next command
      if (.false.) then

      elseif (outs .eq. ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',
     .  ' ''?'' to see menu')
        goto 10

C ... Number of k-points
      elseif (outs(1:4) .eq. 'new ') then

        call words(outs,nw)
        if (nw .ne. 2) goto 98
        call word(outs,2,j1,j2)
        if (allocated(qp)) deallocate(qp)
        if (outs(j1:j2) .eq. 'irr') then
          mxkp = nk1*nk2*nk3
          if (allocated(qfbz)) deallocate(qfbz)
          allocate(qfbz(3,mxkp))
          allocate(gstar_iv(abs(-mxkp-1)))
          if (-mxkp-1<0) gstar_iv(:)=0

          gstar_iv = - 2

          if (allocated(ipq)) deallocate(ipq)
          allocate(ipq(nk1,nk2,nk3))
          allocate(wgt_rv(abs(-mxkp)))
          if (-mxkp<0) wgt_rv(:)=0.0d0

          call info0(20,1,0,' q-points in full BZ ...')
           call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
     .     , 0 , ipq , qfbz , wgt_rv , nq , mxkp , 0 , 0 ) 


          call dpzero ( wgt_rv , mxkp )

          call info2(20,0,0,' Irr. qp ...',0,0)
          allocate(qp(3,nq))
           call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og 
     .     , nsgrp , ipq , qp , wgt_rv , nq , mxkp , gstar_iv , 0 ) 


          lqp = 4
          call rx('irr still in progress')
        else
          j = 0
          j = a2vec(outs(j1:),len(outs(j1:)),j,2,', ',2,-3,1,ix,nq)
          allocate(qp(3,nq))
          lqp = 1
        endif
        call info2(0,0,0,'%4p... new chi0:  %i k-points',nq,0)
        nmag = 0
        havechi0 = 0
        lnsave = .true.


C --- Read chi0 ---
      elseif (outs(1:5) .eq. 'read ') then

        if (nq .le. 0) then
          call info0(0,0,0,'%6p... "new" required before read')
          goto 98
        endif

        call words(outs,nw)
        if (nw .lt. 2) goto 98
        call word(outs,2,j1,j2)

C   --- Kotani style, full matrix ---
c run  fe --pr29 '--chimedit~new 1~read tkpb'
        if (outs(j1:j2) .eq. 'tkpb') then
c          call info0(0,0,0,'%6p... not ready for tkpb')
          do  iq = 1, nq
            fn = 'ChiPM0000.fmat'
            write(fn(6:9),'(i4)') iq       ! idummy
            do  j = 6, 9
              if (fn(j:j) .eq. ' ') fn(j:j) = '0'
            enddo
            ifi = fopng(trim(fn),-1,5)
            rewind ifi
c             open(333, file='ChiPM0001.fmat',form='unformatted' )
            read(ifi) nbloch,natom,nmbas, iqxini,iqxend,nw_i,nw2
            write(6,308) nbloch,natom,nmbas, iqxini,iqxend, nw_i,nw2
 308        format (i4,6i8)
            if (natom.ne.nbas) call rx0('natom ne nbas')

C           First qp: allocate and assign

            nx=200
            if (iq .eq. 1) then
              if (allocated(momsite)) deallocate(momsite)
              if (allocated(mmnorm )) deallocate(mmnorm)
              allocate(imbas(nmbas), momsite(nmbas), mmnorm(nmbas)
     &        ,nblocha(1:natom),svec(1:nbloch,1:nmbas)
     &        ,zzr(nbloch,1) )
              allocate( zxq(nbloch,nbloch,nq,nx) )
              allocate( emesh(nx))
            endif

            read(ifi) imbas(1:nmbas),momsite(1:nmbas),mmnorm(1:nmbas)
            read(ifi) nblocha(1:natom),svec(1:nbloch,1:nmbas)
            read(ifi) zzr(1:nbloch,1)

            do iw=1,nx
              read(ifi)  qpi,emesh(iw), zxq(1:nbloch,1:nbloch,iq,iw)
            enddo
            qp(:,iq) = qpi

C===read in matrix X0
          enddo

          call stonerpb(nq,nx,nmbas,nbloch,qp,momsite,mmnorm,emesh,zxq)


          goto 98


C   --- Kotani style, rigid spin approximation ---
        elseif (outs(j1:j2) .eq. 'tkrs') then

C     ... For each k-point, do
          do  iq = 1, nq
C           Make Kotani-style file name
            fn = 'ChiPM0000.nlfc.mat'
            write(fn(6:9),'(i4)') iq       ! idummy
            do  j = 6, 9
              if (fn(j:j) .eq. ' ') fn(j:j) = '0'
            enddo

C           Open file, read number and list of magnetic sites
            ifi = fopng(trim(fn),-1,1)
            rewind ifi
            read(ifi,*) nmagf       ! number of magnetic atoms, 'nmbas
            call info2(0,0,0,'%4p... reading file '//trim(fn)//
     .      ', %i magnetic sites',nmagf,0)
            if (nmag .gt. 100) call rx('increase size of magat')
            read(ifi,*) magatf(1:nmagf)  ! magnetic atom index  'imbas

C           Sanity check
            if (nmagf .lt. nmag) then
              call info2(0,0,0,
     .        '%8pabort: %i magnetic sites sought but file '//
     .        'contains only %i',nmag,nmagf)
              call fclr(' ',ifi)
              goto 10
            endif

C           If magnetic sites not specified, take from file
            if (iq .eq. 1) then
              nmag = nmagf
              magat(1:nmagf) = magatf(1:nmagf)      ! magat- first q point
            endif
C           Sanity check for subsequent qp
            if (isanrg(nmag,nmagf,nmagf,'        abort: ',
     .      'nmag',.false.)) goto 10
            do  j = 1, max(nmag,nmagf)
              if (magat(j) .ne. magatf(j)) then
                call info0(0,0,0,
     .          '%8pabort: magnetic site list does not match file:')
                print 345, ' sought:',magat(1:nmag)
                print 345, ' file:',  magat(1:nmagf)
  345           format(a10,100i4)
                call fclr(' ',ifi)
                goto 10
              endif
            enddo

C           Get energy mesh ; put into wk
            allocate(wk(100000))
C           Skip next 3 lines

            read(ifi,*) xx; read(ifi,*) xx; read(ifi,*) xx   !neglect the first 3 lines
            j = 0
            do while (.true.)
              read(ifi,*,end=30,err=30) qpi,wk(j+1)   ! readin q and omega
              j = j + 1
            enddo
   30       continue
            if (iq .eq. 1) then
              if(sum(abs(qpi(:))).ne.0d0) call rx('1st /= Gamma pnt')
            endif

C           First qp: allocate and assign
            if (iq .eq. 1) then
              nmbas=nmag
              ne = j
              allocate(emesh(ne))
              emesh(1:ne) = wk(1:ne)
              allocate( momsite(nmagf), mmnorm(nmagf) )
              allocate( eiqrm(nmagf,nq),
     .        eiqrm_r(nmagf,nq),eiqrm_i(nmagf,nq) )
              allocate(chi0(nmagf,nmagf,nq,ne)
     .        ,chi0_r(nmagf,nmagf,nq,ne), chi0_i(nmagf,nmagf,nq,ne) )

c              allocate( ,e1(nmagf),e2(nmagf))
            endif
C           Subsequent qp: assignments and sanity checks
            if (ne .ne. j) then
              call info5(0,0,0,'%8pabort, qp %i:  expected %i '//
     .        'energy points from but read %i',iq,ne,j,0,0)
              goto 10
            endif
            call daxpy(ne,-1d0,emesh,1,wk,1)
            if (ddot(ne,wk,1,wk,1) .gt. 1d-10) then
              call info2(0,0,0,'%8pabort, qp %i:  energy '//
     .        'mesh does not match first qp',iq,0)
              goto 10
            endif
            deallocate(wk)
            qp(:,iq) = qpi
            call info2(0,0,0,'%8pread %i energy points, qp=%3;11,6D',
     .      ne,qpi)

C           Read chi0
            rewind ifi
            read(ifi,*) j; read(ifi,*) j
            read(ifi,*) momsite(1:nmagf)
            read(ifi,*) mmnorm(1:nmagf)
            read(ifi,*) (eiqrm_r(iy,iq),eiqrm_i(iy,iq),iy=1,nmbas)
            eiqrm(:,iq) = dcmplx(eiqrm_r(:,iq),eiqrm_i(:,iq) ) ! <e(iqr)|m>.
            do j = 1, ne
              read(ifi,*,end=40,err=40) qpi,xx , ( (chi0_r(iy,iz,iq,j)
     &        ,chi0_i(iy,iz,iq,j),iy=1,nmbas),iz=1,nmbas)
              chi0(:,:,iq,j)= dcmplx(chi0_r(:,:,iq,j),chi0_i(:,:,iq,j))
            enddo

 40         continue

C           Cleanup for this qp
            call fclr(' ',ifi)
          enddo   ! end of iq loop

c     check output
c          write(6,"(255i5)") nmagf
c          write(6,"(255i5)") magat(1:nmagf)
c          write(6,"(255d23.15)") momsite(1:nmagf)
c          write(6,"(255d23.15)")  mmnorm(1:nmagf)
c          write(6,*)
c          do  iq = 1, nq
c             write(6,"(255d23.15)")  eiqrm(1:nmagf,iq)
c          enddo
c          write(6,*)
c          do  iq = 1, 1
c             do j =1,ne
c                write(6,"(255d23.15)") emesh(j), chi0(:,:,iq,j)
c             enddo
c          enddo

          do i =1,nmagf
            if (abs(momsite(i)).lt.1d-3) call rx('non-magnetic site ?')
          enddo

          call stonerrsa(nq,ne,nmagf,qp,momsite,mmnorm,eiqrm,emesh,chi0)

          do 45 iq = 1, nq
            do 45 j =1,ne
              do 45 i1 = 1,nmbas
                do 45 i2 = 1,nmbas
                  chi0(i1,i2,iq,j)=              !D_{i,j}=<ebar_i|X0|ebar_j>
     .            chi0(i1,i2,iq,j)*momsite(i1)*momsite(i2)
     .            /mmnorm(i1) /mmnorm(i2)
 45       continue


c omg



c     ./lmf lsmo56 --chimedit
c
C         chi0 has been read: cleanup
          havechi0 = 1

        else
          goto 98

        endif

C ... Specify list of magnetic sites
      elseif (outs(1:9) .eq. 'magsites ') then

        if (havechi0 .eq. 0) then
          call info0(0,0,0,'%6p... "read" required before magsites')
          goto 98
        endif

        call words(outs,nw)
        if (nw .ne. 2) goto 98
        call word(outs,2,j1,j2)
        call mkils0(outs(j1:j2),nmag,j)
        if (nmag .le. 0) then
          call info0(0,0,0,'%6p... Bad or null list : '//outs(j1:j2))
          nmag = 0
          goto 98
        endif
        if (nmag .gt. 100) call rx('increase size of magat')
        call mkilst(outs(j1:j2),nmag,magat)
        call ilst2a(magat,nmag,strn)
        call info2(0,0,0,'%3p... %i magnetic site%-1j%?#n==1##s#:  '//
     .  trim(strn),nmag,0)
C       Sanity check
        do  i = 1, nmag
          j = iinear(nmagf,magat(i),magatf,1)
          if (magat(i) .ne. magatf(j)) then
            call ilst2a(magatf,nmagf,strn)
            call info(0,0,0,'%7pabort, site %i is not among chi0 '//
     .      'list: '//trim(strn)//' .. restore chi0 list',magat(i),0)
            nmag = nmagf
            magat(1:nmagf) = magatf(1:nmagf)
            goto 10
          endif
        enddo
        call rx('magsites still in progress')

C ... show
      elseif (outs(1:5) .eq. 'show ') then
        if (havechi0 .eq. 0) then
          call info2(0,0,0,' ... no chi0 read, '//
     .    '%?#n==0#no k-points specified#'//
     .    '%-1jwaiting to read %i k-points',
     .    nq,0)
        elseif (havechi0 .eq. 1) then
          call info5(0,0,0,' ... chi0 read, '//
     .    '%i site%-1j%?#n==1##s#:  nq=%i  ne=%i  emax=%;4d Ry',
     .    nmag,nq,ne,emesh(ne),0)
        else
          call rx('not ready for show')
        endif

C ... Save
      elseif (outs(1:5) .eq. 'save ' .or. outs(1:6) .eq. 'savea ') then
        call rx('not ready for save')
C        lbin = outs(1:5) .eq. 'save '
C        lnsave = .false.

C ... abort
      elseif (outs(1:2) .eq. 'a ') then
        call rx0('aborting chi editor ... no file written')

C ... quit
      elseif (outs(1:2) .eq. 'q '. or. outs(1:5) .eq. 'quit ') then
        if (lnsave .and. havechi0 .gt. 0) then
          print '('' chipm file not saved ... really quit?'')'
          read(*,'(a150)') outs
          call locase(outs)
          if (.not. (outs(1:1) .eq. 'y' .or. outs(1:1) .eq. 'q'))
     .    goto 10
        endif
        call rx0('exit chi editor')

C ... help
      elseif (outs .eq. '?') then
        print 310
        print 311
  310   format(
     .  ' Select one of these options:'/
     .  t4,'new nk|irr',t21,
     .  'New chi0: specify number of k-points to read.'/t21,
     .  'Optional irr => irreducible points'/
     .  t4,'magsites list',t21,'specify magnetic sites in basis'//
     .  t4,'read tk|tkrs',t21,
     .  'read chi0.  new must be input first.'/t21,
     .  'tk   => Kotani style, full matrix'/t21,
     .  'tkrs => Kotani style, rigid spin approximation')

  311   format(/
     .  t4,'show',t21, 'Show summary information')

      else
        print '(1x,''unrecognized option: '',a)', trim(outs)

      endif
      goto 10

   98 call info0(0,0,0,' chimedit:  improper usage of '//trim(outs)//
     .' ... nothing done')
      goto 10

      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)

      end subroutine chimedit


