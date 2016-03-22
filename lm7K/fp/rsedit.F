c$$$ commented out on 22mar2015

c---> takao modified this by hand after Kino's del_pack.
c   python getarg.py < file.F > x
c   mv x file.F

c$$$      subroutine rsedit(sopts,mode,sctrl,ssite,sspec,slat,spot,sbz,
c$$$     .     nbas,nat,nspec)
c$$$
c$$$      use m_struc_def           !Cgetarg
c$$$      use m_globalvariables
c$$$C     - Restart file editor
c$$$C     ----------------------------------------------------------------------
c$$$C     i Inputs/Outputs
c$$$C     i   sopts :command options performed automatically, before reading
c$$$C     i         :from standard input
c$$$C     i   mode  :0 ASA I/O
c$$$C     i         :1 FP  I/O, mesh density
c$$$C     i         :-1 (input only) read nbas,nat,nspec from file
c$$$C     i   sctrl :struct for program flow parameters; see routine uctrl
c$$$C     i     Elts read: zbak lrs
c$$$C     i     Stored:    lrs
c$$$C     i     Passed to: iors
c$$$C     i   ssite :struct for site-specific information; see routine usite
c$$$C     i     Elts read:
c$$$C     i     Stored:
c$$$C     i     Passed to: dcopy
c$$$C     i   sspec :struct for species-specific information; see routine uspec
c$$$C     i     Elts read:
c$$$C     i     Stored:
c$$$C     i     Passed to: dcopy
c$$$C     i   slat  :struct for lattice information; see routine ulat
c$$$C     i     Elts read: nabc
c$$$C     i     Stored:
c$$$C     i     Passed to: dcopy rdovfa iors
c$$$C     i   spot  :struct for information about the potential; see routine upot
c$$$C     i     Elts read: nlml nlma osmrho osmpot oorhat
c$$$C     i     Stored:    osmrho
c$$$C     i     Passed to: dcopy rdovfa iors
c$$$C     i   sbz   :struct for the Brillouin Zone; see routine ubz
c$$$C     i     Elts read:
c$$$C     i     Stored:
c$$$C     i     Passed to: iors
c$$$C     i   nbas  :size of basis
c$$$C     i   nat   :number atoms in basis with augmentation sites
c$$$C     i         :Note: if nat<nbas, there is a requirement that
c$$$C     i         :lmxa>-1 for nat sites, and
c$$$C     i         :and lmxa=-1 for nbas-nat sites
c$$$C     o Outputs
c$$$C     o   rsedit never returns.
c$$$C     o   rst file can be written.
c$$$C     r Remarks
c$$$C     r   The density consists of a smooth part (smrho) plus
c$$$C     r   nbas atom-centered densities inside the MT spheres.
c$$$C     r   Their sum is the full charge density.
c$$$C     r   The local density is represented as the difference of the
c$$$C     r   two valence components in orhoat, plus the core density.
c$$$C     r   Density in the MT spheres:
c$$$C     r      mesh parameters rmt,nr,a;
c$$$C     r      total density rho (times r**2) to lmxl;
c$$$C     r      a spherical potential v0 defining the wave functions within rmt
c$$$C     r      pnu and idmod to lmxa
c$$$C     r   Smooth density
c$$$C     r      real part of complex*16 array smrho contains the density
c$$$C     r      k1,k2,k3 are the physical dimensions of the array
c$$$C     r      n1,n2,n3 are the dimensions of the mesh.
c$$$C     l Local variables
c$$$C     u Updates
c$$$C     u   19 Dec 08 New exch option
c$$$C     u   01 Jul 08 First created
c$$$C     ----------------------------------------------------------------------
c$$$      implicit none
c$$$C     ... Passed parameters
c$$$      character sopts*(*)
c$$$      integer mode,nbas,nat,nspec,ifi,n0
c$$$      parameter (n0=10)
c$$$      type(s_ctrl)::sctrl
c$$$      type(s_site)::ssite(*),ssite1(nbas), ssite2(nbas)
c$$$c     ki      type(s_spec)::sspec(*),sspec1(nspec),sspec2(nspec)
c$$$      type(s_spec)::sspec(nspec),sspec1(nspec),sspec2(nspec)
c$$$      type(s_lat)::slat,slat2
c$$$      type(s_pot)::spot,spot2
c$$$      type(s_bz)::sbz
c$$$
c$$$C     ... Local parameters
c$$$      integer nglob,fopna,fopng,iors,igets,a2vec
c$$$      integer i,j,k,l,im,m,n,j1,j2,js1,js2,nit(2),nspeck(2),isync,nw,
c$$$     .     lexch
c$$$      logical lnsave,lbin,lsopts,have2s,havers,lfail
c$$$      integer nsp,ngabc(3),n1,n2,n3,k1,k2,k3,lrssav,nvl,nchan
c$$$      integer ngabc2(3)
c$$$      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
c$$$c     integer ospec1,ospec2,osite1,osite2,oslat2,ospot2,
c$$$c     .  oorh2,osmrh2,osmpt2,osmrho,osmpot,oorhat,oqmom
c$$$c     kino Dec.26.2011:             type(s_rv1),pointer :: sv_p_oorhat(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        complex(8),pointer :: zv_p_osmrh2(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        complex(8),pointer :: zv_p_osmrho(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        complex(8),pointer :: zv_p_osmpt2(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        complex(8),pointer :: zv_p_osmpot(:) =>NULL()
c$$$
c$$$      integer ,allocatable :: orh2_iv(:)
c$$$      real(8) ,allocatable :: qmom_rv(:)
c$$$
c$$$      integer nlst(2),ix(10),nbask(2),natk(2),isw
c$$$      integer ,allocatable:: ilst(:,:)
c$$$      double precision qbg,dgets,fac(2)
c$$$      character dc*1, fn*120, fileid*68, outs*150, strn*120
c$$$
c$$$      integer :: i_copy_size
c$$$C     ... Heap
c$$$      integer ::iwdummy
c$$$
c$$$C     ... data statements
c$$$C     data vec0 /0d0,0d0,0d0/
c$$$
c$$$C     Defaults
c$$$      fn = 'rst1'
c$$$      havers = .false.
c$$$      have2s = .false.
c$$$C     # spec and # sites in spec 2nd atom density
c$$$      nbask(1) = nbas
c$$$      natk(1) = nat
c$$$      nspeck(1) = nspec
c$$$      nbask(2) = nbas
c$$$      nspeck(2) = nspec
c$$$      natk(2) = nat
c$$$C     true if rs data available, but not saved
c$$$      lnsave = .false.
c$$$C     hangenglob      nsp  = nglob('nsp')
c$$$      nsp  = globalvariables%nsp
c$$$      qbg = sctrl%zbak(1)
c$$$
c$$$
c$$$      nvl=spot%nlml
c$$$      nchan=spot%nlma
c$$$
c$$$
c$$$C     Allocate local species structures; copy original to them
c$$$c     call defrr(ospec1,nspec*nint(sspec(1)))
c$$$c     call defrr(ospec2,nspec*nint(sspec(1)))
c$$$c     call dcopy(nspec*nint(sspec(1)),sspec,1,spec1),1)
c$$$c     call dcopy(nspec*nint(sspec(1)),sspec,1,spec2),1)
c$$$      sspec1(1:nspec)=sspec(1:nspec)
c$$$      sspec2(1:nspec)=sspec(1:nspec)
c$$$C     Allocate local site structures; copy original to them
c$$$c     call defrr(osite1,nbas*nint(ssite(1)))
c$$$c     call defrr(osite2,nbas*nint(ssite(1)))
c$$$c     call dcopy(nbas*nint(ssite(1)),ssite,1,w(osite1),1)
c$$$c     call dcopy(nbas*nint(ssite(1)),ssite,1,w(osite2),1)
c$$$      ssite1(1:nbas)=ssite(1:nbas)
c$$$      ssite2(1:nbas)=ssite(1:nbas)
c$$$C     Allocate second lat structure; copy original to it
c$$$c     call defrr(oslat2,nint(slat(1)))
c$$$c     call dcopy(nint(slat(1)),slat,1,w(oslat2),1)
c$$$      slat2=slat
c$$$C     Allocate second pot structure; copy original to it
c$$$c     call defrr(ospot2,nint(spot(1)))
c$$$c     call dcopy(nint(spot(1)),spot,1,w(ospot2),1)
c$$$      spot2=spot
c$$$
c$$$C     Get pointers for density
c$$$
c$$$ckino Dec.28.2011:        zv_p_osmrho => spot%zv_p_osmrho
c$$$
c$$$ckino Dec.28.2011:        zv_p_osmpot => spot%zv_p_osmpot
c$$$
c$$$c     kino Dec.26.2011:         sv_p_oorhat => spot%sv_p_oorhat
c$$$
c$$$
c$$$C     call upack('pot osmrho oorhat',spot,oorh2,osmrh2,0,0,0)
c$$$C     Allocate space for second density
c$$$C     call upack('lat nabc',slat,ngabc,0,0,0,0)
c$$$C     call fftz30(n1,n2,n3,k1,k2,k3)
c$$$C     call defcc(osmrh2, -k1*k2*k3*nsp)
c$$$C     call pack1('pot osmrho',w(ospot2),osmrh2)
c$$$
c$$$      allocate(qmom_rv(nvl))
c$$$
c$$$
c$$$      dc = sopts(1:1)
c$$$      if (dc .ne. ' ') then
c$$$        print 301
c$$$ 301    format(//' Entering the restart file editor. ',
c$$$     .        'Parsing command-line options ...')
c$$$        lsopts = .true.
c$$$        js2 = 0
c$$$      else
c$$$        print 302
c$$$ 302    format(//' Welcome to the restart file editor.  ',
c$$$     .        'Enter ''?'' to see options.')
c$$$        lsopts = .false.
c$$$      endif
c$$$      if (mode .ne. 1) call rx('rsedit not ready for mode ne 1')
c$$$
c$$$C     ... Return here to resume parsing for arguments
c$$$ 10   continue
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
c$$$C     306 format(' Failed to parse string ',a,' ... try again.')
c$$$ 100  continue
c$$$c     #if PRTNOCR
c$$$c     print '(/'' Option : '',$)'
c$$$c     #else
c$$$      print '(/'' Option : '')'
c$$$c     #endif
c$$$      outs = ' '
c$$$      if (lsopts) then
c$$$        print '(a)', trim(sopts(js1:js2))
c$$$        outs = sopts(js1:js2)
c$$$      else
c$$$        read(*,'(a150)') outs
c$$$      endif
c$$$      call locase(outs)
c$$$
c$$$C     ... Parse and execute the next command
c$$$c     if (.false.) then
c$$$c
c$$$c     else
c$$$
c$$$      if(outs .eq. ' ') then
c$$$        print 304
c$$$ 304    format(' Enter ''q'' to exit, ''a'' to abort',
c$$$     .        ' ''?'' to see menu')
c$$$        goto 10
c$$$
c$$$C     ... Read 2nd density from restart file
c$$$      elseif (outs(1:6).eq.'readb ' .or. outs(1:6).eq.'reada ' .or.
c$$$     .        outs(1:5).eq.'read ' .or. outs(1:7).eq.'readfa ') then
c$$$        if (outs(1:7) .eq. 'readfa ') then
c$$$C     Make space for smooth density
c$$$
c$$$          i_copy_size=size(slat2%nabc)
c$$$          call i8icopy(i_copy_size,slat2%nabc,1,ngabc2,1)
c$$$
c$$$          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
c$$$ckino Dec.28.2011:              allocate(zv_p_osmrh2(abs(-k1*k2*k3*nsp)))
c$$$ckino Dec.28.2011:              if (associated(spot2%zv_p_osmrho))  deallocate(spot2%zv_p_osmrho)
c$$$          if (allocated(spot2%zv_a_osmrho))  deallocate(spot2%zv_a_osmrho)
c$$$ckino Dec.28.2011:              allocate(spot2%zv_p_osmrho(abs(-k1*k2*k3*nsp)))
c$$$          if (allocated(spot2%zv_a_osmrho)) deallocate(spot2%zv_a_osmrho)
c$$$          allocate(spot2%zv_a_osmrho(abs(-k1*k2*k3*nsp)))
c$$$ckino Dec.28.2011:              if (-k1*k2*k3*nsp<0) zv_p_osmrh2(:)=0.0d0
c$$$ckino Dec.28.2011:              if (-k1*k2*k3*nsp<0) spot2%zv_p_osmrho(:)=0.0d0
c$$$          if (-k1*k2*k3*nsp<0) spot2%zv_a_osmrho(:)=0.0d0
c$$$
c$$$
c$$$ckino Dec.28.2011:              spot2%zv_p_osmrho => zv_p_osmrh2
c$$$
c$$$
c$$$          allocate(orh2_iv(3*nbask(2)))
c$$$
c$$$          call rdovfa ( nbask ( 2 ) , nspeck ( 2 ) , ssite2 , sspec2 , 
c$$$     .           slat2 , spot2 , qbg , orh2_iv )
c$$$
c$$$        else
c$$$          call word(outs,2,j1,j2)
c$$$          if (outs(1:6).eq.'reada ') then
c$$$            if (j2 .lt. j1) then
c$$$              fn = 'rsta'
c$$$              ifi = fopna('rsta',-1,1)
c$$$            else
c$$$              fn = outs(j1:j2)
c$$$              ifi = fopng(outs(j1:j2),-1,1)
c$$$            endif
c$$$            lbin = .false.
c$$$          else
c$$$            if (j2 .lt. j1) then
c$$$              fn = 'rst'
c$$$              ifi = fopna('rst',-1,4+1)
c$$$            else
c$$$              fn = outs(j1:j2)
c$$$              ifi = fopng(outs(j1:j2),-1,4+1)
c$$$            endif
c$$$            lbin = .true.
c$$$          endif
c$$$          lrssav = int(sctrl%lrs)
c$$$
c$$$
c$$$          sctrl%lrs=32*0+64+128*0
c$$$
c$$$C     Get the number of atoms and species
c$$$          k = iors ( - 1 , sctrl , iwdummy , iwdummy , slat2 , spot2 , 
c$$$     .           sbz , fileid , nbask ( 2 ) , natk ( 2 ) , nspeck ( 2 ) , orh2_iv 
c$$$     .           , iwdummy , nit ( 2 ) , lbin , ifi )
c$$$
c$$$
c$$$          call info5(10,0,0,' rsedit: '//
c$$$     .           'read 2nd density from file "'//trim(fn)//'" ...  '//
c$$$     .           'nbas = %i  nat = %i  nspec = %i',nbask(2),natk(1),nspeck(2),
c$$$     .           0,0)
c$$$C     Allocate species and site structures, fill with available info
c$$$c     call defrr(ospec2,nspeck(2)*nint(sspec(1)))
c$$$c     call defrr(osite2,nbask(2)*nint(ssite(1)))
c$$$c     call uspec(w(ospec2),-1,0,1,0,0,0)
c$$$c     call usite(w(osite2),-1,0,1,0,0,0)
c$$$          do i=1,nspec
c$$$            call uspec_init(sspec2(i))
c$$$          enddo
c$$$          do i=1,nbas
c$$$            call usite_init(ssite2(i))
c$$$          enddo
c$$$C     call shstru('spec',w(ospec2),1,nspeck(2))
c$$$C     call shstru('site',w(osite2),1,nbask(2))
c$$$          allocate(orh2_iv(3*nbask(2)))
c$$$
c$$$          call pshpr(1)
c$$$          k = iors ( - 2 , sctrl , ssite2 , sspec2 , slat2 , spot2 , sbz 
c$$$     .           , fileid , nbask ( 2 ) , natk ( 2 ) , nspeck ( 2 ) , orh2_iv 
c$$$     .           , iwdummy , nit ( 2 ) , lbin , ifi )
c$$$
c$$$
c$$$C     Check for unassigned species => reduce basis
c$$$          k = 0
c$$$          do  i = 1, nbask(2)
c$$$
c$$$            j=ssite2(i)%spec
c$$$
c$$$            if (j .eq. 0) exit
c$$$            k = i
c$$$          enddo
c$$$          if (k .ne. nbask(2)) then
c$$$            call info2(10,0,0,'%9f(warning) only %i sites had '//
c$$$     .              'associated species ... reducing nbas(2)',k,0)
c$$$          endif
c$$$          nbask(2) = k
c$$$          call poppr
c$$$C     Make space for smooth density
c$$$
c$$$          i_copy_size=size(slat2%nabc)
c$$$          call i8icopy(i_copy_size,slat2%nabc,1,ngabc2,1)
c$$$
c$$$          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
c$$$ckino Dec.28.2011:              allocate(zv_p_osmrh2(abs(-k1*k2*k3*nsp)))
c$$$ckino Dec.28.2011:              if (associated(spot2%zv_p_osmrho))  deallocate(spot2%zv_p_osmrho)
c$$$          if (allocated(spot2%zv_a_osmrho))  deallocate(spot2%zv_a_osmrho)
c$$$ckino Dec.28.2011:              allocate(spot2%zv_p_osmrho(abs(-k1*k2*k3*nsp)))
c$$$          if (allocated(spot2%zv_a_osmrho)) deallocate(spot2%zv_a_osmrho)
c$$$          allocate(spot2%zv_a_osmrho(abs(-k1*k2*k3*nsp)))
c$$$ckino Dec.28.2011:              if (-k1*k2*k3*nsp<0) zv_p_osmrh2(:)=0.0d0
c$$$ckino Dec.28.2011:              if (-k1*k2*k3*nsp<0) spot2%zv_p_osmrho(:)=0.0d0
c$$$          if (-k1*k2*k3*nsp<0) spot2%zv_a_osmrho(:)=0.0d0
c$$$
c$$$ckino Dec.28.2011:              allocate(zv_p_osmpt2(abs(-k1*k2*k3*nsp)))
c$$$ckino Dec.28.2011:              allocate(spot2%zv_p_osmpot(abs(-k1*k2*k3*nsp)))
c$$$          if (allocated(spot2%zv_a_osmpot)) deallocate(spot2%zv_a_osmpot)
c$$$          allocate(spot2%zv_a_osmpot(abs(-k1*k2*k3*nsp)))
c$$$ckino Dec.28.2011:              if (-k1*k2*k3*nsp<0) zv_p_osmpt2(:)=0.0d0
c$$$ckino Dec.28.2011:              if (-k1*k2*k3*nsp<0) spot2%zv_p_osmpot(:)=0.0d0
c$$$          if (-k1*k2*k3*nsp<0) spot2%zv_a_osmpot(:)=0.0d0
c$$$
c$$$
c$$$ckino Dec.28.2011:              spot2%zv_p_osmrho => zv_p_osmrh2
c$$$
c$$$ckino Dec.28.2011:              spot2%zv_p_osmpot => zv_p_osmpt2
c$$$
c$$$
c$$$C     Pack volume
c$$$C     call shstru('spec',w(ospec2),1,nspeck(2))
c$$$C     call shstru('site',w(osite2),1,nbask(2))
c$$$          k = iors ( 1 , sctrl , ssite2 , sspec2 , slat2 , spot2 , sbz 
c$$$     .           , fileid , nbask ( 2 ) , natk ( 2 ) , nspeck ( 2 ) , orh2_iv 
c$$$     .           , iwdummy , nit ( 2 ) , lbin , ifi )
c$$$
c$$$
c$$$          call fclose(ifi)
c$$$        endif
c$$$
c$$$        have2s = .true.
c$$$
c$$$        sctrl%lrs=lrssav
c$$$
c$$$
c$$$        call info0(1,1,0,' Second Density:')
c$$$        call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
c$$$        call prsed2 ( nbask ( 2 ) , nsp , ssite2 , sspec2 , orh2_iv ,
c$$$ckino Dec.28.2011:       .        slat2 , zv_p_osmrh2 , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .        slat2 , spot2%zv_p_osmrho , k1 , k2 , k3 ) 
c$$$     .        slat2 , spot2%zv_a_osmrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$
c$$$        goto 10
c$$$
c$$$C     ... Read restart file
c$$$      elseif (outs(1:4).eq.'rsb ' .or. outs(1:4).eq.'rsa ' .or.
c$$$     .        outs(1:3).eq.'rs' .or. outs(1:5).eq.'rsfa ') then
c$$$        if (outs(1:5) .eq. 'rsfa ') then
c$$$          call rdovfa ( nbask ( 1 ) , nspeck ( 1 ) , ssite1 , sspec1 ,
c$$$c     kino Dec.26.2011:       .     slat , spot , qbg , sv_p_oorhat )
c$$$ckino Dec.28.2011:       .           slat , spot , qbg , spot%sv_p_oorhat ) 
c$$$     .           slat , spot , qbg , spot%sv_a_oorhat )
c$$$
c$$$          fileid = 'rsedit overlap fa densities'
c$$$        else
c$$$          call word(outs,2,j1,j2)
c$$$          if (outs(1:4) .eq. 'rsa ') then
c$$$            if (j2 .lt. j1) then
c$$$              fn = 'rsta'
c$$$            else
c$$$              fn = outs(j1:j2)
c$$$            endif
c$$$            ifi = fopna(fn,-1,1)
c$$$            lbin = .false.
c$$$          else
c$$$            if (j2 .lt. j1) then
c$$$              fn = 'rst'
c$$$            else
c$$$              fn = outs(j1:j2)
c$$$            endif
c$$$            ifi = fopna(fn,-1,4+1)
c$$$            lbin = .true.
c$$$          endif
c$$$          call info0(10,0,-1,
c$$$     .           ' rsedit: read restart file "'//trim(fn)//'"')
c$$$          k = iors ( 1 , sctrl , ssite1 , sspec1 , slat , spot , sbz ,
c$$$c     kino Dec.26.2011:       .      fileid , nbask ( 1 ) , natk ( 1 ) , nspeck ( 1 ) , sv_p_oorhat
c$$$ckino Dec.28.2011:       .           fileid , nbask ( 1 ) , natk ( 1 ) , nspeck ( 1 ) , spot%sv_p_oorhat 
c$$$     .           fileid , nbask ( 1 ) , natk ( 1 ) , nspeck ( 1 ) , spot%sv_a_oorhat 
c$$$     .           , iwdummy , nit , lbin , ifi )
c$$$
c$$$
c$$$C     call shstru('spec',w(ospec2),1,nspeck(2))
c$$$C     call shstru('site',w(osite2),1,nbask(2))
c$$$          call fclose(ifi)
c$$$        endif
c$$$        havers = .true.
c$$$        lnsave = .true.
c$$$        call info2(1,1,0,' Density from '//
c$$$     .        '%?#n==0#rs file:#overlapping free atoms:#',
c$$$     .        isw(outs(1:5).eq.'rsfa '),0)
c$$$
c$$$        i_copy_size=size(slat%nabc)
c$$$        call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
c$$$
c$$$        call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
c$$$c     kino Dec.26.2011:            call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , sv_p_oorhat
c$$$ckino Dec.28.2011:           call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_p_oorhat
c$$$        call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_a_oorhat
c$$$ckino Dec.28.2011:       .        , slat , zv_p_osmrho , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .        , slat , spot%zv_p_osmrho , k1 , k2 , k3 ) 
c$$$     .        , slat , spot%zv_a_osmrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$
c$$$C     ... write atomic densities to file
c$$$      elseif (outs(1:7) .eq. 'wrhoat ') then
c$$$
c$$$        if (havers) then
c$$$c     kino Dec.26.2011:             call prsed3 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , sv_p_oorhat
c$$$ckino Dec.28.2011:              call prsed3 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_p_oorhat
c$$$          call prsed3 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_a_oorhat 
c$$$     .           )
c$$$
c$$$        else
c$$$          call info0(0,0,0,' no rst file read ... nothing writtten')
c$$$        endif
c$$$
c$$$C     ... show
c$$$      elseif (outs(1:5) .eq. 'show ') then
c$$$        if (havers) then
c$$$          call info2(1,1,0,' Density from rs file:  '//
c$$$     .           '%i sites, %i species',nbask(1),nspeck(1))
c$$$C     call upack1('lat nabc',slat,ngabc)
c$$$          call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
c$$$c     kino Dec.26.2011:              call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , sv_p_oorhat
c$$$ckino Dec.28.2011:              call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_p_oorhat
c$$$          call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_a_oorhat
c$$$ckino Dec.28.2011:       .           , slat , zv_p_osmrho , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .           , slat , spot%zv_p_osmrho , k1 , k2 , k3 ) 
c$$$     .           , slat , spot%zv_a_osmrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$        endif
c$$$
c$$$        if (have2s) then
c$$$          call info2(1,1,0,' Second Density:  '//
c$$$     .           '%i sites, %i species',nbask(2),nspeck(2))
c$$$
c$$$          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
c$$$          call prsed2 ( nbask ( 2 ) , nsp , ssite2 , sspec2 , orh2_iv ,
c$$$ckino Dec.28.2011:       .           slat2 , zv_p_osmrh2 , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .           slat2 , spot2%zv_p_osmrho , k1 , k2 , k3 ) 
c$$$     .           slat2 , spot2%zv_a_osmrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$        endif
c$$$
c$$$C     ... Spin averaging of densities, or  Exchange l=0 parts of
c$$$C     set object 1|2 list1 list2
c$$$C     exch 1|2 site1 site2  or exch 1|2 flip site1 site2
c$$$      elseif (outs(1:4) .eq. 'set ' .or. outs(1:5) .eq. 'exch ') then
c$$$
c$$$        lexch = 0
c$$$        if (outs(1:5) .eq. 'exch ') lexch = 1
c$$$        call words(outs,nw)
c$$$        if (nw .lt. 4) goto 98
c$$$
c$$$C     Switches for exchange
c$$$        if (lexch .ne. 0) then
c$$$
c$$$C     Determine whether first or 2nd density
c$$$          call word(outs,2,j1,j2)
c$$$          if (outs(j1:j2) .eq. '1') then
c$$$            lexch = 1
c$$$          elseif (outs(j1:j2) .eq. '2') then
c$$$            lexch = 2
c$$$          else
c$$$            call info0(0,0,0,' 2rd argument must be 1 or 2 in '//
c$$$     .              trim(outs))
c$$$            goto 98
c$$$          endif
c$$$
c$$$C     Check for flip
c$$$          call word(outs,3,j1,j2)
c$$$          if (outs(j1:j2) .eq. 'flip') then
c$$$            if (nw .ne. 5) goto 98
c$$$! Must be 5 arguments in this case
c$$$            isync = 4
c$$$            nw = 3
c$$$          else
c$$$            isync = 2
c$$$            nw = 2
c$$$          endif
c$$$          im = 0
c$$$
c$$$C     Get 1st and 2nd site
c$$$          call numsyv(n)
c$$$          do  k = 1, 2
c$$$            call word(outs,nw+k,j1,j2)
c$$$            call lodsyv('nbas',0,dble(nbask(lexch)),i)
c$$$            i = 0
c$$$            i = a2vec(outs(j1:j2),len(outs(j1:j2)),i,2,' ',1,1,1,ix,
c$$$     .              nlst(k))
c$$$            if (i .le. 0 .or. nlst(k) .gt. nbask(lexch)) then
c$$$              call info0(0,0,0,
c$$$     .                 '%6p... invalid site index : '//outs(j1:j2))
c$$$              goto 98
c$$$            endif
c$$$          enddo
c$$$          call clrsyv(n)
c$$$
c$$$C     Carry out operation on 1st or second density
c$$$          if (lexch .eq. 1) then
c$$$            ix(1) = nbask(1); ix(2) = nbask(1)
c$$$            if (havers) then
c$$$c     ki                call prsed1 ( im * 10 + isync , 1 , nlst , ix , nsp , ssite1
c$$$              call prsed1 (nbas,nbas, nspec
c$$$     .                 ,  im * 10 + isync , 1 , nlst , ix , nsp , ssite1
c$$$c     kino Dec.26.2011:       .          , sspec1 , sv_p_oorhat , ssite1 , sspec1 , sv_p_oorhat , zv_p_osmrho
c$$$ckino Dec.28.2011:       .                 , sspec1 , spot%sv_p_oorhat , ssite1 , sspec1 , spot%sv_p_oorhat , zv_p_osmrho
c$$$ckino Dec.28.2011:       .                 , sspec1 , spot%sv_p_oorhat , ssite1 , sspec1 , spot%sv_p_oorhat , spot%zv_p_osmrho
c$$$ckino Dec.28.2011:       .                 , sspec1 , spot%sv_p_oorhat , ssite1 , sspec1 , spot%sv_p_oorhat , spot%zv_a_osmrho 
c$$$     .                 , sspec1 , spot%sv_a_oorhat , ssite1 , sspec1 , spot%sv_a_oorhat , spot%zv_a_osmrho
c$$$ckino Dec.28.2011:       .                 , ngabc ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , zv_p_osmrho , ngabc
c$$$ckino Dec.28.2011:       .                 , ngabc ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , spot%zv_p_osmrho , ngabc 
c$$$     .                 , ngabc ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , spot%zv_a_osmrho , ngabc 
c$$$     .                 ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , fac , lfail )
c$$$
c$$$
c$$$              if (lfail) then
c$$$                call info0(0,0,0,' Densities failed to synchronize')
c$$$                goto 98
c$$$              endif
c$$$            else
c$$$              call info0(0,0,0,' no rst data ... cannot do: '//trim(outs))
c$$$              goto 10
c$$$            endif
c$$$          else
c$$$            ix(1) = nbask(2); ix(2) = nbask(2)
c$$$            if (have2s) then
c$$$c     ki               call prsed1 ( im * 10 + isync , 1 , nlst , ix , nsp , ssite2
c$$$              call prsed1 (nbas,nbas,nspec, im * 10 + isync , 1 , nlst , ix , nsp , ssite2
c$$$ckino Dec.28.2011:       .                 , sspec2 , orh2_iv , ssite2 , sspec2 , orh2_iv , zv_p_osmrh2
c$$$ckino Dec.28.2011:       .                 , sspec2 , orh2_iv , ssite2 , sspec2 , orh2_iv , spot2%zv_p_osmrho 
c$$$     .                 , sspec2 , orh2_iv , ssite2 , sspec2 , orh2_iv , spot2%zv_a_osmrho
c$$$ckino Dec.28.2011:       .                 , ngabc2 ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , zv_p_osmrh2 ,
c$$$ckino Dec.28.2011:       .                 , ngabc2 ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , spot2%zv_p_osmrho , 
c$$$     .                 , ngabc2 ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , spot2%zv_a_osmrho , 
c$$$     .                 ngabc2 ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , fac , lfail )
c$$$
c$$$
c$$$            else
c$$$              if (lfail) then
c$$$                call info0(0,0,0,' Densities failed to synchronize')
c$$$                goto 98
c$$$              endif
c$$$              goto 10
c$$$            endif
c$$$          endif
c$$$
c$$$C     Switches for set
c$$$        else
c$$$          call word(outs,2,j1,j2)
c$$$          im = -1
c$$$          if (outs(j1:j2) .eq. 'all') im = 4+2+1
c$$$          if (outs(j1:j2) .eq. 'n') im = 4
c$$$          if (outs(j1:j2) .eq. 'v') im = 2
c$$$          if (outs(j1:j2) .eq. 'pnu') im = 1
c$$$          if (outs(j1:j2) .eq. 'ntrue') im = 0
c$$$          if (im .eq. -1) then
c$$$            call info0(0,0,0,
c$$$     .              '     ... mode "'//outs(j1:j2)//'" not recognized')
c$$$            goto 98
c$$$          endif
c$$$
c$$$C     Fourth argument same for each density
c$$$          call word(outs,4,j1,j2)
c$$$          if (outs(j1:j2) .eq. 'zers') then
c$$$            k = im*10+0
c$$$          elseif (outs(j1:j2) .eq. 'zerq') then
c$$$            k = im*10+1
c$$$          elseif (outs(j1:j2) .eq. 'flip') then
c$$$            k = im*10+2
c$$$          else
c$$$            call info0(0,0,0,' failed to recognize 4th argument in: '//
c$$$     .              trim(outs))
c$$$            goto 98
c$$$          endif
c$$$
c$$$C     Carry out operation on 1st or second density
c$$$          call word(outs,3,j1,j2)
c$$$          if (outs(j1:j2) .eq. '1') then
c$$$            if (havers) then
c$$$              call prsed4 ( k , nbask ( 1 ) , nspeck ( 1 ) , nsp , ssite1 ,
c$$$c     kino Dec.26.2011:       .          sspec1 , sv_p_oorhat , slat , zv_p_osmrho , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .                 sspec1 , spot%sv_p_oorhat , slat , zv_p_osmrho , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .                 sspec1 , spot%sv_p_oorhat , slat , spot%zv_p_osmrho , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .                 sspec1 , spot%sv_p_oorhat , slat , spot%zv_a_osmrho , k1 , k2 , k3 ) 
c$$$     .                 sspec1 , spot%sv_a_oorhat , slat , spot%zv_a_osmrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$            else
c$$$              call info0(0,0,0,' no rst data ... cannot do: '//trim(outs))
c$$$              goto 10
c$$$            endif
c$$$          elseif (outs(j1:j2) .eq. '2') then
c$$$            if (have2s) then
c$$$              call prsed4 ( k , nbask ( 2 ) , nspeck ( 2 ) , nsp , ssite2 ,
c$$$ckino Dec.28.2011:       .                 sspec2 , orh2_iv , slat2 , zv_p_osmrh2 , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .                 sspec2 , orh2_iv , slat2 , spot2%zv_p_osmrho , k1 , k2 , k3 ) 
c$$$     .                 sspec2 , orh2_iv , slat2 , spot2%zv_a_osmrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$            else
c$$$              call info0(0,0,0,' no 2nd density... cannot do: '//
c$$$     .                 trim(outs))
c$$$              goto 10
c$$$            endif
c$$$          else
c$$$            call info0(0,0,0,' failed to recognize 3rd argument in '//
c$$$     .              trim(outs))
c$$$            goto 98
c$$$          endif
c$$$        endif
c$$$
c$$$        lnsave = .true.
c$$$
c$$$C     ... Linear combination of densities
c$$$      elseif (outs(1:4) .eq. 'add ') then
c$$$        if (.not. (havers .and. have2s)) goto 97
c$$$
c$$$        call words(outs,nw)
c$$$        if (nw .lt. 3) goto 98
c$$$        call word(outs,2,j1,j2)
c$$$        im = -1
c$$$        if (outs(j1:j2) .eq. 'all') im = 4+2+1
c$$$        if (outs(j1:j2) .eq. 'n') im = 4
c$$$        if (outs(j1:j2) .eq. 'v') im = 2
c$$$        if (outs(j1:j2) .eq. 'pnu') im = 1
c$$$        if (outs(j1:j2) .eq. 'ntrue') im = 0
c$$$        if (im .eq. -1) then
c$$$          call info0(0,0,0,
c$$$     .           '     ... mode "'//outs(j1:j2)//'" not recognized')
c$$$          goto 98
c$$$        endif
c$$$        call word(outs,3,j1,j2)
c$$$        if (outs(j1:j2) .eq. 'sync') then
c$$$          isync = 0
c$$$C     Require 3 words or 5 words
c$$$          if (im .ge. 4 .and. nw .ne. 3 .or.
c$$$     .           im .lt. 4 .and. nw .ne. 5) goto 98
c$$$          nw = 3
c$$$        else
c$$$          isync = 1
c$$$C     Require 4 words or 6 words
c$$$          if (im .ge. 4 .and. nw .ne. 4 .or.
c$$$     .           im .lt. 4 .and. nw .ne. 6) goto 98
c$$$          j = 0
c$$$          j = a2vec(outs(j1:),len(outs(j1:)),j,4,', ',2,-3,2,ix,fac)
c$$$          if (j .ne. 2) then
c$$$            call word(outs,4,j,j2)
c$$$            call info0(0,0,0,
c$$$     .              '     ... Could not parse arguments : '//outs(j1:j2))
c$$$            goto 98
c$$$          endif
c$$$          nw = 4
c$$$        endif
c$$$        call numsyv(n)
c$$$        if (im .lt. 4) then
c$$$          do  k = 1, 2
c$$$            call word(outs,nw+k,j1,j2)
c$$$            call lodsyv('nbas',0,dble(nbask(1)),i)
c$$$            if (k .eq. 2) call lodsyv('nbas',0,dble(nbask(2)),i)
c$$$C     call shosyv(0,0,0,6)
c$$$            call mkils0(outs(j1:j2),nlst(k),j)
c$$$            if (nlst(k) .le. 0) then
c$$$              call info0(0,0,0,
c$$$     .                 '%6p... Bad or null list : '//outs(j1:j2))
c$$$              goto 98
c$$$            endif
c$$$            if (k .eq. 1) then
c$$$              allocate(ilst(nlst(k),2))
c$$$              call mkilst(outs(j1:j2),nlst(k),ilst)
c$$$              call imxmn(nlst(k),ilst,1,l,m)
c$$$            elseif (k .eq. 2) then
c$$$              if (nlst(2) .eq. 1) then
c$$$                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
c$$$                call ivset(ilst(1,2),2,nlst(1),ilst(1,2))
c$$$              elseif (nlst(2) .ne. nlst(1)) then
c$$$                call info2(0,0,0,'%5p... 2nd list: '//outs(j1:j2)//
c$$$     .                    ' has %i entries when %i are needed',nlst(2),
c$$$     .                    nlst(1))
c$$$                goto 98
c$$$              else
c$$$                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
c$$$              endif
c$$$C     Check that elements in list are within range
c$$$              call imxmn(nlst(k),ilst(1,k),1,l,m)
c$$$            endif
c$$$            if (l .lt. 1 .or.  k.eq.1 .and. m .gt. nbask(1)
c$$$     .              .or.  k.eq.2 .and. m .gt. nbask(2)) then
c$$$              call info0(0,0,0,
c$$$     .                 '     ... Entry in list out of range: '//outs(j1:j2))
c$$$              goto 98
c$$$            endif
c$$$          enddo
c$$$        elseif (im .ge. 4) then
c$$$          if (nbask(1) .ne. nbask(2)) then
c$$$            call info0(0,0,0,
c$$$     .              '%10fadda requires 1-1 correspondence between sites'//
c$$$     .              ' in 1st and 2nd densities')
c$$$            goto 98
c$$$          endif
c$$$          nlst(1) = nbask(1)
c$$$          nlst(2) = nbask(1)
c$$$          allocate(ilst(nlst(1),2))
c$$$          strn = ' '
c$$$          call awrit1('1:%i',strn,len(strn),0,nbask(1))
c$$$          call mkilst(strn,nlst(1),ilst)
c$$$          call icopy(nbask(1),ilst,1,ilst(1,2),1)
c$$$        endif
c$$$        call clrsyv(n)
c$$$
c$$$C     call upack1('lat nabc',slat,ngabc)
c$$$C     call upack1('lat nabc',w(oslat2),ngabc2)
c$$$c     ki           call prsed1 ( im * 10 + isync , nlst , ilst , nbask , nsp , ssite1
c$$$        call prsed1 (nbas,nbas,nspec, im * 10 + isync , nlst , ilst , nbask , nsp , ssite1
c$$$c     kino Dec.26.2011:       .     , sspec1 , sv_p_oorhat , ssite2 , sspec2 , orh2_iv , zv_p_osmrho
c$$$ckino Dec.28.2011:       .        , sspec1 , spot%sv_p_oorhat , ssite2 , sspec2 , orh2_iv , zv_p_osmrho
c$$$ckino Dec.28.2011:       .        , sspec1 , spot%sv_p_oorhat , ssite2 , sspec2 , orh2_iv , spot%zv_p_osmrho
c$$$ckino Dec.28.2011:       .        , sspec1 , spot%sv_p_oorhat , ssite2 , sspec2 , orh2_iv , spot%zv_a_osmrho 
c$$$     .        , sspec1 , spot%sv_a_oorhat , ssite2 , sspec2 , orh2_iv , spot%zv_a_osmrho
c$$$ckino Dec.28.2011:       .        , ngabc ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , zv_p_osmrh2 , ngabc2
c$$$ckino Dec.28.2011:       .        , ngabc ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , spot2%zv_p_osmrho , ngabc2 
c$$$     .        , ngabc ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , spot2%zv_a_osmrho , ngabc2 
c$$$     .        ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , fac , lfail )
c$$$
c$$$
c$$$
c$$$
c$$$        if (lfail) then
c$$$          call info0(0,0,0,' Densities failed to synchronize')
c$$$        elseif (isync .eq. 0) then
c$$$          call info0(0,0,0,' Densities synchronized')
c$$$        elseif (isync .eq. 1) then
c$$$C     call upack1('lat nabc',slat,ngabc)
c$$$          call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
c$$$          call info0(1,1,0,' New density from rs file:')
c$$$c     kino Dec.26.2011:              call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , sv_p_oorhat
c$$$ckino Dec.28.2011:              call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_p_oorhat
c$$$          call prsed2 ( nbask ( 1 ) , nsp , ssite1 , sspec1 , spot%sv_a_oorhat
c$$$ckino Dec.28.2011:       .           , slat , zv_p_osmrho , k1 , k2 , k3 )
c$$$ckino Dec.28.2011:       .           , slat , spot%zv_p_osmrho , k1 , k2 , k3 ) 
c$$$     .           , slat , spot%zv_a_osmrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$        endif
c$$$
c$$$        deallocate(ilst)
c$$$        lnsave = .true.
c$$$        goto 10
c$$$
c$$$C     ... Save
c$$$      elseif (outs(1:5) .eq. 'save ' .or. outs(1:6) .eq. 'savea ') then
c$$$        lbin = outs(1:5) .eq. 'save '
c$$$        call word(outs,2,j1,j2)
c$$$        if (j2 .ge. j1) fn = outs(j1:j2)
c$$$        if (lbin) then
c$$$          if (j2 .lt. j1) fn = 'rst'
c$$$          ifi = fopna(fn,-1,4)
c$$$        else
c$$$          if (j2 .lt. j1) fn = 'rsta'
c$$$          ifi = fopna(fn,-1,0)
c$$$        endif
c$$$        call info0(10,0,-1,
c$$$     .        ' rsedit: write restart file "'//trim(fn)//'"')
c$$$        k = iors ( 1 , sctrl , ssite1 , sspec1 , slat , spot , sbz ,
c$$$c     kino Dec.26.2011:       .    fileid , nbask ( 1 ) , natk ( 1 ) , nspeck ( 1 ) , sv_p_oorhat
c$$$ckino Dec.28.2011:       .        fileid , nbask ( 1 ) , natk ( 1 ) , nspeck ( 1 ) , spot%sv_p_oorhat 
c$$$     .        fileid , nbask ( 1 ) , natk ( 1 ) , nspeck ( 1 ) , spot%sv_a_oorhat 
c$$$     .        , iwdummy , nit , lbin , - ifi )
c$$$
c$$$
c$$$
c$$$        lnsave = .false.
c$$$
c$$$C     ... abort
c$$$      elseif (outs(1:2) .eq. 'a ') then
c$$$        call rx0('aborting rs editor ... no file written')
c$$$
c$$$C     ... quit
c$$$      elseif (outs(1:2) .eq. 'q '. or. outs(1:5) .eq. 'quit ') then
c$$$        if (lnsave) then
c$$$          print '('' rs file not saved ... really quit?'')'
c$$$          read(*,'(a150)') outs
c$$$          call locase(outs)
c$$$          if (.not. (outs(1:1) .eq. 'y' .or. outs(1:1) .eq. 'q'))
c$$$     .           goto 10
c$$$        endif
c$$$        call rx0('exit rs editor')
c$$$
c$$$C     ... help
c$$$      elseif (outs .eq. '?') then
c$$$        print 310
c$$$        print 311
c$$$        print 312
c$$$        print 313
c$$$        print 314
c$$$ 310    format(
c$$$     .        ' Select one of these options:'/
c$$$     .        t4,'rs [fn]',t15,
c$$$     .        'read density and associated parameters from binary restart ',
c$$$     .        'file.'/t15,'Use "fn" for file name, ',
c$$$     .        'or "rst.ext" if "fn" is not supplied'/
c$$$     .        t4,'rsa [fn]',t15,
c$$$     .        'read density and associated parameters from ascii restart ',
c$$$     .        'file.'/t15,'Use "fn" for file name, ',
c$$$     .        'or "rst.ext" if "fn" is not supplied'/
c$$$     .        t4,'rsb [fn]',t15, 'same as "rs"'/
c$$$     .        t4,'rsfa [fn]',t15,
c$$$     .        'generate density and associated parameters from atm ',
c$$$     .        'file.'/t15,'Use "fn" for file name, ',
c$$$     .        'or "atm.ext" if "fn" is not supplied')
c$$$ 311    format(/
c$$$     .        t4,'read [fn]',t15,
c$$$     .        'read 2nd density from binary restart file.  Use "fn" for ',
c$$$     .        'file name,'/t15,'or "rst.ext" if "fn" is not supplied'/
c$$$     .        t4,'reada [fn]',t15,
c$$$     .        'read 2nd density from ascii restart file.  Use "fn" for ',
c$$$     .        'file name,'/t15,'or "rst.ext" if "fn" is not supplied'/
c$$$     .        t4,'readb [fn]',t15, 'same as "read"'//
c$$$     .        t4,'wrhoat',t15, 'Write atomic densities to rhoMT.{1..nbas}'//
c$$$     .        t4,'show',t15, 'Show summary information about densities'//
c$$$     .        t4,'...The following affect the charge or magnetic moment',
c$$$     .        ' of a density.'/
c$$$     .        t4,'set all|n       1|2  zers|zerq|flip'/
c$$$     .        t4,'set v|ntrue|pnu 1|2  zers|zerq|flip'/
c$$$     .        t14,'* The first argument specifies which objects are affected.'
c$$$     .        /
c$$$     .        t16,'"all" applies to all quantities, "n" to the ',
c$$$     .        'three density types,'/
c$$$     .        t16,'"v" to the spherical potential, ',
c$$$     .        '"ntrue" to the true local density,'/
c$$$     .        t16,'"pnu" to the linearization parameters'/
c$$$     .        t14,'* The second argument refers to the first or second ',
c$$$     .        'density'/
c$$$     .        t14,'* Third argument = zers: zeros spin part of density:   ',
c$$$     .        'n+ - n- = 0'/
c$$$     .        t14,'  Third argument = zerq: zeros charge part of density: ',
c$$$     .        'n+ + n- = 0'/
c$$$     .        t14,'  Third argument = flip: exchange n+, n- ')
c$$$
c$$$ 312    format(/
c$$$     .        t4,'...The following put a linear combination of the ',
c$$$     .        'two densities into the first, or'/
c$$$     .        t4,'checks whether the two densities are well enough ',
c$$$     .        'synchronized to make it possible'/
c$$$     .        t4,'add all|n',t22,'sync'/
c$$$     .        t4,'add v|ntrue|pnu',t22,'sync lst1 lst2'/
c$$$     .        t4,'add all|n',t22,'fac1 fac2'/
c$$$     .        t4,'add v|ntrue|pnu',t22,'fac1 fac2 lst1 lst2'/
c$$$     .        t14,'* The first argument''s meaning is the same as the ',
c$$$     .        'first argument to "set"'/
c$$$     .        t14,'* Modifier "sync" checks whether the densities are ',
c$$$     .        'compatible, making no changes.'/
c$$$     .        t14,'* If "sync" is missing, the objects specified will be'/
c$$$     .        t16,'replaced by fac1*[1st type] + fac2*[second type].'/
c$$$     .        t14,'* For the site-only cases "v" and "ntrue," ',
c$$$     .        'site list "lst1" is required.'/
c$$$     .        t16,'"lst2" is a matching site list for the 2nd density, OR'/
c$$$     .        t16,'"lst2" is single site (same site added to all sites)'/)
c$$$
c$$$ 313    format(
c$$$     .        t4,'...The following exchanges the l=0 parts of two site',
c$$$     .        ' densities '/
c$$$     .        t4,'exch 1|2',t22,'[flip] site1 site2'/
c$$$     .        )
c$$$
c$$$ 314    format(
c$$$     .        t4,'save  [fn]',t15,'saves restart data in ',
c$$$     .        'binary restart file (name="rst" unless fn supplied)'/
c$$$     .        t4,'savea [fn]',t15,'saves restart data in ',
c$$$     .        'ascii restart file (name="rsta" unless fn supplied)'/
c$$$     .        t4,'q',t15,'to quit the editor'/
c$$$     .        t4,'a',t15,'to abort')
c$$$
c$$$      else
c$$$        print '(1x,''unrecognized option: '',a)', trim(outs)
c$$$
c$$$      endif
c$$$      goto 10
c$$$
c$$$ 97   call info0(0,0,0,'%10pBoth densities must be supplied before'
c$$$     .     //'invoking this command ')
c$$$ 98   call info0(0,0,0,' rsedit:  improper usage of '//trim(outs)//
c$$$     .     ' ... nothing done')
c$$$      goto 10
c$$$
c$$$      if (allocated(orh2_iv)) deallocate(orh2_iv)
c$$$      if (allocated(qmom_rv)) deallocate(qmom_rv)
c$$$
c$$$      end subroutine rsedit
c$$$
c$$$
c$$$c     ki      subroutine prsed1(mode,nlst,ilst,nbask,nsp,site1,spec1,orh1,
c$$$      subroutine prsed1 ( nsite1 , nsite2 , nspec , mode , nlst , ilst 
c$$$     .     , nbask , nsp , site1 , spec1 , sv_p_orh1 , site2 , spec2 , sv_p_orh2 
c$$$     .     , smrho , k11 , k21 , k31 , smrh2 , k12 , k22 , k32 , fac , lfail 
c$$$     .     )
c$$$
c$$$
c$$$
c$$$      use m_struc_def           !Cgetarg
c$$$      use m_globalvariables
c$$$C     - Add fac*2nd density to density
c$$$C     ----------------------------------------------------------------------
c$$$C     i Inputs
c$$$C     i   mode  :0 check synchronization
c$$$C     i         :1 add the density
c$$$C     i         :2 swap l=0 part of 2 local densities.  nlst should be 1
c$$$C     i         :4 like 2, but exchange spin in the swap
c$$$C     i         :10s digit:
c$$$C     i         :0 work with true density only
c$$$C     i         :1 work with pnu
c$$$C     i         :2 work with v0,v1
c$$$C     i         :4 work with n1,n2,n0
c$$$C     i         :Any combination is allowed
c$$$C     i   nlst  :number of site densities to which 2nd density is added
c$$$C     i   ilst  :list of indices to site densities
c$$$C     i         :ilst(:,1) is list of sites to be modified
c$$$C     i         :ilst(:,2) is list from which 2nd density is to be found
c$$$C     i   nsp   :2 for spin-polarized case, otherwise 1
c$$$C     i
c$$$C     i   site1 :struct for site-specific information; see routine usite
c$$$C     i   spec1 :struct for species-specific information; see routine uspec
c$$$C     i   orh1  :vector of offsets to arrays containing local site densities
c$$$C     i   site2 :Analog of site1 for 2nd density
c$$$C     i   spec2 :Analog of spec1 for 2nd density
c$$$C     i   orh2  :Analog of orh1 for 2nd density
c$$$C     i   fac   :Add fac * 2nd density into density
c$$$C     o Outputs
c$$$C     l Local variables
c$$$C     l         :
c$$$C     r Remarks
c$$$C     r
c$$$C     u Updates
c$$$C     u   28 Jun 08
c$$$C     ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C     ... Passed parameters
c$$$      integer mode,nlst,ilst(nlst,2),nbask(2),nsp,
c$$$     .     k11,k21,k31,k12,k22,k32
c$$$      type(s_rv1) :: sv_p_orh2(3,*)
c$$$
c$$$      type(s_rv1) :: sv_p_orh1(3,*)
c$$$
c$$$      real(8):: fac(2)
c$$$c     ki
c$$$      integer:: nsite1,nsite2,nspec
c$$$      type(s_site)::site1(nsite1)
c$$$c     ki      type(s_spec)::spec1(*)
c$$$      type(s_spec)::spec1(nspec)
c$$$      type(s_site)::site2(nsite2)
c$$$c     ki      type(s_spec)::spec2(*)
c$$$      type(s_spec)::spec2(nspec)
c$$$
c$$$      double complex smrho(k11,k21,k31,nsp),smrh2(k12,k22,k32,nsp)
c$$$      logical lfail
c$$$C     ... Local parameters
c$$$      integer il,ib,jb,is,js,lmxa,lmxl,lmxa2,lmxl2,nr1,nr2,nrmx
c$$$      integer nglob,stdo,ipr,nlml1,nlml2,isync,im,n0
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_orhoc2(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_orhoc1(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_ov01(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_ov11(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_ov02(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_ov12(:) =>NULL()
c$$$
c$$$      integer ,allocatable:: flgspc(:)
c$$$      parameter (nrmx=1501,n0=10)
c$$$      double precision rmt1,rmt2,a1,a2,rofi(nrmx),pnu(n0,2),pnu2(n0,2)
c$$$      character spid1*8,spid2*8
c$$$      logical lfails
c$$$C     ... Heap
c$$$
c$$$C     hangenglob      stdo = nglob('stdo')
c$$$      stdo = globalvariables%stdo
c$$$      call getpr(ipr)
c$$$      isync = mod(mode,10)
c$$$      im  = mod(mode/10,10)
c$$$
c$$$      lfails = .false.
c$$$      if (im .ge. 4) then
c$$$        call info0(0,1,-1,' Checking whether sm density meshes '//
c$$$     .        'synchronize ...')
c$$$        if (k11 .ne. k12 .or. k21 .ne. k22 .or. k31 .ne. k32) then
c$$$          write(stdo,212) k11,k21,k31,k12,k22,k32
c$$$ 212      format(' mesh mismatch:'/
c$$$     .           ' first  density n1 n2 n3 = ',3i5/
c$$$     .           ' second density n1 n2 n3 = ',3i5)
c$$$          lfails = .true.
c$$$          if (isync .eq. 1) call rx('sm density mesh mismatch')
c$$$        else
c$$$          call info0(0,0,0,' ok')
c$$$        endif
c$$$        if (isync .eq. 1) then
c$$$          call info2(0,1,-1,' ... Replace smooth density n0 by linear '
c$$$     .           //'combination %d*[n0] + %d*[2nd n0]',fac(1),fac(2))
c$$$          call lcrho(k11*k21*k31*2,nsp,1,1,fac(1),fac(2),smrho,smrh2)
c$$$        endif
c$$$      endif
c$$$
c$$$      if (isync .ge. 2) then
c$$$        call info5(0,1,0,' ... Swap '//
c$$$     .        '%?#n==4#and spin flip ##'//
c$$$     .        '%?#n==0# true site# site #%-1j'//
c$$$     .        '%?#n==0# density n1##%-1j'//
c$$$     .        '%?#n>=4# densities##'//
c$$$     .        '%?#n>=2# potentials##'//
c$$$     .        '%?#n>=1# pnu##',
c$$$     .        isync,im,mod(im,4),mod(im,2),0)
c$$$        if (isync .eq. 3 .or. isync .gt. 4 .or. nlst .gt. 1) then
c$$$          call info0(0,1,0,' ... illegal parameters (aborting)')
c$$$          return
c$$$        endif
c$$$      elseif (isync .eq. 1) then
c$$$        call info5(0,1,0,' ... Replace '//
c$$$     .        '%?#n==0# true site# site #%-1j'//
c$$$     .        '%?#n==0# density n1##%-1j'//
c$$$     .        '%?#n>=4# densities##'//
c$$$     .        '%?#n>=2# potentials##'//
c$$$     .        '%?#n>=1# pnu##'//
c$$$     .        '  by linear combination'//
c$$$     .        ' %d*[1st] + %d*[2nd]',im,
c$$$     .        mod(im,4),mod(im,2),fac(1),fac(2))
c$$$      elseif (isync .eq. 0) then
c$$$        call info0(0,1,0,' Checking whether radial meshes '//
c$$$     .        'synchronize ...')
c$$$      endif
c$$$      call info2(0,0,0,' site 1 : %n:1,3i',nlst,ilst(1,1))
c$$$      call info2(0,0,0,' site 2 : %n:1,3i',nlst,ilst(1,2))
c$$$
c$$$ 200  format(/' ... Improper synchronization between radial meshes'/
c$$$     .     '  site1',8x,'site2',8x,'nr1  nr2',3x,'rmt1   rmt2',
c$$$     .     4x,' a1    a2')
c$$$
c$$$      lfail = .false.
c$$$      allocate(flgspc(nbask(1)))
c$$$      call iinit(flgspc,nbask(1))
c$$$      do  il = 1, nlst
c$$$
c$$$        ib = ilst(il,1)
c$$$        if (ib .gt. nbask(1))
c$$$     .        call rxi('illegal site index for density:',ib)
c$$$        jb = ilst(il,2)
c$$$        if (jb .gt. nbask(2))
c$$$     .        call rxi('illegal index for 2nd density:',jb)
c$$$
c$$$
c$$$        is=site1(ib)%spec
c$$$        i_copy_size=size(site1(ib)%pnu)
c$$$        call dcopy(i_copy_size,site1(ib)%pnu,1,pnu,1)
c$$$ckino Dec.28.2011:           rv_p_ov01 => site1(ib)%rv_p_ov0
c$$$
c$$$ckino Dec.28.2011:           rv_p_ov11 => site1(ib)%rv_p_ov1
c$$$
c$$$
c$$$
c$$$        js=site2(jb)%spec
c$$$        i_copy_size=size(site2(jb)%pnu)
c$$$        call dcopy(i_copy_size,site2(jb)%pnu,1,pnu2,1)
c$$$ckino Dec.28.2011:           rv_p_ov02 => site2(jb)%rv_p_ov0
c$$$
c$$$ckino Dec.28.2011:           rv_p_ov12 => site2(jb)%rv_p_ov1
c$$$
c$$$
c$$$
c$$$        lmxa=spec1(is)%lmxa
c$$$        lmxl=spec1(is)%lmxl
c$$$
c$$$
c$$$        lmxa2=spec2(js)%lmxa
c$$$        lmxl2=spec2(js)%lmxl
c$$$
c$$$        do i_spacks=is,is
c$$$          call spacks_copy('u',spec1(i_spacks)%name,is,is,spid1,i_spacks)
c$$$        enddo
c$$$
c$$$        do i_spacks=js,js
c$$$          call spacks_copy('u',spec2(i_spacks)%name,js,js,spid2,i_spacks)
c$$$        enddo
c$$$
c$$$
c$$$        if (lmxa .eq. -1) goto 10
c$$$
c$$$
c$$$        a1=spec1(is)%a
c$$$        nr1=spec1(is)%nr
c$$$        rmt1=spec1(is)%rmt
c$$$
c$$$
c$$$        a2=spec2(js)%a
c$$$        nr2=spec2(js)%nr
c$$$        rmt2=spec2(js)%rmt
c$$$
c$$$
c$$$C     Radial meshes must match
c$$$        if (abs(a1-a2).gt.1d-6 .or. nr1.ne.nr2 .or.
c$$$     .        abs(rmt1-rmt2).gt.1d-6) then
c$$$
c$$$          if (.not. lfail) write(stdo,200)
c$$$          lfail = .true.
c$$$          write(stdo,202) ib,spid1,jb,spid2,
c$$$     .           nr1,nr2,rmt1,rmt2,a1,a2
c$$$ 202      format(i4,1x,a,i4,1x,a,2i5,1x,2f7.4,1x,2f6.3)
c$$$          if (isync .eq. 1) call rx('site density mismatch')
c$$$
c$$$        endif
c$$$        if (isync .eq. 0) goto 10
c$$$
c$$$        nlml1 = (lmxl+1)**2
c$$$        nlml2 = (lmxl2+1)**2
c$$$        call radmsh(rmt1,a1,nr1,rofi)
c$$$
c$$$C     True local density: modified for im=0 or im>3
c$$$        if (im .lt. 1 .or. im .ge. 4) then
c$$$          if (isync .lt. 2) then
c$$$            call lcrho ( nr1 , nsp , nlml1 , nlml2 , fac ( 1 ) , fac ( 2 
c$$$     .              ) , sv_p_orh1 ( 1 , ib ) %v , sv_p_orh2( 1 , jb )%v )
c$$$
c$$$
c$$$          else
c$$$            call swrho ( isync / 4 , nr1 , nsp , 1 , nlml1 , nlml2 , sv_p_orh1 
c$$$     .              ( 1 , ib ) %v , sv_p_orh2( 1 , jb )%v )
c$$$
c$$$
c$$$          endif
c$$$        endif
c$$$C     Smooth local density and site core density: modified if im>=4
c$$$        if (im .ge. 4) then
c$$$          if (isync .lt. 2) then
c$$$            call lcrho ( nr1 , nsp , nlml1 , nlml2 , fac ( 1 ) , fac ( 2 
c$$$     .              ) , sv_p_orh1 ( 2 , ib ) %v , sv_p_orh2( 2 , jb )%v )
c$$$
c$$$
c$$$            call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , sv_p_orh1 
c$$$     .              ( 3 , ib ) %v , sv_p_orh2( 3 , jb )%v )
c$$$
c$$$
c$$$          else
c$$$            call rx('not ready for this swap')
c$$$          endif
c$$$        endif
c$$$C     Potentials v0,v1: modified if 2's bit im set
c$$$        if (mod(im,4) .ge. 2) then
c$$$          if (isync .lt. 2) then
c$$$ckino Dec.28.2011:                 call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , rv_p_ov01
c$$$ckino Dec.28.2011:                 call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , site1(ib)%rv_p_ov0
c$$$            call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , site1(ib)%rv_a_ov0
c$$$ckino Dec.28.2011:       .              , rv_p_ov02 )
c$$$ckino Dec.28.2011:       .              , site2(jb)%rv_p_ov0 ) 
c$$$     .              , site2(jb)%rv_a_ov0 )
c$$$
c$$$
c$$$ckino Dec.28.2011:                 call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , rv_p_ov11
c$$$ckino Dec.28.2011:                 call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , site1(ib)%rv_p_ov1
c$$$            call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , site1(ib)%rv_a_ov1
c$$$ckino Dec.28.2011:       .              , rv_p_ov12 )
c$$$ckino Dec.28.2011:       .              , site2(jb)%rv_p_ov1 ) 
c$$$     .              , site2(jb)%rv_a_ov1 )
c$$$
c$$$
c$$$          else
c$$$            call rx('not ready for this swap')
c$$$          endif
c$$$        endif
c$$$C     Pnu modified if 1's bit im set
c$$$        if (mod(im,2) .ge. 1) then
c$$$          if (isync .lt. 2) then
c$$$            call lcrho(n0,nsp,1,1,fac(1),fac(2),pnu,pnu2)
c$$$
c$$$            i_copy_size=size(site1(ib)%pnu)
c$$$            call dcopy(i_copy_size,pnu,1,site1(ib)%pnu,1)
c$$$
c$$$          else
c$$$            call rx('not ready for this swap')
c$$$          endif
c$$$        endif
c$$$
c$$$C     Species core density modified if im ge 4
c$$$        if (im .ge. 4 .and. flgspc(is) .eq. 0) then
c$$$
c$$$ckino Dec.28.2011:              rv_p_orhoc1 => spec1(is)%rv_p_orhoc
c$$$
c$$$
c$$$
c$$$ckino Dec.28.2011:              rv_p_orhoc2 => spec2(js)%rv_p_orhoc
c$$$
c$$$
c$$$          if (isync .lt. 2) then
c$$$ckino Dec.28.2011:                 call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , rv_p_orhoc1
c$$$ckino Dec.28.2011:                 call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , spec1(is)%rv_p_orhoc
c$$$            call lcrho ( nr1 , nsp , 1 , 1 , fac ( 1 ) , fac ( 2 ) , spec1(is)%rv_a_orhoc
c$$$ckino Dec.28.2011:       .              , rv_p_orhoc2 )
c$$$ckino Dec.28.2011:       .              , spec2(js)%rv_p_orhoc ) 
c$$$     .              , spec2(js)%rv_a_orhoc )
c$$$
c$$$
c$$$          else
c$$$            call rx('not ready for this swap')
c$$$          endif
c$$$          flgspc(is) = 1
c$$$        endif
c$$$
c$$$ 10     continue
c$$$      enddo
c$$$
c$$$      if (lfail) then
c$$$        call info0(0,0,0,' ... Radial meshes failed to synchronize')
c$$$      elseif (isync .eq. 0) then
c$$$        call info0(0,0,0,' ... Radial meshes synchronized')
c$$$C     elseif (isync .eq. 1) then
c$$$C     call info0(0,0,0,' ... Linear combination successful')
c$$$      endif
c$$$
c$$$      lfail = lfail .or. lfails
c$$$
c$$$      end subroutine prsed1
c$$$
c$$$
c$$$      subroutine prsed2 ( nbas , nsp , ssite , sspec , sv_p_orhoat 
c$$$     .     , slat , smrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$      use m_struc_def           !Cgetarg
c$$$      use m_globalvariables
c$$$C     - Print out density information
c$$$C     ----------------------------------------------------------------------
c$$$C     i Inputs
c$$$C     i   nspec
c$$$C     i   sspec :struct for species-specific information; see routine uspec
c$$$C     i   nspec2:number of species for second density
c$$$C     i   sspec2
c$$$C     o Outputs
c$$$C     l Local variables
c$$$C     l         :
c$$$C     r Remarks
c$$$C     r
c$$$C     u Updates
c$$$C     u   28 Jun 08
c$$$C     ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C     ... Passed parameters
c$$$      integer nbas,nsp,k1,k2,k3
c$$$      type(s_rv1) :: sv_p_orhoat(3,nbas)
c$$$
c$$$      type(s_site)::ssite(*)
c$$$      type(s_spec)::sspec(*)
c$$$      type(s_lat)::slat
c$$$
c$$$      double complex smrho(k1,k2,k3,nsp)
c$$$C     ... Local parameters
c$$$      integer nrmx,nlmx
c$$$      parameter (nrmx=1501, nlmx=64)
c$$$      integer ib,is,nr,lmxa,lmxl,lmxb,nlml,ipr,nglob,stdo,n1,n2,n3
c$$$      double precision z,qc,rg,a,rmt,rsma
c$$$      double precision pi,srfpi,y0,rofi(nrmx),rwgt(nrmx),qv1,qv2,am1,am2
c$$$      double precision ddot,smq,sum2,alat,plat(3,3),qlat(3,3),vol,smag
c$$$      character spid*8
c$$$C     double complex smpot(k1,k2,k3,2)
c$$$C     ... Heap
c$$$
c$$$C     hangenglob      stdo = nglob('stdo')
c$$$      stdo = globalvariables%stdo
c$$$      call getpr(ipr)
c$$$      pi = 4d0*datan(1d0)
c$$$      srfpi = dsqrt(4d0*pi)
c$$$      y0 = 1d0/srfpi
c$$$      n1 = k1
c$$$      n2 = k2
c$$$      n3 = k3
c$$$
c$$$      alat=slat%alat
c$$$      i_copy_size=size(slat%plat)
c$$$      call dcopy(i_copy_size,slat%plat,1,plat,1)
c$$$
c$$$      call dinv33(plat,1,qlat,vol)
c$$$      vol = dabs(vol)*(alat**3)
c$$$      if (nsp .eq. 2) then
c$$$        call daxpy(k1*k2*k3*2,1d0,smrho(1,1,1,2),1,smrho,1)
c$$$      endif
c$$$      call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smq,sum2)
c$$$      smag = 0
c$$$      if (nsp .eq. 2) then
c$$$        call daxpy(k1*k2*k3*2,-1d0,smrho(1,1,1,2),1,smrho,1)
c$$$        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smag,sum2)
c$$$        smag = 2*smag - smq
c$$$      endif
c$$$
c$$$C     if (nsp .eq. 1) write(stdo,100)
c$$$C     if (nsp .eq. 2) write(stdo,201)
c$$$
c$$$      call info5(0,0,-1,' Smooth density: vol=%,2;2d  n123 = %i %i %i',
c$$$     .     vol,k1,k2,k3,0)
c$$$      call info5(0,0,0,'  qsm = %,6;6d%?#n==2#  mag. mom = %,6;6d##',
c$$$     .     smq,nsp,smag,0,0)
c$$$
c$$$      if (nsp .eq. 1) write(stdo,200)
c$$$      if (nsp .eq. 2) write(stdo,201)
c$$$ 200  format(/' Site parameters and charges:'/
c$$$     .     '   site',10x,'z     rmt     nr   a  lmxl   qtrue',6x,'q2')
c$$$ 201  format(/' Site parameters and charges:'/
c$$$     .     '   site',10x,'z     rmt     nr   a  lmxl   qtrue',
c$$$     .     5x,'mtrue',5x,'q2        m2')
c$$$
c$$$      do  ib = 1, nbas
c$$$
c$$$
c$$$        is=ssite(ib)%spec
c$$$
c$$$C     call upack('site spec pnu ov0 ov1',ssite,ib,is,pnu,ov0,ov01)
c$$$C     call upack2('site pz',ssite,ib,pnz)
c$$$
c$$$        z=sspec(is)%z
c$$$        qc=sspec(is)%qc
c$$$        rg=sspec(is)%rg
c$$$
c$$$        do i_spacks=is,is
c$$$          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c$$$        enddo
c$$$
c$$$
c$$$        a=sspec(is)%a
c$$$        nr=sspec(is)%nr
c$$$        rmt=sspec(is)%rmt
c$$$
c$$$
c$$$        rsma=sspec(is)%rsma
c$$$        lmxa=sspec(is)%lmxa
c$$$        lmxl=sspec(is)%lmxl
c$$$        lmxb=sspec(is)%lmxb
c$$$
c$$$        if (lmxa .eq. -1) goto 10
c$$$        nlml = (lmxl+1)**2
c$$$
c$$$        call radmsh(rmt,a,nr,rofi)
c$$$        call radwgt(rmt,a,nr,rwgt)
c$$$
c$$$        call splrho ( 0 , nsp , nr , nlml , sv_p_orhoat( 1 , ib )%v , 
c$$$     .        sv_p_orhoat( 2 , ib )%v , sv_p_orhoat( 3 , ib )%v )
c$$$
c$$$        qv1 = srfpi * ddot ( nr , rwgt , 1 , sv_p_orhoat( 1 , ib )%v 
c$$$     .        , 1 )
c$$$
c$$$        qv2 = srfpi * ddot ( nr , rwgt , 1 , sv_p_orhoat( 2 , ib )%v 
c$$$     .        , 1 )
c$$$
c$$$        call splrho ( 1 , nsp , nr , nlml , sv_p_orhoat( 1 , ib )%v , 
c$$$     .        sv_p_orhoat( 2 , ib )%v , sv_p_orhoat( 3 , ib )%v )
c$$$
c$$$        am1 = 2 * srfpi * ddot ( nr , rwgt , 1 , sv_p_orhoat( 1 , ib )%v 
c$$$     .        , 1 ) - qv1
c$$$
c$$$        am2 = 2 * srfpi * ddot ( nr , rwgt , 1 , sv_p_orhoat( 2 , ib )%v 
c$$$     .        , 1 ) - qv2
c$$$
c$$$
c$$$        if (nsp .eq. 1) then
c$$$          write(stdo,202) ib,spid,z,rmt,nr,a,lmxl,qv1,am1
c$$$        else
c$$$          write(stdo,202) ib,spid,z,rmt,nr,a,lmxl,qv1,am1,qv2,am2
c$$$        endif
c$$$ 202    format(i4,1x,a,f6.1,f9.5,i5,f6.3,i3,4f10.5)
c$$$
c$$$ 10     continue
c$$$      enddo
c$$$
c$$$
c$$$      end subroutine prsed2
c$$$
c$$$
c$$$      subroutine prsed3 ( nbas , nsp , ssite , sspec , sv_p_orhoat 
c$$$     .     )
c$$$
c$$$
c$$$      use m_struc_def           !Cgetarg
c$$$
c$$$C     - Read/write atomic densities
c$$$C     implicit none
c$$$C     ... Passed parameters
c$$$      integer nbas,nsp
c$$$      type(s_rv1) :: sv_p_orhoat(3,nbas)
c$$$
c$$$      type(s_site)::ssite(*)
c$$$      type(s_spec)::sspec(*)
c$$$
c$$$C     ... Local parameters
c$$$      integer nrmx,nlmx
c$$$      parameter (nrmx=1501, nlmx=64)
c$$$      integer ib,is,nr,lmxa,lmxl,nlml
c$$$      double precision a,rmt
c$$$      double precision rofi(nrmx)
c$$$C     ... Heap
c$$$
c$$$
c$$$      do  ib = 1, nbas
c$$$
c$$$
c$$$        is=ssite(ib)%spec
c$$$
c$$$
c$$$        lmxa=sspec(is)%lmxa
c$$$        lmxl=sspec(is)%lmxl
c$$$
c$$$        if (lmxa .eq. -1) goto 10
c$$$
c$$$        a=sspec(is)%a
c$$$        nr=sspec(is)%nr
c$$$        rmt=sspec(is)%rmt
c$$$
c$$$        nlml = (lmxl+1)**2
c$$$        call radmsh(rmt,a,nr,rofi)
c$$$
c$$$        call wrhomt ( 'rhoMT.' , 'density' , ib , sv_p_orhoat( 1 , ib )%v 
c$$$     .        , rofi , nr , nlml , nsp )
c$$$
c$$$
c$$$ 10     continue
c$$$      enddo
c$$$      end subroutine prsed3
c$$$
c$$$
c$$$      subroutine prsed4 ( mode , nbas , nspec , nsp , ssite , sspec 
c$$$     .     , sv_p_orhoat , slat , smrho , k1 , k2 , k3 )
c$$$
c$$$
c$$$      use m_struc_def           !Cgetarg
c$$$      use m_globalvariables
c$$$C     - Remove spin or charge part of atomic densities
c$$$C     ----------------------------------------------------------------------
c$$$C     i Inputs
c$$$C     i   mode  :0 remove spin part, 1 remove charge part, 2, flip spins
c$$$C     i         :10s digit:
c$$$C     i         :0 work with true density only
c$$$C     i         :1 work with pnu
c$$$C     i         :2 work with v0,v1
c$$$C     i         :4 work with n1,n2,n0
c$$$C     i         :Any combination is allowed
c$$$C     i   nbas  :number of atoms
c$$$C     i   sspec :struct for species-specific information; see routine uspec
c$$$C     i   nspec2:number of species for second density
c$$$C     i   sspec2
c$$$C     o Outputs
c$$$C     l Local variables
c$$$C     l         :
c$$$C     r Remarks
c$$$C     r
c$$$C     u Updates
c$$$C     u   28 Jun 08
c$$$C     ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C     ... Passed parameters
c$$$      integer mode,nbas,nspec,nsp,k1,k2,k3
c$$$      type(s_rv1) :: sv_p_orhoat(3,nbas)
c$$$
c$$$c     ki      type(s_site)::ssite(*)
c$$$      type(s_site)::ssite(nbas)
c$$$c     ki      type(s_spec)::sspec(*)
c$$$      type(s_spec)::sspec(nspec)
c$$$      type(s_lat)::slat
c$$$
c$$$      double complex smrho(k1,k2,k3,nsp)
c$$$C     ... Local parameters
c$$$      integer nrmx,nlmx,n0
c$$$      parameter (nrmx=1501, nlmx=64, n0=10)
c$$$      integer ib,is,nr,lmxa,lmxl,lmxb,nlml,ipr,nglob,stdo,im,mode0
c$$$      double precision z,qc,rg,a,rmt,rsma
c$$$      double precision pi,srfpi,y0,rofi(nrmx),rwgt(nrmx),pnu(n0,2)
c$$$c     kino Dec.26.2011:           real(8),pointer :: rv_p_orhoca(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_ov1(:) =>NULL()
c$$$
c$$$ckino Dec.28.2011:        real(8),pointer :: rv_p_ov0(:) =>NULL()
c$$$
c$$$      character spid*8
c$$$C     ... Heap
c$$$      integer ::iwdummy
c$$$
c$$$
c$$$C     hangenglob      stdo = nglob('stdo')
c$$$      stdo = globalvariables%stdo
c$$$      call getpr(ipr)
c$$$      pi = 4d0*datan(1d0)
c$$$      srfpi = dsqrt(4d0*pi)
c$$$      y0 = 1d0/srfpi
c$$$      mode0 = mod(mode,10)
c$$$      im  = mod(mode/10,10)
c$$$      if (mode0 .eq. 0 .and. nsp .eq. 1) return
c$$$      if (mode0 .eq. 2 .and. nsp .eq. 1) return
c$$$
c$$$      if (im .ge. 4) then
c$$$C     Replace smrho-up,smrho-down with smrho, smrho(up-down)
c$$$        call splrho(30,nsp,k1*k2*k3*2,1,smrho,smrho,smrho)
c$$$        call prsed5(mode0,smrho,k1*k2*k3*2,1,nsp)
c$$$C     Replace smrho, smrho(up-down) with smrho-up,smrho-down
c$$$        call splrho(31,nsp,k1*k2*k3*2,1,smrho,smrho,smrho)
c$$$      endif
c$$$
c$$$C     if (nsp .eq. 1) write(stdo,200)
c$$$C     if (nsp .eq. 2) write(stdo,201)
c$$$C     200 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
c$$$C     .  6x,'q2')
c$$$C     201 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
c$$$C     .  5x,'mtrue',5x,'q2        m2')
c$$$
c$$$      do  ib = 1, nbas
c$$$
c$$$
c$$$        is=ssite(ib)%spec
c$$$
c$$$
c$$$        is=ssite(ib)%spec
c$$$        i_copy_size=size(ssite(ib)%pnu)
c$$$        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
c$$$ckino Dec.28.2011:           rv_p_ov0 => ssite(ib)%rv_p_ov0
c$$$
c$$$ckino Dec.28.2011:           rv_p_ov1 => ssite(ib)%rv_p_ov1
c$$$
c$$$
c$$$C     call upack2('site pz',ssite,ib,pnz)
c$$$
c$$$        z=sspec(is)%z
c$$$        qc=sspec(is)%qc
c$$$        rg=sspec(is)%rg
c$$$
c$$$        do i_spacks=is,is
c$$$          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c$$$        enddo
c$$$
c$$$
c$$$        a=sspec(is)%a
c$$$        nr=sspec(is)%nr
c$$$        rmt=sspec(is)%rmt
c$$$
c$$$
c$$$        rsma=sspec(is)%rsma
c$$$        lmxa=sspec(is)%lmxa
c$$$        lmxl=sspec(is)%lmxl
c$$$        lmxb=sspec(is)%lmxb
c$$$
c$$$        if (lmxa .eq. -1) goto 10
c$$$        nlml = (lmxl+1)**2
c$$$
c$$$        call radmsh(rmt,a,nr,rofi)
c$$$        call radwgt(rmt,a,nr,rwgt)
c$$$
c$$$C     Replace rho-up,rho-down with rho, up-down
c$$$        call splrho ( 0 , nsp , nr , nlml , sv_p_orhoat( 1 , ib )%v , 
c$$$     .        sv_p_orhoat( 2 , ib )%v , sv_p_orhoat( 3 , ib )%v )
c$$$
c$$$
c$$$C     True local density: modified for im=0 or im>3
c$$$        if (im .lt. 1 .or. im .ge. 4) then
c$$$          call prsed5 ( mode0 , sv_p_orhoat( 1 , ib )%v , nr , nlml , nsp 
c$$$     .           )
c$$$
c$$$        endif
c$$$        if (im .ge. 4) then
c$$$          call prsed5 ( mode0 , sv_p_orhoat( 2 , ib )%v , nr , nlml , nsp 
c$$$     .           )
c$$$
c$$$          call prsed5 ( mode0 , sv_p_orhoat( 3 , ib )%v , nr , 1 , nsp 
c$$$     .           )
c$$$
c$$$        endif
c$$$
c$$$C     Restore rho-up,rho-down from rho, up-down
c$$$        call splrho ( 1 , nsp , nr , nlml , sv_p_orhoat( 1 , ib )%v , 
c$$$     .        sv_p_orhoat( 2 , ib )%v , sv_p_orhoat( 3 , ib )%v )
c$$$
c$$$
c$$$C     Potentials v0,v1: modified if 2's bit im set
c$$$        if (mod(im,4) .ge. 2) then
c$$$C     call prmx('v0',w(ov0),nr,nr,2)
c$$$ckino Dec.28.2011:              call splrho ( 30 , nsp , nr , 1 , rv_p_ov0 , iwdummy , iwdummy
c$$$ckino Dec.28.2011:              call splrho ( 30 , nsp , nr , 1 , ssite(ib)%rv_p_ov0 , iwdummy , iwdummy
c$$$          call splrho ( 30 , nsp , nr , 1 , ssite(ib)%rv_a_ov0 , iwdummy , iwdummy 
c$$$     .           )
c$$$
c$$$
c$$$ckino Dec.28.2011:              call prsed5 ( mode0 , rv_p_ov0 , nr , 1 , nsp )
c$$$ckino Dec.28.2011:              call prsed5 ( mode0 , ssite(ib)%rv_p_ov0 , nr , 1 , nsp )
c$$$          call prsed5 ( mode0 , ssite(ib)%rv_a_ov0 , nr , 1 , nsp )
c$$$
c$$$ckino Dec.28.2011:              call splrho ( 31 , nsp , nr , 1 , rv_p_ov0 , iwdummy , iwdummy
c$$$ckino Dec.28.2011:              call splrho ( 31 , nsp , nr , 1 , ssite(ib)%rv_p_ov0 , iwdummy , iwdummy
c$$$          call splrho ( 31 , nsp , nr , 1 , ssite(ib)%rv_a_ov0 , iwdummy , iwdummy 
c$$$     .           )
c$$$
c$$$
c$$$C     call prmx('v0',w(ov0),nr,nr,2)
c$$$ckino Dec.28.2011:              call splrho ( 30 , nsp , nr , 1 , rv_p_ov1 , iwdummy , iwdummy
c$$$ckino Dec.28.2011:              call splrho ( 30 , nsp , nr , 1 , ssite(ib)%rv_p_ov1 , iwdummy , iwdummy
c$$$          call splrho ( 30 , nsp , nr , 1 , ssite(ib)%rv_a_ov1 , iwdummy , iwdummy 
c$$$     .           )
c$$$
c$$$
c$$$ckino Dec.28.2011:              call prsed5 ( mode0 , rv_p_ov1 , nr , 1 , nsp )
c$$$ckino Dec.28.2011:              call prsed5 ( mode0 , ssite(ib)%rv_p_ov1 , nr , 1 , nsp )
c$$$          call prsed5 ( mode0 , ssite(ib)%rv_a_ov1 , nr , 1 , nsp )
c$$$
c$$$ckino Dec.28.2011:              call splrho ( 31 , nsp , nr , 1 , rv_p_ov1 , iwdummy , iwdummy
c$$$ckino Dec.28.2011:              call splrho ( 31 , nsp , nr , 1 , ssite(ib)%rv_p_ov1 , iwdummy , iwdummy
c$$$          call splrho ( 31 , nsp , nr , 1 , ssite(ib)%rv_a_ov1 , iwdummy , iwdummy 
c$$$     .           )
c$$$
c$$$
c$$$        endif
c$$$
c$$$C     pnu v0,v1: modified if 1's bit im set
c$$$        if (mod(im,2) .ge. 1) then
c$$$C     call prmx('pnu',pnu,n0,n0,2)
c$$$          call splrho ( 30 , nsp , n0 , 1 , pnu , iwdummy , iwdummy )
c$$$
c$$$          call prsed5(mode0,pnu,n0,1,nsp)
c$$$          call splrho ( 31 , nsp , n0 , 1 , pnu , iwdummy , iwdummy )
c$$$
c$$$
c$$$          i_copy_size=size(ssite(ib)%pnu)
c$$$          call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
c$$$
c$$$C     call prmx('pnu',pnu,n0,n0,2)
c$$$        endif
c$$$
c$$$ 10     continue
c$$$      enddo
c$$$
c$$$C     Species core density modified if im ge 4
c$$$      if (im .ge. 4) then
c$$$        do  is = 1, nspec
c$$$
c$$$
c$$$c     kino Dec.26.2011:             rv_p_orhoca => sspec(is)%rv_p_orhoc
c$$$
c$$$          nr=sspec(is)%nr
c$$$
c$$$c     kino Dec.26.2011:              call splrho ( 30 , nsp , nr , 1 , rv_p_orhoca , iwdummy , iwdummy
c$$$ckino Dec.28.2011:              call splrho ( 30 , nsp , nr , 1 , sspec(is)%rv_p_orhoc , iwdummy , iwdummy
c$$$          call splrho ( 30 , nsp , nr , 1 , sspec(is)%rv_a_orhoc , iwdummy , iwdummy 
c$$$     .           )
c$$$
c$$$
c$$$c     kino Dec.26.2011:             call prsed5 ( mode0 , rv_p_orhoca , nr , 1 , nsp )
c$$$ckino Dec.28.2011:              call prsed5 ( mode0 , sspec(is)%rv_p_orhoc , nr , 1 , nsp )
c$$$          call prsed5 ( mode0 , sspec(is)%rv_a_orhoc , nr , 1 , nsp )
c$$$
c$$$c     kino Dec.26.2011:              call splrho ( 31 , nsp , nr , 1 , rv_p_orhoca , iwdummy , iwdummy
c$$$ckino Dec.28.2011:              call splrho ( 31 , nsp , nr , 1 , sspec(is)%rv_p_orhoc , iwdummy , iwdummy
c$$$          call splrho ( 31 , nsp , nr , 1 , sspec(is)%rv_a_orhoc , iwdummy , iwdummy 
c$$$     .           )
c$$$
c$$$
c$$$
c$$$        enddo
c$$$      endif
c$$$
c$$$      call info0(1,1,0,' New density:')
c$$$      call prsed2 ( nbas , nsp , ssite , sspec , sv_p_orhoat , slat 
c$$$     .     , smrho , k1 , k2 , k3 )
c$$$
c$$$      end subroutine prsed4
c$$$
c$$$      subroutine prsed5(mode,rho,nr,nlml,nsp)
c$$$
c$$$C     - Zero out spin part (mode=0) or charge part (mode=1) of rho
c$$$C     implicit none
c$$$      integer mode,nr,nlml,nsp
c$$$      double precision rho(nr,nlml,nsp)
c$$$
c$$$      if (mode .eq. 2) then
c$$$        call dscal(nr*nlml,-1d0,rho(1,1,2),1)
c$$$        return
c$$$      endif
c$$$      if (mode .eq. 0) then
c$$$        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
c$$$      endif
c$$$      call dpzero(rho,nr*nlml)
c$$$      if (mode .eq. 0) then
c$$$        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
c$$$      endif
c$$$      end subroutine prsed5
c$$$
c$$$
