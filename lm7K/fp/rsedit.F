
c---> takao modified this by hand after Kino's del_pack.
c   python getarg.py < file.F > x
c   mv x file.F

      subroutine rsedit(sopts,mode,sctrl,ssite,sspec,slat,spot,sbz,
     .nbas,nat,nspec)

      use m_struc_def  !Cgetarg
C- Restart file editor
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Ci   sopts :command options performed automatically, before reading
Ci         :from standard input
Ci   mode  :0 ASA I/O
Ci         :1 FP  I/O, mesh density
Ci         :-1 (input only) read nbas,nat,nspec from file
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: zbak lrs
Ci     Stored:    lrs
Ci     Passed to: iors
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
Ci     Passed to: dcopy rdovfa iors
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: nlml nlma osmrho osmpot oorhat
Ci     Stored:    osmrho
Ci     Passed to: dcopy rdovfa iors
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read:
Ci     Stored:
Ci     Passed to: iors
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Co Outputs
Co   rsedit never returns.
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
Cu   19 Dec 08 New exch option
Cu   01 Jul 08 First created
C  ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character sopts*(*)
      integer mode,nbas,nat,nspec,ifi,n0
      parameter (n0=10)
      type(s_ctrl)::sctrl
      type(s_site)::ssite(*),ssite1(nbas), ssite2(nbas)
      type(s_spec)::sspec(*),sspec1(nspec),sspec2(nspec)
      type(s_lat)::slat,slat2
      type(s_pot)::spot,spot2
      type(s_bz)::sbz

C ... Local parameters
      integer nglob,fopna,fopng,iors,igets,a2vec
      integer i,j,k,l,im,m,n,j1,j2,js1,js2,nit(2),nspeck(2),isync,nw,
     .lexch
      logical lnsave,lbin,lsopts,have2s,havers,lfail
      integer nsp,ngabc(3),n1,n2,n3,k1,k2,k3,lrssav,nvl,nchan
      integer ngabc2(3)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
c      integer ospec1,ospec2,osite1,osite2,oslat2,ospot2,
c     .  oorh2,osmrh2,osmpt2,osmrho,osmpot,oorhat,oqmom
      integer:: osmrh2 , osmpt2 , osmrho , osmpot , oorhat
      integer ,allocatable :: orh2_iv(:)
      real(8) ,allocatable :: qmom_rv(:)

      integer nlst(2),ix(10),nbask(2),natk(2),isw
      integer ,allocatable:: ilst(:,:)
      double precision qbg,dgets,fac(2)
      character dc*1, fn*120, fileid*68, outs*150, strn*120

      integer :: i_copy_size
C ... Heap
      integer w(1)
      common /w/ w
C ... data statements
C     data vec0 /0d0,0d0,0d0/

C     Defaults
      fn = 'rst1'
      havers = .false.
      have2s = .false.
C     # spec and # sites in spec 2nd atom density
      nbask(1) = nbas
      natk(1) = nat
      nspeck(1) = nspec
      nbask(2) = nbas
      nspeck(2) = nspec
      natk(2) = nat
C     true if rs data available, but not saved
      lnsave = .false.
      nsp  = nglob('nsp')
      qbg = sctrl%zbak(1)


      nvl=spot%nlml
      nchan=spot%nlma


C     Allocate local species structures; copy original to them
c      call defrr(ospec1,nspec*nint(sspec(1)))
c      call defrr(ospec2,nspec*nint(sspec(1)))
c      call dcopy(nspec*nint(sspec(1)),sspec,1,spec1),1)
c      call dcopy(nspec*nint(sspec(1)),sspec,1,spec2),1)
      sspec1=sspec(1:nspec)
      sspec2=sspec(1:nspec)
C     Allocate local site structures; copy original to them
c      call defrr(osite1,nbas*nint(ssite(1)))
c      call defrr(osite2,nbas*nint(ssite(1)))
c      call dcopy(nbas*nint(ssite(1)),ssite,1,w(osite1),1)
c      call dcopy(nbas*nint(ssite(1)),ssite,1,w(osite2),1)
      ssite1=ssite(1:nbas)
      ssite2=ssite(1:nbas)
C     Allocate second lat structure; copy original to it
c      call defrr(oslat2,nint(slat(1)))
c      call dcopy(nint(slat(1)),slat,1,w(oslat2),1)
      slat2=slat
C     Allocate second pot structure; copy original to it
c      call defrr(ospot2,nint(spot(1)))
c      call dcopy(nint(spot(1)),spot,1,w(ospot2),1)
      spot2=spot

C     Get pointers for density

      osmrho=spot%osmrho
      osmpot=spot%osmpot
      oorhat=spot%oorhat

C     call upack('pot osmrho oorhat',spot,oorh2,osmrh2,0,0,0)
C     Allocate space for second density
C      call upack('lat nabc',slat,ngabc,0,0,0,0)
C      call fftz30(n1,n2,n3,k1,k2,k3)
C      call defcc(osmrh2, -k1*k2*k3*nsp)
C      call pack1('pot osmrho',w(ospot2),osmrh2)

      allocate(qmom_rv(nvl))
Ctakao_ZeroClear_NotRequiered        if (nvl<0) qmom_rv(:)=0.0d0


      dc = sopts(1:1)
      if (dc .ne. ' ') then
        print 301
  301   format(//' Entering the restart file editor. ',
     .  'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the restart file editor.  ',
     .  'Enter ''?'' to see options.')
        lsopts = .false.
      endif
      if (mode .ne. 1) call rx('rsedit not ready for mode ne 1')

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
c#if PRTNOCR
c      print '(/'' Option : '',$)'
c#else
      print '(/'' Option : '')'
c#endif
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif
      call locase(outs)

C ... Parse and execute the next command
c      if (.false.) then
c
c      else

      if(outs .eq. ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',
     .  ' ''?'' to see menu')
        goto 10

C ... Read 2nd density from restart file
      elseif (outs(1:6).eq.'readb ' .or. outs(1:6).eq.'reada ' .or.
     .outs(1:5).eq.'read ' .or. outs(1:7).eq.'readfa ') then
        if (outs(1:7) .eq. 'readfa ') then
C         Make space for smooth density

          i_copy_size=size(slat2%nabc)
          call i8icopy(i_copy_size,slat2%nabc,1,ngabc2,1)

          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
          call defcc(osmrh2, -k1*k2*k3*nsp)

          spot2%osmrho=osmrh2

          allocate(orh2_iv(3*nbask(2)))
Ctakao_ZeroClear_NotRequiered            if (3*nbask(2)<0) orh2_iv(:)=0

          call rdovfa ( nbask ( 2 ) , nspeck ( 2 ) , ssite2 , sspec2 , 
     .    slat2 , spot2 , qbg , orh2_iv )

        else
          call word(outs,2,j1,j2)
          if (outs(1:6).eq.'reada ') then
            if (j2 .lt. j1) then
              fn = 'rsta'
              ifi = fopna('rsta',-1,1)
            else
              fn = outs(j1:j2)
              ifi = fopng(outs(j1:j2),-1,1)
            endif
            lbin = .false.
          else
            if (j2 .lt. j1) then
              fn = 'rst'
              ifi = fopna('rst',-1,4+1)
            else
              fn = outs(j1:j2)
              ifi = fopng(outs(j1:j2),-1,4+1)
            endif
            lbin = .true.
          endif
          lrssav = int(sctrl%lrs)


          sctrl%lrs=32*0+64+128*0

C       Get the number of atoms and species
          k = iors ( - 1 , sctrl , w , w , slat2 , spot2 , sbz , fileid 
     .    , nbask ( 2 ) , natk ( 2 ) , nspeck ( 2 ) , orh2_iv , w , 
     .    nit ( 2 ) , lbin , ifi )

          call info5(10,0,0,' rsedit: '//
     .    'read 2nd density from file "'//trim(fn)//'" ...  '//
     .    'nbas = %i  nat = %i  nspec = %i',nbask(2),natk(1),nspeck(2),
     .    0,0)
C       Allocate species and site structures, fill with available info
c        call defrr(ospec2,nspeck(2)*nint(sspec(1)))
c        call defrr(osite2,nbask(2)*nint(ssite(1)))
c        call uspec(w(ospec2),-1,0,1,0,0,0)
c        call usite(w(osite2),-1,0,1,0,0,0)
          do i=1,nspec
            call uspec_init(sspec2(i))
          enddo
          do i=1,nbas
            call usite_init(ssite2(i))
          enddo
C       call shstru('spec',w(ospec2),1,nspeck(2))
C       call shstru('site',w(osite2),1,nbask(2))
          allocate(orh2_iv(3*nbask(2)))
Ctakao_ZeroClear_NotRequiered          if (3*nbask(2)<0) orh2_iv(:)=0

          call pshpr(1)
          k = iors ( - 2 , sctrl , ssite2 , sspec2 , slat2 , spot2 , sbz 
     .    , fileid , nbask ( 2 ) , natk ( 2 ) , nspeck ( 2 ) , orh2_iv 
     .    , w , nit ( 2 ) , lbin , ifi )

C       Check for unassigned species => reduce basis
          k = 0
          do  i = 1, nbask(2)

            j=ssite2(i)%spec

            if (j .eq. 0) exit
            k = i
          enddo
          if (k .ne. nbask(2)) then
            call info2(10,0,0,'%9f(warning) only %i sites had '//
     .      'associated species ... reducing nbas(2)',k,0)
          endif
          nbask(2) = k
          call poppr
C       Make space for smooth density

          i_copy_size=size(slat2%nabc)
          call i8icopy(i_copy_size,slat2%nabc,1,ngabc2,1)

          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
          call defcc(osmrh2, -k1*k2*k3*nsp)
          call defcc(osmpt2, -k1*k2*k3*nsp)

          spot2%osmrho=osmrh2
          spot2%osmpot=osmpt2

C       Pack volume
C       call shstru('spec',w(ospec2),1,nspeck(2))
C       call shstru('site',w(osite2),1,nbask(2))
          k = iors ( 1 , sctrl , ssite2 , sspec2 , slat2 , spot2 , sbz 
     .    , fileid , nbask ( 2 ) , natk ( 2 ) , nspeck ( 2 ) , orh2_iv 
     .    , w , nit ( 2 ) , lbin , ifi )

          call fclose(ifi)
        endif

        have2s = .true.

        sctrl%lrs=lrssav


        call info0(1,1,0,' Second Density:')
        call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
        call prsed2 ( nbask ( 2 ) , nsp , ssite2 , sspec2 , orh2_iv 
     .  , slat2 , w ( osmrh2 ) , k1 , k2 , k3 )


        goto 10

C ... Read restart file
      elseif (outs(1:4).eq.'rsb ' .or. outs(1:4).eq.'rsa ' .or.
     .outs(1:3).eq.'rs' .or. outs(1:5).eq.'rsfa ') then
        if (outs(1:5) .eq. 'rsfa ') then
          call rdovfa(nbask(1),nspeck(1),ssite1,sspec1,slat,spot,
     .    qbg,w(oorhat))
          fileid = 'rsedit overlap fa densities'
        else
          call word(outs,2,j1,j2)
          if (outs(1:4) .eq. 'rsa ') then
            if (j2 .lt. j1) then
              fn = 'rsta'
            else
              fn = outs(j1:j2)
            endif
            ifi = fopna(fn,-1,1)
            lbin = .false.
          else
            if (j2 .lt. j1) then
              fn = 'rst'
            else
              fn = outs(j1:j2)
            endif
            ifi = fopna(fn,-1,4+1)
            lbin = .true.
          endif
          call info0(10,0,-1,
     .    ' rsedit: read restart file "'//trim(fn)//'"')
          k = iors(1,sctrl,ssite1,sspec1,slat,spot,sbz,fileid,
     .    nbask(1),natk(1),nspeck(1),w(oorhat),w,nit,lbin,ifi)
C       call shstru('spec',w(ospec2),1,nspeck(2))
C       call shstru('site',w(osite2),1,nbask(2))
          call fclose(ifi)
        endif
        havers = .true.
        lnsave = .true.
        call info2(1,1,0,' Density from '//
     .  '%?#n==0#rs file:#overlapping free atoms:#',
     .  isw(outs(1:5).eq.'rsfa '),0)

        i_copy_size=size(slat%nabc)
        call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)

        call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
        call prsed2(nbask(1),nsp,ssite1,sspec1,w(oorhat),
     .  slat,w(osmrho),k1,k2,k3)

C ... write atomic densities to file
      elseif (outs(1:7) .eq. 'wrhoat ') then

        if (havers) then
          call prsed3(nbask(1),nsp,ssite1,sspec1,w(oorhat))
        else
          call info0(0,0,0,' no rst file read ... nothing writtten')
        endif

C ... show
      elseif (outs(1:5) .eq. 'show ') then
        if (havers) then
          call info2(1,1,0,' Density from rs file:  '//
     .    '%i sites, %i species',nbask(1),nspeck(1))
C          call upack1('lat nabc',slat,ngabc)
          call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
          call prsed2(nbask(1),nsp,ssite1,sspec1,w(oorhat),
     .    slat,w(osmrho),k1,k2,k3)
        endif

        if (have2s) then
          call info2(1,1,0,' Second Density:  '//
     .    '%i sites, %i species',nbask(2),nspeck(2))

          call fftz30(ngabc2(1),ngabc2(2),ngabc2(3),k1,k2,k3)
          call prsed2 ( nbask ( 2 ) , nsp , ssite2 , sspec2 , orh2_iv 
     .    , slat2 , w ( osmrh2 ) , k1 , k2 , k3 )

        endif

C ... Spin averaging of densities, or  Exchange l=0 parts of
C     set object 1|2 list1 list2
C     exch 1|2 site1 site2  or exch 1|2 flip site1 site2
      elseif (outs(1:4) .eq. 'set ' .or. outs(1:5) .eq. 'exch ') then

        lexch = 0
        if (outs(1:5) .eq. 'exch ') lexch = 1
        call words(outs,nw)
        if (nw .lt. 4) goto 98

C     Switches for exchange
        if (lexch .ne. 0) then

C       Determine whether first or 2nd density
          call word(outs,2,j1,j2)
          if (outs(j1:j2) .eq. '1') then
            lexch = 1
          elseif (outs(j1:j2) .eq. '2') then
            lexch = 2
          else
            call info0(0,0,0,' 2rd argument must be 1 or 2 in '//
     .      trim(outs))
            goto 98
          endif

C       Check for flip
          call word(outs,3,j1,j2)
          if (outs(j1:j2) .eq. 'flip') then
            if (nw .ne. 5) goto 98
            ! Must be 5 arguments in this case
            isync = 4
            nw = 3
          else
            isync = 2
            nw = 2
          endif
          im = 0

C       Get 1st and 2nd site
          call numsyv(n)
          do  k = 1, 2
            call word(outs,nw+k,j1,j2)
            call lodsyv('nbas',0,dble(nbask(lexch)),i)
            i = 0
            i = a2vec(outs(j1:j2),len(outs(j1:j2)),i,2,' ',1,1,1,ix,
     .      nlst(k))
            if (i .le. 0 .or. nlst(k) .gt. nbask(lexch)) then
              call info0(0,0,0,
     .        '%6p... invalid site index : '//outs(j1:j2))
              goto 98
            endif
          enddo
          call clrsyv(n)

C       Carry out operation on 1st or second density
          if (lexch .eq. 1) then
            ix(1) = nbask(1); ix(2) = nbask(1)
            if (havers) then
              call prsed1(im*10+isync,1,nlst,ix,nsp,
     .        ssite1,sspec1,w(oorhat),
     .        ssite1,sspec1,w(oorhat),
     .        w(osmrho),ngabc(1),ngabc(2),ngabc(3),
     .        w(osmrho),ngabc(1),ngabc(2),ngabc(3),
     .        fac,lfail)
              if (lfail) then
                call info0(0,0,0,' Densities failed to synchronize')
                goto 98
              endif
            else
              call info0(0,0,0,' no rst data ... cannot do: '//trim(outs))
              goto 10
            endif
          else
            ix(1) = nbask(2); ix(2) = nbask(2)
            if (have2s) then
              call prsed1 ( im * 10 + isync , 1 , nlst , ix , nsp , ssite2 
     .        , sspec2 , orh2_iv , ssite2 , sspec2 , orh2_iv , w ( osmrh2 
     .        ) , ngabc2 ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , w ( osmrh2 ) 
     .        , ngabc2 ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , fac , lfail )

            else
              if (lfail) then
                call info0(0,0,0,' Densities failed to synchronize')
                goto 98
              endif
              goto 10
            endif
          endif

C     Switches for set
        else
          call word(outs,2,j1,j2)
          im = -1
          if (outs(j1:j2) .eq. 'all') im = 4+2+1
          if (outs(j1:j2) .eq. 'n') im = 4
          if (outs(j1:j2) .eq. 'v') im = 2
          if (outs(j1:j2) .eq. 'pnu') im = 1
          if (outs(j1:j2) .eq. 'ntrue') im = 0
          if (im .eq. -1) then
            call info0(0,0,0,
     .      '     ... mode "'//outs(j1:j2)//'" not recognized')
            goto 98
          endif

C       Fourth argument same for each density
          call word(outs,4,j1,j2)
          if (outs(j1:j2) .eq. 'zers') then
            k = im*10+0
          elseif (outs(j1:j2) .eq. 'zerq') then
            k = im*10+1
          elseif (outs(j1:j2) .eq. 'flip') then
            k = im*10+2
          else
            call info0(0,0,0,' failed to recognize 4th argument in: '//
     .      trim(outs))
            goto 98
          endif

C       Carry out operation on 1st or second density
          call word(outs,3,j1,j2)
          if (outs(j1:j2) .eq. '1') then
            if (havers) then
              call prsed4(k,nbask(1),nspeck(1),nsp,ssite1,sspec1,
     .        w(oorhat),slat,w(osmrho),k1,k2,k3)
            else
              call info0(0,0,0,' no rst data ... cannot do: '//trim(outs))
              goto 10
            endif
          elseif (outs(j1:j2) .eq. '2') then
            if (have2s) then
              call prsed4 ( k , nbask ( 2 ) , nspeck ( 2 ) , nsp , ssite2 , 
     .        sspec2 , orh2_iv , slat2 , w ( osmrh2 ) , k1 , k2 , k3 )

            else
              call info0(0,0,0,' no 2nd density... cannot do: '//
     .        trim(outs))
              goto 10
            endif
          else
            call info0(0,0,0,' failed to recognize 3rd argument in '//
     .      trim(outs))
            goto 98
          endif
        endif

        lnsave = .true.

C ... Linear combination of densities
      elseif (outs(1:4) .eq. 'add ') then
        if (.not. (havers .and. have2s)) goto 97

        call words(outs,nw)
        if (nw .lt. 3) goto 98
        call word(outs,2,j1,j2)
        im = -1
        if (outs(j1:j2) .eq. 'all') im = 4+2+1
        if (outs(j1:j2) .eq. 'n') im = 4
        if (outs(j1:j2) .eq. 'v') im = 2
        if (outs(j1:j2) .eq. 'pnu') im = 1
        if (outs(j1:j2) .eq. 'ntrue') im = 0
        if (im .eq. -1) then
          call info0(0,0,0,
     .    '     ... mode "'//outs(j1:j2)//'" not recognized')
          goto 98
        endif
        call word(outs,3,j1,j2)
        if (outs(j1:j2) .eq. 'sync') then
          isync = 0
C         Require 3 words or 5 words
          if (im .ge. 4 .and. nw .ne. 3 .or.
     .    im .lt. 4 .and. nw .ne. 5) goto 98
          nw = 3
        else
          isync = 1
C         Require 4 words or 6 words
          if (im .ge. 4 .and. nw .ne. 4 .or.
     .    im .lt. 4 .and. nw .ne. 6) goto 98
          j = 0
          j = a2vec(outs(j1:),len(outs(j1:)),j,4,', ',2,-3,2,ix,fac)
          if (j .ne. 2) then
            call word(outs,4,j,j2)
            call info0(0,0,0,
     .      '     ... Could not parse arguments : '//outs(j1:j2))
            goto 98
          endif
          nw = 4
        endif
        call numsyv(n)
        if (im .lt. 4) then
          do  k = 1, 2
            call word(outs,nw+k,j1,j2)
            call lodsyv('nbas',0,dble(nbask(1)),i)
            if (k .eq. 2) call lodsyv('nbas',0,dble(nbask(2)),i)
C           call shosyv(0,0,0,6)
            call mkils0(outs(j1:j2),nlst(k),j)
            if (nlst(k) .le. 0) then
              call info0(0,0,0,
     .        '%6p... Bad or null list : '//outs(j1:j2))
              goto 98
            endif
            if (k .eq. 1) then
              allocate(ilst(nlst(k),2))
              call mkilst(outs(j1:j2),nlst(k),ilst)
              call imxmn(nlst(k),ilst,1,l,m)
            elseif (k .eq. 2) then
              if (nlst(2) .eq. 1) then
                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
                call ivset(ilst(1,2),2,nlst(1),ilst(1,2))
              elseif (nlst(2) .ne. nlst(1)) then
                call info2(0,0,0,'%5p... 2nd list: '//outs(j1:j2)//
     .          ' has %i entries when %i are needed',nlst(2),
     .          nlst(1))
                goto 98
              else
                call mkilst(outs(j1:j2),nlst(2),ilst(1,2))
              endif
C             Check that elements in list are within range
              call imxmn(nlst(k),ilst(1,k),1,l,m)
            endif
            if (l .lt. 1 .or.  k.eq.1 .and. m .gt. nbask(1)
     .      .or.  k.eq.2 .and. m .gt. nbask(2)) then
              call info0(0,0,0,
     .        '     ... Entry in list out of range: '//outs(j1:j2))
              goto 98
            endif
          enddo
        elseif (im .ge. 4) then
          if (nbask(1) .ne. nbask(2)) then
            call info0(0,0,0,
     .      '%10fadda requires 1-1 correspondence between sites'//
     .      ' in 1st and 2nd densities')
            goto 98
          endif
          nlst(1) = nbask(1)
          nlst(2) = nbask(1)
          allocate(ilst(nlst(1),2))
          strn = ' '
          call awrit1('1:%i',strn,len(strn),0,nbask(1))
          call mkilst(strn,nlst(1),ilst)
          call icopy(nbask(1),ilst,1,ilst(1,2),1)
        endif
        call clrsyv(n)

C        call upack1('lat nabc',slat,ngabc)
C        call upack1('lat nabc',w(oslat2),ngabc2)
        call prsed1 ( im * 10 + isync , nlst , ilst , nbask , nsp , ssite1 
     .  , sspec1 , w ( oorhat ) , ssite2 , sspec2 , orh2_iv , w ( 
     .  osmrho ) , ngabc ( 1 ) , ngabc ( 2 ) , ngabc ( 3 ) , w ( osmrh2 
     .  ) , ngabc2 ( 1 ) , ngabc2 ( 2 ) , ngabc2 ( 3 ) , fac , lfail 
     .  )

        if (lfail) then
          call info0(0,0,0,' Densities failed to synchronize')
        elseif (isync .eq. 0) then
          call info0(0,0,0,' Densities synchronized')
        elseif (isync .eq. 1) then
C         call upack1('lat nabc',slat,ngabc)
          call fftz30(ngabc(1),ngabc(2),ngabc(3),k1,k2,k3)
          call info0(1,1,0,' New density from rs file:')
          call prsed2(nbask(1),nsp,ssite1,sspec1,w(oorhat),
     .    slat,w(osmrho),k1,k2,k3)
        endif

        deallocate(ilst)
        lnsave = .true.
        goto 10

C ... Save
      elseif (outs(1:5) .eq. 'save ' .or. outs(1:6) .eq. 'savea ') then
        lbin = outs(1:5) .eq. 'save '
        call word(outs,2,j1,j2)
        if (j2 .ge. j1) fn = outs(j1:j2)
        if (lbin) then
          if (j2 .lt. j1) fn = 'rst'
          ifi = fopna(fn,-1,4)
        else
          if (j2 .lt. j1) fn = 'rsta'
          ifi = fopna(fn,-1,0)
        endif
        call info0(10,0,-1,
     .  ' rsedit: write restart file "'//trim(fn)//'"')
        k = iors(1,sctrl,ssite1,sspec1,slat,spot,sbz,fileid,
     .  nbask(1),natk(1),nspeck(1),w(oorhat),w,nit,lbin,-ifi)

        lnsave = .false.

C ... abort
      elseif (outs(1:2) .eq. 'a ') then
        call rx0('aborting rs editor ... no file written')

C ... quit
      elseif (outs(1:2) .eq. 'q '. or. outs(1:5) .eq. 'quit ') then
        if (lnsave) then
          print '('' rs file not saved ... really quit?'')'
          read(*,'(a150)') outs
          call locase(outs)
          if (.not. (outs(1:1) .eq. 'y' .or. outs(1:1) .eq. 'q'))
     .    goto 10
        endif
        call rx0('exit rs editor')

C ... help
      elseif (outs .eq. '?') then
        print 310
        print 311
        print 312
        print 313
        print 314
  310   format(
     .  ' Select one of these options:'/
     .  t4,'rs [fn]',t15,
     .  'read density and associated parameters from binary restart ',
     .  'file.'/t15,'Use "fn" for file name, ',
     .  'or "rst.ext" if "fn" is not supplied'/
     .  t4,'rsa [fn]',t15,
     .  'read density and associated parameters from ascii restart ',
     .  'file.'/t15,'Use "fn" for file name, ',
     .  'or "rst.ext" if "fn" is not supplied'/
     .  t4,'rsb [fn]',t15, 'same as "rs"'/
     .  t4,'rsfa [fn]',t15,
     .  'generate density and associated parameters from atm ',
     .  'file.'/t15,'Use "fn" for file name, ',
     .  'or "atm.ext" if "fn" is not supplied')
  311   format(/
     .  t4,'read [fn]',t15,
     .  'read 2nd density from binary restart file.  Use "fn" for ',
     .  'file name,'/t15,'or "rst.ext" if "fn" is not supplied'/
     .  t4,'reada [fn]',t15,
     .  'read 2nd density from ascii restart file.  Use "fn" for ',
     .  'file name,'/t15,'or "rst.ext" if "fn" is not supplied'/
     .  t4,'readb [fn]',t15, 'same as "read"'//
     .  t4,'wrhoat',t15, 'Write atomic densities to rhoMT.{1..nbas}'//
     .  t4,'show',t15, 'Show summary information about densities'//
     .  t4,'...The following affect the charge or magnetic moment',
     .  ' of a density.'/
     .  t4,'set all|n       1|2  zers|zerq|flip'/
     .  t4,'set v|ntrue|pnu 1|2  zers|zerq|flip'/
     .  t14,'* The first argument specifies which objects are affected.'
     .  /
     .  t16,'"all" applies to all quantities, "n" to the ',
     .  'three density types,'/
     .  t16,'"v" to the spherical potential, ',
     .  '"ntrue" to the true local density,'/
     .  t16,'"pnu" to the linearization parameters'/
     .  t14,'* The second argument refers to the first or second ',
     .  'density'/
     .  t14,'* Third argument = zers: zeros spin part of density:   ',
     .  'n+ - n- = 0'/
     .  t14,'  Third argument = zerq: zeros charge part of density: ',
     .  'n+ + n- = 0'/
     .  t14,'  Third argument = flip: exchange n+, n- ')

  312   format(/
     .  t4,'...The following put a linear combination of the ',
     .  'two densities into the first, or'/
     .  t4,'checks whether the two densities are well enough ',
     .  'synchronized to make it possible'/
     .  t4,'add all|n',t22,'sync'/
     .  t4,'add v|ntrue|pnu',t22,'sync lst1 lst2'/
     .  t4,'add all|n',t22,'fac1 fac2'/
     .  t4,'add v|ntrue|pnu',t22,'fac1 fac2 lst1 lst2'/
     .  t14,'* The first argument''s meaning is the same as the ',
     .  'first argument to "set"'/
     .  t14,'* Modifier "sync" checks whether the densities are ',
     .  'compatible, making no changes.'/
     .  t14,'* If "sync" is missing, the objects specified will be'/
     .  t16,'replaced by fac1*[1st type] + fac2*[second type].'/
     .  t14,'* For the site-only cases "v" and "ntrue," ',
     .  'site list "lst1" is required.'/
     .  t16,'"lst2" is a matching site list for the 2nd density, OR'/
     .  t16,'"lst2" is single site (same site added to all sites)'/)

  313   format(
     .  t4,'...The following exchanges the l=0 parts of two site',
     .  ' densities '/
     .  t4,'exch 1|2',t22,'[flip] site1 site2'/
     .  )

  314   format(
     .  t4,'save  [fn]',t15,'saves restart data in ',
     .  'binary restart file (name="rst" unless fn supplied)'/
     .  t4,'savea [fn]',t15,'saves restart data in ',
     .  'ascii restart file (name="rsta" unless fn supplied)'/
     .  t4,'q',t15,'to quit the editor'/
     .  t4,'a',t15,'to abort')

      else
        print '(1x,''unrecognized option: '',a)', trim(outs)

      endif
      goto 10

   97 call info0(0,0,0,'%10pBoth densities must be supplied before'
     .//'invoking this command ')
   98 call info0(0,0,0,' rsedit:  improper usage of '//trim(outs)//
     .' ... nothing done')
      goto 10

      if (allocated(orh2_iv)) deallocate(orh2_iv)
      if (allocated(qmom_rv)) deallocate(qmom_rv)

      end subroutine rsedit


      subroutine prsed1(mode,nlst,ilst,nbask,nsp,site1,spec1,orh1,
     .site2,spec2,orh2,
     .smrho,k11,k21,k31,
     .smrh2,k12,k22,k32,
     .fac,lfail)

      use m_struc_def  !Cgetarg

C- Add fac*2nd density to density
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 check synchronization
Ci         :1 add the density
Ci         :2 swap l=0 part of 2 local densities.  nlst should be 1
Ci         :4 like 2, but exchange spin in the swap
Ci         :10s digit:
Ci         :0 work with true density only
Ci         :1 work with pnu
Ci         :2 work with v0,v1
Ci         :4 work with n1,n2,n0
Ci         :Any combination is allowed
Ci   nlst  :number of site densities to which 2nd density is added
Ci   ilst  :list of indices to site densities
Ci         :ilst(:,1) is list of sites to be modified
Ci         :ilst(:,2) is list from which 2nd density is to be found
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci
Ci   site1 :struct for site-specific information; see routine usite
Ci   spec1 :struct for species-specific information; see routine uspec
Ci   orh1  :vector of offsets to arrays containing local site densities
Ci   site2 :Analog of site1 for 2nd density
Ci   spec2 :Analog of spec1 for 2nd density
Ci   orh2  :Analog of orh1 for 2nd density
Ci   fac   :Add fac * 2nd density into density
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Jun 08
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nlst,ilst(nlst,2),nbask(2),nsp,
     .k11,k21,k31,k12,k22,k32
      integer orh1(3,*),orh2(3,*)
      real(8):: fac(2)
      type(s_site)::site1(*)
      type(s_spec)::spec1(*)
      type(s_site)::site2(*)
      type(s_spec)::spec2(*)

      double complex smrho(k11,k21,k31,nsp),smrh2(k12,k22,k32,nsp)
      logical lfail
C ... Local parameters
      integer il,ib,jb,is,js,lmxa,lmxl,lmxa2,lmxl2,nr1,nr2,nrmx
      integer nglob,stdo,ipr,nlml1,nlml2,isync,im,n0
      integer ov01,ov11,ov02,ov12,orhoc1,orhoc2
      integer ,allocatable:: flgspc(:)
      parameter (nrmx=1501,n0=10)
      double precision rmt1,rmt2,a1,a2,rofi(nrmx),pnu(n0,2),pnu2(n0,2)
      character spid1*8,spid2*8
      logical lfails
C ... Heap
      integer w(1)
      common /w/ w

      stdo = nglob('stdo')
      call getpr(ipr)
      isync = mod(mode,10)
      im  = mod(mode/10,10)

      lfails = .false.
      if (im .ge. 4) then
        call info0(0,1,-1,' Checking whether sm density meshes '//
     .  'synchronize ...')
        if (k11 .ne. k12 .or. k21 .ne. k22 .or. k31 .ne. k32) then
          write(stdo,212) k11,k21,k31,k12,k22,k32
  212     format(' mesh mismatch:'/
     .    ' first  density n1 n2 n3 = ',3i5/
     .    ' second density n1 n2 n3 = ',3i5)
          lfails = .true.
          if (isync .eq. 1) call rx('sm density mesh mismatch')
        else
          call info0(0,0,0,' ok')
        endif
        if (isync .eq. 1) then
          call info2(0,1,-1,' ... Replace smooth density n0 by linear '
     .    //'combination %d*[n0] + %d*[2nd n0]',fac(1),fac(2))
          call lcrho(k11*k21*k31*2,nsp,1,1,fac(1),fac(2),smrho,smrh2)
        endif
      endif

      if (isync .ge. 2) then
        call info5(0,1,0,' ... Swap '//
     .  '%?#n==4#and spin flip ##'//
     .  '%?#n==0# true site# site #%-1j'//
     .  '%?#n==0# density n1##%-1j'//
     .  '%?#n>=4# densities##'//
     .  '%?#n>=2# potentials##'//
     .  '%?#n>=1# pnu##',
     .  isync,im,mod(im,4),mod(im,2),0)
        if (isync .eq. 3 .or. isync .gt. 4 .or. nlst .gt. 1) then
          call info0(0,1,0,' ... illegal parameters (aborting)')
          return
        endif
      elseif (isync .eq. 1) then
        call info5(0,1,0,' ... Replace '//
     .  '%?#n==0# true site# site #%-1j'//
     .  '%?#n==0# density n1##%-1j'//
     .  '%?#n>=4# densities##'//
     .  '%?#n>=2# potentials##'//
     .  '%?#n>=1# pnu##'//
     .  '  by linear combination'//
     .  ' %d*[1st] + %d*[2nd]',im,
     .  mod(im,4),mod(im,2),fac(1),fac(2))
      elseif (isync .eq. 0) then
        call info0(0,1,0,' Checking whether radial meshes '//
     .  'synchronize ...')
      endif
      call info2(0,0,0,' site 1 : %n:1,3i',nlst,ilst(1,1))
      call info2(0,0,0,' site 2 : %n:1,3i',nlst,ilst(1,2))

  200 format(/' ... Improper synchronization between radial meshes'/
     .'  site1',8x,'site2',8x,'nr1  nr2',3x,'rmt1   rmt2',
     .4x,' a1    a2')

      lfail = .false.
      allocate(flgspc(nbask(1)))
      call iinit(flgspc,nbask(1))
      do  il = 1, nlst

        ib = ilst(il,1)
        if (ib .gt. nbask(1))
     .  call rxi('illegal site index for density:',ib)
        jb = ilst(il,2)
        if (jb .gt. nbask(2))
     .  call rxi('illegal index for 2nd density:',jb)


        is=site1(ib)%spec
        i_copy_size=size(site1(ib)%pnu)
        call dcopy(i_copy_size,site1(ib)%pnu,1,pnu,1)
        ov01=site1(ib)%ov0
        ov11=site1(ib)%ov1


        js=site2(jb)%spec
        i_copy_size=size(site2(jb)%pnu)
        call dcopy(i_copy_size,site2(jb)%pnu,1,pnu2,1)
        ov02=site2(jb)%ov0
        ov12=site2(jb)%ov1


        lmxa=spec1(is)%lmxa
        lmxl=spec1(is)%lmxl


        lmxa2=spec2(js)%lmxa
        lmxl2=spec2(js)%lmxl

        do i_spacks=is,is
          call spacks_copy('u',spec1(i_spacks)%name,is,is,spid1,i_spacks)
        enddo

        do i_spacks=js,js
          call spacks_copy('u',spec2(i_spacks)%name,js,js,spid2,i_spacks)
        enddo


        if (lmxa .eq. -1) goto 10


        a1=spec1(is)%a
        nr1=spec1(is)%nr
        rmt1=spec1(is)%rmt


        a2=spec2(js)%a
        nr2=spec2(js)%nr
        rmt2=spec2(js)%rmt


C       Radial meshes must match
        if (abs(a1-a2).gt.1d-6 .or. nr1.ne.nr2 .or.
     .  abs(rmt1-rmt2).gt.1d-6) then

          if (.not. lfail) write(stdo,200)
          lfail = .true.
          write(stdo,202) ib,spid1,jb,spid2,
     .    nr1,nr2,rmt1,rmt2,a1,a2
  202     format(i4,1x,a,i4,1x,a,2i5,1x,2f7.4,1x,2f6.3)
          if (isync .eq. 1) call rx('site density mismatch')

        endif
        if (isync .eq. 0) goto 10

        nlml1 = (lmxl+1)**2
        nlml2 = (lmxl2+1)**2
        call radmsh(rmt1,a1,nr1,rofi)

C       True local density: modified for im=0 or im>3
        if (im .lt. 1 .or. im .ge. 4) then
          if (isync .lt. 2) then
            call lcrho(nr1,nsp,nlml1,nlml2,fac(1),fac(2),
     .      w(orh1(1,ib)),w(orh2(1,jb)))
          else
            call swrho(isync/4,nr1,nsp,1,nlml1,nlml2,
     .      w(orh1(1,ib)),w(orh2(1,jb)))
          endif
        endif
C       Smooth local density and site core density: modified if im>=4
        if (im .ge. 4) then
          if (isync .lt. 2) then
            call lcrho(nr1,nsp,nlml1,nlml2,fac(1),fac(2),
     .      w(orh1(2,ib)),w(orh2(2,jb)))
            call lcrho(nr1,nsp,1,1,fac(1),fac(2),
     .      w(orh1(3,ib)),w(orh2(3,jb)))
          else
            call rx('not ready for this swap')
          endif
        endif
C       Potentials v0,v1: modified if 2's bit im set
        if (mod(im,4) .ge. 2) then
          if (isync .lt. 2) then
            call lcrho(nr1,nsp,1,1,fac(1),fac(2),w(ov01),w(ov02))
            call lcrho(nr1,nsp,1,1,fac(1),fac(2),w(ov11),w(ov12))
          else
            call rx('not ready for this swap')
          endif
        endif
C       Pnu modified if 1's bit im set
        if (mod(im,2) .ge. 1) then
          if (isync .lt. 2) then
            call lcrho(n0,nsp,1,1,fac(1),fac(2),pnu,pnu2)

            i_copy_size=size(site1(ib)%pnu)
            call dcopy(i_copy_size,pnu,1,site1(ib)%pnu,1)

          else
            call rx('not ready for this swap')
          endif
        endif

C       Species core density modified if im ge 4
        if (im .ge. 4 .and. flgspc(is) .eq. 0) then

          orhoc1=spec1(is)%orhoc


          orhoc2=spec2(js)%orhoc

          if (isync .lt. 2) then
            call lcrho(nr1,nsp,1,1,fac(1),fac(2),w(orhoc1),w(orhoc2))
          else
            call rx('not ready for this swap')
          endif
          flgspc(is) = 1
        endif

   10   continue
      enddo

      if (lfail) then
        call info0(0,0,0,' ... Radial meshes failed to synchronize')
      elseif (isync .eq. 0) then
        call info0(0,0,0,' ... Radial meshes synchronized')
C      elseif (isync .eq. 1) then
C        call info0(0,0,0,' ... Linear combination successful')
      endif

      lfail = lfail .or. lfails

      end subroutine prsed1


      subroutine prsed2(nbas,nsp,ssite,sspec,orhoat,slat,smrho,k1,k2,k3)

      use m_struc_def  !Cgetarg

C- Print out density information
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspec
Ci   sspec :struct for species-specific information; see routine uspec
Ci   nspec2:number of species for second density
Ci   sspec2
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Jun 08
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,k1,k2,k3
      integer orhoat(3,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer nrmx,nlmx
      parameter (nrmx=1501, nlmx=64)
      integer ib,is,nr,lmxa,lmxl,lmxb,nlml,ipr,nglob,stdo,n1,n2,n3
      double precision z,qc,rg,a,rmt,rsma
      double precision pi,srfpi,y0,rofi(nrmx),rwgt(nrmx),qv1,qv2,am1,am2
      double precision ddot,smq,sum2,alat,plat(3,3),qlat(3,3),vol,smag
      character spid*8
C     double complex smpot(k1,k2,k3,2)
C ... Heap
      integer w(1)
      common /w/ w

      stdo = nglob('stdo')
      call getpr(ipr)
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0 = 1d0/srfpi
      n1 = k1
      n2 = k2
      n3 = k3

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)

      call dinv33(plat,1,qlat,vol)
      vol = dabs(vol)*(alat**3)
      if (nsp .eq. 2) then
        call daxpy(k1*k2*k3*2,1d0,smrho(1,1,1,2),1,smrho,1)
      endif
      call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smq,sum2)
      smag = 0
      if (nsp .eq. 2) then
        call daxpy(k1*k2*k3*2,-1d0,smrho(1,1,1,2),1,smrho,1)
        call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,smag,sum2)
        smag = 2*smag - smq
      endif

C      if (nsp .eq. 1) write(stdo,100)
C      if (nsp .eq. 2) write(stdo,201)

      call info5(0,0,-1,' Smooth density: vol=%,2;2d  n123 = %i %i %i',
     .vol,k1,k2,k3,0)
      call info5(0,0,0,'  qsm = %,6;6d%?#n==2#  mag. mom = %,6;6d##',
     .smq,nsp,smag,0,0)

      if (nsp .eq. 1) write(stdo,200)
      if (nsp .eq. 2) write(stdo,201)
  200 format(/' Site parameters and charges:'/
     .'   site',10x,'z     rmt     nr   a  lmxl   qtrue',6x,'q2')
  201 format(/' Site parameters and charges:'/
     .'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
     .5x,'mtrue',5x,'q2        m2')

      do  ib = 1, nbas


        is=ssite(ib)%spec

C        call upack('site spec pnu ov0 ov1',ssite,ib,is,pnu,ov0,ov01)
C        call upack2('site pz',ssite,ib,pnz)

        z=sspec(is)%z
        qc=sspec(is)%qc
        rg=sspec(is)%rg

        do i_spacks=is,is
          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
        enddo


        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt


        rsma=sspec(is)%rsma
        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl
        lmxb=sspec(is)%lmxb

        if (lmxa .eq. -1) goto 10
        nlml = (lmxl+1)**2

        call radmsh(rmt,a,nr,rofi)
        call radwgt(rmt,a,nr,rwgt)

        call splrho(0,nsp,nr,nlml,
     .  w(orhoat(1,ib)),w(orhoat(2,ib)),w(orhoat(3,ib)))
        qv1   = srfpi*ddot(nr,rwgt,1,w(orhoat(1,ib)),1)
        qv2   = srfpi*ddot(nr,rwgt,1,w(orhoat(2,ib)),1)
        call splrho(1,nsp,nr,nlml,
     .  w(orhoat(1,ib)),w(orhoat(2,ib)),w(orhoat(3,ib)))
        am1   = 2*srfpi*ddot(nr,rwgt,1,w(orhoat(1,ib)),1) - qv1
        am2   = 2*srfpi*ddot(nr,rwgt,1,w(orhoat(2,ib)),1) - qv2

        if (nsp .eq. 1) then
          write(stdo,202) ib,spid,z,rmt,nr,a,lmxl,qv1,am1
        else
          write(stdo,202) ib,spid,z,rmt,nr,a,lmxl,qv1,am1,qv2,am2
        endif
  202   format(i4,1x,a,f6.1,f9.5,i5,f6.3,i3,4f10.5)

   10   continue
      enddo


      end subroutine prsed2


      subroutine prsed3(nbas,nsp,ssite,sspec,orhoat)

      use m_struc_def  !Cgetarg

C- Read/write atomic densities
C     implicit none
C ... Passed parameters
      integer nbas,nsp
      integer orhoat(3,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C ... Local parameters
      integer nrmx,nlmx
      parameter (nrmx=1501, nlmx=64)
      integer ib,is,nr,lmxa,lmxl,nlml
      double precision a,rmt
      double precision rofi(nrmx)
C ... Heap
      integer w(1)
      common /w/ w


      do  ib = 1, nbas


        is=ssite(ib)%spec


        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl

        if (lmxa .eq. -1) goto 10

        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt

        nlml = (lmxl+1)**2
        call radmsh(rmt,a,nr,rofi)

        call wrhomt('rhoMT.','density',ib,w(orhoat(1,ib)),rofi,nr,
     .  nlml,nsp)

   10   continue
      enddo
      end subroutine prsed3


      subroutine prsed4(mode,nbas,nspec,nsp,ssite,sspec,orhoat,
     .slat,smrho,k1,k2,k3)

      use m_struc_def  !Cgetarg

C- Remove spin or charge part of atomic densities
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 remove spin part, 1 remove charge part, 2, flip spins
Ci         :10s digit:
Ci         :0 work with true density only
Ci         :1 work with pnu
Ci         :2 work with v0,v1
Ci         :4 work with n1,n2,n0
Ci         :Any combination is allowed
Ci   nbas  :number of atoms
Ci   sspec :struct for species-specific information; see routine uspec
Ci   nspec2:number of species for second density
Ci   sspec2
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Jun 08
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,nspec,nsp,k1,k2,k3
      integer orhoat(3,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer nrmx,nlmx,n0
      parameter (nrmx=1501, nlmx=64, n0=10)
      integer ib,is,nr,lmxa,lmxl,lmxb,nlml,ipr,nglob,stdo,im,mode0
      double precision z,qc,rg,a,rmt,rsma
      double precision pi,srfpi,y0,rofi(nrmx),rwgt(nrmx),pnu(n0,2)
      integer ov0,ov1,orhoca
      character spid*8
C ... Heap
      integer w(1)
      common /w/ w

      stdo = nglob('stdo')
      call getpr(ipr)
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0 = 1d0/srfpi
      mode0 = mod(mode,10)
      im  = mod(mode/10,10)
      if (mode0 .eq. 0 .and. nsp .eq. 1) return
      if (mode0 .eq. 2 .and. nsp .eq. 1) return

      if (im .ge. 4) then
C       Replace smrho-up,smrho-down with smrho, smrho(up-down)
        call splrho(30,nsp,k1*k2*k3*2,1,smrho,smrho,smrho)
        call prsed5(mode0,smrho,k1*k2*k3*2,1,nsp)
C       Replace smrho, smrho(up-down) with smrho-up,smrho-down
        call splrho(31,nsp,k1*k2*k3*2,1,smrho,smrho,smrho)
      endif

C      if (nsp .eq. 1) write(stdo,200)
C      if (nsp .eq. 2) write(stdo,201)
C  200 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
C     .  6x,'q2')
C  201 format(/'   site',10x,'z     rmt     nr   a  lmxl   qtrue',
C     .  5x,'mtrue',5x,'q2        m2')

      do  ib = 1, nbas


        is=ssite(ib)%spec


        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
        ov0=ssite(ib)%ov0
        ov1=ssite(ib)%ov1

C        call upack2('site pz',ssite,ib,pnz)

        z=sspec(is)%z
        qc=sspec(is)%qc
        rg=sspec(is)%rg

        do i_spacks=is,is
          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
        enddo


        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt


        rsma=sspec(is)%rsma
        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl
        lmxb=sspec(is)%lmxb

        if (lmxa .eq. -1) goto 10
        nlml = (lmxl+1)**2

        call radmsh(rmt,a,nr,rofi)
        call radwgt(rmt,a,nr,rwgt)

C       Replace rho-up,rho-down with rho, up-down
        call splrho(0,nsp,nr,nlml,
     .  w(orhoat(1,ib)),w(orhoat(2,ib)),w(orhoat(3,ib)))

C       True local density: modified for im=0 or im>3
        if (im .lt. 1 .or. im .ge. 4) then
          call prsed5(mode0,w(orhoat(1,ib)),nr,nlml,nsp)
        endif
        if (im .ge. 4) then
          call prsed5(mode0,w(orhoat(2,ib)),nr,nlml,nsp)
          call prsed5(mode0,w(orhoat(3,ib)),nr,1,nsp)
        endif

C       Restore rho-up,rho-down from rho, up-down
        call splrho(1,nsp,nr,nlml,
     .  w(orhoat(1,ib)),w(orhoat(2,ib)),w(orhoat(3,ib)))

C       Potentials v0,v1: modified if 2's bit im set
        if (mod(im,4) .ge. 2) then
C         call prmx('v0',w(ov0),nr,nr,2)
          call splrho(30,nsp,nr,1,w(ov0),w,w)
          call prsed5(mode0,w(ov0),nr,1,nsp)
          call splrho(31,nsp,nr,1,w(ov0),w,w)
C         call prmx('v0',w(ov0),nr,nr,2)
          call splrho(30,nsp,nr,1,w(ov1),w,w)
          call prsed5(mode0,w(ov1),nr,1,nsp)
          call splrho(31,nsp,nr,1,w(ov1),w,w)
        endif

C       pnu v0,v1: modified if 1's bit im set
        if (mod(im,2) .ge. 1) then
C         call prmx('pnu',pnu,n0,n0,2)
          call splrho(30,nsp,n0,1,pnu,w,w)
          call prsed5(mode0,pnu,n0,1,nsp)
          call splrho(31,nsp,n0,1,pnu,w,w)

          i_copy_size=size(ssite(ib)%pnu)
          call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)

C         call prmx('pnu',pnu,n0,n0,2)
        endif

   10   continue
      enddo

C     Species core density modified if im ge 4
      if (im .ge. 4) then
        do  is = 1, nspec


          orhoca=sspec(is)%orhoc
          nr=sspec(is)%nr

          call splrho(30,nsp,nr,1,w(orhoca),w,w)
          call prsed5(mode0,w(orhoca),nr,1,nsp)
          call splrho(31,nsp,nr,1,w(orhoca),w,w)

        enddo
      endif

      call info0(1,1,0,' New density:')
      call prsed2(nbas,nsp,ssite,sspec,orhoat,slat,smrho,k1,k2,k3)
      end subroutine prsed4

      subroutine prsed5(mode,rho,nr,nlml,nsp)

C- Zero out spin part (mode=0) or charge part (mode=1) of rho
C     implicit none
      integer mode,nr,nlml,nsp
      double precision rho(nr,nlml,nsp)

      if (mode .eq. 2) then
        call dscal(nr*nlml,-1d0,rho(1,1,2),1)
        return
      endif
      if (mode .eq. 0) then
        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
      endif
      call dpzero(rho,nr*nlml)
      if (mode .eq. 0) then
        call dswap(nr*nlml,rho(1,1,1),1,rho(1,1,2),1)
      endif
      end subroutine prsed5


