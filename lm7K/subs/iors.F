Cchp1 cmdvar= real(8),rv_p_,oves
#define FP 1
      integer function iors(mode,sctrl,ssite,sspec,slat,spot,sbz,
     .fid,nbas,nat,nspec,orhoat,qnu,nit,lbin,ifi)

      use m_struc_def  !Cgetarg
      use m_struc_func 
C- I/O for charge density
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Ci   mode  :0 ASA I/O
Ci         :1 FP  I/O, mesh density
Ci         :-1 (input only) read nbas,nat,nspec from file
Ci         :-2 (input only) FP species and site data; poke into strux:
Ci         :    lat->alat plat ngabc
Ci         :    site->spec clabel pos pos0 force vel lmxa
Ci         :    spec->name a nr rmt z lmxa lmxl kmxt p pz qc ctail etail
Ci                    stc orhoc nxi exi chfa rsmfa lfoca rfoca
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: iors lgors mpibc1
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: pos class spec force vel pnu pz ov0 ov1
Ci     Stored:    pos pos0 force vel pnu pz ov0 ov1
Ci                NB: data only stored on file read (ifi>0)
Ci     Passed to: iors pvsms2 mpibc1
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: a nr rmt z lmxa lmxl kmxt p pz lfoca qc orhoc idmod
Ci                rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
Ci                chfa rsmfa
Ci     Stored:    name a nr rmt z lmxa lmxl kmxt p pz qc ctail etail
Ci                stc orhoc nxi exi chfa rsmfa
Ci                NB: data only stored on file read (ifi>0)
Ci     Passed to: iors spacks pvsms2 mpibc1
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat nabc dist
Ci     Stored:    *
Ci     Passed to: iors mpibc1
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: osmrho (FP) oves (ASA)
Ci     Stored:    oves (ASA)
Ci     Passed to: iors mpibc1
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: ef def w
Ci     Stored:    ef def w
Ci                NB: data only stored on file read (ifi>0)
Ci     Passed to: *
Ci   fid   :string containing identification
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lbin  :T file I/O in binary mode
Ci         :F file I/O in ascii mode
Cio Inputs/Outputs
Cio  orhoat:vector of offsets containing site density
Cio        :(FP only)
Cio  qnu   :vector of moments (ASA only --- mode=0)
Cio  spot->smrho:smoothed interstitial density ; see Remarks
Cio        :(FP only)
Cio  nit   :iteration number
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
Cr      NOTE: on input, arrays for rhoat and v0 are allocated here
Cr   Smooth density
Cr      real part of complex*16 array smrho contains the density
Cr      k1,k2,k3 are the physical dimensions of the array
Cr      n1,n2,n3 are the dimensions of the mesh.
Cr   Additional information stored:
Cr      fid: file identifier, a string of length 64 characters or less.
Cr      parameters relating to coordinates and molecular dynamics.
Cr   On input, iors tries to transform data format where needed:
Cr      lmxl incompatible: pad or truncate l-components
Cr      FT mesh changed: map onto new mesh
Cl Local variables
Cl   lrs switches:  0=>use rst file data; 1=>ignore rst file data
Cl   lrs(1) site positions
Cl   lrs(2) starting fermi level
Cl   lrs(3) starting pnu's
Cm MPI
Cm   Master process reads and broadcasts line by line. err= and end= are
Cm   troublesome since the slave processes may hang if the rst file is
Cm   corrupted or incompatible. For now if iors returns < 1 lmfp will
Cm   exit and hope the slave processes follow suit!
Cu Updates
Cu   01 Jul 08 New mode -2
Cu   25 Aug 07 Bug fix, case floating orbitals not positioned at end
Cu   10 Jul 07 Cleaner error exit, MPI
Cu   20 Jun 06 Repackaged MPI
Cu   07 Jul 05 (version update 1.04)
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   04 Feb 05 Spin-split file non-spin-polarized density (FP)
Cu   26 Apr 03 Added MPI parallelization for ASA
Cu   11 Jan 03 Bug fix: corrected calc. vol (might have been<0)
Cu   10 Dec 02 File shears atom positions by shear transformation
Cu             (slat->plat) (file plat)^-1
Cu   19 Feb 02 (version update 1.03)
Cu             File now contains nspec
Cu             File contents backwardly compatible with prior versions.
Cu             New mode (-1)
Cu             Routine's argument list changed.
Cu   15 Feb 02 (ATP) Added MPI parallelization for fp
Cu   15 Jan 02 ascii version now labels site rho1,rho2,..
Cu   27 Aug 01 Extended to local orbitals.
Cu   17 May 01 Added ASA file I/O.  New argument list.
Cu   27 Apr 01 Added lbin switch
Cu   25 Jun 00 spin polarized
Cu   21 Apr 00 Adapted from nfp rw_rs
C  ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lbin
      integer mode,nbas,nat,nspec,nit,ifi,orhoat(3,nbas),n0
      parameter (n0=10)
      type(s_ctrl)::sctrl
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_pot)::spot
      type(s_bz)::sbz

      double precision qnu(3,n0,*)
      character fid*(*)
C ... Local parameters
      integer procid,master,mpipid,nproc
      integer i,i0,i1,i2,i3,i4,ib,ipr,iprint,ic,is,is0,isp,jb,k1,k2,k3,
     .igetss,jfi,k11,k21,k31,kmax,kmax0,kmxv,l,lfoc,lfoc0,lmxa,
     .lmxa0,lmxb,lmxb0,lmxl,lmxl0,lmxr,lmxv,lmxv0,lrel,lrel0,n11,n21,
     .n31,nbas0,nspec0,nglob,nsp,nlml,nlml0,npan,npan0,nr,nr0,nsp0,
     .nxi,nat0,ibaug
Cchp1       integer:: orhoca , ov0 , ov1 , ngabc(3) , n1 , n2 , n3 , lgunit 
Cchp1      ., stdo , isw , osmrho , oves
       integer:: orhoca , ov0 , ov1 , ngabc(3) , n1 , n2 , n3 , lgunit 
     . , stdo , isw , osmrho 
       real(8),pointer :: rv_p_oves(:)

      complex(8) ,allocatable :: h_zv(:)
      real(8) ,allocatable :: rwgt_rv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer idmod(n0),idmoz(n0),lrs(10),nl
      logical isanrg,lgors,lfail,ltmp1,ltmp2,latvec,lshear,cmdopt,mlog
      double precision a,a0,alat,alat0,cof,eh,fac,qc,rfoc,rfoc0,rmt,
     .rmt0,rsma,rsma0,rsmfa,rsmr,rsmr0,rsmv,rsmv0,stc,sum,vfac,vol,
     .vol0,vs,vs1,z,z0
      double precision pnu(n0,2),pnz(n0,2),ql(n0,2*n0),pos(3),
     .pos0(3),force(3),vel(3),plat(9),plat0(9),qlat(9),qlat0(9),
     .exi(n0),hfc(n0,2),vec0(3),wk(100),rh,vrmax(2),pnus(n0,2),
     .pnzs(n0,2),dval
      character spid*8,spid0*8,fid0*68,line*20,msg*23,use*80,ignore*80,
     .msgw*17,datimp*32,usernm*32,hostnm*32,jobid*32,ffmt*32,ifmt*32

C ... Heap
      integer w(1)
      common /w/ w
C ... data statements
      data vec0 /0d0,0d0,0d0/


C     call wkprnt(1)

C ... MPI setup
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,ignore)


      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)

      call dinv33(plat,1,qlat,fac)
      vol = dabs(fac)*alat**3
      ipr    = iprint()
      stdo   = lgunit(1)
      if (mode .eq. 0) then
        vs   = -1.05d0
      elseif (mode .eq. -2 .or. mode .eq. -1 .or. mode .eq. 1) then
        vs   =  1.04d0
      endif
      nsp    = nglob('nsp')
      lrel   = nglob('lrel')
      lrs ( 1 ) = isw ( iand(32,int(sctrl%lrs)) .ne.0 )

      lrs ( 2 ) = isw ( iand(64,int(sctrl%lrs)) .ne.0 )

      lrs ( 3 ) = isw ( iand(128,int(sctrl%lrs)) .ne.0 )

      if (mode .eq. -2) then
        lrs(1) = 0
        lrs(3) = 0
      endif
      msg    = '         File mismatch:'
      msgw   = '         warning:'
      iors   = -1
C     ndyn = igets('move size',smove)
      line = 'header'
      if (mode .eq. 1) then
        call fftz30(n1,n2,n3,k1,k2,k3)

        osmrho=spot%osmrho

      endif
      ffmt = '(5f15.10)'
      ifmt = '(20i5)'
C ... Hardwired for now
      npan = 1
      nl = nglob('nl')

C ... Message
      if (ipr .ge. 10) then
        fid0 = 'read'
        i = 5
        if (ifi .lt. 0) then
          fid0 = 'write'
          i = i+1
        endif
        fid0(i:) = ' restart file ('
        i = i+15
        fid0(i:) = 'ascii'
        if (lbin) fid0(i:) = 'binary'
        i = i+5
        if (lbin) i = i+1
        fid0(i:) = ', asa'
        if (mode .eq. -1) fid0(i:) = ' dimensioning data'
        if (mode .eq. 1) fid0(i:) = ', mesh density'
        i = i+5
        if (mode .eq. -1) i = i+18-5
        if (mode .eq. 1) i = i+14-5
        fid0(i:i) = ')'
        i = i+1
        write(stdo,700) fid0(1:i)
  700   format(/' iors  : ',a)
      endif

C --- Input ---
      if (ifi .gt. 0) then
        jfi = ifi
        if (procid .eq. master) then
          rewind jfi
        endif
        use    = '         use from  restart file:'
        ignore = '         ignore in restart file:'

C ...   Read header and general information
        line = 'header'
C       MPI check to see if at least 1st record can be read
C       Abort with error message if file is missing (lfail = .true.)
        lfail = .false.
        if (nproc .gt. 0) then
          if (procid .eq. master) then
            if (lbin) then
              lfail = .true.
              read(jfi,end=996,err=996) vs1
              lfail = .false.
              rewind jfi
  996         continue
            else
              lfail = .true.
              read(jfi,*,end=997,err=997)
              read(jfi,*,end=997,err=997) vs1
              lfail = .false.
              rewind jfi
  997         continue
            endif
          endif
          call mpibc1(lfail,1,1,mlog,'iors','read error')
          if (lfail) goto 998
        endif
        if (procid .eq. master) then
          if (lbin) then
            read(jfi,end=998,err=998) vs1
            read(jfi) fid0
            read(jfi) datimp,usernm,hostnm
            if (abs(vs1) .le. 1.021) then
              nspec0 = nspec
              read(jfi) nbas0,nsp0,npan0,lrel0
              nat0 = nbas0
            else if (abs(vs1) .le. 1.031) then
              read(jfi) nbas0,nsp0,npan0,lrel0,nspec0
              nat0 = nbas0
            else
              read(jfi) nbas0,nat0,nsp0,npan0,lrel0,nspec0
            endif
            read(jfi) nit
            read(jfi) alat0,vol0,plat0
          else
            read(jfi,*,end=998,err=998)
            read(jfi,*,end=998,err=998) vs1
            read(jfi,'(1x,a68)',end=998,err=998) fid0
            read(jfi,'(1x,a32)',end=998,err=998) datimp
            read(jfi,'(1x,a32)',end=998,err=998) usernm
            read(jfi,'(1x,a32)',end=998,err=998) hostnm
            if (abs(vs1) .le. 1.021) then
              nspec0 = nspec
              read(jfi,*,end=998,err=998) nbas0,nsp0,npan0,lrel0,nit
              nat0 = nbas0
            else if (abs(vs1) .le. 1.031) then
              read(jfi,*,end=998,err=998) nbas0,nsp0,npan0,lrel0,nit,
     .        nspec0
              nat0 = nbas0
            else
              read(jfi,*,end=998,err=998)
     .        nbas0,nat0,nsp0,npan0,lrel0,nit,nspec0
            endif
            read(jfi,*) alat0,vol0
            read(jfi,*) plat0
          endif
          if (mode .eq. -2) then

            slat%alat=alat0
            i_copy_size=size(slat%plat)
            call dcopy(i_copy_size,plat0,1,slat%plat,1)

            alat = alat0
            vol = vol0
            call dcopy(9,plat0,1,plat,1)
            call dinv33(plat,1,qlat,fac)
          endif
          call dinv33(plat0,1,qlat0,fac)
          lshear = .not. latvec(3,1d-6,qlat,plat0)
          fid = fid0
          call strip(fid,i,i1)
          call strip(datimp,i,i2)
          call strip(usernm,i,i3)
          call strip(hostnm,i,i4)
          if (ipr .ge. 40) write(stdo,710) fid(1:i1),
     .    usernm(1:i3),hostnm(1:i4),datimp(1:i2)
  710     format(9x,'id -  ',a
     .    /9x,'written by -  ',a,' on ',a,' at: ',a)

          if (mode .eq. -1) then
            nbas = nbas0
            nat  = nat0
            nspec = nspec0
            return
          endif

C       if (isanrg(nspec0, nspec,nspec,msg,'nspec',.true.)) goto 999

C       Number of real atoms may not increase
          if (nat .gt. nat0) then
            if (isanrg(nat0,  nat,nat,msg,'nat',.true.)) goto 999
          elseif (nat .eq. nat0 .and. nbas .ne. nbas0) then
            call info2(10,0,0,
     .      '%9f(warning) mismatch in nbas ... skipping sites'//
     .      '%N%18f expected nbas=%i but rst file has nbas=%i',
     .      nbas,nbas0)
C       OK to reduce nat (e.g. E.S.); these sites will be skipped
          elseif (nat .lt. nat0) then
            call info0(20,0,0,
     .      '%9f(warning) rst mismatch in nat ... skipping sites')
          endif
C       if (isanrg(nsp0,  nsp,nsp,  msg,'nsp', .true.)) goto 999
          if (nsp0 .lt. nsp) call info0(20,0,0,
     .    '%9f(warning) rst file not spin pol .. splitting spins')
          if (isanrg(npan0, npan,npan,msg,'npan',.true.)) goto 999
          lfail = isanrg(lrel0, lrel,lrel,msgw,'lrel',.false.)
          call fsanrg(abs(vs1),1.00d0,abs(vs),0d0,' ','file''s version',
     .    .true.)
        endif
        call mpibc1(nbas0,1,2,mlog,'iors','nbas0')
        call mpibc1(nat0,1,2,mlog,'iors','nat0')
        call mpibc1(plat,9,4,mlog,'iors','plat')
        call mpibc1(nit,1,2,mlog,'iors','nit')

C   --- Read smooth charge density ---
        if (mode .eq. 0) then
C   ... FP mesh density
        elseif (mode .eq. 1 .or. mode .eq. -2) then
#if FP
          line = 'smoothed density'
          if (procid .eq. master) then
            if (lbin) then
              read(jfi,err=999,end=999) n11,n21,n31
            else
              read(jfi,*,err=999,end=999)
              read(jfi,*,err=999,end=999) n11,n21,n31
            endif

            if (mode .eq. -2) then
              n1 = n11
              n2 = n21
              n3 = n31

              i_copy_size=size(slat%nabc)
              call ii8copy(i_copy_size,ngabc,1,slat%nabc,1)

              call fftz30(n1,n2,n3,k1,k2,k3)
              call defcc (osmrho,-k1*k2*k3*nsp)
              call dpdftr(n1,n2,n3,k1,k2,k3,nsp0,w(osmrho),lbin,jfi)
              call rlse(osmrho)

            elseif (n11.eq.n1 .and. n21.eq.n2 .and. n31.eq.n3) then
c          if (ipr .ge. 10) write(stdo,450) n11,n21,n31,n1,n2,n3
cccccccccccccccccc
c          write(6,*)'yyyyyyyy 222111',n11,n21,n31,n1,n2,n3,k1,k2,k3
c          write(6,*)'yyyyyyyy 222444',nsp0,lbin,jfi,osmrho
cccccccccccccccccc
              call dpdftr(n1,n2,n3,k1,k2,k3,nsp0,w(osmrho),lbin,jfi)
              if (nsp .gt. nsp0) then
C           call zprm3('smrho before split',0,w(osmrho),k1,k2,k3)
cccccccccccccccccc
c          write(6,*)'yyyyyyyy 222333'
cccccccccccccccccc
                i = k1*k2*k3*2
                call dscal(i,0.5d0,w(osmrho),1)
                call dpscop(w(osmrho),w(osmrho),i,1,1+i,1d0)
C           call zprm3('smrho after split',0,w(osmrho),k1,k2,k3*2)
              endif
c   ... or read and remesh
            else
              if (ipr .ge. 10) write(stdo,450) n11,n21,n31,n1,n2,n3
  450         format(9x,'remesh density from  ',i4,'  *',i4,'  *',i4,
     .        '    to  ',i4,'  *',i4,'  *',i4)
              call fftz30(n11,n21,n31,k11,k21,k31)
              allocate(h_zv(abs(-k11*k21*k31*nsp)))
              if (-k11*k21*k31*nsp<0) h_zv(:)=0.0d0
              call dpdftr ( n11 , n21 , n31 , k11 , k21 , k31 , nsp0 , h_zv 
     .        , lbin , jfi )
              if (nsp .gt. nsp0) then
C           call zprm3('smrho before split',0,w(oh),k1,k2,k3)
                i = k11*k21*k31*2
                call dscal ( i , 0.5d0 , h_zv , 1 )

                call dpscop ( h_zv , h_zv , i , 1 , 1 + i , 1d0 )

              endif
              call pshpr(50)
              i = 0
              if (n1 .eq. 2*n11 .and. n2 .eq. 2*n21 .and. n3 .eq. 2*n31) i=3
              call chgmsh ( i , plat , nsp , n11 , n21 , n31 , k11 , k21 , 
     .        k31 , h_zv , n1 , n2 , n3 , k1 , k2 , k3 , w ( osmrho ) )
              call poppr
              if (allocated(h_zv)) deallocate(h_zv)

            endif
            if (mode .eq. -2) goto 115
C   ... If cell volume changed, scale smooth density to maintain charge
            vfac = vol0/vol
            if (dabs(vfac-1d0) .gt. 1d-8) then
              if (ipr .ge. 10) write(stdo,460) vol0,vfac
  460         format(9x,'volume changed from',f8.2,
     .        ' :  scale smooth density by',f8.4)
              call dpcopy(w(osmrho),w(osmrho),1,2*k1*k2*k3*nsp,vfac)
            endif
          endif
          call mpibc1(w(osmrho),2*k1*k2*k3*nsp,4,mlog,'iors','smrho')
  115     continue
#else
          call rx('smooth density not implemented')
#endif
        endif

C   --- Read information related to dynamics ---
C       For compatibility with nfp, read record into wk
        if (procid .eq. master) then
          if (lbin) then
            call dpdump(wk,100,jfi)
          else
            read(jfi,*,err=999,end=999)
            call dfdump(wk,100,jfi)
          endif
        endif
        call mpibc1(wk,3,4,mlog,'iors','wk')
        if (lrs(2) .ne. 0) then
          call awrit0('%a ef window,',ignore,len(ignore),0)
        else
          call awrit0('%a ef window,',use,len(use),0)

          sbz%ef=wk(1)
          sbz%def=wk(2)
          sbz%w=wk(3)

        endif

C   --- Read atomic positions,forces,velocities ---
        line = 'site data'
        if (lrs(1) .ne. 0) then
          call awrit0('%a positions,',ignore,len(ignore),0)
        else
          call awrit0('%a positions,',use,len(use),0)
          if (lshear) call awrit0('%a%7b(sheared),',use,len(use),0)
C         Must be enough sites to read from rst file
          if (nbas0 .lt. nbas) then
            call info0(1,0,0,'%9f oops ... cannot use '//
     .      'file site positions (site mismatch)')
            if (isanrg(nbas0, nbas,nbas,msg,'nbas',.true.)) goto 999
          endif
        endif

        do  ib = 1, nbas0
          if (procid .eq. master) then
            if (lbin) then
              read(jfi,err=999,end=999) jb,pos,force,vel
            else
              read(jfi,'(6x,i10)',err=999,end=999) jb
              read(jfi,*,err=999,end=999) pos
              read(jfi,*,err=999,end=999) force
              read(jfi,*,err=999,end=999) vel
            endif
          endif
          if (mode .eq. -2) then

            i_copy_size=size(ssite(ib)%pos)
            call dcopy(i_copy_size,pos,1,ssite(ib)%pos,1)
            i_copy_size=size(ssite(ib)%pos0)
            call dcopy(i_copy_size,pos,1,ssite(ib)%pos0,1)
            i_copy_size=size(ssite(ib)%force)
            call dcopy(i_copy_size,force,1,ssite(ib)%force,1)
            i_copy_size=size(ssite(ib)%vel)
            call dcopy(i_copy_size,vel,1,ssite(ib)%vel,1)

          endif

          call mpibc1(pos,3,4,mlog,'iors','pos')
          call mpibc1(force,3,4,mlog,'iors','force')
          call mpibc1(vel,3,4,mlog,'iors','vel')
          if (ib .gt. nbas) goto 10

C         rst file positions in pos0
          if (lshear) then
            call dgemm('T','N',3,1,3,1d0,qlat0,3,pos,3,0d0,pos0,3)
            call dgemm('N','N',3,1,3,1d0,plat,3,pos0,3,0d0,pos,3)
          endif
          call dcopy(3,pos,1,pos0,1)
          if (lrs(1) .ne. 0) then

            i_copy_size=size(ssite(ib)%pos)
            call dcopy(i_copy_size,ssite(ib)%pos,1,pos,1)

            call dpzero(vel,3)
          endif

          i_copy_size=size(ssite(ib)%pos)
          call dcopy(i_copy_size,pos,1,ssite(ib)%pos,1)
          i_copy_size=size(ssite(ib)%pos0)
          call dcopy(i_copy_size,pos0,1,ssite(ib)%pos0,1)
          i_copy_size=size(ssite(ib)%force)
          call dcopy(i_copy_size,force,1,ssite(ib)%force,1)
          i_copy_size=size(ssite(ib)%vel)
          call dcopy(i_copy_size,vel,1,ssite(ib)%vel,1)

   10     continue
        enddo

C   --- Read information for local densities ---
        if (lrs(3) .ne. 0) then
          call awrit0('%a pnu,',ignore,len(ignore),0)
        else
          call awrit0('%a pnu,',use,len(use),0)
        endif
        if (ipr .ge. 10) then
          call awrit1('%a%?#(p>32)#%b # *#',use,len(use),-stdo,0)
          call awrit1('%a%?#(p>32)#%b # *#',ignore,len(ignore),-stdo,0)
        endif
        if (procid .eq. master) then
          if (.not. lbin) read(jfi,*,err=999,end=999)
        endif
        ibaug = 0
        do  ib = 1, nbas

          ic=ssite(ib)%class
          is=ssite(ib)%spec

C         mode = -2: get species is, lmxa, spid
          if (mode .eq. -2 .and. lbin) then
C           Possibly a floating orbital.
C           If 1st entry > nspec => probably nr, not is.  Assume floating
            read(jfi) is
            backspace jfi
            if (is .gt. nspec) then
              call info2(10,0,0,
     .        '%9f(warning) species mismatch (%i) for ib=%i ... '//
     .        'ignore and assign to species 0',is,ib)

              ssite(ib)%spec=0

              goto 20
            endif
            read(jfi) is,spid,lmxa,lmxl,nr,rmt,a,z,qc
            backspace jfi
            do i_spacks=is,is
              call spacks_copy('p',sspec(i_spacks)%name,is,is,spid,i_spacks)
            enddo


            ssite(ib)%spec=is

            do i_spacks=ib,ib
              call spacks_copy('p',ssite(i_spacks)%clabel,ib,ib,spid,i_spacks)
            enddo


            sspec(is)%lmxa=lmxa

          elseif (mode .eq. -2) then
            read(jfi,'(19x,a8)') spid
            read(jfi,*)
            read(jfi,*) is,lmxa,lmxl,nr
            do i_spacks=is,is
              call spacks_copy('p',sspec(i_spacks)%name,is,is,spid,i_spacks)
            enddo


            ssite(ib)%spec=is


            sspec(is)%lmxa=lmxa

            backspace jfi
            backspace jfi
            backspace jfi
          endif
C         is = -1 -> spec struc does not have these parameters
C         lskip = .false.
          if (is .ne. -1) then
            do i_spacks=is,is
              call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
            enddo


            a=sspec(is)%a
            nr=sspec(is)%nr
            rmt=sspec(is)%rmt
            z=sspec(is)%z


            lmxa=sspec(is)%lmxa
            lmxl=sspec(is)%lmxl
            kmax=sspec(is)%kmxt

C           lmxa = -1 -> site has no augmentation
            if (lmxa .eq. -1) then
C             if (ib .gt. nat0) goto 20
C             lskip = .true.
              goto 20
            endif
          endif
          ibaug = ibaug+1
          if (procid .eq. master) then
            if (lbin) then
              read(jfi) is0,spid0,lmxa0,lmxl0,nr0,rmt0,a0,z0,qc
              if (ib .gt. nbas) goto 20
            else
              read(jfi,'(19x,a8)') spid0
              read(jfi,*)
              read(jfi,*) is0,lmxa0,lmxl0,nr0
              read(jfi,*)
              read(jfi,*) rmt0,a0,z0,qc
              if (ib .gt. nbas) goto 20
            endif
            if (ipr .ge. 40) then
              if (ib .le. nat) write(stdo,380) ib,is0,spid0
              if (ib .gt. nat) write(stdo,380) ib,is0,spid0, ' (skip)'
            endif
  380       format('   atom',i4,'    species',i4,':',a:a)
c     ... read(but don't use) extra info since record is present
            if (lbin) then
              read(jfi) rsma0,rsmr0,rsmv0,lmxv0,lmxr,lmxb0,kmax0
            else
              read(jfi,*)
              read(jfi,*) rsma0,rsmr0,rsmv0,lmxv0,lmxr,lmxb0,kmax0
            endif
          endif
          call mpibc1(lmxa0,1,2,mlog,'iors','lmxa0')
          call mpibc1(lmxl0,1,2,mlog,'iors','lmxl0')
          call mpibc1(nr0,1,2,mlog,'iors','nr0')
          call mpibc1(a0,1,4,mlog,'iors','a0')
          call mpibc1(qc,1,4,mlog,'iors','qc')
C         is = -1 -> pack spec struc with file parameters
          if (is .eq. -1 .or. mode .eq. -2) then
            is = is0
            do i_spacks=is,is
              call spacks_copy('p',sspec(i_spacks)%name,is,is,spid0,i_spacks)
            enddo

            spid = spid0

            sspec(is)%a=a0
            sspec(is)%nr=nr0
            sspec(is)%rmt=rmt0
            sspec(is)%z=z0


            sspec(is)%lmxa=lmxa0
            sspec(is)%lmxl=lmxl0
            sspec(is)%kmxt=kmax0


            a=sspec(is)%a
            nr=sspec(is)%nr
            rmt=sspec(is)%rmt
            z=sspec(is)%z


            lmxa=sspec(is)%lmxa
            lmxl=sspec(is)%lmxl
            kmax=sspec(is)%kmxt

          endif

C     ... Read augmentation data
          call dpzero(pnu,2*n0)
          call dpzero(pnz,2*n0)
          call dpzero(ql,2*n0)
C         if (lmxa0 .lt. lmxa) call defpq(z,lmxa,nsp,pnu,ql)

          i_copy_size=size(sspec(is)%p)
          call dcopy(i_copy_size,sspec(is)%p,1,pnu,1)
          i_copy_size=size(sspec(is)%pz)
          call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)

          call iinit(idmod,n0)
          call iinit(idmoz,n0)
          if (procid .eq. master) then
            do  isp = 1, nsp0
              if (lbin) then
                read(jfi) (pnu(l+1,isp), l=0,lmxa0)
                read(jfi) (pnz(l+1,isp), l=0,lmxa0)
              else
                if (isp .eq. 1) read(jfi,*)
                read(jfi,*) (pnu(l+1,isp), l=0,lmxa0)
                read(jfi,*) (pnz(l+1,isp), l=0,lmxa0)
              endif
              if (nsp .gt. nsp0) then
                do  l = 0, lmxa0
                  pnu(l+1,2) = pnu(l+1,1)
                  pnz(l+1,2) = pnz(l+1,1)
                enddo
              endif
            enddo
          endif
          do  isp = 1, nsp
            call mpibc1(pnu(1,isp),lmxa0+1,4,mlog,'iors','pnu')
            call mpibc1(pnz(1,isp),lmxa0+1,4,mlog,'iors','pnu')
C       ... For backwards compatibility: prior versions wrote pnu for pnz
            do  l = 0, lmxa0+1
              if (pnu(l+1,isp) .eq. mod(pnz(l+1,isp),10d0))
     .        pnz(l+1,isp) = 0
            enddo
          enddo
C         read(but don't use) site idmod for compatibility with nfp
          if (procid .eq. master) then
            if (lbin) then
              read(jfi) (idmod(l+1), l=0,lmxa0)
              read(jfi) (idmoz(l+1), l=0,lmxa0)
            else
              read(jfi,*) (idmod(l+1), l=0,lmxa0)
              read(jfi,*) (idmoz(l+1), l=0,lmxa0)
            endif
          endif

C     ... Store pnu and pz in site and spec data
C         if (ib .le. nat) then

          i_copy_size=size(sspec(is)%p)
          call dcopy(i_copy_size,sspec(is)%p,1,pnus,1)
          i_copy_size=size(sspec(is)%pz)
          call dcopy(i_copy_size,sspec(is)%pz,1,pnzs,1)

          if  (lrs(3) .eq. 0 .or. mode .eq. -2) then

            i_copy_size=size(ssite(ib)%pnu)
            call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
            i_copy_size=size(ssite(ib)%pz)
            call dcopy(i_copy_size,pnz,1,ssite(ib)%pz,1)


            i_copy_size=size(sspec(is)%p)
            call dcopy(i_copy_size,pnu,1,sspec(is)%p,1)
            i_copy_size=size(sspec(is)%pz)
            call dcopy(i_copy_size,pnz,1,sspec(is)%pz,1)

            if ( mode .eq. -2 ) then

              i_copy_size=size(sspec(is)%p)
              call dcopy(i_copy_size,sspec(is)%p,1,pnus,1)
              i_copy_size=size(sspec(is)%pz)
              call dcopy(i_copy_size,sspec(is)%pz,1,pnzs,1)
            endif

C         endif

C     ... Verify lowest valence pnu compatible with file
            lfail = .false.
            ltmp1 = .false.
            ltmp2 = .false.
            do  i = 1, lmxa+1
              vec0(1) = mod(pnz(i,1),10d0)
              if (vec0(1) .eq. 0) vec0(1) = pnu(i,1)
              vec0(2) = mod(pnzs(i,1),10d0)
              if (vec0(2) .eq. 0) vec0(2) = pnus(i,1)
              ltmp1 = ltmp1 .or. int(pnu(i,1)) .ne. int(pnus(i,1))
              ltmp2 = ltmp2 .or.
     .        int(mod(pnz(i,1),10d0)) .ne. int(mod(pnzs(i,1),10d0))
              lfail = lfail .or. min(int(pnu(i,1)),int(vec0(1))) .ne.
     .        min(int(pnus(i,1)),int(vec0(2)))
            enddo
            if (ltmp1 .and. ipr.ge.20)
     .      write(stdo,203) ib,spid,'file pnu',(pnu(i,1), i=1,lmxa+1)
            if (ltmp1 .and. ipr.ge.20)
     .      write(stdo,204) 'given pnu is',(pnus(i,1), i=1,lmxa+1)
            if (ltmp2 .and. ipr.ge.20)
     .      write(stdo,203) ib,spid,'file pz ',(pnz(i,1), i=1,lmxa+1)
            if (ltmp2 .and. ipr.ge.20)
     .      write(stdo,204) 'given pz  is',(pnzs(i,1), i=1,lmxa+1)
  203       format(9x,'site',i4,':',a,':',a,' is',8f6.2)
  204       format(26x,a,8f6.2)
            if (lfail .and. lrs(3) .eq. 1) then
              call rx('iors: file''s pnu is incompatible with input')
            endif

            nlml0 = (lmxl0+1)**2
            nlml = (lmxl+1)**2
            if (nr .le. 0)   nr = nr0
            if (a .le. 1d-6) a = a0

            if (procid .eq. master) then
C     ... Sanity checks, or inform about changed parameters
C         if (is0 .ne. is) call xxerri('species pointer',is,is0)
              call fsanrg(rmt0,rmt,rmt,1d-3,msg,'rmt',.true.)
              call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.false.)
              call fsanrg(z0,z,z,1d-6,msg,'z',.true.)
              call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
              lfail = isanrg(nr0,nr,nr,msgw,'nr',.false.)
              if (isanrg(lmxl,  0,lmxa,  msg,'lmxl', .false.)) goto 999

#if FP
              if (kmax0.ne.kmax.and.ipr.ge.10)
     .        write(stdo,201) ib,spid,'kmax',kmax0,kmax
#endif
              if (lmxa0.ne.lmxa.and.ipr.ge.10)
     .        write(stdo,201) ib,spid,'lmax',lmxa0,lmxa
  201         format(9x,'site',i4,', species ',a,
     .        ': augmentation ',a,' changed from',i2,' to',i2)
            endif

C         Case read but skip over this site data
          else
            lmxl = lmxl0
            nr = nr0
          endif

C     --- Allocate and read arrays for local density and potential ---
          nlml0 = (lmxl0+1)**2
          nlml = (lmxl+1)**2

          call defrr(ov0,           nr*nsp)
C     ... FP local densities rho1,rho2,rhoc and potentials v0, v1
          if (mode .eq. 1 .or. mode .eq. -2) then
            if (nr .ne. nr0)
     .      call rx('iors not set up to convert radial mesh')
            call defrr(orhoat(1,ib),  nr*nlml*nsp)
            call defrr(orhoat(2,ib),  nr*nlml*nsp)
            call defrr(orhoat(3,ib),  nr*nsp)
            call defrr(ov1,           nr*nsp)

            if (procid .eq. master) then
              if (.not. lbin) read(jfi,*)
              call dpdbyl(w(orhoat(1,ib)),nr0,nlml0,nlml,nsp0,nsp,lbin,
     .        jfi)
              if (.not. lbin) read(jfi,*)
              call dpdbyl(w(orhoat(2,ib)),nr0,nlml0,nlml,nsp0,nsp,lbin,
     .        jfi)
              if (nlml0.gt.nlml.and.ipr.ge.10)
     .        write(stdo,202) ib,spid,'truncate',nlml0,nlml
              if (nlml0.lt.nlml.and.ipr.ge.10)
     .        write(stdo,202) ib,spid,'inflate',nlml0,nlml
  202         format(9x,'site',i4,', species ',a,
     .        ': ',a,' local density from nlm=',i3,' to',i3)
              if (.not. lbin) read(jfi,*)
              call dpdbyl(w(orhoat(3,ib)),nr0,1,1,nsp0,nsp,lbin,jfi)
              if (.not. lbin) read(jfi,*)
              call dpdbyl(w(ov0),nr0,1,1,nsp0,nsp,lbin,jfi)
              if (.not. lbin) read(jfi,*)
              call dpdbyl(w(ov1),nr0,1,1,nsp0,nsp,lbin,jfi)
              if (nsp0 .lt. nsp) then
                call dscal(nr0*2,2d0,w(ov0),1)
                call dscal(nr0*2,2d0,w(ov1),1)
              endif
C           call prmx('iors: v0',w(ov0),nr0,nr0,nsp)
            endif
C           if (ib .gt. nat) then
C              call rlse(ov0)
C              goto 20
C           endif

            if (mode .eq. -2) then
              call rlse(orhoat(1,ib))
              call rlse(ov0)

            else
              call mpibc1(w(orhoat(1,ib)),nr*nlml*nsp,4,mlog,'iors',
     .        'rhoat(1)')
              call mpibc1(w(orhoat(2,ib)),nr*nlml*nsp,4,mlog,'iors',
     .        'rhoat(2)')
              call mpibc1(w(orhoat(3,ib)),nr*nlml*nsp,4,mlog,'iors',
     .        'rhoat(3)')
              call mpibc1(w(ov0),nr*nsp,4,mlog,'iors','v0')
              call mpibc1(w(ov1),nr*nsp,4,mlog,'iors','v1')
            endif
C     ... ASA sphere densities
          else
            if (lbin) call rx('iors not ready for binary ASA input')
C           call dpdbyl(w(ov0),nr0,1,1,nsp0,nsp,lbin,jfi)
            call dpzero(ql,3*n0*nsp)
            if (procid .eq. master) then
              call aiomom(spid0,pnu,qnu(1,1,1+(ib-1)*nsp),
     .        idmod,n0,lmxa,nsp,rh,vrmax,jfi)
              if (vs1 .lt. -1.041d0) then
                read(jfi,43) vrmax(1)
   43           format(12x,f11.7)

Cchp1                 oves=spot%oves
                 rv_p_oves => spot%rv_p_oves 


Cchp1                 call dvset(w(oves),ic,ic,vrmax)
                 call dvset ( rv_p_oves , ic , ic , vrmax ) 


Cchp1                 spot%oves=oves
                 spot%rv_p_oves => rv_p_oves 


              endif
              if ( lrs ( 3 ) .eq.0 ) then

                i_copy_size=size(ssite(ib)%pnu)
                call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
              endif

            endif
          endif
C     ... store data in strucs

          sspec(is)%a=a
          sspec(is)%nr=nr
          sspec(is)%rmt=rmt
          sspec(is)%z=z


          sspec(is)%lmxa=lmxa
          sspec(is)%lmxl=lmxl
          sspec(is)%kmxt=kmax
          sspec(is)%qc=qc


          ssite(ib)%ov0=ov0
          ssite(ib)%ov1=ov1

   20     continue
        enddo
        if (isanrg(ibaug, nat,nat,  msg,'nat', .false.)) goto 999

#if FP
        if ( lshear .and. iand(256,int(sctrl%lrs)) .ne.0 ) then


          i_copy_size=size(slat%dist)
          call dcopy(i_copy_size,slat%dist,1,wk,1)

C         call prmx('wk',wk,3,3,3)
          call dgemm('N','T',3,1,3,1d0,plat,3,qlat0,3,0d0,wk,3)
C         call prmx('wk',wk,3,3,3)
          call pvsms2(ssite,sspec,wk,nbas,nsp,orhoat)
        endif
#endif

C   --- Read data on free-atom core states and fit to fa density ---
        line = 'species data'
        if (procid .eq. master) then
          if (.not. lbin) read(jfi,*,err=999,end=999)
        endif
        do  30  is = 1, nspec

          a=sspec(is)%a
          nr=sspec(is)%nr
          lmxa=sspec(is)%lmxa

          if (lmxa .eq. -1) goto 30
          if (procid .eq. master) then
            if (lbin) then
              read(jfi,err=999,end=999) nr0,a0,qc,cof,eh,stc,lfoc0,rfoc0
            else
              read(jfi,*,err=999,end=999)
              read(jfi,*,err=999,end=999)
              read(jfi,*,err=999,end=999) nr0,lfoc0
              read(jfi,*,err=999,end=999)
              read(jfi,*,err=999,end=999) a0,qc,cof,eh
              read(jfi,*,err=999,end=999)
              read(jfi,*,err=999,end=999) stc,rfoc0
            endif

            sspec(is)%lfoca=lfoc0
            sspec(is)%rfoca=rfoc0

            lfail = isanrg(nr0,nr,nr,msgw,'nr',.false.)
            call fsanrg(a0,a,a,0d-9,msg,'spec a',.true.)
          endif
          call mpibc1(qc,1,4,mlog,'iors','qc')
          call mpibc1(cof,1,4,mlog,'iors','cof')
          call mpibc1(eh,1,4,mlog,'iors','eh')
          call mpibc1(stc,1,4,mlog,'iors','stc')

C     ... FP core densities
          if (mode .eq. 1 .or. mode .eq. -2) then
            call defrr(orhoca, nr*nsp)
            if (procid .eq. master) then
              if (nr .ne. nr0)
     .        call rx('iors not set up to convert core radial mesh')
              if (lbin) then
                call dpdump(w(orhoca),nr*nsp0,jfi)
              else
                call dfdump(w(orhoca),nr*nsp0,jfi)
              endif
C           spin-split core density
              if (nsp .gt. nsp0) then
                i = nr
                call dscal(i,0.5d0,w(orhoca),1)
                call dpscop(w(orhoca),w(orhoca),i,1,1+i,1d0)
              endif
            endif
            if (mode .eq. 1) then
              call mpibc1(w(orhoca),nr*nsp,4,mlog,'iors','rhoca')
            else
              call rlse(orhoca)

              do i_spacks=is,is
                call spacks_copy('p',sspec(i_spacks)%coreh,is,is,' ',i_spacks)
              enddo

              do i_spacks=is,is
                call spacks_copy('p',sspec(i_spacks)%pb1,is,is,' ',i_spacks)
              enddo

              do i_spacks=is,is
                call spacks_copy('p',sspec(i_spacks)%pb2,is,is,' ',i_spacks)
              enddo

            endif
            call dpzero(exi,n0)
            call dpzero(hfc,n0*2)
            if (procid .eq. master) then
              if (lbin) then
                read(jfi,err=999,end=999) rsmfa,nxi
                read(jfi,err=999,end=999)
     .          ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp0)
              else
                read(jfi,*,err=999,end=999)
                read(jfi,*,err=999,end=999) rsmfa,nxi
                read(jfi,*,err=999,end=999) (exi(i),i=1,nxi)
                call dfdump(hfc(1,1),nxi,jfi)
                if (nsp0 .eq. 2) call dfdump(hfc(1,2),nxi,jfi)
              endif
              if (nsp .gt. nsp0) then
                i = n0
                call dscal(i,0.5d0,hfc,1)
                call dpscop(hfc,hfc,i,1,1+i,1d0)
              endif
            endif
            call mpibc1(rsmfa,1,4,mlog,'iors','rsmfa')
            call mpibc1(nxi,1,2,mlog,'iors','nxi')
            call mpibc1(exi,nxi,4,mlog,'iors','exi')
            call mpibc1(hfc,nsp*nxi,4,mlog,'iors','hfc')

C     ... For now, mode 0 (ASA) stores no core data
          elseif (mode .eq. 0) then
            orhoca = 0
            rsmfa = 0
            nxi = 0
            call dpzero(exi,n0)
            call dpzero(hfc,2*n0)
          endif
C         call p_spec_hfc(sspec,is,rsmfa,nxi,exi,hfc)
C         Pack species data

          sspec(is)%ctail=cof
          sspec(is)%etail=eh
          sspec(is)%stc=stc
          sspec(is)%orhoc=orhoca


          sspec(is)%nxi=nxi
          i_copy_size=size(sspec(is)%exi)
          call dcopy(i_copy_size,exi,1,sspec(is)%exi,1)
          i_copy_size=size(sspec(is)%chfa)
          call dcopy(i_copy_size,hfc,1,sspec(is)%chfa,1)
          sspec(is)%rsmfa=rsmfa

   30   continue

c   ... Copy or rescale cores, in case foca was switched on or off
        if (mode .eq. 1) then
          do  ib = 1, nbas
            is = int(ssite(ib)%spec)


            a=sspec(is)%a
            nr=sspec(is)%nr
            rmt=sspec(is)%rmt
            lmxa=sspec(is)%lmxa


            lfoc=sspec(is)%lfoca
            qc=sspec(is)%qc
            orhoca=sspec(is)%orhoc

            if (lmxa .eq. -1) goto 40
            if (lfoc .gt. 0) then
              call dpcopy(w(orhoca),w(orhoat(3,ib)),1,nr*nsp,1d0)
            else
              allocate(rwgt_rv(nr))
Ctakao_ZeroClear_NotRequiered              if (nr<0) rwgt_rv(:)=0.0d0

              call radwgt ( rmt , a , nr , rwgt_rv )

              call radsum ( nr , nr , 1 , nsp , rwgt_rv , w ( orhoat ( 3 
     .        , ib ) ) , sum )

              fac = 1d0
              if (dabs(sum) .gt. 1d-6) fac = qc/sum
              if (dabs(fac-1d0).gt.1d-7 .and. ipr.ge.30)
     .        write(stdo,787) ib,qc,sum,fac
  787         format(' fix core chg: ib=',i4,'  qc,sum,fac=',3f12.6)
              call dpcopy(w(orhoat(3,ib)),w(orhoat(3,ib)),1,nr*nsp,fac)
              if (allocated(rwgt_rv)) deallocate(rwgt_rv)

            endif
   40       continue
          enddo
        endif

C   ... MPI broadcast all structures
C        if (procid .eq. master) then
C          sttime = MPI_WTIME()
C        endif
cki        call mpibc1(ssite,nint(ssite(1)%size)*nbas,4,mlog,'iors','ssite')
        do i_site=1,nbas
        call mpibc1_s_site(ssite(i_site),mlog,'iors','ssite')
        enddo
cki        call mpibc1(sspec,nint(sspec(1)%size)*nspec,4,mlog,'iors','sspec')
        do i_spec=1,nspec
        call mpibc1_s_spec(sspec(i_spec),mlog,'iors','sspec')
        enddo
cki        call mpibc1(sctrl,nint(sctrl%size),4,mlog,'iors','sctrl')
        call mpibc1_s_ctrl(sctrl,mlog,'iors','sctrl')
cki        call mpibc1(slat,nint(slat%size),4,mlog,'iors','slat')
        call mpibc1_s_lat(slat,mlog,'iors','slat')
cki        call mpibc1(spot,nint(spot%size),4,mlog,'iors','spot')
        call mpibc1_s_pot(spot,mlog,'iors','spot')
cki        call mpibc1(sbz,nint(sbz%size),4,mlog,'iors','sbz')
        call mpibc1_s_bz(sbz,mlog,'iors','sbz')
        if (mode .eq. 1) then
        elseif (mode .eq. 0) then
          call mpibc1(qnu,3*n0*nsp*nbas,4,mlog,'iors','qnu')
        endif
C        if (procid .eq. master) then
C          entime = MPI_WTIME()
C          call info2(20,0,0,' iors  : MPI broadcast took %;3d sec',
C     .      (entime-sttime),0)
C        endif

C        call shstru('spec',sspec,1,nspec)
C        call shstru('site',ssite,1,nbas)
C        call rx0('done')


C --- Output ---
      else
        if (procid .ne. master) then
          iors = 0
          return
        endif
        jfi = -ifi
        rewind jfi
C ...   Write header and general information
        fid0 = fid
        call strip(fid0,i0,i1)
        call fextg(datimp)
        jobid = datimp(2:)
        call ftime(datimp)
        hostnm = ' '
        usernm = ' '
        call gtenv('HOST',hostnm)
        call gtenv('USER',usernm)
        call strip(datimp,i,i2)
        call strip(usernm,i,i3)
        call strip(hostnm,i,i4)
        if (ipr .ge. 40) write(stdo,710) fid(1:i1),
     .  usernm(1:i3),hostnm(1:i4),datimp(1:i2)

  721   format('----------------------- ',a,' -----------------------')
        if (lbin) then
          write(jfi) vs
          write(jfi) fid0
          write(jfi) datimp,usernm,hostnm,jobid
          if (abs(vs) .le. 1.021) then
            write(jfi) nbas,nsp,npan,lrel
          else if (abs(vs) .le. 1.031) then
            write(jfi) nbas,nsp,npan,lrel,nspec
          else
            write(jfi) nbas,nat,nsp,npan,lrel,nspec
          endif
          write(jfi) nit
          write(jfi) alat,vol,plat
        else
          write(jfi,721) 'header'
          write(jfi,'(1f15.10,t50,a)') vs,
     .    '! version; id; date; name; host'
          write(jfi,'(1x,a68)') fid0
          write(jfi,'(1x,a32)') datimp
          write(jfi,'(1x,a32)') usernm
          write(jfi,'(1x,a32)') hostnm
C         write(jfi,'(''  nbas nsp  npan lrel nit'')')
          if (abs(vs) .le. 1.021) then
            write(jfi,'(5i5,t50,a)') nbas,nsp,npan,lrel,nit,
     .      '! nbas nsp npan lrel nit'
          else if (abs(vs) .le. 1.031) then
            write(jfi,'(6i5,t50,a)') nbas,nsp,npan,lrel,nit,nspec,
     .      '! nbas nsp npan lrel nit nspec'
          else
            write(jfi,'(7i5,t50,a)') nbas,nat,nsp,npan,lrel,nit,nspec,
     .      '! nbas nat nsp npan lrel nit nspec'
          endif
          write(jfi,'(f15.10,1pe20.10,t50,a)') alat,vol,
     .    '! alat vol; plat'
          write(jfi,'(3f15.10)') plat
        endif

C   --- Write smooth charge density ---
#if FP
        if (lbin) then
          write(jfi) n1,n2,n3
        else
          write(jfi,721) 'smooth density'
          write(jfi,ifmt) n1,n2,n3
        endif
        call dpdftr(n1,n2,n3,k1,k2,k3,nsp,w(osmrho),lbin,-jfi)
#endif

C   --- Write information related to dynamics ---
        call dpzero(wk,100)

        wk(1)=sbz%ef
        wk(2)=sbz%def
        wk(3)=sbz%w

        if (lbin) then
          call dpdump(wk,100,-jfi)
        else
          write(jfi,721) 'geometric data'
          call dfdump(wk,100,-jfi)
        endif
        do  110  ib = 1, nbas

          i_copy_size=size(ssite(ib)%pos)
          call dcopy(i_copy_size,ssite(ib)%pos,1,pos,1)
          i_copy_size=size(ssite(ib)%force)
          call dcopy(i_copy_size,ssite(ib)%force,1,force,1)
          i_copy_size=size(ssite(ib)%vel)
          call dcopy(i_copy_size,ssite(ib)%vel,1,vel,1)

          if (lbin) then
            write(jfi) ib,pos,force,vel
          else
            write(jfi,'('' Site:'',i5)') ib
            write(jfi,345) pos,   'position'
            write(jfi,345) force, 'force'
            write(jfi,345) vel,   'velocity'
  345       format(3f15.10,t50,'! ',a)
          endif
  110   continue

C   --- Write information for local densities ---
        if (ipr .ge. 50) write(stdo,364)
        if (.not. lbin) write(jfi,721) 'site densities'
        do  120  ib = 1, nbas


          ic=ssite(ib)%class
          is=ssite(ib)%spec

          do i_spacks=is,is
            call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
          enddo


          a=sspec(is)%a
          nr=sspec(is)%nr
          rmt=sspec(is)%rmt


          z=sspec(is)%z
          qc=sspec(is)%qc
          i_copy_size=size(sspec(is)%idmod)
          call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1)


          rsma=sspec(is)%rsma
          lmxa=sspec(is)%lmxa
          lmxl=sspec(is)%lmxl
          lmxb=sspec(is)%lmxb


          kmxv=sspec(is)%kmxv
          rsmv=sspec(is)%rsmv
          kmax=sspec(is)%kmxt


          i_copy_size=size(ssite(ib)%pnu)
          call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
          i_copy_size=size(ssite(ib)%pz)
          call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
          ov0=ssite(ib)%ov0
          ov1=ssite(ib)%ov1

          if (lmxa .eq. -1) goto 120

          if (lbin) then
            write(jfi) is,spid,lmxa,lmxl,nr,rmt,a,z,qc
          else
            write(jfi,'('' site'',i5,''   spec: '',a8)') ib, spid
            write(jfi,'('' spec  lmxa lmxl  nr'')')
            write(jfi,ifmt) is,lmxa,lmxl,nr
            write(jfi,'(t7,'' rmt'',t21,''a'',t36,''z'',t51,''qc'')')
            write(jfi,ffmt) rmt,a,z,qc
          endif
C     ... Some extra info... lots of it useless or obsolete
          lmxr = 0
          lmxv = 0
          rsmr = 0
          if (lbin) then
            write(jfi) rsma,rsmr,rsmv,lmxv,lmxr,lmxb,kmax
          else
            write(jfi,'('' Next two lines for nfp compatibility'')')
            write(jfi,ffmt) rsma,rsmr,rsmv
            write(jfi,ifmt) lmxv,lmxr,lmxb,kmax
          endif

C     ... Write augmentation data
          do  122  isp = 1, nsp
            if (lbin) then
              write(jfi) (pnu(l+1,isp), l=0,lmxa)
              write(jfi) (pnz(l+1,isp), l=0,lmxa)
            else
              if (isp .eq. 1) write(jfi,'('' pnu, pz, idmod, idmoz'')')
              write(jfi,ffmt) (pnu(l+1,isp), l=0,lmxa)
              write(jfi,ffmt) (pnz(l+1,isp), l=0,lmxa)
            endif
  122     continue
C         Write for compatibility with nfp
          if (lbin) then
            write(jfi) (idmod(l+1), l=0,lmxa)
            write(jfi) (idmod(l+1), l=0,lmxa)
          else
            write(jfi,ifmt) (idmod(l+1), l=0,lmxa)
            write(jfi,ifmt) (idmod(l+1), l=0,lmxa)
          endif

C     ... Write arrays for local density and potential
          nlml = (lmxl+1)**2
          if (mode .eq. 1) then
            if (.not. lbin) write(jfi,'('' rho1'')')
            call dpdbyl(w(orhoat(1,ib)),nr,nlml,nlml,nsp,nsp,lbin,-jfi)
            if (.not. lbin) write(jfi,'('' rho2'')')
            call dpdbyl(w(orhoat(2,ib)),nr,nlml,nlml,nsp,nsp,lbin,-jfi)
            if (.not. lbin) write(jfi,'('' rhoc'')')
            call dpdbyl(w(orhoat(3,ib)),nr,1,1,nsp,nsp,lbin,-jfi)
            if (.not. lbin) write(jfi,'('' v0'')')
            call dpdbyl(w(ov0),nr,1,1,nsp,nsp,lbin,-jfi)
            if (.not. lbin) write(jfi,'('' v1'')')
            call dpdbyl(w(ov1),nr,1,1,nsp,nsp,lbin,-jfi)
          else
            if (lbin) then
              call rx('iors not ready for this mode')
C             call dpdump(pl,nl*nsp,-jfi)
C             call dpdump(w(orhoat(1,ib)),nl*nsp*3,-jfi)
C             call dpdbyl(w(ov0),nr,1,1,nsp,nsp,lbin,-jfi)
            else
              rh = 0
              vrmax(1) = 0
              vrmax(2) = 0
              call aiomom(spid,pnu,qnu(1,1,1+(ib-1)*nsp),
     .        idmod,n0,lmxa,nsp,rh,vrmax,-jfi)

Cchp1               oves=spot%oves
               rv_p_oves => spot%rv_p_oves 


Cchp1               vrmax(1) = dval(w(oves),ic)
               vrmax ( 1 ) = dval ( rv_p_oves , ic ) 

              write(jfi,143) vrmax(1)
  143         format('  Ves(rmax)=',f11.7)
C             Not yet; for future versions
C             call aiopot(nr,nsp,a,rmt,-99d0,w(ov0),-jfi)
            endif
          endif

          if (ipr .ge. 50) then
            write(stdo,349) ib,spid,lmxa,lmxl,rmt,nr,a,
     .      (pnu(l+1,1),l=0,lmxa)
            if (nsp .eq. 2)  write(stdo,350) (pnu(l+1,2), l=0,lmxa)
C           if (npan .eq. 2) write(stdo,350) (pz(l+1,1), l=0,lmxa)
C            if (npan .eq. 2 .and. nsp .eq. 2)
C     .        write(stdo,350) (pz(l+1,2), l=0,lmxa)
          endif
  120   continue

C   --- Write data on free-atom core states and fit to fa density ---
        if (.not. lbin) write(jfi,721) 'species data'

        do  130  is = 1, nspec

          a=sspec(is)%a
          nr=sspec(is)%nr
          qc=sspec(is)%qc
          lmxa=sspec(is)%lmxa


          lfoc=sspec(is)%lfoca
          rfoc=sspec(is)%rfoca


          cof=sspec(is)%ctail
          eh=sspec(is)%etail
          stc=sspec(is)%stc
          orhoca=sspec(is)%orhoc


          nxi=sspec(is)%nxi
          i_copy_size=size(sspec(is)%exi)
          call dcopy(i_copy_size,sspec(is)%exi,1,exi,1)
          i_copy_size=size(sspec(is)%chfa)
          call dcopy(i_copy_size,sspec(is)%chfa,1,hfc,1)
          rsmfa=sspec(is)%rsmfa

          if (lmxa .eq. -1) goto 130

          if (lbin) then
            write(jfi) nr,a,qc,cof,eh,stc,lfoc,rfoc
          else
            do i_spacks=is,is
              call spacks_copy('u',sspec(i_spacks)%name,is,is,spid0,i_spacks)
            enddo

            write(jfi,'('' spec'',i5,2x,'' label '',a)') is, spid0
            write(jfi,'(''  nr   lfoc'')')
            write(jfi,ifmt) nr,lfoc
            write(jfi,'(t7,''  a '',t21,''qc'',t36,''cof'',t61,''eh'')')
            write(jfi,'(2f15.10,1p,2e25.15)') a,qc,cof,eh
            write(jfi,'(t17,'' stc'',t31,''rfoc'')')
            write(jfi,'(1p,2e25.15)') stc,rfoc
          endif
C     ... For now, ASA stores no core data
          if (mode .eq. 1) then
            if (lbin) then
              call dpdump(w(orhoca),nr*nsp,-jfi)
              write(jfi) rsmfa,nxi
              write(jfi) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp)
            else
              call dfdump(w(orhoca),nr*nsp,-jfi)
              write(jfi,'('' core tail: rsmfa,n,e,cof=''/f15.10,i4)')
     .        rsmfa,nxi
              write(jfi,ffmt) (exi(i),i=1,nxi)
              call dfdump(hfc(1,1),nxi,-jfi)
              if (nsp .eq. 2) call dfdump(hfc(1,2),nxi,-jfi)
            endif
          endif
  130   continue

      endif

  349 format(i11,':',a4,2i2,f9.5,i5,f6.3,1x,8f6.3)
  350 format(41x,8f6.3)
  364 format(/9x,'ib:spc la ll   rmt     nr   a     pnu')

      iors = 0
      return

C ... Error handling
  998 if (ipr .gt. 0)
     .write(stdo,'('' iors  : empty file ... nothing read'')')
      return
  999 continue
      if (ipr .gt. 0)
     .write(stdo,'('' iors  : read failed in: '',a)') line

      end function iors


