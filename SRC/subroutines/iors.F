      module m_iors
      use m_struc_def
      public iors
      type(s_rv1),public,allocatable :: v1pot(:),v0pot(:)
      private
      
      contains
      integer function iors(nit,rwrw,irs3,irs5)
      use m_density,only: osmrho, orhoat !these are allocated
      use m_bndfp,only: m_bndfp_ef_SET,eferm 

      use m_supot,only: lat_nabc
      use m_struc_func,only: mpibc1_s_spec,mpibc1_s_site
      use m_lmfinit,only: lat_alat,nsp,lrel,nl,ssite=>v_ssite,sspec=>v_sspec, nbas,nat,nspec,n0  
      use m_lattic,only: lat_plat 
      use m_ext,only:sname
!!- I/O for charge density to rst or rsta. ssite sspec are readin 
!! read write
!!     smrho, rhoat
!!     ssite: pos, pos0, force, pnu pz ov0,ov1 
!      sspec:
!!           a nr rmt z lmxa lmxl kmxt p pz lfoca qc orhoc idmod
!!           rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
!!           chfa rsmfa
!      
!! ef,def --> goto bndfp via m_bndfp_ef_set      
!!      
Cixxx   fid   :string containing identification
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lbin  :=T file I/O in binary mode
Ci         :xxxF file I/O in ascii mode
!!      
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
Cuxxx   10 Dec 02 File shears atom positions by shear transformation
Cuxxx             (slat->plat) (file plat)^-1
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
      implicit none
      logical:: lbin=.true.
      integer::  irs3,irs5,nit , ifi , i_site,i_spec!,i_copy_size !mode=1 ,
      character*256:: fid=''
C ... Local parameters
      integer procid,master,mpipid,nproc
      integer i,i0,i1,i2,i3,i4,ib,ipr,iprint,ic,is,is0,isp,jb,k1,k2,k3,
     .igetss,jfi,k11,k21,k31,kmax,kmax0,kmxv,l,lfoc,lfoc0,lmxa,
     .lmxa0,lmxb,lmxb0,lmxl,lmxl0,lmxr,lmxv,lmxv0,lrel0,n11,n21,
     .n31,nbas0,nspec0,nlml,nlml0,npan,npan0,nr,nr0,nsp0,
     .nxi,nat0,ibaug
      integer:: ngabc(3) , n1 , n2 , n3 , lgunit , stdo , isw
      complex(8) ,allocatable :: h_zv(:)
      real(8) ,allocatable :: rwgt_rv(:)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer idmod(n0),idmoz(n0) !,lrs(10)
      logical isanrg,lfail,ltmp1,ltmp2,latvec,cmdopt,mlog !,lshear
      double precision a,a0,alat,alat0,cof,eh,fac,qc,rfoc,rfoc0,rmt,
     .rmt0,rsma,rsma0,rsmfa,rsmr,rsmr0,rsmv,rsmv0,stc,sum,vfac,vol,vol0,vs,vs1,z,z0
      double precision pnu(n0,2),pnz(n0,2),ql(n0,2*n0),pos(3),
     .pos0(3),force(3),plat(3,3),plat0(3,3),qlat(3,3),qlat0(3,3), !,vel(3)
     .exi(n0),hfc(n0,2),vec0(3),wk(100),rh,vrmax(2),pnus(n0,2),
     .pnzs(n0,2),dval
      character spid*8,spid0*8,fid0*68,line*20,msg*23,use*80,ignore*80,
     .msgw*17,datimp*32,usernm*32,hostnm*32,jobid*32,ffmt*32,ifmt*32
      integer:: fextg, i_dummy_fextg,ifile_handle
      character(*)::rwrw
      data vec0 /0d0,0d0,0d0/
c     call tcn('iors')
      open(newunit=ifi,file='rst.'//trim(sname),form='unformatted') 
      
C ... MPI setup
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,ignore)
      alat=lat_alat
      plat=lat_plat
      ngabc=lat_nabc
      call dinv33(plat,1,qlat,fac)
      vol = dabs(fac)*alat**3
      ipr    = iprint()
      stdo   = lgunit(1)
      vs   =  2d0 !1.04d0 !version of rst file vs=2d0 at 2022-5-15. Only we support reading vs=1.04 only.
      msg  = '         File mismatch:'
      msgw = '         warning:'
      iors = -1
      line = 'header'
      call fftz30(n1,n2,n3,k1,k2,k3)
      ffmt = '(5f15.10)'
      ifmt = '(20i5)'
      npan = 1 !Hardwired for now
      write(stdo,"(/a)")' iors  : '//trim(rwrw)//' rst restart file (binary mesh density)'
C --- Input ---
      if (trim(rwrw)=='read') then
        jfi = ifi
        use    = '         use from  restart file:'
        ignore = '         ignore in restart file:'
        line = 'header'
! MPI check to see if at least 1st record can be read. Abort with error message if file is missing (lfail = .true.)
        lfail = .false.
c        if (nproc .gt. 0) then
        if (procid .eq. master) then
           lfail = .true.
           read(jfi,end=996,err=996) vs1
           lfail = .false.
           rewind jfi
 996       continue
        endif
        call mpibc1_logical(lfail,1,'iors_read error')
        if (lfail) goto 998
c        endif
        if (procid .eq. master) then
            read(jfi,end=998,err=998) vs1
            read(jfi) fid0
            read(jfi) datimp,usernm,hostnm
c            if (abs(vs1) .le. 1.021) then
c              nspec0 = nspec
c              read(jfi) nbas0,nsp0,npan0,lrel0
c              nat0 = nbas0
c            else if (abs(vs1) .le. 1.031) then
c              read(jfi) nbas0,nsp0,npan0,lrel0,nspec0
c              nat0 = nbas0
c            else
              read(jfi) nbas0,nat0,nsp0,npan0,lrel0,nspec0
c            endif
            read(jfi) nit
            read(jfi) alat0,vol0,plat0
          call dinv33(plat0,1,qlat0,fac)
c          lshear = .not. latvec(3,1d-6,qlat,plat0)
          fid = fid0
          call strip(fid,i,i1)
          call strip(datimp,i,i2)
          call strip(usernm,i,i3)
          call strip(hostnm,i,i4)
          if (ipr .ge. 40) write(stdo,710) fid(1:i1),usernm(1:i3),hostnm(1:i4),datimp(1:i2)
  710     format(9x,'id -  ',a,/9x,'written by -  ',a,' on ',a,' at: ',a)
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
            write(stdo,*)'   (warning) rst mismatch in nat ... skipping sites'
          endif
          if (nsp0 .lt. nsp) write(stdo,*)'   (warning) rst file not spin pol .. splitting spins'
          if (isanrg(npan0, npan,npan,msg,'npan',.true.)) goto 999
          lfail = isanrg(lrel0, lrel,lrel,msgw,'lrel',.false.)
          call fsanrg(abs(vs1),1.00d0,abs(vs),0d0,' ','file''s version',.true.)
        endif
        call mpibc1_int(nbas0,1,'iors_nbas0')
        call mpibc1_int(nat0,1,'iors_nat0')
        call mpibc1_real(plat,9,'iors_plat')
        call mpibc1_int(nit,1,'iors_nit')
C   --- Read smooth charge density ---
        allocate(osmrho(k1*k2*k3*nsp))
        osmrho(:)=0d0
        line = 'smoothed density'
          if (procid .eq. master) then
              read(jfi,err=999,end=999) n11,n21,n31
            if (n11.eq.n1 .and. n21.eq.n2 .and. n31.eq.n3) then
              read(jfi)osmrho ! call dpdftr ( n1 , n2 , n3 , k1 , k2 , k3 , nsp0 , osmrho, lbin , jfi )
              if (nsp .gt. nsp0) then
                i = k1*k2*k3*2
                call dscal ( i , 0.5d0 , osmrho , 1 )
                call dpscop ( osmrho , osmrho , i , 1 , 1 + i , 1d0 )
              endif
            else !... or read and remesh
              if (ipr .ge. 10) write(stdo,450) n11,n21,n31,n1,n2,n3
  450         format(9x,'remesh density from  ',i4,'  *',i4,'  *',i4,
     .        '    to  ',i4,'  *',i4,'  *',i4)
              call fftz30(n11,n21,n31,k11,k21,k31)
              allocate(h_zv(k11*k21*k31*nsp))
              read(jfi)h_zv 
              if (nsp .gt. nsp0) then
                i = k11*k21*k31*2
                call dscal ( i , 0.5d0 , h_zv , 1 )
                call dpscop ( h_zv , h_zv , i , 1 , 1 + i , 1d0 )
              endif
              call pshpr(50)
              i = 0
              if (n1 .eq. 2*n11 .and. n2 .eq. 2*n21 .and. n3 .eq. 2*n31) i=3
              call chgmsh ( i , plat , nsp , n11 , n21 , n31 , k11 , k21 ,
     .         k31 , h_zv , n1 , n2 , n3 , k1 , k2 , k3 , osmrho )
              call poppr
              if (allocated(h_zv)) deallocate(h_zv)
            endif
C   ... If cell volume changed, scale smooth density to maintain charge
            vfac = vol0/vol
            if (dabs(vfac-1d0) .gt. 1d-8) then
              if (ipr .ge. 10) write(stdo,460) vol0,vfac
  460         format(9x,'volume changed from',f8.2,' :  scale smooth density by',f8.4)
              call dpcopy ( osmrho, osmrho, 1 , 2 * k1 * k2 * k3* nsp , vfac )
            endif
          endif
          call mpibc1_complex(osmrho, size(osmrho), 'iors_smrho' )
  115     continue
C   --- Read information related to dynamics ---
C       For compatibility with nfp, read record into wk
        if (procid .eq. master) then
            read(jfi) wk !call dpdump(wk,100,jfi)
        endif
        call mpibc1_real(wk,1,'iors:eferm')
c        if (irs4 .ne. 0) then
c          ignore=trim(ignore)//'ef window,' !binit=T in bndfp==> set eferm=0d0
c        else
        use=trim(use)//'use window,'
        call m_bndfp_ef_SET(wk(1)) !bz_ef00) !,bz_def00)
c        endif
C   --- Read atomic positions,forces,velocities ---
        line = 'site data'
        if (irs3 .ne. 0) then
          ignore=trim(ignore)//'positions,'
        else
          use=trim(use)//' positions,'
c          if (lshear) use=trim(use)//' (sheared)'
C         Must be enough sites to read from rst file
          if (nbas0 .lt. nbas) then
            call info0(1,0,0,'%9f oops ... cannot use file site positions (site mismatch)')
            if (isanrg(nbas0, nbas,nbas,msg,'nbas',.true.)) goto 999
          endif
        endif
        do  ib = 1, nbas0
          if (procid .eq. master) then
              read(jfi,err=999,end=999) jb,pos,force!,vel
          endif
          call mpibc1_real(pos,3,'iors_pos')
          call mpibc1_real(force,3,'iors_force')
c          call mpibc1_real(vel,3,'iors_vel')
          if (ib .gt. nbas) goto 10
C         rst file positions in pos0
c          if (lshear) then
c            call dgemm('T','N',3,1,3,1d0,qlat0,3,pos,3,0d0,pos0,3)
c           call dgemm('N','N',3,1,3,1d0,plat,3,pos0,3,0d0,pos,3)
c          endif
          pos0=pos 
          if (irs3 .ne. 0) then !overwrite pos,vel
            pos=ssite(ib)%pos
c            vel=0d0 
          endif
           ssite(ib)%pos  =pos
           ssite(ib)%pos0 =pos0
           ssite(ib)%force=force
c           ssite(ib)%vel  =vel
   10     continue
        enddo
C   --- Read information for local densities ---
        if (irs5 .ne. 0) then
          ignore=trim(ignore)//' pnu,'
        else
          use=trim(use)//' pnu,'
        endif
        if (ipr .ge. 10) then
           write(stdo,*)trim(use)
           write(stdo,*)trim(ignore)
        endif
        allocate(orhoat(3,nbas), v1pot(nbas),v0pot(nbas))
        ibaug = 0
        do  ib = 1, nbas
          ic=ssite(ib)%class
          is=ssite(ib)%spec !     is = -1 -> spec struc does not have these parameters
          if (is /= -1) then
            spid=sspec(is)%name
            a=sspec(is)%a
            nr=sspec(is)%nr
            rmt=sspec(is)%rmt
            z=sspec(is)%z
            lmxa=sspec(is)%lmxa
            lmxl=sspec(is)%lmxl
            kmax=sspec(is)%kmxt
            if (lmxa .eq. -1) goto 20
          endif
          ibaug = ibaug+1
          if (procid== master) then
              read(jfi) is0,spid0,lmxa0,lmxl0,nr0,rmt0,a0,z0,qc
              if (ib .gt. nbas) goto 20
            if (ipr .ge. 40) then
              if (ib .le. nat) write(stdo,380) ib,is0,spid0
              if (ib .gt. nat) write(stdo,380) ib,is0,spid0, ' (skip)'
            endif
  380       format('   atom',i4,'    species',i4,':',a:a)
c     ... read(but don't use) extra info since record is present
              read(jfi) rsma0,rsmr0,rsmv0,lmxv0,lmxr,lmxb0,kmax0
          endif
          call mpibc1_int(lmxa0,1,'iors_lmxa0')
          call mpibc1_int(lmxl0,1,'iors_lmxl0')
          call mpibc1_int(nr0,1,'iors_nr0')
          call mpibc1_real(a0,1,'iors_a0')
          call mpibc1_real(qc,1,'iors_qc')
          if (is == -1 ) call rx('iors: need check for is==-1')
          pnu=0d0
          pnz=0d0
          ql=0d0
          pnu=sspec(is)%p
          pnz=sspec(is)%pz
          idmod=0
          idmoz=0
          if (procid .eq. master) then
            do  isp = 1, nsp0
                read(jfi) (pnu(l+1,isp), l=0,lmxa0)
                read(jfi) (pnz(l+1,isp), l=0,lmxa0)
              if (nsp .gt. nsp0) then
                do  l = 0, lmxa0
                  pnu(l+1,2) = pnu(l+1,1)
                  pnz(l+1,2) = pnz(l+1,1)
                enddo
              endif
            enddo
          endif
          do  isp = 1, nsp
            call mpibc1_real(pnu(1,isp),lmxa0+1,'iors_pnu')
            call mpibc1_real(pnz(1,isp),lmxa0+1,'iors_pnu')
C       ... For backwards compatibility: prior versions wrote pnu for pnz
            do  l = 0, lmxa0+1
              if (pnu(l+1,isp) .eq. mod(pnz(l+1,isp),10d0))
     .        pnz(l+1,isp) = 0
            enddo
          enddo
          if (procid .eq. master) then
              read(jfi) (idmod(l+1), l=0,lmxa0)
              read(jfi) (idmoz(l+1), l=0,lmxa0)
          endif
          pnus =sspec(is)%p
          pnzs =sspec(is)%pz
          if  (irs5 .eq. 0 ) then 
            ssite(ib)%pnu=pnu
            ssite(ib)%pz=pnz
            sspec(is)%p=pnu
            sspec(is)%pz=pnz
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
            if (ltmp1 .and. ipr.ge.20) write(stdo,203) ib,spid,'file pnu',(pnu(i,1), i=1,lmxa+1)
            if (ltmp1 .and. ipr.ge.20) write(stdo,204) 'given pnu is',(pnus(i,1), i=1,lmxa+1)
            if (ltmp2 .and. ipr.ge.20) write(stdo,203) ib,spid,'file pz ',(pnz(i,1), i=1,lmxa+1)
            if (ltmp2 .and. ipr.ge.20) write(stdo,204) 'given pz  is',(pnzs(i,1), i=1,lmxa+1)
  203       format(9x,'site',i4,':',a,':',a,' is',8f6.2)
  204       format(26x,a,8f6.2)
            if (lfail .and. irs5 .eq. 1) then
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
              if (kmax0.ne.kmax.and.ipr.ge.10)
     .        write(stdo,201) ib,spid,'kmax',kmax0,kmax
              if (lmxa0.ne.lmxa.and.ipr.ge.10)
     .        write(stdo,201) ib,spid,'lmax',lmxa0,lmxa
  201         format(9x,'site',i4,', species ',a,
     .        ': augmentation ',a,' changed from',i2,' to',i2)
            endif
          else
            lmxl = lmxl0
            nr = nr0
         endif
!!         
          nlml0 = (lmxl0+1)**2
          nlml  = (lmxl+1)**2
          if (nr .ne. nr0) call rx('iors not set up to convert radial mesh')
          allocate(orhoat(1,ib)%v(nr*nlml*nsp)) !FP local densities rho1,rho2,rhoc and potentials v0, v1
          allocate(orhoat(2,ib)%v(nr*nlml*nsp))
          allocate(orhoat(3,ib)%v(nr*nsp))
          allocate(ssite(ib)%rv_a_ov0(nr*nsp))
          allocate(ssite(ib)%rv_a_ov1(nr*nsp))
          allocate(v0pot(ib)%v(nr*nsp)) 
          allocate(v1pot(ib)%v(nr*nsp))
          if (procid .eq. master) then
             print *,'nnnnnnnnread',ib,nr*nlml*nsp,nr,nlml,nsp
             call readrho(ifi,nr,nlml0,nsp0,nlml,nsp,orhoat(1,ib)%v)
             call readrho(ifi,nr,nlml0,nsp0,nlml,nsp,orhoat(2,ib)%v)
             call readrho(ifi,nr,1,nsp0,1,nsp,orhoat(3,ib)%v)
             call readrhos(ifi,nr,nsp0,nsp,ssite(ib)%rv_a_ov0)
             call readrhos(ifi,nr,nsp0,nsp,ssite(ib)%rv_a_ov1)
             if(nlml0.gt.nlml.and.ipr.ge.10) write(stdo,202) ib,spid,'truncate',nlml0,nlml
             if(nlml0.lt.nlml.and.ipr.ge.10) write(stdo,202) ib,spid,'inflate',nlml0,nlml
 202         format(9x,'site',i4,', species ',a,': ',a,' local density from nlm=',i3,' to',i3)
             v0pot(ib)%v = ssite(ib)%rv_a_ov0
             v1pot(ib)%v = ssite(ib)%rv_a_ov1
          endif
          call mpibc1_real( orhoat(1,ib)%v, size(orhoat(1,ib)%v), 'iors_rhoat(1)' )
          call mpibc1_real( orhoat(2,ib)%v, size(orhoat(2,ib)%v), 'iors_rhoat(2)' )
          call mpibc1_real( orhoat(3,ib)%v, size(orhoat(3,ib)%v), 'iors_rhoat(3)' )
          call mpibc1_real( ssite(ib)%rv_a_ov0 , size(ssite(ib)%rv_a_ov0) , 'iors_v0' )
          call mpibc1_real( ssite(ib)%rv_a_ov1 , size(ssite(ib)%rv_a_ov1) , 'iors_v1' )
          sspec(is)%a=a
          sspec(is)%nr=nr
          sspec(is)%rmt=rmt
          sspec(is)%z=z
          sspec(is)%lmxa=lmxa
          sspec(is)%lmxl=lmxl
          sspec(is)%kmxt=kmax
          sspec(is)%qc=qc
   20     continue
        enddo
        if (isanrg(ibaug, nat,nat,  msg,'nat', .false.)) goto 999
c        if ( lshear.and.irsrot) then 
c          wk(1:9)=lat_dist !rhoat rotation by lat_dist matrix
c          call dgemm('N','T',3,1,3,1d0,plat,3,qlat0,3,0d0,wk,3)
c          call pvsms2 ( ssite , sspec , wk , nbas , nsp , orhoat )
c        endif
C   --- Read data on free-atom core states and fit to fa density ---
        line = 'species data'
        do  30  is = 1, nspec
           a   =sspec(is)%a
           nr  =sspec(is)%nr
           lmxa=sspec(is)%lmxa
           if (lmxa .eq. -1) goto 30
           if (procid .eq. master) then
                 read(jfi,err=999,end=999) nr0,a0,qc,cof,eh,stc,lfoc0,rfoc0
              sspec(is)%lfoca=lfoc0
              sspec(is)%rfoca=rfoc0
              lfail = isanrg(nr0,nr,nr,msgw,'nr',.false.)
              call fsanrg(a0,a,a,0d-9,msg,'spec a',.true.)
           endif
           call mpibc1_real(qc,1, 'iors_qc')
           call mpibc1_real(cof,1,'iors_cof')
           call mpibc1_real(eh,1, 'iors_eh')
           call mpibc1_real(stc,1,'iors_stc')
C     ... FP core densities
           if (allocated(sspec(is)%rv_a_orhoc)) deallocate(sspec(is)%rv_a_orhoc)
           allocate(sspec(is)%rv_a_orhoc(abs(nr*nsp)))
           if (nr*nsp<0) sspec(is)%rv_a_orhoc(:)=0.0d0
           if (procid .eq. master) then
              if (nr .ne. nr0) call rx('iors not set up to convert core radial mesh')
                 read(jfi) sspec(is)%rv_a_orhoc !, nr * nsp0 , jfi )
              if (nsp .gt. nsp0) then !spin-split core density
                 i = nr
                 call dscal ( i , 0.5d0 , sspec(is)%rv_a_orhoc , 1 )
                 call dpscop ( sspec(is)%rv_a_orhoc , sspec(is)%rv_a_orhoc , i , 1 , 1 + i , 1d0  )
              endif
           endif
           call mpibc1_real(sspec(is)%rv_a_orhoc, size(sspec(is)%rv_a_orhoc), 'iors_rhoca'  )
           call dpzero(exi,n0)
           call dpzero(hfc,n0*2)
           if (procid .eq. master) then
                 read(jfi,err=999,end=999) rsmfa,nxi
                 read(jfi,err=999,end=999) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp0)
              if (nsp .gt. nsp0) then
                 i = n0
                 call dscal(i,0.5d0,hfc,1)
                 call dpscop(hfc,hfc,i,1,1+i,1d0)
              endif
           endif
           call mpibc1_real(rsmfa,1,'iors_rsmfa')
           call mpibc1_int(nxi,1,'iors_nxi')
           call mpibc1_real(exi,nxi,'iors_exi')
           call mpibc1_real(hfc,nsp*nxi,'iors_hfc')
           sspec(is)%ctail=cof
           sspec(is)%etail=eh
           sspec(is)%stc=stc
           sspec(is)%nxi=nxi
           sspec(is)%exi=exi
           sspec(is)%chfa=hfc
           sspec(is)%rsmfa=rsmfa
 30     continue

c   ... Copy or rescale cores, in case foca was switched on or off
          do  ib = 1, nbas
            is = int(ssite(ib)%spec)
            a=sspec(is)%a
            nr=sspec(is)%nr
            rmt=sspec(is)%rmt
            lmxa=sspec(is)%lmxa
            lfoc=sspec(is)%lfoca
            qc=sspec(is)%qc
            if (lmxa .eq. -1) goto 40
            if (lfoc .gt. 0) then
              call dpcopy ( sspec(is)%rv_a_orhoc , orhoat( 3 , ib )%v , 1 , nr * nsp , 1d0 )
            else
              allocate(rwgt_rv(nr))
              call radwgt ( rmt , a , nr , rwgt_rv )
              call radsum ( nr , nr , 1 , nsp , rwgt_rv , orhoat( 3 , ib )%v , sum )
              fac = 1d0
              if (dabs(sum) .gt. 1d-6) fac = qc/sum
              if (dabs(fac-1d0).gt.1d-7 .and. ipr.ge.30)
     .        write(stdo,787) ib,qc,sum,fac
  787         format(' fix core chg: ib=',i4,'  qc,sum,fac=',3f12.6)
              call dpcopy ( orhoat( 3 , ib )%v , orhoat( 3 , ib )%v, 1 , nr * nsp , fac )
              if (allocated(rwgt_rv)) deallocate(rwgt_rv)
            endif
   40       continue
          enddo
        do i_site=1,nbas
          call mpibc1_s_site(ssite(i_site),'iors_ssite')
        enddo
        do i_spec=1,nspec
          call mpibc1_s_spec(sspec(i_spec),'iors_sspec')
        enddo
C --- Output ---
      else
        if (procid .ne. master) then
          iors = 0
          return
        endif
        jfi = ifi
c        rewind jfi
        fid0 = fid
        call strip(fid0,i0,i1)
        jobid = sname !datimp(2:)
        call ftime(datimp)
        hostnm = ' '
        usernm = ' '
        call get_environment_variable('HOST',hostnm)
        call get_environment_variable('USER',usernm)
        call strip(datimp,i,i2)
        call strip(usernm,i,i3)
        call strip(hostnm,i,i4)
        if (ipr .ge. 40) write(stdo,710) fid(1:i1), usernm(1:i3),hostnm(1:i4),datimp(1:i2)
  721   format('----------------------- ',a,' -----------------------')
          write(jfi) vs
          write(jfi) fid0
          write(jfi) datimp,usernm,hostnm,jobid
c          if (abs(vs) .le. 1.021) then
c            write(jfi) nbas,nsp,npan,lrel
c          else if (abs(vs) .le. 1.031) then
c            write(jfi) nbas,nsp,npan,lrel,nspec
c          else
            write(jfi) nbas,nat,nsp,npan,lrel,nspec
c          endif
          write(jfi) nit
          write(jfi) alat,vol,plat

C   --- Write smooth charge density ---
          write(jfi) n1,n2,n3
          write(jfi) osmrho!, lbin , -jfi )

C   --- Write information related to dynamics ---
        wk=1d99 !call dpzero(wk,100)
        wk(1)= eferm !sbz%ef !dummy ! we use wk(1) only wk(2:100) are dummy
          write(jfi) wk !call dpdump(wk,100,-jfi)
        do  110  ib = 1, nbas
          pos  =ssite(ib)%pos
          force=ssite(ib)%force
c          vel  =ssite(ib)%vel
            write(jfi) ib,pos,force!,vel
  110   continue
C   --- Write information for local densities ---
        if (ipr .ge. 50) write(stdo,364)
        do  120  ib = 1, nbas
          ic=ssite(ib)%class
          is=ssite(ib)%spec
          spid=sspec(is)%name
          a=sspec(is)%a
          nr=sspec(is)%nr
          rmt=sspec(is)%rmt
          z=sspec(is)%z
          qc=sspec(is)%qc
          idmod=sspec(is)%idmod
          rsma=sspec(is)%rsma
          lmxa=sspec(is)%lmxa
          lmxl=sspec(is)%lmxl
          lmxb=sspec(is)%lmxb
          kmxv=sspec(is)%kmxv
          rsmv=sspec(is)%rsmv
          kmax=sspec(is)%kmxt
         pnu=ssite(ib)%pnu
         pnz=ssite(ib)%pz
          if (lmxa .eq. -1) goto 120
            write(jfi) is,spid,lmxa,lmxl,nr,rmt,a,z,qc
C     ... Some extra info... lots of it useless or obsolete
          lmxr = 0
          lmxv = 0
          rsmr = 0
            write(jfi) rsma,rsmr,rsmv,lmxv,lmxr,lmxb,kmax
C     ... Write augmentation data
          do  122  isp = 1, nsp
              write(jfi) (pnu(l+1,isp), l=0,lmxa)
              write(jfi) (pnz(l+1,isp), l=0,lmxa)
  122     continue
C         Write for compatibility with nfp
            write(jfi) (idmod(l+1), l=0,lmxa)
            write(jfi) (idmod(l+1), l=0,lmxa)
C     ... Write arrays for local density and potential
          nlml = (lmxl+1)**2
              print *,'nnnnnnnwrite',ib,nr*nlml*nsp,nr,nlml,nsp
            write(jfi) orhoat( 1 , ib )%v !, nr , nlml , nlml , nsp  , nsp , lbin , -jfi )
            write(jfi) orhoat( 2 , ib )%v !, nr , nlml , nlml , nsp  , nsp , lbin , -jfi )
            write(jfi) orhoat( 3 , ib )%v !, nr , 1 , 1 , nsp , nsp   , lbin , -jfi )
            write(jfi) ssite(ib)%rv_a_ov0 !, nr , 1 , 1 , nsp , nsp , lbin , -jfi  )
            write(jfi) ssite(ib)%rv_a_ov1 !, nr , 1 , 1 , nsp , nsp , lbin , -jfi  )
            if (ipr .ge. 50) then
               write(stdo,349) ib,spid,lmxa,lmxl,rmt,nr,a, (pnu(l+1,1),l=0,lmxa)
               if (nsp .eq. 2)  write(stdo,350) (pnu(l+1,2), l=0,lmxa)
            endif
  120   continue
C   --- Write data on free-atom core states and fit to fa density ---
        do  130  is = 1, nspec
          a  =sspec(is)%a
          nr =sspec(is)%nr
          qc =sspec(is)%qc
          lmxa=sspec(is)%lmxa
          lfoc=sspec(is)%lfoca
          rfoc=sspec(is)%rfoca
          cof =sspec(is)%ctail
          eh  =sspec(is)%etail
          stc =sspec(is)%stc
          nxi =sspec(is)%nxi
          exi=sspec(is)%exi
          hfc=sspec(is)%chfa
          rsmfa = sspec(is)%rsmfa
          if (lmxa .eq. -1) goto 130
            write(jfi) nr,a,qc,cof,eh,stc,lfoc,rfoc
C     ... For now, ASA stores no core data
              write(jfi) sspec(is)%rv_a_orhoc 
              write(jfi) rsmfa,nxi
              write(jfi) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp)
  130   continue
      endif
  349 format(i11,':',a4,2i2,f9.5,i5,f6.3,1x,8f6.3)
  350 format(41x,8f6.3)
  364 format(/9x,'ib:spc la ll   rmt     nr   a     pnu')
      iors = 0
      return
  998 if (ipr .gt. 0) write(stdo,'('' iors  : empty file ... nothing read'')')
      return
  999 continue
      if (ipr .gt. 0) write(stdo,'('' iors  : read failed in: '',a)') line
c      call tcx('iors')
      end function iors
      end module m_iors
