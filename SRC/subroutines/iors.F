      module m_iors
      contains
      integer function iors(fid,nit,lbin,ifi)
      use m_supot,only: zv_a_osmrho, sv_p_orhoat=>sv_a_oorhat

      use m_supot,only: lat_nabc
      use m_struc_func,only:mpibc1_s_spec,mpibc1_s_site
      use m_bndfp,only: m_bndfp_ef_SET,bz_def,bz_ef
      use m_lmfinit,only: irs3,irs4,irs5,irsrot,lat_alat,nsp,lrel,nl,ssite=>v_ssite,sspec=>v_sspec,
     &    nbas,nat,nspec 
      use m_lattic,only: lat_dist, lat_plat
!!- I/O for charge density to rst or rsta. ssite sspec are readin 
!! read write
!!     smrho, rhoat
!!     ssite: pos, pos0, force, vel, pnu pz ov0,ov1 
!      sspec:
!!           a nr rmt z lmxa lmxl kmxt p pz lfoca qc orhoc idmod
!!           rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
!!           chfa rsmfa
!      
!! ef,def --> goto bndfp via m_bndfp_ef_set      
!!      
Ci   fid   :string containing identification
Ci   nbas  :size of basis
Ci   nat   :number atoms in basis with augmentation sites
Ci         :Note: if nat<nbas, there is a requirement that
Ci         :lmxa>-1 for nat sites, and
Ci         :and lmxa=-1 for nbas-nat sites
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lbin  :T file I/O in binary mode
Ci         :F file I/O in ascii mode
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
      implicit none
C ... Passed parameters
      logical lbin
      integer::  nit , ifi , n0,i_site,i_spec,i_copy_size !mode=1 ,
      parameter (n0=10)
c      type(s_rv1) :: sv_p_orhoat(3,nbas)
c      type(s_site)::ssite(nbas)
c      type(s_spec)::sspec(nspec)
      character fid*(*)
C ... Local parameters
      integer procid,master,mpipid,nproc
      integer i,i0,i1,i2,i3,i4,ib,ipr,iprint,ic,is,is0,isp,jb,k1,k2,k3,
     .igetss,jfi,k11,k21,k31,kmax,kmax0,kmxv,l,lfoc,lfoc0,lmxa,
     .lmxa0,lmxb,lmxb0,lmxl,lmxl0,lmxr,lmxv,lmxv0,lrel0,n11,n21,
     .n31,nbas0,nspec0,nlml,nlml0,npan,npan0,nr,nr0,nsp0,
     .nxi,nat0,ibaug
      integer:: ngabc(3) , n1 , n2 , n3 , lgunit , stdo , isw
      complex(8),allocatable :: zv_a_osmrhox(:)
      complex(8) ,allocatable :: h_zv(:)
      real(8) ,allocatable :: rwgt_rv(:)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer idmod(n0),idmoz(n0) !,lrs(10)
      logical isanrg,lfail,ltmp1,ltmp2,latvec,lshear,cmdopt,mlog
      double precision a,a0,alat,alat0,cof,eh,fac,qc,rfoc,rfoc0,rmt,
     .rmt0,rsma,rsma0,rsmfa,rsmr,rsmr0,rsmv,rsmv0,stc,sum,vfac,vol,
     .vol0,vs,vs1,z,z0
      double precision pnu(n0,2),pnz(n0,2),ql(n0,2*n0),pos(3),
     .pos0(3),force(3),vel(3),plat(3,3),plat0(3,3),qlat(3,3),qlat0(3,3),
     .exi(n0),hfc(n0,2),vec0(3),wk(100),rh,vrmax(2),pnus(n0,2),
     .pnzs(n0,2),dval
      character spid*8,spid0*8,fid0*68,line*20,msg*23,use*80,ignore*80,
     .msgw*17,datimp*32,usernm*32,hostnm*32,jobid*32,ffmt*32,ifmt*32
      logical:: l_dummy_aiomom, aiomom
      integer:: fextg, i_dummy_fextg
      data vec0 /0d0,0d0,0d0/
      real(8):: bz_ef00,bz_def00
C ... MPI setup
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,ignore)
      alat=lat_alat
      plat=lat_plat
      ngabc=lat_nabc
      call dinv33(plat,1,qlat,fac)
c      RWiors=''
      vol = dabs(fac)*alat**3
      ipr    = iprint()
      stdo   = lgunit(1)
      vs   =  1.04d0
c      nsp = globalvariables%nsp
c      lrel   = globalvariables%lrel
c      lrs(1) = irs3 !isw ( iand(32,int(sctrl%lrs)) .ne.0 )
c      lrs(2) = irs4 !isw ( iand(64,int(sctrl%lrs)) .ne.0 )
c      lrs(3) = irs5 !isw ( iand(128,int(sctrl%lrs)) .ne.0 )
      msg  = '         File mismatch:'
      msgw = '         warning:'
      iors = -1
C     ndyn = igets('move size',smove)
      line = 'header'
      call fftz30(n1,n2,n3,k1,k2,k3)
      ffmt = '(5f15.10)'
      ifmt = '(20i5)'
C ... Hardwired for now
      npan = 1
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
        fid0(i:) = ', mesh density'
        i = i+5
        i = i+14-5
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
            call info0(20,0,0,'%9f(warning) rst mismatch in nat ... skipping sites')
          endif
C       if (isanrg(nsp0,  nsp,nsp,  msg,'nsp', .true.)) goto 999
          if (nsp0 .lt. nsp) call info0(20,0,0,'%9f(warning) rst file not spin pol .. splitting spins')
          if (isanrg(npan0, npan,npan,msg,'npan',.true.)) goto 999
          lfail = isanrg(lrel0, lrel,lrel,msgw,'lrel',.false.)
          call fsanrg(abs(vs1),1.00d0,abs(vs),0d0,' ','file''s version',.true.)
        endif
        call mpibc1(nbas0,1,2,mlog,'iors','nbas0')
        call mpibc1(nat0,1,2,mlog,'iors','nat0')
        call mpibc1(plat,9,4,mlog,'iors','plat')
        call mpibc1(nit,1,2,mlog,'iors','nit')
C   --- Read smooth charge density ---
          line = 'smoothed density'
          if (procid .eq. master) then
            if (lbin) then
              read(jfi,err=999,end=999) n11,n21,n31
            else
              read(jfi,*,err=999,end=999)
              read(jfi,*,err=999,end=999) n11,n21,n31
            endif
            if (n11.eq.n1 .and. n21.eq.n2 .and. n31.eq.n3) then
              call dpdftr ( n1 , n2 , n3 , k1 , k2 , k3 , nsp0 , zv_a_osmrho, lbin , jfi )
              if (nsp .gt. nsp0) then
                i = k1*k2*k3*2
                call dscal ( i , 0.5d0 , zv_a_osmrho , 1 )
                call dpscop ( zv_a_osmrho , zv_a_osmrho , i , 1 , 1 + i , 1d0 )
              endif
c   ... or read and remesh
            else
              if (ipr .ge. 10) write(stdo,450) n11,n21,n31,n1,n2,n3
  450         format(9x,'remesh density from  ',i4,'  *',i4,'  *',i4,
     .        '    to  ',i4,'  *',i4,'  *',i4)
              call fftz30(n11,n21,n31,k11,k21,k31)
              allocate(h_zv(k11*k21*k31*nsp))
              h_zv(:)=0.0d0
              call dpdftr ( n11 , n21 , n31 , k11 , k21 , k31 , nsp0 , h_zv , lbin , jfi )
              if (nsp .gt. nsp0) then
                i = k11*k21*k31*2
                call dscal ( i , 0.5d0 , h_zv , 1 )
                call dpscop ( h_zv , h_zv , i , 1 , 1 + i , 1d0 )
              endif
              call pshpr(50)
              i = 0
              if (n1 .eq. 2*n11 .and. n2 .eq. 2*n21 .and. n3 .eq. 2*n31) i=3
              call chgmsh ( i , plat , nsp , n11 , n21 , n31 , k11 , k21 ,
     .         k31 , h_zv , n1 , n2 , n3 , k1 , k2 , k3 , zv_a_osmrho )
              call poppr
              if (allocated(h_zv)) deallocate(h_zv)
            endif
C   ... If cell volume changed, scale smooth density to maintain charge
            vfac = vol0/vol
            if (dabs(vfac-1d0) .gt. 1d-8) then
              if (ipr .ge. 10) write(stdo,460) vol0,vfac
  460         format(9x,'volume changed from',f8.2,' :  scale smooth density by',f8.4)
              call dpcopy ( zv_a_osmrho, zv_a_osmrho, 1 , 2 * k1 * k2 * k3* nsp , vfac )
            endif
          endif
          call mpibc1 ( zv_a_osmrho , 2 * k1 * k2 * k3 * nsp , 4 , mlog, 'iors' , 'smrho' )
  115     continue
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
        if (irs4 .ne. 0) then
          call awrit0('%a ef window,',ignore,len(ignore),0)
        else
          call awrit0('%a ef window,',use,len(use),0)
          bz_ef00=wk(1)
          bz_def00=wk(2)
          call m_bndfp_ef_SET(bz_ef00,bz_def00)
c          bz_w= -9999 !dummy wk(3)
        endif
C   --- Read atomic positions,forces,velocities ---
        line = 'site data'
        if (irs3 .ne. 0) then
          call awrit0('%a positions,',ignore,len(ignore),0)
        else
          call awrit0('%a positions,',use,len(use),0)
          if (lshear) call awrit0('%a%7b(sheared),',use,len(use),0)
C         Must be enough sites to read from rst file
          if (nbas0 .lt. nbas) then
            call info0(1,0,0,'%9f oops ... cannot use file site positions (site mismatch)')
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
          if (irs3 .ne. 0) then
            pos=ssite(ib)%pos
            vel=0d0 
          endif
           ssite(ib)%pos  =pos
           ssite(ib)%pos0 =pos0
           ssite(ib)%force=force
           ssite(ib)%vel  =vel
   10     continue
        enddo

C   --- Read information for local densities ---
        if (irs5 .ne. 0) then
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
C         is = -1 -> spec struc does not have these parameters
C         lskip = .false.
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
          if (is == -1 ) call rx('iors: need check for is==-1')
C         is = -1 -> pack spec struc with file parameters
c$$$         if(is==-1) then       
c$$$            is = is0
c$$$c            do i_spacks=is,is
c$$$c              call spacks_copy('p',sspec(i_spacks)%name,is,is,spid0,i_spacks)
c$$$c            enddo
c$$$            sspec(is)%name=spid0
c$$$            spid = spid0
c$$$            sspec(is)%a=a0
c$$$            sspec(is)%nr=nr0
c$$$            sspec(is)%rmt=rmt0
c$$$            sspec(is)%z=z0
c$$$            sspec(is)%lmxa=lmxa0
c$$$            sspec(is)%lmxl=lmxl0
c$$$            sspec(is)%kmxt=kmax0
c$$$            a=sspec(is)%a
c$$$            nr=sspec(is)%nr
c$$$            rmt=sspec(is)%rmt
c$$$            z=sspec(is)%z
c$$$            lmxa=sspec(is)%lmxa
c$$$            lmxl=sspec(is)%lmxl
c$$$            kmax=sspec(is)%kmxt
c$$$          endif
C     ... Read augmentation data
c          call dpzero(pnu,2*n0)
c          call dpzero(pnz,2*n0)
c          call dpzero(ql,2*n0)
          pnu=0d0
          pnz=0d0
          ql=0d0
          i_copy_size=size(sspec(is)%p)
          call dcopy(i_copy_size,sspec(is)%p,1,pnu,1)
          i_copy_size=size(sspec(is)%pz)
          call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)
          idmod=0
          idmoz=0
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
          if  (irs5 .eq. 0 ) then 
            i_copy_size=size(ssite(ib)%pnu)
            call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
            i_copy_size=size(ssite(ib)%pz)
            call dcopy(i_copy_size,pnz,1,ssite(ib)%pz,1)
            i_copy_size=size(sspec(is)%p)
            call dcopy(i_copy_size,pnu,1,sspec(is)%p,1)
            i_copy_size=size(sspec(is)%pz)
            call dcopy(i_copy_size,pnz,1,sspec(is)%pz,1)
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
C         Case read but skip over this site data
          else
            lmxl = lmxl0
            nr = nr0
          endif

C     --- Allocate and read arrays for local density and potential ---
          nlml0 = (lmxl0+1)**2
          nlml = (lmxl+1)**2
          if (allocated(ssite(ib)%rv_a_ov0)) deallocate(ssite(ib)%rv_a_ov0)
          allocate(ssite(ib)%rv_a_ov0(abs(nr*nsp)))
C     ... FP local densities rho1,rho2,rhoc and potentials v0, v1
            if (nr .ne. nr0) call rx('iors not set up to convert radial mesh')
            allocate(sv_p_orhoat(1,ib)%v(abs(nr*nlml*nsp)))
            allocate(sv_p_orhoat(2,ib)%v(abs(nr*nlml*nsp)))
            allocate(sv_p_orhoat(3,ib)%v(abs(nr*nsp)))
            if (allocated(ssite(ib)%rv_a_ov1)) deallocate(ssite(ib)%rv_a_ov1)
            allocate(ssite(ib)%rv_a_ov1(abs(nr*nsp)))
            if (procid .eq. master) then
              if (.not. lbin) read(jfi,*)
              call dpdbyl ( sv_p_orhoat( 1 , ib )%v , nr0 , nlml0 , nlml , nsp0 , nsp , lbin , jfi )
              if (.not. lbin) read(jfi,*)
              call dpdbyl ( sv_p_orhoat( 2 , ib )%v , nr0 , nlml0 , nlml , nsp0 , nsp , lbin , jfi )
              if (nlml0.gt.nlml.and.ipr.ge.10) write(stdo,202) ib,spid,'truncate',nlml0,nlml
              if (nlml0.lt.nlml.and.ipr.ge.10) write(stdo,202) ib,spid,'inflate',nlml0,nlml
  202         format(9x,'site',i4,', species ',a,': ',a,' local density from nlm=',i3,' to',i3)
              if (.not. lbin) read(jfi,*)
              call dpdbyl(sv_p_orhoat(3,ib)%v,nr0,1,1,nsp0, nsp , lbin , jfi )
              if (.not. lbin) read(jfi,*)
              call dpdbyl ( ssite(ib)%rv_a_ov0, nr0, 1, 1, nsp0, nsp , lbin , jfi )
              if (.not. lbin) read(jfi,*)
              call dpdbyl ( ssite(ib)%rv_a_ov1, nr0, 1, 1, nsp0, nsp , lbin , jfi )
              if (nsp0 .lt. nsp) then
                call dscal ( nr0 * 2 , 2d0 , ssite(ib)%rv_a_ov0 , 1 )
                call dscal ( nr0 * 2 , 2d0 , ssite(ib)%rv_a_ov1 , 1 )
              endif
            endif
            call mpibc1( sv_p_orhoat( 1 , ib )%v , nr * nlml * nsp , 4 , mlog , 'iors' , 'rhoat(1)' )
            call mpibc1( sv_p_orhoat( 2 , ib )%v , nr * nlml * nsp , 4 , mlog , 'iors' , 'rhoat(2)' )
            call mpibc1( sv_p_orhoat( 3 , ib )%v , nr  * nsp , 4 ,       mlog , 'iors' , 'rhoat(3)' )
            call mpibc1( ssite(ib)%rv_a_ov0 , nr * nsp , 4 , mlog , 'iors' , 'v0' )
            call mpibc1( ssite(ib)%rv_a_ov1 , nr * nsp , 4 , mlog , 'iors' , 'v1' )
c          endif
C     ... store data in strucs
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
        if ( lshear.and.irsrot) then 
          wk(1:9)=lat_dist !rhoat rotation by lat_dist matrix
          call dgemm('N','T',3,1,3,1d0,plat,3,qlat0,3,0d0,wk,3)
          call pvsms2 ( ssite , sspec , wk , nbas , nsp , sv_p_orhoat )
        endif
C   --- Read data on free-atom core states and fit to fa density ---
        line = 'species data'
        if (procid .eq. master) then
          if (.not. lbin) read(jfi,*,err=999,end=999)
        endif
        do  30  is = 1, nspec
          a   =sspec(is)%a
          nr  =sspec(is)%nr
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
            if (allocated(sspec(is)%rv_a_orhoc)) deallocate(sspec(is)%rv_a_orhoc)
            allocate(sspec(is)%rv_a_orhoc(abs(nr*nsp)))
            if (nr*nsp<0) sspec(is)%rv_a_orhoc(:)=0.0d0
            if (procid .eq. master) then
              if (nr .ne. nr0) call rx('iors not set up to convert core radial mesh')
              if (lbin) then
                call dpdump ( sspec(is)%rv_a_orhoc , nr * nsp0 , jfi )
              else
                call dfdump ( sspec(is)%rv_a_orhoc , nr * nsp0 , jfi )
              endif
C           spin-split core density
              if (nsp .gt. nsp0) then
                i = nr
                call dscal ( i , 0.5d0 , sspec(is)%rv_a_orhoc , 1 )
                call dpscop ( sspec(is)%rv_a_orhoc , sspec(is)%rv_a_orhoc , i , 1 , 1 + i , 1d0  )
              endif
            endif
            call mpibc1 ( sspec(is)%rv_a_orhoc , nr * nsp , 4 , mlog , 'iors' , 'rhoca'  )
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
          sspec(is)%ctail=cof
          sspec(is)%etail=eh
          sspec(is)%stc=stc
          sspec(is)%nxi=nxi
          i_copy_size=size(sspec(is)%exi)
          call dcopy(i_copy_size,exi,1,sspec(is)%exi,1)
          i_copy_size=size(sspec(is)%chfa)
          call dcopy(i_copy_size,hfc,1,sspec(is)%chfa,1)
          sspec(is)%rsmfa=rsmfa
   30   continue

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
              call dpcopy ( sspec(is)%rv_a_orhoc , sv_p_orhoat( 3 , ib )%v , 1 , nr * nsp , 1d0 )
            else
              allocate(rwgt_rv(nr))
              call radwgt ( rmt , a , nr , rwgt_rv )
              call radsum ( nr , nr , 1 , nsp , rwgt_rv , sv_p_orhoat( 3 , ib )%v , sum )
              fac = 1d0
              if (dabs(sum) .gt. 1d-6) fac = qc/sum
              if (dabs(fac-1d0).gt.1d-7 .and. ipr.ge.30)
     .        write(stdo,787) ib,qc,sum,fac
  787         format(' fix core chg: ib=',i4,'  qc,sum,fac=',3f12.6)
              call dpcopy ( sv_p_orhoat( 3 , ib )%v , sv_p_orhoat( 3 , ib )%v, 1 , nr * nsp , fac )
              if (allocated(rwgt_rv)) deallocate(rwgt_rv)
            endif
   40       continue
          enddo
        do i_site=1,nbas
          call mpibc1_s_site(ssite(i_site),mlog,'iors','ssite')
        enddo
        do i_spec=1,nspec
          call mpibc1_s_spec(sspec(i_spec),mlog,'iors','sspec')
        enddo
c        RWiors='read'
C --- Output ---
      else
c        RWiors='write'
        if (procid .ne. master) then
          iors = 0
          return
        endif
        jfi = -ifi
        rewind jfi
C ...   Write header and general information
        fid0 = fid
        call strip(fid0,i0,i1)
        i_dummy_fextg= fextg(datimp)
        jobid = datimp(2:)
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
            write(jfi,'(5i5,t50,a)') nbas,nsp,npan,lrel,nit,'! nbas nsp npan lrel nit'
          else if (abs(vs) .le. 1.031) then
            write(jfi,'(6i5,t50,a)') nbas,nsp,npan,lrel,nit,nspec,'! nbas nsp npan lrel nit nspec'
          else
            write(jfi,'(7i5,t50,a)') nbas,nat,nsp,npan,lrel,nit,nspec,
     .      '! nbas nat nsp npan lrel nit nspec'
          endif
          write(jfi,'(f15.10,1pe20.10,t50,a)') alat,vol,'! alat vol; plat'
          write(jfi,'(3f15.10)') plat
        endif

C   --- Write smooth charge density ---
        if (lbin) then
          write(jfi) n1,n2,n3
        else
          write(jfi,721) 'smooth density'
          write(jfi,ifmt) n1,n2,n3
        endif
        call dpdftr ( n1 , n2 , n3 , k1 , k2 , k3 , nsp , zv_a_osmrho, lbin , - jfi )

C   --- Write information related to dynamics ---
        wk=0d0 !call dpzero(wk,100)
        wk(1)=bz_ef !sbz%ef
        wk(2)=bz_def !sbz%def
        wk(3)=-9999 !dummy bz_w !sbz%w
        if (lbin) then
          call dpdump(wk,100,-jfi)
        else
          write(jfi,721) 'geometric data'
          call dfdump(wk,100,-jfi)
        endif
        do  110  ib = 1, nbas
          pos  =ssite(ib)%pos
          force=ssite(ib)%force
          vel  =ssite(ib)%vel
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
          spid=sspec(is)%name
          a=sspec(is)%a
          nr=sspec(is)%nr
          rmt=sspec(is)%rmt
          z=sspec(is)%z
          qc=sspec(is)%qc
          i_copy_size=size(sspec(is)%idmod)
          call icopy(i_copy_size,sspec(is)%idmod,1,idmod,1)
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
c         pnu=ssite(ib)%pnu
c         pnz=ssite(ib)%pz
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
            if (.not. lbin) write(jfi,'('' rho1'')')
            call dpdbyl ( sv_p_orhoat( 1 , ib )%v , nr , nlml , nlml , nsp  , nsp , lbin , - jfi )
            if (.not. lbin) write(jfi,'('' rho2'')')
            call dpdbyl ( sv_p_orhoat( 2 , ib )%v , nr , nlml , nlml , nsp  , nsp , lbin , - jfi )
            if (.not. lbin) write(jfi,'('' rhoc'')')
            call dpdbyl ( sv_p_orhoat( 3 , ib )%v , nr , 1 , 1 , nsp , nsp   , lbin , - jfi )
            if (.not. lbin) write(jfi,'('' v0'')')
            call dpdbyl ( ssite(ib)%rv_a_ov0 , nr , 1 , 1 , nsp , nsp , lbin , - jfi  )
            if (.not. lbin) write(jfi,'('' v1'')')
            call dpdbyl ( ssite(ib)%rv_a_ov1 , nr , 1 , 1 , nsp , nsp , lbin , - jfi  )
            if (ipr .ge. 50) then
               write(stdo,349) ib,spid,lmxa,lmxl,rmt,nr,a, (pnu(l+1,1),l=0,lmxa)
               if (nsp .eq. 2)  write(stdo,350) (pnu(l+1,2), l=0,lmxa)
            endif
  120   continue

C   --- Write data on free-atom core states and fit to fa density ---
        if (.not. lbin) write(jfi,721) 'species data'
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
          i_copy_size=size(sspec(is)%exi)
          call dcopy(i_copy_size,sspec(is)%exi,1,exi,1)
          i_copy_size=size(sspec(is)%chfa)
          call dcopy(i_copy_size,sspec(is)%chfa,1,hfc,1)
          rsmfa = sspec(is)%rsmfa
          if (lmxa .eq. -1) goto 130
          if (lbin) then
            write(jfi) nr,a,qc,cof,eh,stc,lfoc,rfoc
          else
c            do i_spacks=is,is
c              call spacks_copy('u',sspec(i_spacks)%name,is,is,spid0,i_spacks)
c            enddo
            spid0=sspec(is)%name
            write(jfi,'('' spec'',i5,2x,'' label '',a)') is, spid0
            write(jfi,'(''  nr   lfoc'')')
            write(jfi,ifmt) nr,lfoc
            write(jfi,'(t7,''  a '',t21,''qc'',t36,''cof'',t61,''eh'')')
            write(jfi,'(2f15.10,1p,2e25.15)') a,qc,cof,eh
            write(jfi,'(t17,'' stc'',t31,''rfoc'')')
            write(jfi,'(1p,2e25.15)') stc,rfoc
          endif
C     ... For now, ASA stores no core data
            if (lbin) then
              call dpdump ( sspec(is)%rv_a_orhoc , nr * nsp , - jfi )
              write(jfi) rsmfa,nxi
              write(jfi) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp)
            else
              call dfdump ( sspec(is)%rv_a_orhoc , nr * nsp , - jfi )
              write(jfi,'('' core tail: rsmfa,n,e,cof=''/f15.10,i4)')
     .        rsmfa,nxi
              write(jfi,ffmt) (exi(i),i=1,nxi)
              call dfdump(hfc(1,1),nxi,-jfi)
              if (nsp .eq. 2) call dfdump(hfc(1,2),nxi,-jfi)
            endif
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
      end function iors

      end module m_iors
