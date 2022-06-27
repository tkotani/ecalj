module m_iors
  use m_struc_def
  use m_lgunit,only:stdo
  public iors
  type(s_rv1),public,allocatable :: v1pot(:),v0pot(:)
  private

contains
  integer function iors(nit,rwrw) 
    use m_density,only: osmrho, orhoat !these are allocated
    use m_bndfp,only: m_bndfp_ef_SET,eferm

    use m_supot,only: lat_nabc
    use m_struc_func,only: mpibc1_s_spec,mpibc1_s_site
    use m_lmfinit,only: alat=>lat_alat,nsp,lrel,nl,ssite=>v_ssite,sspec=>v_sspec, nbas,nat,nspec,n0,&
         idmodis=>idmod,slabl,rsma,rsmfa
    use m_lattic,only: plat=>lat_plat,vol=>lat_vol,qlat=>lat_qlat
    use m_ext,only:sname
    use m_ftox
    !!- I/O for charge density to rst or rsta. ssite sspec are readin
    !! read write
    !!     smrho, rhoat
    !!     ssite: pos, pos0, force, pnu pz ov0,ov1
    !      sspec:
    !!           a nr rmt z lmxa lmxl kmxt p pz lfoca qc orhoc idmod
    !!           rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
    !!           chfa rsmfa

    !! ef,def --> goto bndfp via m_bndfp_ef_set
    !!
    ! xxx   fid   :string containing identification
    !i   nbas  :size of basis
    !i   nat   :number atoms in basis with augmentation sites
    !i         :Note: if nat<nbas, there is a requirement that
    !i         :lmxa>-1 for nat sites, and
    !i         :and lmxa=-1 for nbas-nat sites
    !i   ifi   :file logical unit, but >0 for read, <0 for write
    !i   lbin  :=T file I/O in binary mode
    !i         :xxxF file I/O in ascii mode
    !!
    !r Remarks
    !r   The density consists of a smooth part (smrho) plus
    !r   nbas atom-centered densities inside the MT spheres.
    !r   Their sum is the full charge density.
    !r   The local density is represented as the difference of the
    !r   two valence components in orhoat, plus the core density.
    !r   Density in the MT spheres:
    !r      mesh parameters rmt,nr,a;
    !r      total density rho (times r**2) to lmxl;
    !r      a spherical potential v0 defining the wave functions within rmt
    !r      pnu and idmod to lmxa
    !r      NOTE: on input, arrays for rhoat and v0 are allocated here
    !r   Smooth density
    !r      real part of complex*16 array smrho contains the density
    !r      k1,k2,k3 are the physical dimensions of the array
    !r      n1,n2,n3 are the dimensions of the mesh.
    !r   Additional information stored:
    !r      fid: file identifier, a string of length 64 characters or less.
    !r      parameters relating to coordinates and molecular dynamics.
    !r   On input, iors tries to transform data format where needed:
    !r      lmxl incompatible: pad or truncate l-components
    !r      FT mesh changed: map onto new mesh
    !l Local variables
    !l   lrs switches:  0=>use rst file data; 1=>ignore rst file data
    !l   lrs(1) site positions
    !l   lrs(2) starting fermi level
    !l   lrs(3) starting pnu's
    !m MPI
    !m   Master process reads and broadcasts line by line. err= and end= are
    !m   troublesome since the slave processes may hang if the rst file is
    !m   corrupted or incompatible. For now if iors returns < 1 lmfp will
    !m   exit and hope the slave processes follow suit!
    !u Updates
    !u   01 Jul 08 New mode -2
    !u   25 Aug 07 Bug fix, case floating orbitals not positioned at end
    !u   10 Jul 07 Cleaner error exit, MPI
    !u   20 Jun 06 Repackaged MPI
    !u   07 Jul 05 (version update 1.04)
    !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
    !u   04 Feb 05 Spin-split file non-spin-polarized density (FP)
    !u   26 Apr 03 Added MPI parallelization for ASA
    !u   11 Jan 03 Bug fix: corrected calc. vol (might have been<0)
    ! xxx   10 Dec 02 File shears atom positions by shear transformation
    ! xxx             (slat->plat) (file plat)^-1
    !u   19 Feb 02 (version update 1.03)
    !u             File now contains nspec
    !u             File contents backwardly compatible with prior versions.
    !u             New mode (-1)
    !u             Routine's argument list changed.
    !u   15 Feb 02 (ATP) Added MPI parallelization for fp
    !u   15 Jan 02 ascii version now labels site rho1,rho2,..
    !u   27 Aug 01 Extended to local orbitals.
    !u   17 May 01 Added ASA file I/O.  New argument list.
    !u   27 Apr 01 Added lbin switch
    !u   25 Jun 00 spin polarized
    !u   21 Apr 00 Adapted from nfp rw_rs
    !  ----------------------------------------------------------------------
    implicit none
    logical:: lbin=.true.
    integer::  nit , ifi , i_site,i_spec!,i_copy_size !mode=1 ,
    character*256:: fid=''
    ! ... Local parameters
    integer :: procid,master,mpipid,nproc
    integer :: i,i0,i1,i2,i3,i4,ib,ipr,iprint,ic,is,is0,isp,jb,k1,k2,k3, &
         igetss,jfi,k11,k21,k31,kmax,kmax0,l,lfoc,lfoc0,lmxa, & !kmxv,
         lmxa0,lmxb,lmxb0,lmxl,lmxl0,lmxr,lmxv,lmxv0,lrel0,n11,n21, &
         n31,nbas0,nspec0,nlml,nlml0,npan,npan0,nr,nr0,nsp0, &
         nxi,nat0,ibaug
    integer:: ngabc(3) , n1 , n2 , n3 , isw
    complex(8) ,allocatable :: h_zv(:)
    real(8) ,allocatable :: rwgt_rv(:)
    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    integer :: idmod(n0),idmoz(n0) !,lrs(10)
    logical :: isanrg,lfail,ltmp1,ltmp2,latvec,cmdopt,mlog !,lshear
    double precision :: a,a0,alat0,cof,eh,fac,qc,rfoc,rfoc0,rmt, &
         rmt0,rsma0,rsmr,rsmr0,rsmv,rsmv0,stc,sum,vfac,vol0,vs,vs1,z,z0
    double precision :: pnu(n0,2),pnz(n0,2),ql(n0,2*n0),pos(3)=99999, &
         force(3),plat0(3,3),qlat0(3,3), &
         exi(n0),hfc(n0,2),vec0(3),wk(100),rh,vrmax(2),pnus(n0,2), &
         pnzs(n0,2),dval
    character spid*8,spid0*8,fid0*68,line*20,msg*23,use*80,ignore*80, &
         msgw*17,datimp*32,usernm*32,hostnm*32,jobid*32,ffmt*32,ifmt*32
    integer:: fextg, i_dummy_fextg,ifile_handle,n
    character(*)::rwrw
    data vec0 /0d0,0d0,0d0/
    !     call tcn('iors')
    ! ... MPI setup
    nproc  = mpipid(0)
    procid = mpipid(1)
    master = 0
    mlog = cmdopt('--mlog',6,0,ignore)
    ngabc=lat_nabc
    ipr    = iprint()
    !      stdo   = lgunit(1)
    vs   =  2d0 !1.04d0 !version of rst file vs=2d0 at 2022-5-15. Only we support reading vs=1.04 only.
    msg  = '         File mismatch:'
    msgw = '         warning:'
    iors = -1
    line = 'header'
    call fftz30(n1,n2,n3,k1,k2,k3)
    ffmt = '(5f15.10)'
    ifmt = '(20i5)'
    npan = 1 !Hardwired for now
    if(ipr>0) write(stdo,"(/a)")' iors  : '//trim(rwrw)//' rst restart file (binary mesh density)'
    open(newunit=ifi,file='rst.'//trim(sname),form='unformatted')
    ! --- Input ---
    if (trim(rwrw)=='read') then
       jfi = ifi
       use    = '         use from  restart file:'
       ignore = '         ignore in restart file:'
       line = 'header'
       ! MPI check to see if at least 1st record can be read. Abort with error message if file is missing (lfail = .true.)
       lfail = .false.
       if (procid == master) then
          lfail = .true.
          read(jfi,end=996,err=996) vs1
          lfail = .false.
          rewind jfi
996       continue
       endif
       call mpibc1_logical(lfail,1,'iors_read error')
       if (lfail) goto 998
       if (procid == master) then
          read(jfi,end=998,err=998) vs1
          read(jfi) fid0
          read(jfi) datimp,usernm,hostnm
          read(jfi) nbas0,nat0,nsp0,npan0,lrel0,nspec0
          read(jfi) nit
          read(jfi) alat0,vol0
          if (ipr >= 40) write(stdo,710) trim(fid0),trim(usernm),trim(hostnm),trim(datimp)
          if (nat > nat0) then !  Number of real atoms may not increase
             if (isanrg(nat0, nat,nat,msg,'nat', .TRUE. )) goto 999
          elseif (nat==nat0 .AND. nbas/=nbas0) then
             write(stdo,ftox)' (warning) mismatch in nbas ... skipping sites'
             write(stdo,ftox)' expected nbas=',nbas,'but rst file has nbas=',nbas0
             !     OK to reduce nat (e.g. E.S.); these sites will be skipped
          elseif (nat < nat0) then
             write(stdo,*)'   (warning) rst mismatch in nat ... skipping sites'
          endif
          if (nsp0 < nsp) write(stdo,*)'   (warning) rst file not spin pol .. splitting spins'
          if(lrel0/=lrel) write(stdo,ftox) 'warning lrel_ctrl lrelrst=',lrel0,lrel
       endif
710    format(9x,'id -  ',a,/9x,'written by -  ',a,' on ',a,' at: ',a)
       call mpibc1_int(nbas0,1,'iors_nbas0')
       call mpibc1_int(nat0,1,'iors_nat0')
       call mpibc1_real(plat,9,'iors_plat')
       call mpibc1_int(nit,1,'iors_nit')
       !   --- Read smooth charge density ---
       allocate(osmrho(k1*k2*k3*nsp))
       osmrho=0d0
       line = 'smoothed density'
       if (procid == master) then
          read(jfi,err=999,end=999) n11,n21,n31
          if (n11 == n1 .AND. n21 == n2 .AND. n31 == n3) then
             read(jfi)osmrho
             if (nsp > nsp0) then
                n =k1*k2*k3
                osmrho(1:n   )=.5d0*osmrho(1:n)
                osmrho(1+n:n+n)= osmrho(1:n)
             endif
          else                 !... or read and remesh
             if (ipr >= 10) write(stdo,450) n11,n21,n31,n1,n2,n3
450          format(9x,'remesh density from  ',i4,'  *',i4,'  *',i4,'    to  ',i4,'  *',i4,'  *',i4)
             call fftz30(n11,n21,n31,k11,k21,k31)
             allocate(h_zv(k11*k21*k31*nsp))
             read(jfi)h_zv
             if (nsp > nsp0) then
                n  = k11*k21*k31
                h_zv(1:n    ) = .5d0 *h_zv(1:n)
                h_zv(1+n:n+n) = h_zv(1:n)
             endif
             call pshpr(50)
             i = 0
             if (n1 == 2*n11 .AND. n2 == 2*n21 .AND. n3 == 2*n31) i=3
             call chgmsh ( i , plat , nsp , n11 , n21 , n31 , k11 , k21 , &
                  k31 , h_zv , n1 , n2 , n3 , k1 , k2 , k3 , osmrho )
             call poppr
             if (allocated(h_zv)) deallocate(h_zv)
          endif
          !     ... If cell volume changed, scale smooth density to maintain charge
          vfac = vol0/vol
          if (dabs(vfac-1d0) > 1d-8) then
             if (ipr >= 10) write(stdo,460) vol0,vfac
460          format(9x,'volume changed from',f8.2,' :  scale smooth density by',f8.4)
             osmrho=osmrho*vfac
          endif
       endif
       call mpibc1_complex(osmrho, size(osmrho), 'iors_smrho' )
115    continue
       !   --- Read information related to dynamics ---
       if (procid == master) read(jfi) wk
       call mpibc1_real(wk,1,'iors:eferm')
       use=trim(use)//'use window,'
       call m_bndfp_ef_SET(wk(1)) !bz_ef00) !,bz_def00)
       line = 'site data' !Read atomic positions,forces,velocities ---
!       if (irs3 /= 0) then
!          ignore=trim(ignore)//'positions,'
!       else
!          use=trim(use)//' positions,'
!          if (nbas0 < nbas) then
!             call info0(1,0,0,'%9f oops ... cannot use file site positions (site mismatch)')
!             if (isanrg(nbas0, nbas,nbas,msg,'nbas', .TRUE. )) goto 999
!          endif
!       endif
       do ib = 1, nbas0
          if (procid == master) read(jfi,err=999,end=999) jb,pos,force
          !call mpibc1_real(pos,  3,'iors_pos')
          call mpibc1_real(force,3,'iors_force')
          if(ib > nbas) cycle
          !if(irs3 /= 0) pos=ssite(ib)%pos
          !ssite(ib)%pos  =pos
          !ssite(ib)%pos0 =pos
          ssite(ib)%force=force
       enddo
       !   --- Read information for local densities ---
!       if (irs5 /= 0) then
!          ignore=trim(ignore)//' pnu,'
!       else
          use=trim(use)//' pnu,'
!       endif
       if (ipr >= 10) then
          write(stdo,*)trim(use)
          write(stdo,*)trim(ignore)
       endif
       allocate(orhoat(3,nbas), v1pot(nbas),v0pot(nbas))
       ibaug = 0
       do  ib = 1, nbas
          ic=ssite(ib)%class
          is=ssite(ib)%spec !     is = -1 -> spec struc does not have these parameters
          if (is /= -1) then
             spid=slabl(is) !sspec(is)%name
             a=sspec(is)%a
             nr=sspec(is)%nr
             rmt=sspec(is)%rmt
             z=sspec(is)%z
             lmxa=sspec(is)%lmxa
             lmxl=sspec(is)%lmxl
             kmax=sspec(is)%kmxt
             if (lmxa == -1) goto 20
          endif
          ibaug = ibaug+1
          if (procid== master) then
             read(jfi) is0,spid0,lmxa0,lmxl0,nr0,rmt0,a0,z0,qc
             if (ib > nbas) goto 20
             if (ipr >= 40) then
                if (ib <= nat) write(stdo,380) ib,is0,spid0
                if (ib > nat) write(stdo,380) ib,is0,spid0, ' (skip)'
             endif
380          format('   atom',i4,'    species',i4,':',a:a)
             !     ... read(but don't use) extra info since record is present
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
          pnu=ssite(ib)%pnu
          pnz=ssite(ib)%pz
          idmod=0
          idmoz=0
          if (procid == master) then
             do  isp = 1, nsp0
                read(jfi) (pnu(l+1,isp), l=0,lmxa0)
                read(jfi) (pnz(l+1,isp), l=0,lmxa0)
                if (nsp > nsp0) then
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
             !       ... For backwards compatibility: prior versions wrote pnu for pnz
             do  l = 0, lmxa0+1
                if (pnu(l+1,isp) == mod(pnz(l+1,isp),10d0)) &
                     pnz(l+1,isp) = 0
             enddo
          enddo
          if (procid == master) then
             read(jfi) (idmod(l+1), l=0,lmxa0)
             read(jfi) (idmoz(l+1), l=0,lmxa0)
          endif
          !          pnus =sspec(is)%p
          !          pnzs =sspec(is)%pz
!          if  (irs5 == 0 ) then
             ssite(ib)%pnu=pnu
             ssite(ib)%pz=pnz
             !            sspec(is)%p=  pnu !int?
             !            sspec(is)%pz= pnz
             !     ... Verify lowest valence pnu compatible with file
             !            lfail = .false.
             !            ltmp1 = .false.
             !            ltmp2 = .false.
             !            do  i = 1, lmxa+1
             !              vec0(1) = mod(pnz(i,1),10d0)
             !              if (vec0(1) .eq. 0) vec0(1) = pnu(i,1)
             !              vec0(2) = mod(pnzs(i,1),10d0)
             !              if (vec0(2) .eq. 0) vec0(2) = pnus(i,1)
             !              ltmp1 = ltmp1 .or. int(pnu(i,1)) .ne. int(pnus(i,1))
             !              ltmp2 = ltmp2 .or.
             !     .        int(mod(pnz(i,1),10d0)) .ne. int(mod(pnzs(i,1),10d0))
             !              lfail = lfail .or. min(int(pnu(i,1)),int(vec0(1))) .ne.
             !     .        min(int(pnus(i,1)),int(vec0(2)))
             !            enddo
             if (ipr >= 20) write(stdo,203) ib,spid,'file pnu',(pnu(i,1), i=1,lmxa+1)
             !            if (ltmp1 .and. ipr.ge.20) write(stdo,204) 'given pnu is',(pnus(i,1), i=1,lmxa+1)
             if (ipr >= 20) write(stdo,203) ib,spid,'file pz ',(pnz(i,1), i=1,lmxa+1)
             !            if (ltmp2 .and. ipr.ge.20) write(stdo,204) 'given pz  is',(pnzs(i,1), i=1,lmxa+1)
203          format(9x,'site',i4,':',a,':',a,' is',8f6.2)
204          format(26x,a,8f6.2)
!             if (lfail .AND. irs5 == 1) then
!                call rx('iors: file''s pnu is incompatible with input')
!             endif
             nlml0 = (lmxl0+1)**2
             nlml = (lmxl+1)**2
             if (nr <= 0)   nr = nr0
             if (a <= 1d-6) a = a0
             if (procid == master) then
                call fsanrg(rmt0,rmt,rmt,1d-3,msg,'rmt',.true.)
                call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.false.)
                call fsanrg(z0,z,z,1d-6,msg,'z',.true.)
                call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
                lfail = isanrg(nr0,nr,nr,msgw,'nr',.false.)
                if (isanrg(lmxl,  0,lmxa,  msg,'lmxl', .FALSE. )) goto 999
                if (kmax0 /= kmax .AND. ipr >= 10) write(stdo,201) ib,spid,'kmax',kmax0,kmax
                if (lmxa0 /= lmxa .AND. ipr >= 10) write(stdo,201) ib,spid,'lmax',lmxa0,lmxa
201             format(9x,'site',i4,', species ',a,': augmentation ',a,' changed from',i2,' to',i2)
             endif
!          else
!             lmxl = lmxl0
!             nr = nr0
!          endif
          !!
          nlml0 = (lmxl0+1)**2
          nlml  = (lmxl+1)**2
          if (nr /= nr0) call rx('iors not set up to convert radial mesh')
          allocate(orhoat(1,ib)%v(nr*nlml*nsp)) !FP local densities rho1,rho2,rhoc and potentials v0, v1
          allocate(orhoat(2,ib)%v(nr*nlml*nsp))
          allocate(orhoat(3,ib)%v(nr*nsp))
          allocate(ssite(ib)%rv_a_ov0(nr*nsp))
          allocate(ssite(ib)%rv_a_ov1(nr*nsp))
          allocate(v0pot(ib)%v(nr*nsp))
          allocate(v1pot(ib)%v(nr*nsp))
          if (procid == master) then
             !             print *,'nnnnnnnnread',ib,nr*nlml*nsp,nr,nlml,nsp
             call readrho(ifi,nr,nlml0,nsp0,nlml,nsp,orhoat(1,ib)%v)
             call readrho(ifi,nr,nlml0,nsp0,nlml,nsp,orhoat(2,ib)%v)
             call readrho(ifi,nr,1,nsp0,1,nsp,orhoat(3,ib)%v)
             call readrhos(ifi,nr,nsp0,nsp,ssite(ib)%rv_a_ov0)
             call readrhos(ifi,nr,nsp0,nsp,ssite(ib)%rv_a_ov1)
             if(nlml0 > nlml .AND. ipr >= 10) write(stdo,202) ib,spid,'truncate',nlml0,nlml
             if(nlml0 < nlml .AND. ipr >= 10) write(stdo,202) ib,spid,'inflate',nlml0,nlml
202          format(9x,'site',i4,', species ',a,': ',a,' local density from nlm=',i3,' to',i3)
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
20        continue
       enddo
       if (isanrg(ibaug, nat,nat,  msg,'nat', .FALSE. )) goto 999
       !   --- Read data on free-atom core states and fit to fa density ---
       line = 'species data'
       do  30  is = 1, nspec
          a   =sspec(is)%a
          nr  =sspec(is)%nr
          lmxa=sspec(is)%lmxa
          if (lmxa == -1) goto 30
          if (procid == master) then
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
          !     ... FP core densities
          if (allocated(sspec(is)%rv_a_orhoc)) deallocate(sspec(is)%rv_a_orhoc)
          allocate(sspec(is)%rv_a_orhoc(abs(nr*nsp)))
          if (nr*nsp<0) sspec(is)%rv_a_orhoc(:)=0.0d0
          if (procid == master) then
             if (nr /= nr0) call rx('iors not set up to convert core radial mesh')
             read(jfi) sspec(is)%rv_a_orhoc !, nr * nsp0 , jfi )
             if (nsp > nsp0) then !spin-split core density
                i = nr
                call dscal ( i , 0.5d0 , sspec(is)%rv_a_orhoc , 1 )
                call dpscop ( sspec(is)%rv_a_orhoc , sspec(is)%rv_a_orhoc , i , 1 , 1 + i , 1d0  )
             endif
          endif
          call mpibc1_real(sspec(is)%rv_a_orhoc, size(sspec(is)%rv_a_orhoc), 'iors_rhoca'  )
          call dpzero(exi,n0)
          call dpzero(hfc,n0*2)
          if (procid == master) then
             read(jfi,err=999,end=999) rsmfa(is),nxi
             read(jfi,err=999,end=999) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp0)
             if (nsp > nsp0) then
                i = n0
                call dscal(i,0.5d0,hfc,1)
                call dpscop(hfc,hfc,i,1,1+i,1d0)
             endif
          endif
          call mpibc1_real(rsmfa(is),1,'iors_rsmfa')
          call mpibc1_int(nxi,1,'iors_nxi')
          call mpibc1_real(exi,nxi,'iors_exi')
          call mpibc1_real(hfc,nsp*nxi,'iors_hfc')
          sspec(is)%ctail=cof
          sspec(is)%etail=eh
          sspec(is)%stc=stc
          sspec(is)%nxi=nxi
          sspec(is)%exi=exi
          sspec(is)%chfa=hfc
!          sspec(is)%rsmfa=rsmfa
30     enddo

       !   ... Copy or rescale cores, in case foca was switched on or off
       do  ib = 1, nbas
          is = int(ssite(ib)%spec)
          a=sspec(is)%a
          nr=sspec(is)%nr
          rmt=sspec(is)%rmt
          lmxa=sspec(is)%lmxa
          lfoc=sspec(is)%lfoca
          qc=sspec(is)%qc
          if (lmxa == -1) goto 40
          if (lfoc > 0) then
             call dpcopy ( sspec(is)%rv_a_orhoc , orhoat( 3 , ib )%v , 1 , nr * nsp , 1d0 )
          else
             allocate(rwgt_rv(nr))
             call radwgt ( rmt , a , nr , rwgt_rv )
             call radsum ( nr , nr , 1 , nsp , rwgt_rv , orhoat( 3 , ib )%v , sum )
             fac = 1d0
             if (dabs(sum) > 1d-6) fac = qc/sum
             if (dabs(fac-1d0) > 1d-7 .AND. ipr >= 30) &
                  write(stdo,787) ib,qc,sum,fac
787          format(' fix core chg: ib=',i4,'  qc,sum,fac=',3f12.6)
             call dpcopy ( orhoat( 3 , ib )%v , orhoat( 3 , ib )%v, 1 , nr * nsp , fac )
             if (allocated(rwgt_rv)) deallocate(rwgt_rv)
          endif
40        continue
       enddo
       do i_site=1,nbas
          call mpibc1_s_site(ssite(i_site),'iors_ssite')
       enddo
       do i_spec=1,nspec
          call mpibc1_s_spec(sspec(i_spec),'iors_sspec')
       enddo
       ! --- Output for master---
    else
       if (procid /= master) then
          iors = 0
          call rx('iors: something wrong duplicated writing by cores?')
       endif
       jfi = ifi
       !        rewind jfi
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
       if (ipr >= 40) write(stdo,710) fid(1:i1), usernm(1:i3),hostnm(1:i4),datimp(1:i2)
721    format('----------------------- ',a,' -----------------------')
       write(jfi) vs
       write(jfi) fid0
       write(jfi) datimp,usernm,hostnm,jobid
       write(jfi) nbas,nat,nsp,npan,lrel,nspec
       write(jfi) nit
       write(jfi) alat,vol,plat
       write(jfi) n1,n2,n3
       write(jfi) osmrho
       !   --- Write information related to dynamics ---
       wk=1d99 !call dpzero(wk,100)
       wk(1)= eferm !sbz%ef !dummy ! we use wk(1) only wk(2:100) are dummy
       write(jfi) wk !call dpdump(wk,100,-jfi)
       do  110  ib = 1, nbas
          !pos  =ssite(ib)%pos
          force=ssite(ib)%force     !          vel  =ssite(ib)%vel
          write(jfi) ib,pos,force   !,vel
110    enddo
       !   --- Write information for local densities ---
       if (ipr >= 50) write(stdo,364)
       do  120  ib = 1, nbas
          ic=ssite(ib)%class
          is=ssite(ib)%spec
          spid=slabl(is) !sspec(is)%name
          a=sspec(is)%a
          nr=sspec(is)%nr
          rmt=sspec(is)%rmt
          z=sspec(is)%z
          qc=sspec(is)%qc
          idmod=idmodis(:,is) !sspec(is)%idmod
          !rsma=sspec(is)%rsma
          lmxa=sspec(is)%lmxa
          lmxl=sspec(is)%lmxl
          lmxb=sspec(is)%lmxb
          !kmxv=sspec(is)%kmxv
          rsmv=sspec(is)%rsmv
          kmax=sspec(is)%kmxt
          pnu=ssite(ib)%pnu
          pnz=ssite(ib)%pz
          if (lmxa == -1) goto 120
          write(jfi) is,spid,lmxa,lmxl,nr,rmt,a,z,qc ! Some extra info. lots of it useless or obsolete
          lmxr = 0
          lmxv = 0
          rsmr = 0
          write(jfi) rsma(is),rsmr,rsmv,lmxv,lmxr,lmxb,kmax !  ... Write augmentation data
          do  122  isp = 1, nsp
             write(jfi) (pnu(l+1,isp), l=0,lmxa)
             write(jfi) (pnz(l+1,isp), l=0,lmxa)
122       enddo
          write(jfi) (idmod(l+1), l=0,lmxa) !         Write for compatibility with nfp
          write(jfi) (idmod(l+1), l=0,lmxa)
          !     ... Write arrays for local density and potential
          nlml = (lmxl+1)**2
          !print *,'nnnnnnnwrite',ib,nr*nlml*nsp,nr,nlml,nsp
          write(jfi) orhoat( 1 , ib )%v !, nr , nlml , nlml , nsp  , nsp , lbin , -jfi )
          write(jfi) orhoat( 2 , ib )%v !, nr , nlml , nlml , nsp  , nsp , lbin , -jfi )
          write(jfi) orhoat( 3 , ib )%v !, nr , 1 , 1 , nsp , nsp   , lbin , -jfi )
          write(jfi) ssite(ib)%rv_a_ov0 !, nr , 1 , 1 , nsp , nsp , lbin , -jfi  )
          write(jfi) ssite(ib)%rv_a_ov1 !, nr , 1 , 1 , nsp , nsp , lbin , -jfi  )
          if (ipr >= 50) then
             write(stdo,349) ib,spid,lmxa,lmxl,rmt,nr,a, (pnu(l+1,1),l=0,lmxa)
             if (nsp == 2)  write(stdo,350) (pnu(l+1,2), l=0,lmxa)
          endif
120    enddo
       !   --- Write data on free-atom core states and fit to fa density ---
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
          !rsmfa = sspec(is)%rsmfa
          if (lmxa == -1) goto 130
          write(jfi) nr,a,qc,cof,eh,stc,lfoc,rfoc
          write(jfi) sspec(is)%rv_a_orhoc
          write(jfi) rsmfa(is),nxi
          write(jfi) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp)
130    enddo
    endif
349 format(i11,':',a4,2i2,f9.5,i5,f6.3,1x,8f6.3)
350 format(41x,8f6.3)
364 format(/9x,'ib:spc la ll   rmt     nr   a     pnu')
    iors = 0
    close(ifi)
    return
998 continue
    if (ipr > 0) write(stdo,'('' iors  : empty file ... nothing read'')')
    close(ifi)
    return
999 continue
    if (ipr > 0) write(stdo,'('' iors  : read failed in: '',a)') line
    close(ifi)
    !      call tcx('iors')
  end function iors

  subroutine readrho(ifi,nr,nlm0,nsp0,nlm,nsp,aout)!size controled read
    integer :: nr,nlm,nlm0,nsp0,nsp,ifi
    real(8):: aread(nr,nlm0,nsp0)
    real(8):: aout (nr,nlm,nsp)
    integer :: jfi,isp,ilm,i,nlmx
    double precision :: xx
    nlmx = min(nlm,nlm0)
    read(ifi) aread
    if(nsp0==2 .AND. nsp==1) then
       aout(1:nr,1:nlmx,1) = aread(1:nr,1:nlmx,1)+aread(1:nr,1:nlmx,2)
    elseif(nsp0==1 .AND. nsp==2) then
       aout(1:nr,1:nlmx,1)=.5d0*aread(1:nr,1:nlmx,1)
       aout(1:nr,1:nlmx,2)=.5d0*aread(1:nr,1:nlmx,1)
    elseif(nsp0==nsp) then
       aout(1:nr,1:nlmx,1:nsp)=aread(1:nr,1:nlmx,1:nsp)
    else
       call rx('readrho nsp wrong')
    endif
    aout(1:nr,nlmx+1:nlm0,1:nsp)=0d0
  end subroutine readrho

  subroutine readrhos(ifi,nr,nsp0,nsp,aout)
    integer :: nr,nsp0,nsp,ifi
    real(8):: aread(nr,nsp0)
    real(8):: aout (nr,nsp)
    read(ifi) aread
    aout(1:nr,1:nsp0)=aread(1:nr,1:nsp0)
    if(nsp==2 .AND. nsp0==1) aout(1:nr,2)=aread(1:nr,1)
  end subroutine readrhos

end module m_iors
