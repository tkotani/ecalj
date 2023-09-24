!> I/O rst file. charge density and sspec are readin
module m_iors
  use m_struc_def
  use m_lgunit,only:stdo
  public iors
  private
contains
  integer function iors(nit,rwrw) 
    use m_density,only: osmrho, orhoat,v1pot,v0pot,pnuall,pnzall,eferm !Main I/O. these are allocated. In addition sspc is written
    use m_supot,only: n1,n2,n3
    use m_struc_func,only: mpibc1_s_spec 
    use m_lmfinit,only: alat=>lat_alat,nsp,lrel,nl,ispec,sspec=>v_sspec, nbas,nat,nspec,n0, idmodis=>idmod,slabl,readpnu
    use m_lattic,only: plat=>lat_plat,vol=>lat_vol,qlat=>lat_qlat
    use m_ext,only:sname
    use m_ftox
    use m_chgmsh,only:chgmsh
    !! I/O data
    !!     smrho, rhoat
    !      sspec:
    !!           a nr rmt z lmxa lmxl kmxt p pz lfoca qc orhoc idmod
    !!           rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
    !!           chfa rsmfa
    !!
    !i   nbas  :size of basis
    !i   nat   :number atoms in basis with augmentation sites
    !i         :Note: if nat<nbas, there is a requirement that
    !i         :lmxa>-1 for nat sites, and
    !i         :and lmxa=-1 for nbas-nat sites
    !!
    !r Remarks
    !r   The density consists of a smooth part (smrho) plus nbas atom-centered densities inside the MT spheres.
    !r   Their sum is the full charge density. 
    !r   The local density is represented as the difference of the two valence components in orhoat, plus the core density.
    !r   Density in the MT spheres:
    !r      mesh parameters rmt,nr,a;
    !r      total density rho (times r**2) to lmxl;
    !r      a spherical potential v0 defining the wave functions within rmt
    !r      pnu and idmod to lmxa
    !r      NOTE: on input, arrays for rhoat and v0 are allocated here
    !r   Smooth density
    !r      real part of complex*16 array smrho contains the density
    !r      n1,n2,n3 are the dimensions of the mesh.
    !r   Additional information stored:
    !r      fid: file identifier, a string of length 64 characters or less.
    !r      parameters relating to coordinates and molecular dynamics.
    !r   On input, iors tries to transform data format where needed:
    !r      lmxl incompatible: pad or truncate l-components
    !r      FT mesh changed: map onto new mesh
    implicit none
    integer::  nit , ifi , i_site,i_spec!,i_copy_size !mode=1 ,
    character*256:: fid=''
    integer :: procid,master,mpipid,nproc
    integer :: i,i0,i1,i2,i3,i4,ib,ipr,iprint,ic,is,is0,isp,jb, &
         igetss,jfi,kmax,kmax0,l,lmxa, & !kmxv,lfoc,lfoc0
         lmxa0,lmxb,lmxb0,lmxl,lmxl0,lmxr,lmxv,lmxv0,lrel0,n11,n21, &
         n31,nbas0,nspec0,nlml,nlml0,npan,npan0,nr,nr0,nsp0, &
         nxi,nat0,ibaug
    integer:: isw
    complex(8) ,allocatable :: h_zv(:)
    real(8) ,allocatable :: rwgt_rv(:)
    integer :: idmod(n0),idmoz(n0) 
    logical :: isanrg,lfail,ltmp1,ltmp2,latvec,skiprstpnu!,cmdopt0 !,lshear
    double precision :: a,a0,alat0,cof,eh,fac,qc,rmt, & !rfoc,rfoc0
         rmt0,rsma0,rsmv0,stc,sum,vfac,vol0,vs,vs1,z,z0
    real(8),pointer:: pnu(:,:),pnz(:,:)
    real(8):: ql(n0,2*n0),pos(3)=99999, &
         forcexxx(3)=999,plat0(3,3),qlat0(3,3), &
         exi(n0),hfc(n0,2),vec0(3),wk(100),rh,vrmax(2),pnus(n0,2), &
         pnzs(n0,2),dval,rsmfa,eferm0
    character spid*8,spid0*8,fid0*68,line*20,msg*23,use*80,ignore*80, &
         msgw*17,datimp*24,usernm*32,hostnm*32,jobid*32,ffmt*32,ifmt*32
    integer:: fextg, i_dummy_fextg,n
    character(*)::rwrw
    data vec0 /0d0,0d0,0d0/
    nproc  = mpipid(0)
    procid = mpipid(1)
    master = 0
    ipr    = iprint()
    vs   =  2d0 !1.04d0 !version of rst file vs=2d0 at 2022-5-15. Only we support reading vs=1.04 only.
    msg  = '         File mismatch:'
    msgw = '         warning:'
    iors = -1
    line = 'header'
    ffmt = '(5f15.10)'
    ifmt = '(20i5)'
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
          read(jfi) !fid0
          read(jfi) datimp,usernm,hostnm
          read(jfi) nbas0,nat0,nsp0 !,npan0,lrel0,nspec0
          read(jfi) nit
          read(jfi) alat0,vol0
          if(ipr >= 40) write(stdo,710) trim(usernm),trim(hostnm),trim(datimp)
          if(nbas/=nbas0) then
             write(stdo,ftox)' (warning) mismatch in nbas ... skipping sites'
             write(stdo,ftox)' expected nbas=',nbas,'but rst file has nbas=',nbas0
          endif
          if (nsp0 < nsp) write(stdo,*)'   (warning) rst file not spin pol .. splitting spins'
       endif
710    format(/9x,'written by -  ',a,' on ',a,' at: ',a)
       call mpibc1_int(nbas0,1,'iors_nbas0')
       call mpibc1_int(nat0,1,'iors_nat0')
       call mpibc1_real(plat,9,'iors_plat')
       call mpibc1_int(nit,1,'iors_nit')
       !   --- Read smooth charge density ---
       allocate(osmrho(n1*n2*n3,nsp))
       osmrho=0d0
       line = 'smoothed density'
       if (procid == master) then
          read(jfi,err=999,end=999) n11,n21,n31
          if (n11 == n1 .AND. n21 == n2 .AND. n31 == n3) then
             n =n1*n2*n3
             read(jfi) osmrho(1:n1*n2*n3,1:nsp0)
             if (nsp > nsp0) then
                osmrho(1:n,1)=.5d0*osmrho(1:n,1)
                osmrho(1:n,2)= osmrho(1:n,1)
             endif
          else                 !... or read and remesh
             if (ipr >= 10) write(stdo,450) n11,n21,n31,n1,n2,n3
450          format(9x,'remesh density from  ',i4,'  *',i4,'  *',i4,'    to  ',i4,'  *',i4,'  *',i4)
             allocate(h_zv(n11*n21*n31*nsp))
             read(jfi)h_zv(1:n11*n21*n31*nsp0)
             if (nsp > nsp0) then
                n  = n11*n21*n31
                h_zv(1:n    ) = .5d0 *h_zv(1:n)
                h_zv(1+n:n+n) = h_zv(1:n)
             endif
             call pshpr(50)
             i = 0
             if (n1 == 2*n11 .AND. n2 == 2*n21 .AND. n3 == 2*n31) i=3
             call chgmsh ( i , plat , nsp , n11 , n21 , n31 , n11,n21,n31 , h_zv , n1 , n2 , n3 , n1 , n2 , n3 , osmrho )
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
       if (procid == master) read(jfi) eferm0
       call mpibc1_real(eferm0,1,'iors:eferm')
       use=trim(use)//'use window,'
       eferm=eferm0
       line = 'site data' !Read atomic positions,forcexxxs,velocities ---
       do ib = 1, nbas0 
          if (procid == master) read(jfi)
       enddo
       !   --- Read information for local densities ---
       use=trim(use)//' pnu,'
       if (ipr >= 10) then
          write(stdo,*)trim(use)
          write(stdo,*)trim(ignore)
       endif
       allocate(orhoat(3,nbas), v1pot(nbas),v0pot(nbas))
       ibaug = 0
       do  ib = 1, nbas
          is=ispec(ib) !  is = -1 -> spec struc does not have these parameters
          if (is /= -1) then
             spid=slabl(is)
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
             read(jfi) lmxv0,lmxr,lmxb0,kmax0
          endif
          call mpibc1_int(lmxa0,1,'iors_lmxa0')
          call mpibc1_int(lmxl0,1,'iors_lmxl0')
          call mpibc1_int(nr0,1,'iors_nr0')
          call mpibc1_real(a0,1,'iors_a0')
          call mpibc1_real(qc,1,'iors_qc')
          if (is == -1 ) call rx('iors: need check for is==-1')
          if(readpnu) then 
             read(jfi)
             read(jfi)
          else   
             pnu=>pnuall(:,:,ib)
             pnz=>pnzall(:,:,ib)
             if (procid == master) then
                read(jfi) ((pnu(l+1,isp), l=0,lmxa0),isp=1,nsp0)
                read(jfi) ((pnz(l+1,isp), l=0,lmxa0),isp=1,nsp0)
                do  isp = 1, nsp0
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
                   if (pnu(l+1,isp) == mod(pnz(l+1,isp),10d0)) pnz(l+1,isp) = 0
                enddo
             enddo
             pnuall(:,1:nsp,ib)=pnu(:,1:nsp)
             pnzall(:,1:nsp,ib)=pnz(:,1:nsp)
             if (ipr >= 20) write(stdo,203) ib,spid,'file pnu',(pnu(i,1), i=1,lmxa+1)
             if (ipr >= 20) write(stdo,203) ib,spid,'file pz ',(pnz(i,1), i=1,lmxa+1)
          endif
       
          idmod=0
          idmoz=0
          if (procid == master) then
             read(jfi) (idmod(l+1), l=0,lmxa0)
             read(jfi) (idmoz(l+1), l=0,lmxa0)
          endif
203       format(9x,'site',i4,':',a,':',a,' is',8f6.2)
204       format(26x,a,8f6.2)
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
201          format(9x,'site',i4,', species ',a,': augmentation ',a,' changed from',i2,' to',i2)
          endif
          nlml0 = (lmxl0+1)**2
          nlml  = (lmxl+1)**2
          if (nr /= nr0) call rx('iors not set up to convert radial mesh')
          allocate(orhoat(1,ib)%v(nr*nlml*nsp)) !FP local densities rho1,rho2,rhoc and potentials v0, v1
          allocate(orhoat(2,ib)%v(nr*nlml*nsp))
          allocate(orhoat(3,ib)%v(nr*nsp))
          allocate(v0pot(ib)%v(nr*nsp))
          allocate(v1pot(ib)%v(nr*nsp))
          if (procid == master) then
             call readrho(ifi,nr,nlml0,nsp0,nlml,nsp,orhoat(1,ib)%v)
             call readrho(ifi,nr,nlml0,nsp0,nlml,nsp,orhoat(2,ib)%v)
             call readrho(ifi,nr,1,nsp0,1,nsp,orhoat(3,ib)%v)
             call readrhos(ifi,nr,nsp0,nsp,v0pot(ib)%v) !ssite(ib)%rv_a_ov0)
             call readrhos(ifi,nr,nsp0,nsp,v1pot(ib)%v)  !ssite(ib)%rv_a_ov1)
             if(nlml0 > nlml .AND. ipr >= 10) write(stdo,202) ib,spid,'truncate',nlml0,nlml
             if(nlml0 < nlml .AND. ipr >= 10) write(stdo,202) ib,spid,'inflate',nlml0,nlml
202          format(9x,'site',i4,', species ',a,': ',a,' local density from nlm=',i3,' to',i3)
          endif
          call mpibc1_real( orhoat(1,ib)%v, size(orhoat(1,ib)%v), 'iors_rhoat(1)' )
          call mpibc1_real( orhoat(2,ib)%v, size(orhoat(2,ib)%v), 'iors_rhoat(2)' )
          call mpibc1_real( orhoat(3,ib)%v, size(orhoat(3,ib)%v), 'iors_rhoat(3)' )
          call mpibc1_real( v0pot(ib)%v,size(v0pot(ib)%v) , 'iors_v0' )
          call mpibc1_real( v1pot(ib)%v,size(v1pot(ib)%v) , 'iors_v1' )
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
             read(jfi,err=999,end=999) nr0,a0,qc,cof,eh,stc !,lfoc0 !,rfoc0
             lfail = isanrg(nr0,nr,nr,msgw,'nr',.false.)
             call fsanrg(a0,a,a,0d-9,msg,'spec a',.true.)
          endif
          call mpibc1_real(qc,1, 'iors_qc')
          call mpibc1_real(cof,1,'iors_cof')
          call mpibc1_real(eh,1, 'iors_eh')
          call mpibc1_real(stc,1,'iors_stc')
          !     ... FP core densities
          if (allocated(sspec(is)%rv_a_orhoc)) deallocate(sspec(is)%rv_a_orhoc)
          allocate(sspec(is)%rv_a_orhoc(nr*nsp))
          if (procid == master) then
             if (nr /= nr0) call rx('iors not set up to convert core radial mesh')
             read(jfi) sspec(is)%rv_a_orhoc(1:nr*nsp0) !, nr * nsp0 , jfi )
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
             read(jfi,err=999,end=999) rsmfa,nxi
             read(jfi,err=999,end=999) ((exi(i),hfc(i,isp),i=1,nxi),isp=1,nsp0)
             if (nsp > nsp0) then
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
          sspec(is)%qc=qc
30     enddo
       !   ... Copy or rescale cores, in case foca was switched on or off
       do  ib = 1, nbas
          is = ispec(ib)
          a=sspec(is)%a
          nr=sspec(is)%nr
          rmt=sspec(is)%rmt
          lmxa=sspec(is)%lmxa
          qc=sspec(is)%qc
          if (lmxa == -1) goto 40
          if (sspec(is)%lfoca > 0) then
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
       do i_spec=1,nspec
          call mpibc1_s_spec(sspec(i_spec),'iors_sspec')
       enddo
    endif
!=======================================================================    
    if(rwrw=='write') then ! --- Output for master---
       if (procid /= master) then
          iors = 0
          call rx('iors: something wrong duplicated writing by cores?')
       endif
       jfi = ifi
!       fid0 = fid
       jobid = sname 
       call ftime(datimp)
       hostnm = ' '
       usernm = ' '
       call get_environment_variable('HOST',hostnm)
       call get_environment_variable('USER',usernm)
       if (ipr >= 40) write(stdo,710) trim(usernm),trim(hostnm),trim(datimp) !trim(fid), 
721    format('----------------------- ',a,' -----------------------')
       write(jfi) vs
       write(jfi) !fid0
       write(jfi) datimp,usernm,hostnm,jobid
       write(jfi) nbas,0,nsp,0,0,0 !lrel,nspec
       write(jfi) nit
       write(jfi) alat,vol,plat
       write(jfi) n1,n2,n3
       write(jfi) osmrho
       write(jfi) eferm 
       do ib = 1, nbas
          write(jfi)!forcexxx !ssite(ib)%force
       enddo
       if (ipr >= 50) write(stdo,364)
       do  120  ib = 1, nbas
          is=ispec(ib)   
          spid=slabl(is) 
          a=sspec(is)%a
          nr=sspec(is)%nr
          rmt=sspec(is)%rmt
          z=sspec(is)%z
          qc=sspec(is)%qc
          idmod=idmodis(:,is)
          lmxa=sspec(is)%lmxa
          lmxl=sspec(is)%lmxl
          lmxb=sspec(is)%lmxb
          kmax=sspec(is)%kmxt
          pnu=>pnuall(:,1:nsp,ib)
          pnz=>pnzall(:,1:nsp,ib)
          if (lmxa == -1) cycle
          write(jfi) is,spid,lmxa,lmxl,nr,rmt,a,z,qc ! Some extra info. lots of it useless or obsolete
          lmxr = 0
          lmxv = 0
          write(jfi) lmxv,lmxr,lmxb,kmax !  ... Write augmentation data
          write(jfi) ((pnu(l+1,isp), l=0,lmxa),isp=1,nsp)
          write(jfi) ((pnz(l+1,isp), l=0,lmxa),isp=1,nsp)
          write(jfi) (idmod(l+1), l=0,lmxa) !         Write for compatibility with nfp
          write(jfi) (idmod(l+1), l=0,lmxa)
          nlml = (lmxl+1)**2 !     ... Write arrays for local density and potential
          write(jfi) orhoat( 1 , ib )%v 
          write(jfi) orhoat( 2 , ib )%v 
          write(jfi) orhoat( 3 , ib )%v 
          write(jfi) v0pot(ib)%v !ssite(ib)%rv_a_ov0 
          write(jfi) v1pot(ib)%v !ssite(ib)%rv_a_ov1 
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
          cof =sspec(is)%ctail
          eh  =sspec(is)%etail
          stc =sspec(is)%stc
          nxi =sspec(is)%nxi
          exi=sspec(is)%exi
          hfc=sspec(is)%chfa
          rsmfa = sspec(is)%rsmfa
          if (lmxa == -1) cycle
          write(jfi) nr,a,qc,cof,eh,stc ,0,0 !lfoc,rfoc
          write(jfi) sspec(is)%rv_a_orhoc
          write(jfi) rsmfa,nxi
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
