Cdelw1 var ['oqwk', '3', 'oqwk', 'nlml', 'real(8)', 'rv_w_', 'symvvl']
Cdelw1 var ['osym', '3', 'osym', 'nlml*nlml*nrclas', 'real(8)', 'rv_w_', 'symvvl']
Cdelw1 var ['og', '-1', 'undef', '0', 'undef', 'x', 'symvvl']
Cdelw1 var ['oag', '-1', 'undef', '0', 'undef', 'x', 'symvvl']
Cdelw1 undel [['undef', 'wref'], 'og', 'symvvl']
Cdelw1 undel [['undef', 'wref'], 'oag', 'symvvl']
Cgetarg...info...           structure ['mshvmt', 'lat', 'slat']
Cgetarg...info...           structure ['mshvmt', 'site', 'ssite']
Cgetarg...info...           structure ['mshvmt', 'spec', 'sspec']
Cgetarg...info...           structure ['symvvl', 'lat', 'slat']
Cgetarg...info...           structure ['symvvl', 'site', 'ssite']
Cgetarg...info...           structure ['symvvl', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['mshvmt', 'm_struc_def']
Cgetarg...info...           use_to_add ['symvvl', 'm_struc_def']
#define F90 1
      subroutine mshvmt(nbas,ssite,sspec,slat,ng,gv,
     .  kv,cv,k1,k2,k3,smpot,vval)
       
       use m_struc_def  !Cgetarg

C- Makes potential at MT surfaces given potential on a uniform mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    class spec pos
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: rmt lmxl
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat nabc nsgrp osymgr oag
Ci     Stored:    *
Ci     Passed to: *
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   cv    :work array
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   smpot :estat potential 
Co Outputs
Co   vval  :coffs to YL expansion of es potential at MT boundary
Co         :for each site, computed from the mesh density.
Cr Remarks
Cr   A PW exp(i.q.r) has a one-center expansion at radius r
Cr      sum_L C_L Y_L(r) where C_L = 4 pi i^l j_l(|rq|) Y_L(q)
Cr   Routine symvvl symmetrizes the vval generated here.
Cb Bugs
Cb   Possible to make ves for sites with lmxl=-1, which tells
Cb   value of ves at point.  However, vval doesn't have the 
Cb   space allocated.  So skip for now
Cu Updates
Cu   01 Jul 05 handle sites with lmxl=-1 
Cu   22 Aug 01 Newly created.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,nbas,ng,kv(ng,3)
Cgetarg       double precision ssite(1),sspec(1),gv(ng,3),slat(1),
Cgetarg      .  vval(1)
       real(8):: gv(ng,3) , vval(1) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex smpot(k1,k2,k3),cv(ng)
C ... Local parameters
      integer i,ib,is,lmxx,nlmx,iv0,lmxl,nlm,ngabc(3),
     .  n1,n2,n3,m,ilm,l,ipr
      double precision alat,pi,tpiba,tau(3),rmt,fac,plat(3,3)
C     double precision q(3),qlat(3,3)
      double complex vvali,fprli
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      parameter (lmxx=6, nlmx=(lmxx+1)**2)
Cki#error have #if F90 directive
#if F90
      double precision,allocatable:: phil(:,:),yl(:,:)
      double precision,allocatable:: gv2(:,:),agv(:),cgp(:),sgp(:)
#else
      double precision phil(ng,0:lmxx),yl(ng,nlmx)
      double precision gv2(ng,3),agv(ng),cgp(ng),sgp(ng)
#endif

      call tcn('mshvmt')
Cki#error have #if F90 directive
#if F90
      allocate(phil(ng,0:lmxx),yl(ng,nlmx))
      allocate(gv2(ng,3),agv(ng),cgp(ng),sgp(ng))
#endif
      call getpr(ipr)
      pi = 4d0*datan(1d0)
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg       call upack('lat alat plat nabc',slat,alat,plat,ngabc,0,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 

      tpiba = 2*pi/alat
      call gvgetf(ng,1,kv,k1,k2,k3,smpot,cv)

C --- YL(G)*G**l, agv=|G| for each g ---
      call dpcopy(gv,gv2,1,3*ng,tpiba)
      call ropyln(ng,gv2(1,1),gv2(1,2),gv2(1,3),lmxx,ng,yl,agv)
      do  i = 1, ng
        agv(i) = sqrt(agv(i))
      enddo

C --- For each ib in nbas, do ---
      iv0 = 0
      do  ib = 1, nbas
Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg         call upack('site spec pos',ssite,ib,is,tau,0,0)
         
         is=ssite(ib)%spec
         i_copy_size=size(ssite(ib)%pos) 
         call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1) 

Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg         call upack('spec rmt lmxl',sspec,is,rmt,lmxl,0,0)
         
         rmt=sspec(is)%rmt
         lmxl=sspec(is)%lmxl

        if (lmxl .eq. -1) goto 10
        nlm = (lmxl+1)**2
        if (nlm .gt. nlmx) call rxi('mshvmt: increase nlmx to',nlm)
C       Add a negligibly small amount to rmt to handle case rmt=0
        rmt = rmt+1d-32

C   --- j_l(|rmt*q|)/rmt**l for each G and l=0..lmax ---
C       Does not evolve correctly in the correct large r limit
        call ropbes(agv,rmt**2,lmxl,cgp,sgp,phil,ng)
C Patch for now        
C       do  i = 1, ng
C         call besslr(agv(i)**2*rmt**2,0,0,0,phil(i,0),yl)
C       enddo

C   ... Phases exp(-i.G.tau), fast version
C       call suphs0(plat,ng,gv,gv2)
C       call dinv33(plat,1,qlat,fac)
C       call dpzero(q,3)
C       call suphas(q,tau,ng,gv2,n1,n2,n3,qlat,cgp,sgp)
C   ... Phases calculated straightforwardly.  Fast enough not to matter.
        call dscal(3,alat,tau,1)
        do  i = 1, ng
          fac = -(tau(1)*gv2(i,1)+tau(2)*gv2(i,2)+tau(3)*gv2(i,3))
          cgp(i) = dcos(fac)
          sgp(i) = dsin(fac)
        enddo

C   --- Sum_G 4*pi*(i*rmt)**l j_l(|rmt*G|)/(rmt*G)**l YL(G) G**l ---
C       call dpzero(vval(iv0+1),nlm)
        ilm = 0
        fprli = 4*pi
        do  l  = 0, lmxl
          do  m = -l, l
            ilm = ilm+1
            vvali = 0
            do  i = 2, ng
              vvali = vvali + (phil(i,l)*yl(i,ilm))*
     .                        (cv(i)*dcmplx(cgp(i),-sgp(i)))
            enddo
            vval(ilm+iv0) = fprli*vvali
          enddo
          fprli = fprli*(0d0,1d0)*rmt
        enddo

C   ... Printout
C        if (ipr .gt. 0) then
C          do  ilm = 1, nlm
C            if (ilm .eq. 1) then
C              write(stdo,650) ib,ilm,vval(ilm+iv0)
C            elseif (dabs(vval(ilm+iv0)) .gt. 1d-6) then
C              write(stdo,651)    ilm,vval(ilm+iv0)
C            endif
C  650              format(i4,i6,2f12.6)
C  651                     format(4x,i6,f12.6)
C          enddo
C        endif

        iv0 = iv0 + nlm
   10   continue
      enddo
Cki#error have #if F90 directive
#if F90
      deallocate(phil,yl)
      deallocate(gv2,agv,cgp,sgp)
#endif
      call tcx('mshvmt')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine mshvmt 


      subroutine symvvl(nbas,ssite,sspec,slat,vval,vrmt)
       
       use m_struc_def  !Cgetarg

C- Symmetrizes the potential at the MT boundary.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl rmt
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat nsgrp osymgr oag alat nabc
Co Outputs
Cio  vval  :On input,  unsymmetrized potential
Cio        :On output, elements of potential for sites in the same
Cio        :class are symmetrized.
Co Outputs
Co   vrmt  :spherical average of potential (i.e. Y0*vval(l=0)) returned
Co         :for each site.
Cr Remarks
Cr   This routine symmetrizes any vector of the same structure as vval.
Cb Bugs
Cb   Possible to make ves for sites with lmxl=-1, which tells
Cb   value of ves at point.  However, vval doesn't have the 
Cb   space allocated.  So skip for now
Cu Updates
Cu   01 Jul 05 handle sites with lmxl=-1 
Cu   23 Aug 01 Newly created.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas
Cgetarg       double precision ssite(1),sspec(1),slat(1),vval(1),vrmt(nbas)
       real(8):: vval(1) , vrmt(nbas) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

C ... Local parameters
      integer ic,ib,ilm,mxint,nclass,ipa(nbas),nrclas,stdo,iv0
      integer ipc(nbas),ips(nbas),lmxl(nbas)
      double precision pos(3,nbas),posc(3,nbas),plat(3,3),pi,y0
Cdelw1 do not change  [['undef', 'wref'], 'og', 'symvvl']
Cdelw1 do not change  [['undef', 'wref'], 'oag', 'symvvl']
Cdelw1       integer igetss,nlml,lgunit,ipr,jpr,ngrp,og,oag,oqwk,osym,nn,iclbas
       integer:: igetss , nlml , lgunit , ipr , jpr , ngrp , og , oag 
     . , nn , iclbas 
       real(8) ,allocatable :: rv_w_oqwk(:)
       real(8) ,allocatable :: rv_w_osym(:)

C ... Heap
      integer w(1)
      common /w/ w

      stdo = lgunit(1)
      call getpr(ipr)
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           integer :: ngrp , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: og , integer(8) :: slat%osymgr 1
Cgetarg...info...           integer :: oag , integer(8) :: slat%oag 1
Cgetarg       call upack('lat plat nsgrp osymgr oag',slat,plat,ngrp,og,oag,0)
       
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       ngrp=slat%nsgrp
Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 254
Cdelw1 do not change og because of ['undef', 'wref']
       og=slat%osymgr
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 255
Cdelw1 do not change oag because of ['undef', 'wref']
       oag=slat%oag

Cgetarg       call spackv(10,'site class',ssite,1,nbas,ipc)
       i_copy_size=1; 
       do i_spackv=1,nbas 
       call spackv_array_copy_i8_i('u',ssite(i_spackv)%class,i_copy_size,i_spackv+1-1,ipc)
       enddo

Cgetarg       call spackv(10,'site spec', ssite,1,nbas,ips)
       i_copy_size=1; 
       do i_spackv=1,nbas 
       call spackv_array_copy_i8_i('u',ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,ips)
       enddo

Cgetarg       call spackv(10,'site pos',ssite,1,nbas,pos)
       i_copy_size=size(ssite(1)%pos)
       do i_spackv=1,nbas 
       call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,pos)
       enddo

      nclass = mxint(nbas,ipc)
      do  ib = 1, nbas
Cgetarg         lmxl(ib) = igetss('spec lmxl',ips(ib),sspec)
         lmxl ( ib ) = int(sspec(ips(ib))%lmxl) 

      enddo

      do  ic = 1, nclass

C   ... Make nrclas,ipa,posc
        call psymr0(lmxl,ic,nbas,ipc,pos,posc,ipa,nrclas)
        if (nrclas .gt. 0) then
        ib = iclbas(ic,ipc)
        if (lmxl(ib) .gt. -1) then
        nlml = (lmxl(ib)+1)**2
        if (ipr .ge. 50) write(stdo,800) ic,nrclas,nlml
  800   format(' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)

Cdelw1         call defrr(oqwk,nlml)
         allocate(rv_w_oqwk(nlml))
Ctakao_ZeroClear_NotRequiered          if (nlml<0) rv_w_oqwk(:)=0.0d0

Cdelw1         call defrr(osym,nlml*nlml*nrclas)
         allocate(rv_w_osym(nlml*nlml*nrclas))
Ctakao_ZeroClear_NotRequiered          if (nlml*nlml*nrclas<0) rv_w_osym(:)=0.0d0

Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 295
Cdelw1 do not change og because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 295
Cdelw1 do not change oag because of ['undef', 'wref']
Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 295
Cdelw1 do not change og because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 295
Cdelw1 do not change oag because of ['undef', 'wref']
Cdelw1         call symqmp(nrclas,nlml,nlml,plat,posc,ngrp,w(og),w(oag),
Cdelw1      .    w(oqwk),ipa,w(osym),vval,nn)
         call symqmp ( nrclas , nlml , nlml , plat , posc , ngrp , w ( 
     .   og ) , w ( oag ) , rv_w_oqwk , ipa , rv_w_osym , vval , nn ) 

Cdelw1 w_varlist remains: oqwk osym 
Cdelw1 w_varlistundel: [og oag]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: [oqwk osym]
        if (allocated(rv_w_osym)) deallocate(rv_w_osym)
        if (allocated(rv_w_oqwk)) deallocate(rv_w_oqwk)

       endif
       endif
      enddo

C ... Extract vrmt = l=0 term for each site, and printout
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      if (ipr .ge. 45) write(stdo,221)
  221 format(/' site class  ilm      vval',6x,'ves(rmax)')
      iv0 = 0
      do  ib = 1, nbas
        if (lmxl(ib) .eq. -1) goto 10

        nlml = (lmxl(ib)+1)**2
        vrmt(ib) = vval(1+iv0)*y0

C   ... Printout
        ic = ipc(ib)
        jpr = 0
        if (ipr .gt. 60) jpr = 2
        if (ib .eq. iclbas(ic,ipc)) then
          if (ipr .ge. 45) jpr = 1
          if (ipr .ge. 50) jpr = 2
        endif
        if (jpr .gt. 0) then
          do  ilm = 1, nlml
            if (ilm .eq. 1) then
              write(stdo,650) ib,ic,ilm,vval(ilm+iv0),vrmt(ib)
            elseif (dabs(vval(ilm+iv0)) .gt. 1d-6  .and. jpr.gt.1) then
              write(stdo,651)    ilm,vval(ilm+iv0)
            endif
  650       format(i4,2i6,2f12.6)
  651       format(10x,i6,f12.6)
          enddo
        endif

        iv0 = iv0 + nlml
   10  continue
      enddo

Cgetarg       end

       end subroutine symvvl 


