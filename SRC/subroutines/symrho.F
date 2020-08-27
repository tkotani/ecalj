      subroutine symrho ( ssite , sspec ,  lf , smrho , sv_p_orhoat !slat ,
     .     , qbyl , hbyl , f )
      use m_work,only: iv_a_oistab , rv_a_osymgr, rv_a_oag, rv_a_ogv , iv_a_oips0 , zv_a_obgv ,iv_a_okv
      use m_struc_def  !Cgetarg
      use m_lmfinit,only: globalvariables
c      use m_globalvariables
      use m_lmfinit,only:lat_nabc
C- Symmetrize charge density and related quantities
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    class pos
Ci     Passed to: spackv
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl lmxa nr
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy plat qlat oistab nsgrp osymgr oag
Ci     Stored:    *
Ci     Passed to: *
Ci   lf    :>0 symmetrize forces
Ci Inputs/Outputs
Cio  smrho :smooth density
Cio        :Symmetrized on output
Cio  orhoat:vector of offsets containing site density
Cio        :Symmetrized on output
Cio  qbyl  :site- and l-decomposed charges
Cio        :Symmetrized on output
Cio  hbyl  :site- and l-decomposed one-electron energies
Cio        :Symmetrized on output
Cio  f     :forces
Cio        :Symmetrized on output
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   19 Jun 00 Packaged from nfp symrat.f and symsmr.f
C ----------------------------------------------------------------------
      implicit none
      integer:: lf
      type(s_rv1) :: sv_p_orhoat(3,*)
      real(8):: f(*) , qbyl(*) , hbyl(*)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
c      type(s_lat)::slat
      double complex smrho(*)
      integer nsp,nbas,ngabc(3),n1,n2,n3,k1,k2,k3,nglob,iprint,i_copy_size
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      call tcn('symrho')
      call info(30,1,0,' Symmetrize density..',0,0)
      nbas = globalvariables%nbas
      nsp = globalvariables%nsp
c      i_copy_size=size(lat_nabc)
c     call icopy(i_copy_size,lat_nabc,1,ngabc,1)
      ngabc=lat_nabc
      call fftz30(n1,n2,n3,k1,k2,k3)
      call symrat ( ssite , sspec , nbas , nsp , lf , sv_p_orhoat !slat , 
     . , qbyl , hbyl , f )
      if ( iprint ( ) .gt. 50 ) call prrhat ( nbas , ssite , sspec 
     . , sv_p_orhoat )
      call symsmr(nsp,k1,k2,k3,smrho) !slat,
      call tcx('symrho')
      end subroutine symrho

      subroutine symrat ( ssite , sspec , nbas , nsp , lf , !slat , 
     . sv_p_orhoat , qbyl , hbyl , f )
      use m_struc_def  !Cgetarg
      use m_work,only: iv_a_oistab , rv_a_osymgr, rv_a_oag
      use m_lmfinit,only: lat_plat,lat_qlat,lat_nsgrp
C- Symmetrize the atomic charge densities and the forces.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    class pos
Ci     Passed to: spackv
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl lmxa nr
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy plat qlat oistab nsgrp osymgr oag
Ci     Stored:    *
Ci     Passed to: *
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lf    :>0 symmetrize forces
Ci Inputs/Outputs
Cio  orhoat:vector of offsets containing site density
Cio        :Symmetrized on output
Cio  qbyl  :site- and l-decomposed charges
Cio        :Symmetrized on output
Cio  hbyl  :site- and l-decomposed one-electron energies
Cio        :Symmetrized on output
Cio  f     :forces
Cio        :Symmetrized on output
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer:: nbas , nsp , n0 , lf
      type(s_rv1) :: sv_p_orhoat(3,nbas)
      parameter (n0=10)
      real(8):: f(3,nbas) , qbyl(n0,nsp,nbas) , hbyl(n0,nsp,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
c      type(s_lat)::slat
      integer:: stdo , ib0 , ic , ipr , iprint , is , lgunit , lmxa 
     .       , lmxl , nclass , ngrp , nlml , nlmx , nr , nrclas , igetss , 
     .       mxint , ival,i_copy_size,i_spackv
      integer ,allocatable :: ipa_iv(:)
      integer ,allocatable :: ipc_iv(:)
      integer ,allocatable :: ips_iv(:)
      real(8) ,allocatable :: pos_rv(:)
      real(8) ,allocatable :: pos0_rv(:,:)
      real(8) ,allocatable :: rho_rv(:)
      real(8) ,allocatable :: sym_rv(:)
      double precision plat(3,3),qlat(3,3)

      call tcn('symrat')
      stdo = lgunit(1)
      ipr = iprint()
c      i_copy_size=size(lat_plat)
c      call dcopy(i_copy_size,lat_plat,1,plat,1)
c      i_copy_size=size(lat_qlat)
c     call dcopy(i_copy_size,lat_qlat,1,qlat,1)
      plat=lat_plat
      qlat=lat_qlat
      ngrp=lat_nsgrp
      allocate(ips_iv(nbas))
      allocate(ipc_iv(nbas))
      allocate(pos0_rv(3,nbas))
C ... Separate sites into symmetry classes
C      call symcls(nbas,ssite,sspec,slat,ngrp,w(og),w(oag),
C     .  w(oips),w(opos0),nclass,w(oipc))
C      if (ipr .ge. 40) write(stdo,300) nclass
C  300 format(/' symrat: number of symmetry classes is',i3)

c      i_copy_size=1;
      do i_spackv=1,nbas
c        call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %class 
c     .  , i_copy_size , i_spackv + 1 - 1 , ipc_iv )
        ipc_iv(i_spackv)= ssite( i_spackv)%class
      enddo
c      i_copy_size=size(ssite(1)%pos)
      do i_spackv=1,nbas
c        call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
c     .  , i_copy_size , i_spackv + 1 - 1 , pos0_rv )
        pos0_rv(:,i_spackv)= ssite(i_spackv)%pos 
      enddo
      nclass = mxint ( nbas , ipc_iv )
C --- Start loop over classes ---
      allocate(ipa_iv(nbas))
      allocate(pos_rv(3*nbas))
      do  ic = 1, nclass
        call psymr0 ( - 2 , ic , nbas , ipc_iv , pos0_rv , pos_rv 
     .  , ipa_iv , nrclas )
        if (nrclas .gt. 0) then
          ib0 = ival ( ipa_iv , 1 )
          is = int(ssite(ib0)%spec)
          lmxl=sspec(is)%lmxl
          lmxa=sspec(is)%lmxa
          nr=sspec(is)%nr
          nlml = (lmxl+1)**2
          if (ipr .ge. 40) write(stdo,800) ic,nrclas,nlml
  800     format(/' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)

C   ... Make the projectors; make to at least to l=1 for forces
          nlmx = max0(nlml,4)
          allocate(sym_rv(nlmx*nlmx*nrclas))
          call symprj ( nrclas , nlmx , ngrp , nbas , iv_a_oistab , rv_a_osymgr
     .       , rv_a_oag , plat , qlat , pos_rv , sym_rv )
C   ... Apply the projectors to rhoat
          if (lmxl .gt. -1) then
            allocate(rho_rv(nr*nlml*nsp))
            call psymr1 ( nrclas , ipa_iv , nr , nlml , nsp , nlmx , sym_rv 
     .       , rho_rv , sv_p_orhoat , 1 )
            call psymr1 ( nrclas , ipa_iv , nr , nlml , nsp , nlmx , sym_rv 
     .       , rho_rv , sv_p_orhoat , 2 )
C   ... Symmetrize site charges and eigval sum
            call psymrq ( nrclas , nsp , ipa_iv , lmxa , qbyl , hbyl )
          endif
C   ... Symmetrize the forces
          if ( lf .ne. 0 ) call psymrf ( nrclas , ipa_iv , nlmx , sym_rv 
     .    , f )
          if (allocated(rho_rv)) deallocate(rho_rv)
          if (allocated(sym_rv)) deallocate(sym_rv)
        endif
      enddo
      if (allocated(pos_rv)) deallocate(pos_rv)
      if (allocated(ipa_iv)) deallocate(ipa_iv)
      if (allocated(pos0_rv)) deallocate(pos0_rv)
      if (allocated(ipc_iv)) deallocate(ipc_iv)
      if (allocated(ips_iv)) deallocate(ips_iv)
      call tcx('symrat')
      end subroutine symrat

      subroutine psymrf(nrclas,ipa,nlmx,s,f)
C- Symmetrize forces
      implicit none
C ... Passed parameters
      integer nrclas,nlmx,ipa(nrclas)
      double precision s(nlmx,nlmx,nrclas),f(3,1)
C ... Local parameters
      integer ia,ib
      double precision x(3)
      x(1) = 0d0
      x(2) = 0d0
      x(3) = 0d0
      do  ia = 1, nrclas
        ib = ipa(ia)
        x(1)= x(1)+s(4,4,ia)*f(1,ib)+s(4,2,ia)*f(2,ib)+s(4,3,ia)*f(3,ib)
        x(2)= x(2)+s(2,4,ia)*f(1,ib)+s(2,2,ia)*f(2,ib)+s(2,3,ia)*f(3,ib)
        x(3)= x(3)+s(3,4,ia)*f(1,ib)+s(3,2,ia)*f(2,ib)+s(3,3,ia)*f(3,ib)
      enddo
      do  ia = 1, nrclas
        ib = ipa(ia)
        f(1,ib) = (s(4,4,ia)*x(1)+s(2,4,ia)*x(2)+s(3,4,ia)*x(3))*nrclas
        f(2,ib) = (s(4,2,ia)*x(1)+s(2,2,ia)*x(2)+s(3,2,ia)*x(3))*nrclas
        f(3,ib) = (s(4,3,ia)*x(1)+s(2,3,ia)*x(2)+s(3,3,ia)*x(3))*nrclas
      enddo
      end subroutine psymrf

      subroutine psymrq(nrclas,nsp,ipa,lmxa,qbyl,hbyl)
C- Symmetrize l-decomposed site charges and eval sums
      implicit none
C ... Passed parameters
      integer nrclas,nsp,lmxa,ipa(nrclas),n0
      parameter (n0=10)
      double precision qbyl(n0,nsp,1),hbyl(n0,nsp,1)
C ... Local parameters
      integer stdo,ia,ib,iprint,l,lgunit,isp,icopy_size
      double precision qsum(n0,2),hsum(n0,2),fac
      stdo = lgunit(1)
      call dpzero(qsum,2*n0)
      call dpzero(hsum,2*n0)
      fac = 1d0/nrclas
      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
          do  l = 0, lmxa
            qsum(l+1,isp) = qsum(l+1,isp) + qbyl(l+1,isp,ib)*fac
            hsum(l+1,isp) = hsum(l+1,isp) + hbyl(l+1,isp,ib)*fac
          enddo
        enddo
      enddo
      if (iprint() .ge. 35) then
        write(stdo,"(a)")' psymrq: symmetrized l-charge qbyl, hbyl and ebar     '
        write(stdo,"(a)")'   ebar=hbyl/qbyl= center of gravity of occpied states'
        write(stdo,770) (ipa(ia),ia = 1,nrclas)
c  770   format(' symmetrized qbyl,hbyl for class containing atom site ib=',20i3)
  770   format(' atom site ib=',255i3)
        if (nsp .eq. 1) write(stdo,780)
     .  (l,qsum(l+1,1),hsum(l+1,1), hsum(l+1,1)/qsum(l+1,1), l=0,lmxa)
        if (nsp .eq. 2) write(stdo,781)
     .  (l,(qsum(l+1,isp),hsum(l+1,isp),hsum(l+1,isp)/qsum(l+1,isp),isp=1,nsp), l=0,lmxa)
  780   format(' l=',i5,3f12.6)
  781   format(' l=',i5,3f12.6,'   spin 2',3f12.6)
      endif
      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
          do  l = 0, lmxa
            qbyl(l+1,isp,ib) = qsum(l+1,isp)
            hbyl(l+1,isp,ib) = hsum(l+1,isp)
          enddo
        enddo
      enddo
      end subroutine psymrq

      subroutine psymr1 ( nrclas , ipa , nr , nlml , nsp , nlmx , sym 
     . , rho , sv_p_orhoat , icmp )
      use m_struc_def, only: s_rv1
C- Symmetrize density for one class of atoms
      implicit none
C ... Passed parameters
      integer nrclas,nsp
      integer:: ipa(nrclas) , nlmx , nr , nlml , icmp
      type(s_rv1) :: sv_p_orhoat(3,1)

      double precision sym(nlmx,nlmx,nrclas),rho(nr,nlml,nsp)
C ... Local parameters
      integer stdo,ia,ib,iprint,lgunit,nn
      double precision wgt
C ... Accumulate symmetrized true density on first site
      stdo = lgunit(1)
      call dpzero(rho, nr*nlml*nsp)
      do  ia = 1, nrclas
        ib = ipa(ia)
        call pxsmr1 ( 1d0 , nr , nlml , nsp , sym ( 1 , 1 , ia ) , sv_p_orhoat( icmp , ib )%v 
     .   , rho , nn )
      enddo
C ... Copy to all sites in class
      wgt = nrclas
      do  ia = 1, nrclas
        ib = ipa(ia)
        call dpzero ( sv_p_orhoat( icmp , ib )%v , nr * nlml * nsp )
        call pysmr1 ( wgt , nr , nlml , nsp , sym ( 1 , 1 , ia ) , rho 
     .   , sv_p_orhoat( icmp , ib )%v , nn )
      enddo
      if (iprint() .ge. 40) write(stdo,100) nn,nlml*nlml
  100 format(' psymr: did',i5,'  of',i5)
      end subroutine psymr1


      subroutine symsmr(nsp,k1,k2,k3,smrho) !slat,
      use m_work,only:  iv_a_oips0, iv_a_okv,rv_a_ogv,zv_a_obgv
      use m_struc_def           !Cgetarg
      use m_lmfinit,only: lat_nabc,lat_ng,lat_nsgrp
C- Symmetrize the smooth charge density
      implicit none
C ... Passed parameters
      integer nsp,k1,k2,k3
c      type(s_lat)::slat
      double complex smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer:: n1 , n2 , n3 , ng , ngrp , ngabc(3) , isp
      complex(8) ,allocatable :: csym_zv(:)
      complex(8) ,allocatable :: cv_zv(:)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      call tcn('symsmr')
c      i_copy_size=size(lat_nabc)
c      call icopy(i_copy_size,lat_nabc,1,ngabc,1)
      ngabc=lat_nabc
      ng=lat_ng
      ngrp=lat_nsgrp
      if (ngrp .gt. 1) then
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3)
        allocate(cv_zv(ng))
        allocate(csym_zv(ng))
        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
C       call zprm3('smrho before poke',smrho,k1,k2,k3)
ccccccccccccccccccccc
      write(6,*)" aaa222xxx",sum(abs(smrho(:,:,:,1)))
ccccccccccccccccccccc
        do  isp = 1, nsp
          call gvgetf ( ng , 1 , iv_a_okv , k1 , k2 , k3 , smrho ( 1 , 
     .     1 , 1 , isp ) , cv_zv )
ccccccccccccccccccccc
      write(6,*)" aaa333xxx",sum(abs(smrho(:,:,:,1)))
ccccccccccccccccccccc
          call gvsym ( ng , rv_a_ogv , iv_a_oips0 , zv_a_obgv , cv_zv , 
     .       csym_zv )
c          call dpadd ( csym_zv , cv_zv , 1 , ng * 2 , - 1d0 )
ccccccccccccccccccccc
      write(6,*)" aaa444xxx",ng,k1,k2,k3,sum(abs(cv_zv)),sum(abs(csym_zv))
ccccccccccccccccccccc
           csym_zv=csym_zv - cv_zv
          call gvaddf ( ng , iv_a_okv , k1 , k2 , k3 , csym_zv , smrho 
     .     ( 1 , 1 , 1 , isp ) )
        enddo
ccccccccccccccccccccc
      write(6,*)" aaa555xxx",sum(abs(smrho(:,:,:,1)))
ccccccccccccccccccccc
        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
        if (allocated(csym_zv)) deallocate(csym_zv)
        if (allocated(cv_zv)) deallocate(cv_zv)
C ... Force density to be real and positive
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3)
C        do  10  i23 = 1, k2*k3
C        do  10  i1  = 1, k1
C   10   smrho(i1,i23,1) = dble(smrho(i1,i23,1))
      else
        call info(30,1,1,' Smooth density not symmetrized (ngrp=1)',0,0)
      endif
      call tcx('symsmr')
      end subroutine symsmr

      subroutine rhopos(smrho,k1,k2,k3,n1,n2,n3)
C- Make smrho real and positive
      implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3
      double complex smrho(k1,k2,k3)
C ... Local parameters
      integer stdo,lgunit,i1,i2,i3,nneg
      double precision rmin,xx
      stdo = lgunit(1)
      nneg = 0
      rmin = 999
      do  10  i3 = 1, n3
      do  10  i2 = 1, n2
      do  10  i1 = 1, n1
        xx = dble(smrho(i1,i2,i3))
        rmin = min(rmin,xx)
        if (xx .lt. 0) then
          nneg = nneg+1
          xx = 1d-8
        endif
        smrho(i1,i2,i3) = xx
   10 continue
      if (nneg .gt. 0) write(stdo,333) nneg,rmin
  333 format(' rhopos (warning): mesh density negative at',i6,
     .' points.  min=',f13.8)
      end subroutine rhopos

      subroutine rhoqm(smrho,k1,k2,k3,n1,n2,n3,nsp,vol,qsum)
C- Return charge, magnetic moment of smooth density
C ----------------------------------------------------------------------
Ci Inputs
Ci   smrho :smooth density on uniform mesh
Ci   k1..k3:
Ci   n1..n3:
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   vol   :cell volume
Co Outputs
Co   qsum  :qsum(1) = smrho(+) + smrho(-)
Co         :qsum(2) = smrho(+) - smrho(-) (for nsp=2 only)
Cl Local variables
Cl         :
Cr Remarks
Cr   Input smrho is assumed to be (rho1, rho2)
Cr   If instead smrho=(rho1+rho2,rho1-rho2) => qsum(1,2) = q+amom, q-amom
Cu Updates
Cu   13 Dec 08 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3,nsp
      double complex smrho(k1,k2,k3,nsp)
      double precision vol,qsum(2)
C ... Local parameters
      integer i,i1,i2,i3
      double precision sumi,q1,fac
      qsum(1) = 0
      qsum(2) = 0
      fac = vol/(n1*n2*n3)
      q1 = 0
      do  i = 1, nsp
        sumi = 0
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              sumi = sumi + dble(smrho(i1,i2,i3,i))
            enddo
          enddo
        enddo
        if (i .eq. 2) qsum(2) = qsum(2) + q1-sumi
        q1 = sumi
        qsum(1) = qsum(1) + sumi
      enddo
      qsum(1) = fac*qsum(1)
      qsum(2) = fac*qsum(2)
C     write(*,333) qsum
C 333 format(' rhoqm : istl charge, moment = ',2f13.7)
      end subroutine rhoqm


