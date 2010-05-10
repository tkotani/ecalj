Cdelw1 var ['oips', '0', 'oips', 'nbas', 'integer', 'iv_w_', 'symrat']
Cdelw1 var ['oipc', '0', 'oipc', 'nbas', 'integer', 'iv_w_', 'symrat']
Cdelw1 var ['opos0', '0', 'opos0', '3*nbas', 'real(8)', 'rv_w_', 'symrat']
Cdelw1 var ['oipa', '0', 'oipa', 'nbas', 'integer', 'iv_w_', 'symrat']
Cdelw1 var ['opos', '0', 'opos', '3*nbas', 'real(8)', 'rv_w_', 'symrat']
Cdelw1 var ['osym', '2', 'osym', 'nlmx*nlmx*nrclas', 'real(8)', 'rv_w_', 'symrat']
Cdelw1 var ['oistab', '-1', 'undef', '0', 'undef', 'x', 'symrat']
Cdelw1 var ['og', '-1', 'undef', '0', 'undef', 'x', 'symrat']
Cdelw1 var ['oag', '-1', 'undef', '0', 'undef', 'x', 'symrat']
Cdelw1 var ['orho', '3', 'orho', 'nr*nlml*nsp', 'real(8)', 'rv_w_', 'symrat']
Cdelw1 var ['ocv', '1', 'ocv', 'ng', 'complex(8)', 'zv_w_', 'symsmr']
Cdelw1 var ['ocsym', '1', 'ocsym', 'ng', 'complex(8)', 'zv_w_', 'symsmr']
Cdelw1 var ['okv', '-1', 'undef', '0', 'undef', 'x', 'symsmr']
Cdelw1 var ['ogv', '-1', 'undef', '0', 'undef', 'x', 'symsmr']
Cdelw1 var ['oips0', '-1', 'undef', '0', 'undef', 'x', 'symsmr']
Cdelw1 var ['obgv', '-1', 'undef', '0', 'undef', 'x', 'symsmr']
Cdelw1 undel [['undef', 'wref'], 'oistab', 'symrat']
Cdelw1 undel [['undef', 'wref'], 'og', 'symrat']
Cdelw1 undel [['undef', 'wref'], 'oag', 'symrat']
Cdelw1 undel [['undef', 'wref'], 'okv', 'symsmr']
Cdelw1 undel [['undef', 'wref'], 'ogv', 'symsmr']
Cdelw1 undel [['undef', 'wref'], 'oips0', 'symsmr']
Cdelw1 undel [['undef', 'wref'], 'obgv', 'symsmr']
Cdelw1 undel [['call'], 'orhoat', 'psymr1']
Cgetarg...info...           structure ['symrho', 'lat', 'slat']
Cgetarg...info...           structure ['symrat', 'lat', 'slat']
Cgetarg...info...           structure ['symrat', 'site', 'ssite']
Cgetarg...info...           structure ['symrat', 'spec', 'sspec']
Cgetarg...info...           structure ['symsmr', 'lat', 'slat']
Cgetarg...info...           use_to_add ['symrat', 'm_struc_def']
Cgetarg...info...           use_to_add ['symrho', 'm_struc_def']
Cgetarg...info...           use_to_add ['symsmr', 'm_struc_def']
      subroutine symrho(ssite,sspec,slat,lf,smrho,orhoat,qbyl,hbyl,f)
       
       use m_struc_def  !Cgetarg

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
C     implicit none
      integer lf,orhoat(3,*)
Cgetarg       double precision ssite(1),sspec(1),slat(1),f(*),qbyl(*),hbyl(*)
       real(8):: f(*) , qbyl(*) , hbyl(*) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex smrho(*)
      integer nsp,nbas,ngabc(3),n1,n2,n3,k1,k2,k3,nglob,iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))

      call tcn('symrho')
      call info(30,1,0,' Symmetrize density..',0,0)
      nbas = nglob('nbas')
      nsp = nglob('nsp')
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg       call upack1('lat nabc',slat,ngabc)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 

      call fftz30(n1,n2,n3,k1,k2,k3)

      call symrat(ssite,sspec,slat,nbas,nsp,lf,orhoat,qbyl,hbyl,f)
      if (iprint() .gt. 50) call prrhat(nbas,ssite,sspec,orhoat)
      call symsmr(slat,nsp,k1,k2,k3,smrho)
      call tcx('symrho')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine symrho 


      subroutine symrat(ssite,sspec,slat,nbas,nsp,lf,orhoat,qbyl,hbyl,f)
       
       use m_struc_def  !Cgetarg

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
C     implicit none
C ... Passed parameters
      integer nbas,nsp,orhoat(3,nbas),n0,lf
      parameter (n0=10)
Cgetarg       double precision ssite(1),sspec(1),slat(1),
Cgetarg      .  f(3,nbas),qbyl(n0,nsp,nbas),hbyl(n0,nsp,nbas)
       real(8):: f(3,nbas) , qbyl(n0,nsp,nbas) , hbyl(n0,nsp,nbas) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

C ... Local parameters
Cdelw1 do not change  [['undef', 'wref'], 'oag', 'symrat']
Cdelw1 do not change  [['undef', 'wref'], 'og', 'symrat']
Cdelw1 do not change  [['undef', 'wref'], 'oistab', 'symrat']
Cdelw1       integer stdo,ib0,ic,ipr,iprint,is,lgunit,lmxa,lmxl,nclass,
Cdelw1      .  ngrp,nlml,nlmx,nr,nrclas,oag,ocg,ocy,og,oidxcg,oipa,oipc,oips,
Cdelw1      .  ojcg,opos,opos0,orho,osym,igetss,mxint,oistab,ival
       integer:: stdo , ib0 , ic , ipr , iprint , is , lgunit , lmxa 
     . , lmxl , nclass , ngrp , nlml , nlmx , nr , nrclas , oag , ocg 
     . , ocy , og , oidxcg , ojcg , igetss , mxint , oistab , ival 
       integer ,allocatable :: iv_w_oipa(:)
       integer ,allocatable :: iv_w_oipc(:)
       integer ,allocatable :: iv_w_oips(:)
       real(8) ,allocatable :: rv_w_opos(:)
       real(8) ,allocatable :: rv_w_opos0(:)
       real(8) ,allocatable :: rv_w_orho(:)
       real(8) ,allocatable :: rv_w_osym(:)

      double precision plat(9),qlat(9)
C ... Heap
      integer w(1)
      common /w/ w

      call tcn('symrat')
      stdo = lgunit(1)
      ipr = iprint()

Cgetarg...info...           integer :: ocg , integer(8) :: slat%ocg 1
Cgetarg...info...           integer :: ojcg , integer(8) :: slat%ojcg 1
Cgetarg...info...           integer :: oidxcg , integer(8) :: slat%oidxcg 1
Cgetarg...info...           integer :: ocy , integer(8) :: slat%ocy 1
Cgetarg       call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
       
       ocg=slat%ocg
       ojcg=slat%ojcg
       oidxcg=slat%oidxcg
       ocy=slat%ocy

Cgetarg...info...           real(8) :: plat(9) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: qlat(9) , real(8) :: slat%qlat 9
Cgetarg...info...           integer :: oistab , integer(8) :: slat%oistab 1
Cgetarg       call upack('lat plat qlat oistab',slat,plat,qlat,oistab,0,0)
       
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       i_copy_size=size(slat%qlat) 
       call dcopy(i_copy_size,slat%qlat,1,qlat,1) 
Cdelw1 warning(1) , probably  oistab  is not defined yet at linenumber= 156
Cdelw1 do not change oistab because of ['undef', 'wref']
       oistab=slat%oistab

Cgetarg...info...           integer :: ngrp , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: og , integer(8) :: slat%osymgr 1
Cgetarg...info...           integer :: oag , integer(8) :: slat%oag 1
Cgetarg       call upack('lat nsgrp osymgr oag',slat,ngrp,og,oag,0,0)
       
       ngrp=slat%nsgrp
Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 164
Cdelw1 do not change og because of ['undef', 'wref']
       og=slat%osymgr
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 165
Cdelw1 do not change oag because of ['undef', 'wref']
       oag=slat%oag


Cdelw1       call defi (oips, nbas)
       allocate(iv_w_oips(nbas))
       if (nbas<0) iv_w_oips(:)=0

Cdelw1       call defi (oipc, nbas)
       allocate(iv_w_oipc(nbas))
       if (nbas<0) iv_w_oipc(:)=0

Cdelw1       call defrr(opos0,3*nbas)
       allocate(rv_w_opos0(3*nbas))
       if (3*nbas<0) rv_w_opos0(:)=0.0d0


C ... Separate sites into symmetry classes
C      call symcls(nbas,ssite,sspec,slat,ngrp,w(og),w(oag),
C     .  w(oips),w(opos0),nclass,w(oipc))
C      if (ipr .ge. 40) write(stdo,300) nclass
C  300 format(/' symrat: number of symmetry classes is',i3)

Cgetarg       call spackv(10,'site class',ssite,1,nbas,w(oipc))
       i_copy_size=1; 
       do i_spackv=1,nbas 
Cdelw1        call spackv_array_copy_i8_i('u',ssite(i_spackv)%class,i_copy_size,i_spackv+1-1,w(oipc))
        call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %class 
     .  , i_copy_size , i_spackv + 1 - 1 , iv_w_oipc ) 

       enddo

Cgetarg       call spackv(10,'site pos',ssite,1,nbas,w(opos0))
       i_copy_size=size(ssite(1)%pos)
       do i_spackv=1,nbas 
Cdelw1        call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos0))
        call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .  , i_copy_size , i_spackv + 1 - 1 , rv_w_opos0 ) 

       enddo

Cdelw1       nclass = mxint(nbas,w(oipc))
       nclass = mxint ( nbas , iv_w_oipc ) 


C --- Start loop over classes ---
Cdelw1       call defi  (oipa,    nbas)
       allocate(iv_w_oipa(nbas))
       if (nbas<0) iv_w_oipa(:)=0

Cdelw1       call defrr (opos,    3*nbas)
       allocate(rv_w_opos(3*nbas))
       if (3*nbas<0) rv_w_opos(:)=0.0d0


      do  ic = 1, nclass
Cdelw1         call psymr0(-2,ic,nbas,w(oipc),w(opos0),w(opos),w(oipa),nrclas)
         call psymr0 ( - 2 , ic , nbas , iv_w_oipc , rv_w_opos0 , rv_w_opos 
     .   , iv_w_oipa , nrclas ) 

        if (nrclas .gt. 0) then
Cdelw1         ib0 = ival(w(oipa),1)
         ib0 = ival ( iv_w_oipa , 1 ) 

Cgetarg         is = igetss('site spec',ib0,ssite)
         is = int(ssite(ib0)%spec) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg         call upack('spec lmxl lmxa nr',sspec,is,lmxl,lmxa,nr,0)
         
         lmxl=sspec(is)%lmxl
         lmxa=sspec(is)%lmxa
         nr=sspec(is)%nr

        nlml = (lmxl+1)**2
        if (ipr .ge. 40) write(stdo,800) ic,nrclas,nlml
  800   format(/' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)

C   ... Make the projectors; make to at least to l=1 for forces
        nlmx = max0(nlml,4)
Cdelw1         call defrr (osym,   nlmx*nlmx*nrclas)
         allocate(rv_w_osym(nlmx*nlmx*nrclas))
         if (nlmx*nlmx*nrclas<0) rv_w_osym(:)=0.0d0

Cdelw1 warning(1) , probably  oistab  is not defined yet at linenumber= 219
Cdelw1 do not change oistab because of ['undef', 'wref']
Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 219
Cdelw1 do not change og because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 219
Cdelw1 do not change oag because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oistab  is not defined yet at linenumber= 219
Cdelw1 do not change oistab because of ['undef', 'wref']
Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 219
Cdelw1 do not change og because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oag  is not defined yet at linenumber= 219
Cdelw1 do not change oag because of ['undef', 'wref']
Cdelw1         call symprj(nrclas,nlmx,ngrp,nbas,w(oistab),w(og),w(oag),plat,
Cdelw1      .    qlat,w(opos),w(osym))
         call symprj ( nrclas , nlmx , ngrp , nbas , w ( oistab ) , w 
     .   ( og ) , w ( oag ) , plat , qlat , rv_w_opos , rv_w_osym ) 


C   ... Apply the projectors to rhoat
        if (lmxl .gt. -1) then
Cdelw1         call defrr(orho,  nr*nlml*nsp)
         allocate(rv_w_orho(nr*nlml*nsp))
         if (nr*nlml*nsp<0) rv_w_orho(:)=0.0d0

Cdelw1         call psymr1(nrclas,w(oipa),nr,nlml,nsp,nlmx,w(osym),w(orho),
Cdelw1      .    orhoat,1)
         call psymr1 ( nrclas , iv_w_oipa , nr , nlml , nsp , nlmx , rv_w_osym 
     .   , rv_w_orho , orhoat , 1 ) 

Cdelw1         call psymr1(nrclas,w(oipa),nr,nlml,nsp,nlmx,w(osym),w(orho),
Cdelw1      .    orhoat,2)
         call psymr1 ( nrclas , iv_w_oipa , nr , nlml , nsp , nlmx , rv_w_osym 
     .   , rv_w_orho , orhoat , 2 ) 


C   ... Symmetrize site charges and eigval sum
Cdelw1         call psymrq(nrclas,nsp,w(oipa),lmxa,qbyl,hbyl)
         call psymrq ( nrclas , nsp , iv_w_oipa , lmxa , qbyl , hbyl ) 

        endif

C   ... Symmetrize the forces
Cdelw1         if (lf .ne. 0) call psymrf(nrclas,w(oipa),nlmx,w(osym),f)
         if ( lf .ne. 0 ) call psymrf ( nrclas , iv_w_oipa , nlmx , rv_w_osym 
     .   , f ) 


Cdelw1 rlse name= osym old_list= oips oipc opos0 oipa opos osym orho 
Cdelw1 rlse name= osym new_list= oips oipc opos0 oipa opos 
Cdelw1         call rlse (osym)
         if (allocated(rv_w_orho)) deallocate(rv_w_orho)
         if (allocated(rv_w_osym)) deallocate(rv_w_osym)

        endif
      enddo

Cdelw1 rlse name= oips old_list= oips oipc opos0 oipa opos 
Cdelw1 rlse name= oips new_list= (None)
Cdelw1       call rlse (oips)
       if (allocated(rv_w_opos)) deallocate(rv_w_opos)
       if (allocated(iv_w_oipa)) deallocate(iv_w_oipa)
       if (allocated(rv_w_opos0)) deallocate(rv_w_opos0)
       if (allocated(iv_w_oipc)) deallocate(iv_w_oipc)
       if (allocated(iv_w_oips)) deallocate(iv_w_oips)


      call tcx('symrat')

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oistab og oag]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine symrat 


      subroutine psymrf(nrclas,ipa,nlmx,s,f)

C- Symmetrize forces
C     implicit none
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
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine psymrf 


      subroutine psymrq(nrclas,nsp,ipa,lmxa,qbyl,hbyl)

C- Symmetrize l-decomposed site charges and eval sums
C     implicit none
C ... Passed parameters
      integer nrclas,nsp,lmxa,ipa(nrclas),n0
      parameter (n0=10)
      double precision qbyl(n0,nsp,1),hbyl(n0,nsp,1)
C ... Local parameters
      integer stdo,ia,ib,iprint,l,lgunit,isp
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

      if (iprint() .ge. 40) then
        write(stdo,770) (ipa(ia),ia = 1,nrclas)
  770   format(' symmetrized qbyl,hbyl for class containing ib=',20i3)
        if (nsp .eq. 1) write(stdo,780)
     .    (l,qsum(l+1,1),hsum(l+1,1), l=0,lmxa)
        if (nsp .eq. 2) write(stdo,781)
     .    (l,(qsum(l+1,isp),hsum(l+1,isp),isp=1,nsp), l=0,lmxa)
  780   format(i5,2f12.6)
  781   format(i5,2f12.6,'   spin 2',2f12.6)
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

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine psymrq 


      subroutine psymr1(nrclas,ipa,nr,nlml,nsp,nlmx,sym,rho,orhoat,icmp)

C- Symmetrize density for one class of atoms
C     implicit none
C ... Passed parameters
      integer nrclas,nsp
Cdelw1 do not change  [['call'], 'orhoat', 'psymr1']
      integer ipa(nrclas),orhoat(3,1),nlmx,nr,nlml,icmp
      double precision sym(nlmx,nlmx,nrclas),rho(nr,nlml,nsp)
C ... Local parameters
      integer stdo,ia,ib,iprint,lgunit,nn
      double precision wgt
C ... Heap
      integer w(1)
      common /w/ w

C ... Accumulate symmetrized true density on first site
      stdo = lgunit(1)
      call dpzero(rho, nr*nlml*nsp)
      do  ia = 1, nrclas
        ib = ipa(ia)
        call pxsmr1(1d0,nr,nlml,nsp,sym(1,1,ia),w(orhoat(icmp,ib)),rho,
     .    nn)
      enddo

C ... Copy to all sites in class
      wgt = nrclas
      do  ia = 1, nrclas
      ib = ipa(ia)
      call dpzero(w(orhoat(icmp,ib)), nr*nlml*nsp)
      call pysmr1(wgt,nr,nlml,nsp,sym(1,1,ia),rho,w(orhoat(icmp,ib)),nn)
      enddo

      if (iprint() .ge. 40) write(stdo,100) nn,nlml*nlml
  100 format(' psymr: did',i5,'  of',i5)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orhoat]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine psymr1 


      subroutine symsmr(slat,nsp,k1,k2,k3,smrho)
       
       use m_struc_def  !Cgetarg

C- Symmetrize the smooth charge density
C     implicit none
C ... Passed parameters
      integer nsp,k1,k2,k3
Cgetarg       double precision slat(1)
       type(s_lat)::slat

      double complex smrho(k1,k2,k3,nsp)
C ... Local parameters
Cdelw1 do not change  [['undef', 'wref'], 'obgv', 'symsmr']
Cdelw1 do not change  [['undef', 'wref'], 'ogv', 'symsmr']
Cdelw1 do not change  [['undef', 'wref'], 'oips0', 'symsmr']
Cdelw1 do not change  [['undef', 'wref'], 'okv', 'symsmr']
Cdelw1       integer n1,n2,n3,ng,ngrp,obgv,ocsym,ocv,ogv,oips0,okv,ngabc(3),
Cdelw1      .  isp
       integer:: n1 , n2 , n3 , ng , ngrp , obgv , ogv , oips0 , okv 
     . , ngabc(3) , isp 
       complex(8) ,allocatable :: zv_w_ocsym(:)
       complex(8) ,allocatable :: zv_w_ocv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C ... Heap
      integer w(1)
      common /w/ w

      call tcn('symsmr')
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           integer :: ng , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           integer :: okv , integer(8) :: slat%okv 1
Cgetarg       call upack('lat nabc ng ogv okv',slat,ngabc,ng,ogv,okv,0)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
       ng=slat%ng
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 399
Cdelw1 do not change ogv because of ['undef', 'wref']
       ogv=slat%ogv
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 400
Cdelw1 do not change okv because of ['undef', 'wref']
       okv=slat%okv

Cgetarg...info...           integer :: ngrp , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: oips0 , integer(8) :: slat%oips0 1
Cgetarg...info...           integer :: obgv , integer(8) :: slat%obgv 1
Cgetarg       call upack('lat nsgrp oips0 obgv',slat,ngrp,oips0,obgv,0,0)
       
       ngrp=slat%nsgrp
Cdelw1 warning(1) , probably  oips0  is not defined yet at linenumber= 408
Cdelw1 do not change oips0 because of ['undef', 'wref']
       oips0=slat%oips0
Cdelw1 warning(1) , probably  obgv  is not defined yet at linenumber= 409
Cdelw1 do not change obgv because of ['undef', 'wref']
       obgv=slat%obgv

      if (ngrp .gt. 1) then
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3)
Cdelw1         call defcc(ocv,     ng)
         allocate(zv_w_ocv(ng))
         if (ng<0) zv_w_ocv(:)=0.0d0

Cdelw1         call defcc(ocsym,   ng)
         allocate(zv_w_ocsym(ng))
         if (ng<0) zv_w_ocsym(:)=0.0d0

        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
C       call zprm3('smrho before poke',smrho,k1,k2,k3)
        do  isp = 1, nsp
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 418
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 418
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1           call gvgetf(ng,1,w(okv),k1,k2,k3,smrho(1,1,1,isp),w(ocv))
           call gvgetf ( ng , 1 , w ( okv ) , k1 , k2 , k3 , smrho ( 1 , 
     .     1 , 1 , isp ) , zv_w_ocv ) 

Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 419
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oips0  is not defined yet at linenumber= 419
Cdelw1 do not change oips0 because of ['undef', 'wref']
Cdelw1 warning(1) , probably  obgv  is not defined yet at linenumber= 419
Cdelw1 do not change obgv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 419
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oips0  is not defined yet at linenumber= 419
Cdelw1 do not change oips0 because of ['undef', 'wref']
Cdelw1 warning(1) , probably  obgv  is not defined yet at linenumber= 419
Cdelw1 do not change obgv because of ['undef', 'wref']
Cdelw1           call gvsym(ng,w(ogv),w(oips0),w(obgv),w(ocv),w(ocsym))
           call gvsym ( ng , w ( ogv ) , w ( oips0 ) , w ( obgv ) , zv_w_ocv 
     .     , zv_w_ocsym ) 

Cdelw1           call dpadd(w(ocsym),w(ocv),1,ng*2,-1d0)
           call dpadd ( zv_w_ocsym , zv_w_ocv , 1 , ng * 2 , - 1d0 ) 

Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 421
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 421
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1           call gvaddf(ng,w(okv),k1,k2,k3,w(ocsym),smrho(1,1,1,isp))
           call gvaddf ( ng , w ( okv ) , k1 , k2 , k3 , zv_w_ocsym , smrho 
     .     ( 1 , 1 , 1 , isp ) ) 

C         call zprm3('smrho after poke',smrho,k1,k2,k3)
        enddo
        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
Cdelw1 rlse name= ocv old_list= ocv ocsym 
Cdelw1 rlse name= ocv new_list= (None)
Cdelw1         call rlse(ocv)
         if (allocated(zv_w_ocsym)) deallocate(zv_w_ocsym)
         if (allocated(zv_w_ocv)) deallocate(zv_w_ocv)


C ... Force density to be real and positive
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3)
C        do  10  i23 = 1, k2*k3
C        do  10  i1  = 1, k1
C   10   smrho(i1,i23,1) = dble(smrho(i1,i23,1))

      else
        call info(30,1,1,' Smooth density not symmetrized (ngrp=1)',0,0)
      endif

      call tcx('symsmr')

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [okv ogv oips0 obgv]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine symsmr 


      subroutine rhopos(smrho,k1,k2,k3,n1,n2,n3)

C- Make smrho real and positive
C     implicit none
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
     .  ' points.  min=',f13.8)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
C     implicit none
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

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine rhoqm 


