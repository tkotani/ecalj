!! -----------------------------------------
subroutine  tetwt5x_dtet4(npm, ncc, &
     q, eband1,eband2, &
     qbas,ginv,efermi, &
     ntetf, nqbzw, nband,nqbz, &
     nctot,ecore,  idtetf,qbzw,ib1bz, &
     job, &
     iwgt,nbnb,  &        ! & job=0
     demin, demax,   &    ! & job=0  real(4) --->real(8) 2022june for simplicity
     frhis,nwhis, &
     nbnbx,ibjb,nhwtot, ihw,nhw,jhw,& ! & job=1
     whw,             &               ! & job=1
     iq,isp1,isp2,nqibz,eibzmode,nwgt, &! &  2016jun
     nbmx,ebmx,mtet, &
     wan)                 !
  use m_keyvalue,only: getkeyvalue
  use m_tetrakbt,only: tetrakbt_init, tetrakbt
  !> Obtain weights (imaginary part) for Dielectric function by tetrahedron method.
  !! ------------------------------------------------------------------------------
  !! \param  q       = q-vector in x(q,iw)  2*pi*q(1:3)/alat is the true q.
  !! \param  qbas    = base reciprocal lattice vectors
  !! \param  ginv    = inverse of qbas
  !! \param  efermi  = Fermi level in Rydberg
  !! \param  eband1 = eigenvalues for   k
  !! \param  eband2 = eigenvalues for q+k
  !! \param  ecore  = eigenvalues for core
  !!              ,where k runs in the Full 1st BZ, not only in IBZ.
  !! \param  idtetf,qbzw,ib1bz = outputs from TETFBZF,
  !!         which should be called before calling tetwt4.
  !! \param nctot   = the number of core
  !! \param nband   = total number of bands
  !! \param nqbz    = number of k-points in the 1st BZ
  !! \param frhis(1:nwhis+1) histogram bin.
  !!        The ihis-th bin is [frhis(ihis),frhis(ihis+1)].
  !! \param  nbnbmx
  !! \param  nhwtot  the size od whw
  !! \param  ibjb(nctot+nband,nband,nqbz);ibjb index for given ib jb k.
  !! \param  ihw(ibjb,kx):  omega index, specify the division of the histogram.
  !! \param  nhw(ibjb,kx): the number of histogram data
  !! \param  jhw(ibjb,kx): pointer to whw
  !! \remark
  !! + This version tetwt5x_dtet4 Feb2006 works for timereversal=off (npm=2 case)
  !!   - When job=1, we get
  !!   - wtthis(ihis) = \int_hislow^hisup d\omega \times
  !!   -                \int d^3k f(e(k)) (1-f(e(q+k))) \delta (omg- e(q+k) + e(k) )
  !!   - wtthis is stored into whw. See below for indexing.
  !! + job=0 is to get  wgt,nbnb, demin,demax.
  !!   - They are just for the allocation of arrays and set up required indexes.
  !!   - wgt is now used only for counting nonzero pairs of ib(occupied) jb(unoccupied)
  !!   - demin demax is maximum and minimum possible values of the excitation energies.
  !! + job=1
  !!   - whw(jhw(ibjb,kx)) \to whw(jhw(ibjb,kx)+nhw(ibjb),kx)-1 ) where ibjb(ib,jb,kx)
  !!   - histogram weights for given ib,jb,kx for histogram divisions
  !!   - from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.

  !-----------------------------------
  !r ek   = eigenvalues at k-points in the 1st BZ
  !r ekq  = eigenvalues at k+q, k   in the 1st BZ
  !r
  !r See J.Rath&A.J.Freeman PRB11(6) p.2109(1975).
  !r
  !r subroutine lindtet6 is the main part to calculate
  !r the microcell integral of wtthis(ihis) = \int_lower(ihis)^upper(ihis) d\omega \int d^3k f(e(k)) (1-f(e(q+k))) \delta (omg- e(q+k) + e(k) )
  !r where e(q+k) is unoccupied, e(k) is occupied,
  !r and f(E) denote the Fermi distribution funciton.
  !r
  !r
  ! r The numbering of band index is not unique due to degeneracy.
  ! r It affects how to choose tetrahedron.  As a result, it affect on
  ! r the integration weight, and it might break the crystal symmetry.
  ! r So I add symmetrization at the last in this routine so as to recover the symmetry.
  !r
  !-----------------------------c
  ! takao kotani Aug 2000.
  ! Sergey Faleev 2002
  ! takao kotaniMar 2002
  ! takao nov2003; mtet mode.
  !     now only for mtet=(/1,1,1/) normal mode, or
  !              for mtet=(/2,2,2/) mode (each tetrahedron is devided into 2x2x2 tetrahedrons.

  ! takao dec2003; matrix_linear() mode
  ! okumura Jan2019
  ! wan: skip ebmx cutoff for magnon calculation
  implicit none
  integer:: dummy4doxygen
  integer:: npm,jpm,ibxmx,jbxmx,jbx,nrankc1,nrankc2,nnn1,nnn2,nnni,nnnj,ncc
  !---in out -------------------------------
  integer :: nband,nqbz,nctot,ntetf,nqbzw, &
       idtetf(0:3,ntetf),ib1bz(nqbzw), nbnb(nqbz,npm)
  real(8) :: q(3), &
       eband1(nband,nqbz), eband2(nband,nqbz), &
       qbas(3,3),ginv(3,3),efermi, &
       efermia, efermib, dmua, dmub, &
       ecore(nctot), qbzw(3,nqbzw) !(n1+1)*(n2+1)*(n3+1))
  complex(8):: omg
  logical :: iwgt(nband+nctot,nband+ncc,nqbz,npm), exb
  !----------------------------------
  integer:: itet,ic, ib,jb
  real(8) :: qk(3),qkm(3),qbz(3) 
  !- For tetrahedra
  integer:: kk(0:3),kq(0:3),kr(0:3),i ,j
  real(8)   ::  det33 , kvec(3,0:3), ea(0:3), eb(0:3) ,x(0:3),am(3,3)
  complex(8):: wtt(0:3,3)
  integer:: noccx_kxx, noccx_k, noccx_kq !, noccx1
  !      real(8),parameter:: eps=0d0 !1d-12 ! cutoff check to determine cancellation.
  logical ::prt=.false.
  real(8),target :: ek_(nband+nctot,0:3), ekq_(nband+nctot,0:3)
  integer:: irnk1, nrank1, irnk2, nrank2,nibib,kx, &
       ires((nband+nctot)**2), iof1,iof2, &
       ini1(nband+nctot),ied1(nband+nctot), ixi1,ixi2,ixe1,ixe2, &
       ini2(nband+nctot),ied2(nband+nctot)
  real(8):: ekqxx(nband+nctot),summ,voltot,volt
  complex(8)::  wmean,a1,a2,wgt2
  integer:: ik,ibx,nbnc,nb1 !kqxx(nqbz),
  logical :: ipr=.false.
  real(8),parameter:: pii=3.1415926535897932d0
  real(8),pointer:: eocc(:),eunocc(:)
  integer:: idim,ivec
  real(8)::  x_(0:3),cut
  integer:: job
  real(8) :: demax(nband+nctot,nband+ncc,nqbz,npm),demax_,demaxx &
       ,demin(nband+nctot,nband+ncc,nqbz,npm),demin_,deminn
  !      real(8) :: demax(nband+nctot,nband,nqbz),demax_,demaxx
  !     &          ,demin(nband+nctot,nband,nqbz),demin_,deminn
  integer::ixx
  ! job=1
  integer:: nbnbx,nhwtot,nwhis,inihis,ihis,ikx,ibib, ikk &
       ,ioff,isum,ini,ied,jini,iini,nnn,    ntetmx,ntetmin
  integer:: ihw(nbnbx,nqbz,npm),  & ! omega pointer
       nhw(nbnbx,nqbz,npm),  &  !number of data
       jhw(nbnbx,nqbz,npm),  &  !histo-weight pointer for whw(*)
       ibjb(nctot+nband,nband+ncc,nqbz,npm)
  real(8) :: whw(nhwtot)    & ! histo-weight
       , frhis(nwhis+1), wtthis2(nwhis,0:3), wtthis(nwhis),piofvoltot
  logical ::chkwrt=.false.,wxx
  real(8),parameter:: pi=3.1415926535897932d0
  ! devided tet
  integer:: nmtet,nqbzwm,nqbzm,ntetfm
  integer,allocatable:: idtetfm(:,:,:),ib1bzm(:) !,index_qbzm(:,:,:)
  real(8),allocatable:: qbzm(:,:),qbzwm(:,:)
  real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)
  integer::kkm(0:3),kvecm(1:3, 0:3),ifeig,ifmtet,nqnumm,kqxxm,kqxxm2 &
       ,n,nsp,nband_x,kp,im !,n_index_qbzm
  real(8),allocatable:: ekzz1(:,:),ekzz2(:,:),eigtet(:,:,:),wtet(:,:,:)
  real(8)::qdummy(3)
  integer::isp1,isp2,iq,nqibz,mtet(3),  iqx,ispx, ix1,iy1,iz1,verbose,iqindx
  !      logical ::mtett !,readgwinput
  real(8):: kkv(3,0:3)
  logical :: matrix_linear
  real(8):: scissors_x0, ddw
  integer :: nbandmx_kq,nbandmx_k
  !      real(8):: fermi_dS(3,3) !not yet used...
  logical:: eibzmode, usetetrakbt
  integer:: nwgt(nqbz)
  real(8):: ebmx,  voltet
  integer:: nbmx,ifile_handle
  logical,optional:: wan
!!! tetrakbt
  real(8):: temperature     ![K], temporally
  logical:: interbandonly=.false.,intrabandonly=.false.,cmdopt0
  !---------------------------------------------------------------------
  call getkeyvalue("GWinput","tetrakbt",usetetrakbt,default=.false.)
  write(6,"(' tetwt5: job efermi usetetrakbt:=',i2,d13.6,l)") job,efermi,usetetrakbt
  if (usetetrakbt) then
     call tetrakbt_init()
  endif
  !      if(verbose()>=150) chkwrt=.true.
  voltot = abs(det33(qbas))
  piofvoltot = pi/4d0/voltot
  !      fermi_dS=0d0
  if(job==0) iwgt = .FALSE. 

  !$$$!! 2016jun
  !$$$      call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband )
  !$$$      call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
  !$$$
  !$$$!! multitet=T ==> micro tetrahedron method (divided-tetrahedron). Not used so much now...
  !$$$      mtet=(/1,1,1/)
  !$$$      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
  !      if(sum(abs(mtet))/=3) then
  !        mtett=.true.
  !        if(sum(abs(mtet))<3)
  !     &  print *, ' we use devided-tetrahedron scheme mtet=',mtet
  !      else
  !        mtett=.false.
  !      endif
  !$$$      if(mtett) then
  !$$$        ifmtet=ifile_handle()
  !$$$        open (ifmtet, file='mtet',form='unformatted')
  !$$$        read(ifmtet) nmtet,nqbzwm,nqbzm,ntetfm !,n_index_qbzm
  !$$$        allocate(
  !$$$     &       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
  !$$$        ! Index for tetrahedron;    qbzmw(idtetfm) gives extended q vector.
  !$$$     &       ib1bzm(nqbzwm), qbzm(3,nqbzm) !,index_qbzm(n_index_qbzm,n_index_qbzm,n_index_qbzm)
  !$$$        ! qbzm(1:3,ib1bz(iq)) gives q vector within the 1st bz.
  !$$$     &       , wtet(0:3,nmtet,ntetf) )
  !$$$        read(ifmtet) idtetfm,ib1bzm,qbzm,qbzwm,wtet  !,index_qbzm
  !$$$        close(ifmtet)
  !$$$        ifeig=501
  !$$$        open (ifeig,file='eigmtet',form='unformatted')
  !$$$        read(ifeig) nband_x,nqnumm,nsp
  !$$$        print *,'readin eigmtet ', nband_x,nqnumm,nsp
  !$$$        if(nband_x /=nband ) then
  !$$$          call rx( 'tetwt5: nband_x /=nband')
  !$$$        endif
  !$$$        allocate( eigtet(nband,nqnumm,nsp) )
  !$$$        do iqx=1,nqnumm
  !$$$          do ispx=1,nsp
  !$$$            read (ifeig) eigtet(1:nband,iqx,ispx)
  !$$$            if(verbose()>150) write(6,"('iq=',i3,'  eig(1:5)=',5f10.4)")iqx,eigtet(1:5,iqx,ispx)
  !$$$          enddo
  !$$$        enddo
  !$$$        close(ifeig)
  !$$$      else
  nmtet = 1
  nqbzwm= nqbzw
  nqbzm = nqbz
  allocate( idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm), &
       ib1bzm(nqbzwm),           qbzm(3,nqbzm) , wtet(0:3,nmtet,ntetf) )
  idtetfm(:,1,:)=idtetf
  qbzwm =qbzw
  ib1bzm=ib1bz
  wtet=0.25d0
  !$$$      endif

  !!
  nbnc = nband+nctot
  nb1  = nband+1
  allocate( ekxx1(nband+nctot,nqbz),ekxx2(nband+nctot,nqbz))
  do kx = 1, nqbz
     ekxx1( 1:nband, kx) = eband1(1:nband,kx)
     ekxx2( 1:nband, kx) = eband2(1:nband,kx)
     ekxx1( nband+1: nband+nctot, kx) = ecore(1:nctot)
     ekxx2( nband+1: nband+nctot, kx) = ecore(1:nctot)
  enddo
  !! Read eigenvalues at q and q+k ---------------------------------------
  !!  ekzz1 for k
  !!  ekzz2 for q+k.
  allocate( ekzz1(nband+nctot,nqbzm),ekzz2(nband+nctot,nqbzm))
  !$$$      if(mtett) then
  !$$$        print *, ' mtett mode nqbzm nqibz=',nqbzm,nqibz
  !$$$        do kx = 1, nqbzm
  !$$$          ekzz1( 1:nband, kx)              = eigtet(1:nband,kx,isp1)
  !$$$          ekzz1( nband+1: nband+nctot, kx) = ecore(1:nctot)
  !$$$          ekzz2( nband+1: nband+nctot, kx) = ecore(1:nctot)
  !$$$        enddo
  !$$$        if(iq<=nqibz) then
  !$$$          do kx = 1, nqbzm
  !$$$            kqxxm = iqindx(q(1:3)+qbzm(1:3,kx),ginv,qbzm,nqbzm)
  !$$$            ekzz2(1:nband, kx) = eigtet(1:nband,kqxxm,isp2) !ekzz1(1:nband, kqxxm)
  !$$$          enddo
  !$$$        else
  !$$$          do kx = 1, nqbzm
  !$$$            kp  = nqbzm *(iq - nqibz) + kx
  !$$$            ekzz2(1:nband, kx)= eigtet(1:nband,kp,isp2)
  !$$$          enddo
  !$$$        endif
  !$$$        deallocate( eigtet)!,index_qbzm)
  !$$$      else
  ekzz1=ekxx1
  ekzz2=ekxx2
  !$$$      endif

  !! Add eigenvalue shift. scissors_x0() is defined in switch.F
  if(scissors_x0()/=0d0) then
     call addsciss(scissors_x0(),efermi,(nband+nctot)*nqbzm,    ekzz1)
     call addsciss(scissors_x0(),efermi,(nband+nctot)*nqbzm,    ekzz2)
  endif

  !! Check
  volt = 0d0
  do itet = 1, ntetf
     do im   = 1, nmtet
        kvec(1:3,0:3) = qbzwm (1:3, idtetfm(0:3,im,itet) )
        do i = 1,3
           kvec(1:3,i) = kvec(1:3,i) - kvec(1:3,0)
        enddo
        volt = volt + abs(det33(kvec(1:3,1:3))/6d0)
        !        write(6,"('itet im vol=',2i5,d13.5)") itet,im,abs(det33(kvec(1:3,1:3))/6d0)
     enddo
  enddo
  !      if(abs(volt-voltot)>1d-10) call rx( ' tetwt: abs(volt-voltot)>1d-10')

  if(job==0) then
     demin=  1d10
     demax= - 1d10
  endif
  !     call getkeyvalue("GWinput","defermia",dmua,default=0d0)
  !     call getkeyvalue("GWinput","defermib",dmub,default=0d0)
  !     call getkeyvalue("GWinput","exband",exb,default=.false.)
  !     efermia = efermi + dmua
  !     efermib = efermi + dmub
  !      write(6,*) "defermia, defermib:", dmua,dmub
  !      write(6,*) " efermia,  efermib:", efermia,efermib
  efermia = efermi
  efermib = efermi
  interbandonly=cmdopt0('--interbandonly')
  intrabandonly=cmdopt0('--intrabandonly')
  !! === Loop over tetrahedron ===
  do 1000 itet = 1, ntetf !;
     kk (0:3) = ib1bz( idtetf(0:3,itet) )     !  k
     if(eibzmode) then
        if(sum(nwgt(kk(0:3)))==0) cycle
     endif
     kkv(1:3, 0:3) = qbzw (1:3, idtetf(0:3,itet) )
     do 1100 im = 1,nmtet
        kkm (0:3)       = ib1bzm( idtetfm(0:3,im,itet) ) !  k   in micro-tet
        kvec(1:3, 0:3) = qbzwm ( 1:3, idtetfm(0:3,im,itet) )
        do i = 1,3
           am(1:3,i) = kvec(1:3,i-1) - kvec(1:3,3)
        enddo
        voltet = abs(det33(am)/6d0)
        ek_ ( 1:nband+nctot, 0:3) = ekzz1( 1:nband+nctot, kkm(0:3)) ! k
        ekq_( 1:nband+nctot, 0:3) = ekzz2( 1:nband+nctot, kkm(0:3)) ! k+q
        !          noccx_k = noccx1 ( ek_(1:nband, 0:3),4,nband, efermi)!the highest number of occupied states
        !          noccx_kq= noccx1 (ekq_(1:nband, 0:3),4,nband, efermi)
        !     the highest number of occupied states
        noccx_k = maxval(count( ek_(1:nband, 0:3)<efermi,dim=1) )
        noccx_kq= maxval(count( ekq_(1:nband,0:3)<efermi,dim=1) )
        nbandmx_k   = nband
        nbandmx_kq  = nband
        !! exclude wannier model: (okumura, 2017/06/13)
        !! ebmx band cutoff 2016Jun
        if (present(wan) .AND. wan) then
           continue
        else
           do i=1,nbmx
              if( maxval(ek_(i, 0:3)) >ebmx) then
                 nbandmx_k = i-1
                 exit
              endif
           enddo
           do i=1,nbmx
              if( maxval(ekq_(i, 0:3)) >ebmx) then
                 nbandmx_kq = i-1
                 exit
              endif
           enddo
        endif
        !!
        do jpm = 1,npm
           if(jpm==1) then
              ibxmx = noccx_k + nctot
              jbxmx = nbandmx_kq
           else
              ibxmx = nbandmx_k
              jbxmx = noccx_kq + nctot
           endif
           do ibx  = 1, ibxmx !noccx_k + nctot  !   occupied
              do jbx  = 1, jbxmx !nband             ! unoccupied
                 if(ibx<=noccx_k .OR. jpm==2  ) then
                    ib = ibx
                 else
                    ib = ibx - noccx_k + nband
                 endif
                 if(jbx<=noccx_kq .OR. jpm==1  ) then
                    jb = jbx
                 else
                    jb = jbx - noccx_kq + nband
                 endif
                 !! --intraband mode  takao@nov2021
                 if(intrabandonly .AND. (ib/=jb) .AND. job==1) cycle
                 if(interbandonly .AND. (ib==jb) .AND. job==1) cycle
                 !! This mechanism treat ek_ and ekq_ as occpied or unoccupied.
                 if(jpm==1) then
                    eocc   => ek_ (ib,0:3)
                    eunocc => ekq_(jb,0:3)
                 else
                    eunocc => ek_ (ib,0:3)
                    eocc   => ekq_(jb,0:3)
                 endif
                 if( minval(eocc) <= efermia .AND.  maxval(eunocc) >= efermib ) then
                    continue
                 else
                    cycle
                 endif
                 ! f( maxval(eunocc(:)-eocc(:)) <dmub-dmua ) cycle ! this makes a bit effective.
                 if( maxval(eunocc(:)-eocc(:)) <0d0 ) cycle ! this makes a bit effective.
                 if(job==0 ) then  !takao
                    iwgt(ib,jb,kk(0:3),jpm)= .true.
                    x(0:3) = .5d0*(eocc-eunocc) ! + omg !Denominator. unit in Hartree.
                    demax_ =  maxval(-x(0:3)) ! in Hartree
                    demin_ =  minval(-x(0:3)) ! in Hartree
                    do ixx=0,3
                       demax(ib,jb,kk(ixx),jpm) = max(demax_, demax(ib,jb,kk(ixx),jpm))
                       demin(ib,jb,kk(ixx),jpm) = min(demin_, demin(ib,jb,kk(ixx),jpm))
                    enddo
                    cycle
                 endif
                 x(0:3) = .5d0*(eocc-eunocc) ! + omg !Denominator. unit in Hartree.
                 wtthis2 = 0d0
                 if(chkwrt) then
                    write(6,"('### Goto lindtet6: itet ib jb Ef=',i8,2i5,d12.4,' ###')" ) itet,ib,jb,efermi
                    write(6,"('  eocc  - Ef= ',4f10.3)") eocc  -efermi
                    write(6,"('  eunocc- Ef= ',4f10.3)") eunocc-efermi
                    write(6,"('  -x(a.u.)= ',4f10.3)") -x(0:3)
                 endif
                 if(usetetrakbt) then
                    call tetrakbt(eunocc-efermib, eocc-efermia, x, voltet, frhis,nwhis,   wtthis2)
                 else
                    call lindtet6( kkv,kvec,eocc,eunocc, x, efermia,efermib,frhis,nwhis,  wtthis2)
                 endif
                 if(chkwrt) then
                    write(6,"('  sum nwhis wtthis2=',i5,d13.5)") nwhis,sum(wtthis2)
                    write(6,"('  === ihis [range(a.u.)] wtthis2(0:3) ===')")
                    do ihis = 1,nwhis
                       ddw = frhis(ihis+1)-frhis(ihis)
                       if(chkwrt .AND. sum(abs(wtthis2(ihis,:)))/=0d0) then
                          write(6,"('ttt',3i4,f10.5,4d12.4)") &
                               itet,ib,jb,(frhis(ihis)+frhis(ihis+1))/2d0,wtthis2(ihis,:)/ddw
                       endif
                    enddo
                 endif
                 do ikx = 0,3
                    ibib = ibjb(ib,jb,kk(ikx),jpm)
                    jini = jhw(ibib,kk(ikx),jpm)
                    iini = ihw(ibib,kk(ikx),jpm)
                    nnn =  nhw(ibib,kk(ikx),jpm)
                    if(matrix_linear()) then
                       whw(jini:jini+nnn-1)=whw(jini:jini+nnn-1) + wtthis2(iini:iini+nnn-1,ikx) *piofvoltot
                    else
                       whw(jini:jini+nnn-1)=whw(jini:jini+nnn-1) + wtthis2(iini:iini+nnn-1,0) *piofvoltot &
                            * 4*wtet(ikx,im,itet) ! piofvoltot= pi/voltot/4
                    endif
                 enddo
              enddo
           enddo
        enddo
1100 enddo
1000 enddo
  deallocate(idtetfm, qbzwm,ib1bzm, qbzm)

  !! === Symmetrization of wgt and whw   ===
  !! NOTE: We just enforce the same weight for degenerated bands.
  !!       In the previous do loop 1000, we judged band connectivity just by the band index.
  !!       If they are degenerated, we have some unbiguity.
  !!       (which band connedt to which in the microtetrahedron).
  do kx = 1, nqbz     ! ipr = .false.; if(ipr)  print *,' kx =',kx
     if(eibzmode) then
        if(nwgt(kx)==0) cycle
     endif
     call chkdgn( ekxx1(:,kx), nband, nrank1, ini1,ied1,0 ,ipr)
     call chkdgn( ekxx2(:,kx), nband, nrank2, ini2,ied2,0 ,ipr)
     nrankc1 = 0
     if(nctot/=0) then
        call chkdgn(ecore, nctot, &
             nrankc1, ini1(nrank1+1), ied1(nrank1+1), nband,ipr)
     endif
     nrankc2 = 0
     if(nctot/=0) then
        call chkdgn(ecore, nctot, &
             nrankc2, ini2(nrank2+1), ied2(nrank2+1), nband,ipr)
     endif
     !       if(ipr) print *,' kx nrank =',kx,nrank1,nrank2,nrankc
     do jpm=1,npm
        if(jpm==1) then
           nnn1=nrankc1
           nnn2=0
        else
           nnn1=0
           nnn2=nrankc2
        endif
        do irnk1 = 1, nrank1 + nnn1
           do irnk2 = 1, nrank2 + nnn2
              ixi1  = ini1(irnk1);  ixe1 = ied1(irnk1)
              ixi2  = ini2(irnk2);  ixe2 = ied2(irnk2)
              if(job==0) then
                 wxx=.false.
                 if( count(iwgt(ixi1:ixe1,ixi2:ixe2,kx,jpm))>0 ) wxx= .TRUE. 
                 iwgt(ixi1:ixe1,ixi2:ixe2, kx,jpm ) =  wxx
                 demaxx = maxval(demax(ixi1:ixe1,ixi2:ixe2, kx,jpm ) )
                 demax(ixi1:ixe1,ixi2:ixe2, kx,jpm ) = demaxx
                 deminn = minval(demin(ixi1:ixe1,ixi2:ixe2, kx,jpm ) )
                 demin(ixi1:ixe1,ixi2:ixe2, kx,jpm ) = deminn
              else
                 !! --- skip for no data section
                 isum=0
                 do ib = ixi1,ixe1
                    do jb = ixi2,ixe2
                       if(ibjb(ib,jb,kx,jpm)/=0) isum = isum + 1
                    enddo
                 enddo
                 if(isum/=0 .AND. isum /= (ixe1-ixi1+1)*(ixe2-ixi2+1) ) then
                    ! top2rx 2013.08.09 kino                  stop 'tetwt5: isum srrange---1'
                    call rx( 'tetwt5: isum srrange---1')
                 endif
                 if(isum==0) cycle
                 !! --- Get symmetrized whw
                 wtthis=0d0
                 do ib = ixi1,ixe1
                    do jb = ixi2,ixe2
                       ibib = ibjb(ib,jb,kx,jpm)
                       ini = ihw(ibib,kx,jpm)
                       ied = ihw(ibib,kx,jpm)+ nhw(ibib,kx,jpm)-1
                       ioff= jhw(ibib,kx,jpm)
                       wtthis(ini:ied) = &
                            wtthis(ini:ied) + whw(ioff:ioff+nhw(ibib,kx,jpm)-1)
                    enddo
                 enddo
                 wtthis = wtthis/((ixe1-ixi1+1)*(ixe2-ixi2+1) )
                 do ib = ixi1,ixe1
                    do jb = ixi2,ixe2
                       ibib = ibjb(ib,jb,kx,jpm)
                       ini = ihw(ibib,kx,jpm)
                       ied = ihw(ibib,kx,jpm)+ nhw(ibib,kx,jpm)-1
                       ioff= jhw(ibib,kx,jpm)
                       whw(ioff:ioff+nhw(ibib,kx,jpm)-1) = wtthis(ini:ied)
                    enddo
                 enddo
              endif
           enddo
        enddo
     enddo
1120 continue
  enddo ! end of kx loop

  if(job==0) then
     do jpm =1, npm
        do ik  =1, nqbz
           nbnb(ik,jpm) = 0
           if(jpm==1) then
              nnni=nctot
              nnnj=0
           else
              nnni=0
              nnnj=nctot
           endif
           do i =1, nband+nnni
              do j =1, nband+nnnj
                 !        if(abs(wgt(i,j,ik) )>0d0 .and. j>nband) stop 'tetwt5: bug ?' ! bug checker
                 !        if(abs(wgt(i,j,ik) )>eps ) then
                 if(iwgt(i,j,ik,jpm)) then
                    nbnb(ik,jpm) = nbnb(ik,jpm)+1
                    !c          write(56,"(3i4,2d13.6)") i,j,ik, wgt(i,j,ik)
                    !          if(abs(wgt(i,j,ik) )>0.1d0 )
                    !     &      write(6,"(' k i ik=',3i4,'  wgt= ',2d13.6)")
                    !     &      i,j,ik, wgt(i,j,ik)
                 endif
              enddo
           enddo !;  write(6,*)' ik=',ik,' num of nonzero wgt=',nbnb(ik)
        enddo
     enddo
     !       write(6,*)' max num of nonzero wgt(k)=',maxval(nbnb)
     !       write(6,*)' tot num of nonzero wgt   =',sum(nbnb)
     !       write(6,*)' sum of wgt   =',sum(wgt)
     !       write(6,*)' tetwt4: end '; call cputid  (0)
     write(6,"('  tetwt5_dtet3: maxval(nbnb(1:nqbz,1:npm)) sum(nbnb)=',2i8 &
          ,' count(iwgt) =',i10)") maxval(nbnb),sum(nbnb),count(iwgt) !,sum(wgt)
  endif
  if (allocated(idtetfm)) deallocate(idtetfm)
  if (allocated(ib1bzm)) deallocate(ib1bzm)
  if (allocated(qbzm)) deallocate(qbzm)
  if (allocated(qbzwm)) deallocate(qbzwm)
  if (allocated(ekxx1)) deallocate(ekxx1)
  if (allocated(ekxx2)) deallocate(ekxx2)
  if (allocated(ekzz1)) deallocate(ekzz1)
  if (allocated(ekzz2)) deallocate(ekzz2)
  if (allocated(eigtet)) deallocate(eigtet)
  if (allocated(wtet)) deallocate(wtet)
  if (allocated(idtetfm)) deallocate(idtetfm)
  if (allocated(ib1bzm)) deallocate(ib1bzm)
  if (allocated(qbzm)) deallocate(qbzm)
  if (allocated(qbzwm)) deallocate(qbzwm)
  if (allocated(ekxx1)) deallocate(ekxx1)
  if (allocated(ekxx2)) deallocate(ekxx2)
  if (allocated(ekzz1)) deallocate(ekzz1)
  if (allocated(ekzz2)) deallocate(ekzz2)
  if (allocated(eigtet)) deallocate(eigtet)
  if (allocated(wtet)) deallocate(wtet)
end subroutine tetwt5x_dtet4
!------------------------------------------------------------------------



!----------------------------------------------------------------
subroutine hisrange(frhis,nwhis, demin,demax, &
     ihw,nhw)
  !- determine the pointer ihw and number of pointer.
  implicit none
  integer:: nwhis,ihw,nhw,ihis
  !      real(8):: frhis(nwhis+1), demin,demax
  real(8):: frhis(nwhis+1)
  real(8):: demin,demax
  !  Range for each division is from frhis(ihis-1) to frhis(ihis).
  !      print *,' hisrange nwhis=',nwhis,demin,demax
  ihw=-9999;nhw=-9999
  do ihis = 1,nwhis
     if (demin < frhis(ihis+1) ) then
        ihw = ihis
        exit
     endif
  enddo
  do ihis = ihw, nwhis
     if (demax < frhis(ihis+1) ) then
        nhw = ihis - ihw +1
        exit
     endif
  enddo
  if(ihw==-9999 .OR. nhw==-9999) then
     print *,' error info=',nwhis, demin,demax,ihw,nhw
     do ihis = 1, nwhis
        write(6,*) ihis, frhis(ihis+1)
     enddo
     ! top2rx 2013.08.09 kino        stop ' hisrange: wrong'
     call rx( ' hisrange: wrong')
  endif
end subroutine hisrange


!------------------------------------------------------------------------
subroutine lindtet6(kkv,kvec, ea, eb, x, efermia, efermib, frhis, nwhis,  wtthis )
  !     m  fermi_ds)
  !- Calculate the imaginary part of \int dk1 dk2 dk3 f(ea)(1-f(eb))/(\omega+x ), that is,
  !  the microcell integral of wtthis(ihis) = \int_lower(ihis)^upper(ihis) d\omega \int d^3k f(ea(k)) (1-f(eb(q+k))) \delta (omg + x(k) )
  ! f(E) denote the Fermi distribution funciton. Only for T=0.
  ! Tetrahedon is specified by values at 4 corners; kvec(1:3, 1:4), ea(1:4), eb(1:4), x(1:4)
  !r This code is based on J.Rath&A.J.Freeman PRB11(6) p.2109(1975).
  !r The \sum_ihis wthis(ihis) = the total volume of the microcell = \int d^3 k
  ! akao/Feb/2002 -----------------------------------------------------------------
  !! ds is the accumulation variable for the fermi surface.
  implicit none
  integer :: ieaord(1:4),i,isig,n,itmp,ix
  real(8)    ::  kvec(1:3, 1:4), x(1:4), ea(1:4), &
       kk(3,1:8),xx(1:8),ee(1:4), am(3,3), eb(1:4),ebf(1:8), vcell,efermi,etest ,  kkv(3,4), efermia, efermib
  integer:: nwhis
  real(8):: frhis(nwhis+1),wtthis(nwhis,4),det33
  logical:: chkwrt=.false.
  !      real(8):: fermi_ds(3,3)

  !! this block is equivalent to sortea(ea,ieaord,n,isig) which gives ieaord (ordering index of ea)
  ! ifdef EXPAND_SORTEA
  n=4
  !      isig = 1
  ! option noparallel
  do i = 1,n
     ieaord(i) = i
  enddo
  ! option noparallel
  do ix= 2,n
     ! option noparallel
     do i=ix,2,-1
        if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
           itmp = ieaord(i-1)
           ieaord(i-1) = ieaord(i)
           ieaord(i) = itmp
           !            isig= -isig
           cycle
        endif
        exit
     enddo
  enddo
  ieaord(1:4) = ieaord(4:1:-1)
  ! else
  !      call sortea( ea,ieaord, 4 ,isig); ieaord(1:4) = ieaord(4:1:-1)
  ! endif
  ! the order E_4,E_3,E_2,E_1 This suits for Rath&Freeman.
  kk(1:3,1:4) = kvec(1:3,ieaord(1:4))
  ! 4 corners  denoted by kvec(:,1:4), ee(1:4), and xx(1:4)
  ee (1:4)    = ea (ieaord(1:4)) - efermia
  xx (1:4)    = x  (ieaord(1:4))
  ebf(1:4)    = eb (ieaord(1:4)) - efermib

  ! ccccccccccc
  !      write(6,"(' i iea=',2i4)") (i,ieaord(i),i=1,4)
  !      write(6,"(' i=',i3,' xx ee ebf =',3f10.5,'  kk=',3f10.5)")
  !     &     (i,xx(i),ee(i),ebf(i),kk(1:3,i),i=4,1,-1)
  ! ccccccccccc


  if( 0d0<=ee(4) ) then
     if(chkwrt) write(6,*) 'lindtet6: fig 000'
     !        wttx = (0d0,0d0)
  elseif( ee(4) < 0d0 .AND. 0d0<= ee(3) ) then   !!! Fig 1.
     if(chkwrt) write(6,*) 'lindtet6: fig 1 xxx'
     call midk3(kk,ee,xx,ebf, 4,2,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1 !K1 is on the like k4---k2.
     call midk3(kk,ee,xx,ebf, 4,1,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
     call midk3(kk,ee,xx,ebf, 4,3,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3

     !! Corners of occupied states at the Fermi energy are kk(1:3,i+4) i=1,2,3.
     !! This makes a Fermi surface
     !        if(sum(abs(ea-eb))<1d-6) then
     !          call fermids(kk(:,2+4),kk(:,1+4),kk(:,3+4),fermi_dS)
     !        endif

     !! still developing
     call  inttetra6(kkv,kk,xx,ebf,(/4,1+4,2+4,3+4/), & ! &  k4,K1,K2,K3
     frhis,  nwhis, &
          wtthis )
  elseif( ee(3) < 0d0 .AND. 0d0<= ee(2) ) then   !!! Fig 2.
     if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx'
     call midk3(kk,ee,xx,ebf, 4,2,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1
     call midk3(kk,ee,xx,ebf, 4,1,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
     call midk3(kk,ee,xx,ebf, 3,1,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3
     call midk3(kk,ee,xx,ebf, 3,2,  kk(1,4+4),xx(4+4),ebf(4+4)) !K4

     !        if(sum(abs(ea-eb))<1d-6) then
     !          call fermids(kk(:,1+4),kk(:,2+4), kk(:,3+4),fermi_dS) !K1-K2-K3
     !          call fermids(kk(:,1+4),kk(:,2+4), kk(:,4+4),fermi_dS) !K1-K2-K4
     !        endif

     if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx2'
     call  inttetra6(kkv,kk,xx,ebf,(/4, 3,1+4,2+4/), & ! &  k4,k3,K1,K2
          frhis,  nwhis, &
          wtthis )
     if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx3'
     call  inttetra6(kkv,kk,xx,ebf,(/3,2+4,3+4,1+4/), & ! &  k3,K2,K3,K1
          frhis,  nwhis, &
          wtthis )
     if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx4'
     call  inttetra6(kkv,kk,xx,ebf,(/3,1+4,3+4,4+4/), & ! &  k3,K1,K3,K4
          frhis,  nwhis, &
          wtthis )
     if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx5'
     !-----------
  elseif( ee(2) < 0d0 .AND. 0d0<= ee(1) ) then   !!! Fig 3.
     if(chkwrt) write(6,*) 'lindtet6: fig 3 xxx'
     call midk3(kk,ee,xx,ebf, 1,4,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1
     call midk3(kk,ee,xx,ebf, 1,2,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
     call midk3(kk,ee,xx,ebf, 1,3,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3

     !        if(sum(abs(ea-eb))<1d-6) then
     !          call fermids(kk(:,1+4),kk(:,2+4), kk(:,3+4),fermi_dS) !K1-K2-K3
     !        endif

     call  inttetra6(kkv,kk,xx,ebf,(/3, 4,3+4, 2/), & ! &  k3,k4,K3,k2
          frhis,  nwhis, &
          wtthis )
     call  inttetra6(kkv,kk,xx,ebf,(/4,1+4,2+4,3+4/),&  ! &  k4,K1,K2,K3
          frhis,  nwhis, &
          wtthis )
     call  inttetra6(kkv,kk,xx,ebf,(/4, 2,2+4,3+4/), & ! &  k4,k2,K2,K3
          frhis,  nwhis, &
          wtthis )
  else
     if(chkwrt) write(6,*) 'lindtet6: fig 4 xxx'
     call  inttetra6(kkv,kk,xx,ebf, (/1,2,3,4/), & ! &  k1,k2,k3,k4
          frhis, nwhis, &
          wtthis )
  endif
  if(chkwrt) write(6,*) 'end of lindtet6',wtthis
  !      stop 'yyyyyyyyyyyyyyyyyyyyyy'
end subroutine lindtet6


!-----------------------------------------------------
subroutine inttetra6(kkv,kk_,xx_,ebf,itetx, &
     frhis, nwhis, &
     wtthis )
  ! calculate tetrahedron integral Eq.(16).
  ! the four corners and denoted by itetx.
  !i kk (k1-k4) and xx (value of denominator)
  !i Kx (K1-K4) and xkx
  !r The four corners are selected from 8 points. It is specified by itetx.
  !r wtthis is accumulated
  implicit none
  integer:: itetx(4),ix,i,ieaord(4),isig,n,itmp
  real(8):: kk_(3,1:8),xx_(1:8), kk(3,1:8),xx(1:8),ebf(1:4),ebfin(1:8) ,ee(4)  ,kin(3,4), xin(4)
  integer:: nwhis
  real(8):: frhis(nwhis+1),wtthis(nwhis,4)
  logical :: chkwrt=.false.

  real(8)   ::  kkv(3,4)
  ! ---
  if(chkwrt) then
     print *,' inttetra6: ---------------------------'
     !        print *,' ebf  =', ebf(1:4)
     !        print *,' ebfKx=', ebfKx(1:4)
  endif

  !--- kk xin ebfin ---
  !      do i = 1,4
  !        if( itetx(i) < 10 ) then   ! kk (k1-k4 in Paper)
  !          ix = itetx(i)
  !          xin(    i) = xx_(    ix)
  !          ebfin(i)   = ebf(ix)
  !          kin(1:3,i) = kk_(1:3,ix)
  !        else                       ! Kx (K1-K4 in Paper)
  !          ix = itetx(i)/10
  !          xin(    i) = xKx_(    ix)
  !          ebfin(i)   = ebfKx(ix)
  !          kin(1:3,i) = Kx_ (1:3,ix)
  !        endif
  !      enddo

  !--- kk xin ebfin ---
  xin(    1:4) = xx_(    itetx(1:4))
  kin(1:3,1:4) = kk_(1:3,itetx(1:4))
  ebfin(  1:4) = ebf(    itetx(1:4))

  if(chkwrt) then
     write(6,"(' i=',i3,' xin ebfin =',2f10.5,'  kin=',3f10.5)") &
          (i,xin(i),ebfin(i),kin(1:3,i),i=4,1,-1)
  endif

  !... ee decendent order ee(4)>ee(3)>ee(2)>ee(1)
! #ifdef EXPAND_SORTEA
!   n = 4
!   isig = 1
!   ! option noparallel
!   do i = 1,n
!      ieaord(i) = i
!   enddo
!   ! option noparallel
!   do ix= 2,n
!      ! option noparallel
!      do i=ix,2,-1
!         if( ebfin(ieaord(i-1)) >ebfin(ieaord(i) ) ) then
!            itmp = ieaord(i-1)
!            ieaord(i-1) = ieaord(i)
!            ieaord(i) = itmp
!            isig= -isig
!            cycle
!         endif
!         exit
!      enddo
!   enddo

! #else
  call sortea( ebfin,ieaord, 4 ,isig)  ! the order E_4,E_3,E_2,E_1 This suits for Rath&Freeman.
  !#endif
  kk(1:3,1:4) = kin  (1:3,ieaord(1:4))   ! 4 corners  denoted by kkv(:,1:4), ee(1:4), and xx(1:4)
  ee(    1:4) = ebfin(    ieaord(1:4))
  xx(    1:4) = xin  (    ieaord(1:4))

  if(chkwrt)write(6,"(' i=',i3,' xx ee =',2f10.5,'  kk=',3f10.5)") &
       (i,xx(i),ee(i),kk(1:3,i),i=4,1,-1)

  if( 0d0>=ee(4) ) then
     !        wttx = (0d0,0d0)
  elseif( ee(4) > 0d0 .AND. 0d0>= ee(3) ) then   !!! Fig 1.
     if(chkwrt) write(6,*) 'inttetra5: fig 1'
     call midk(kk,ee,xx, 4,2,  kk(1,1+4),xx(1+4)) !K1 -> Kx(:,1), x(K1) -> xkx(1). K1 is on the like k4---k2.
     call midk(kk,ee,xx, 4,1,  kk(1,2+4),xx(2+4)) !K2
     call midk(kk,ee,xx, 4,3,  kk(1,3+4),xx(3+4)) !K3
     call inttetrac6(kkv,kk,xx, (/4,1+4,2+4,3+4/), & ! &  k4,K1,K2,K3
          frhis, nwhis, &
          wtthis )
  elseif( ee(3) > 0d0 .AND. 0d0>= ee(2) ) then   !!! Fig 2.
     if(chkwrt) write(6,*) 'inttetra5: fig 2'
     call midk(kk,ee,xx, 4,2,  kk(1,1+4),xx(1+4)) !K1
     call midk(kk,ee,xx, 4,1,  kk(1,2+4),xx(2+4)) !K2
     call midk(kk,ee,xx, 3,1,  kk(1,3+4),xx(3+4)) !K3
     call midk(kk,ee,xx, 3,2,  kk(1,4+4),xx(4+4)) !K4
     call inttetrac6(kkv,kk,xx, (/4, 3,1+4,2+4/),  &! &  k4,k3,K1,K2
          frhis, nwhis, &
          wtthis )
     call inttetrac6(kkv,kk,xx, (/3,2+4,3+4,1+4/), & ! &  k3,K2,K3,K1
          frhis, nwhis, &
          wtthis )
     call inttetrac6(kkv,kk,xx, (/3,1+4,3+4,4+4/), & ! &  k3,K1,K3,K4
          frhis, nwhis, &
          wtthis )
  elseif( ee(2) > 0d0 .AND. 0d0>= ee(1) ) then   !!! Fig 3.
     if(chkwrt) write(6,*) 'inttetra5: fig 3'
     call midk(kk,ee,xx, 1,4,  kk(1,1+4),xx(1+4)) !K1
     call midk(kk,ee,xx, 1,2,  kk(1,2+4),xx(2+4)) !K2
     call midk(kk,ee,xx, 1,3,  kk(1,3+4),xx(3+4)) !K3
     call inttetrac6(kkv,kk,xx, (/3, 4,3+4, 2/),  &! &  k3,k4,K3,k2
          frhis, nwhis, &
          wtthis )
     call inttetrac6(kkv,kk,xx, (/4,1+4,2+4,3+4/), & ! &  k4,K1,K2,K3
          frhis, nwhis, &
          wtthis )
     call inttetrac6(kkv,kk,xx, (/4, 2,2+4,3+4/), & ! &  k4,k2,K2,K3
          frhis, nwhis, &
          wtthis )
  else
     if(chkwrt) write(6,*) 'inttetra5: fig 4'
     call inttetrac6(kkv,kk,xx, (/1,2,3,4/),  &! &  k1,k2,k3,k4
          frhis, nwhis, &
          wtthis )
  endif
end subroutine inttetra6

!-----------------------------------
subroutine inttetrac6(kkv,kk,xx, itetx, &
     frhis, nwhis, &
     wtthis )
  !- Calculate tetrahedron integral Eq.(16).
  ! The four corners and denoted by itetx(1:4).
  !i kk (k1-k4) and xx (value of denominator)
  !i Kx (K1-K4) and xkx
  !r The four corners are selected from 8 points. It is specified by itetx(1:4).
  !r wtthis is accumulating.
  implicit none
  integer:: itetx(4),ix,i
  real(8) ::  kk(3,1:8),xx(1:8), am(3,3) &
       ,kin(3,4), xin(4), det33, intttvv,  intttv,intvv,work !Kx(3,1:4),xKx(1:4),
  integer:: nwhis
  real(8):: frhis(nwhis+1),wtthis(nwhis,4),voltet,kkv4bm(3),bm(3,3)
  logical:: chkwrt=.false.,matrix_linear
  real(8) ::  kkv(3,4), kkvkin(4,4)

  if(chkwrt) print *, ' inttetrac6: === '
  xin(    1:4) = xx(    itetx(1:4))
  kin(1:3,1:4) = kk(1:3,itetx(1:4))
  do i = 1,3
     am(1:3,i) = kin(1:3,i) - kin(1:3,4)
  enddo
  voltet = abs(det33(am)/6d0) ! \omega (volume of tetrahedra) = abs(det33(am)/6d0) See Eq. (17).

  !--- kkvkin: kin is decomplosed into kkv.!!!!!!!!!!
  !   e.g. kin (:,1) =  \sum_i kkv(:,i) * kkvkin(i,1)
  if(matrix_linear()) then
     do i = 1,3
        am(1:3,i) = kkv(1:3,i) - kkv(1:3,4)
     enddo
     call minv33tp(am, bm)
     kkv4bm(1)= sum( kkv(:,4)*bm(:,1) )
     kkv4bm(2)= sum( kkv(:,4)*bm(:,2) )
     kkv4bm(3)= sum( kkv(:,4)*bm(:,3) )
     do i=1,4
        kkvkin(1,i)= sum( kin(:,i)*bm(:,1) ) - kkv4bm(1)
        kkvkin(2,i)= sum( kin(:,i)*bm(:,2) ) - kkv4bm(2)
        kkvkin(3,i)= sum( kin(:,i)*bm(:,3) ) - kkv4bm(3)
        kkvkin(4,i)= 1d0 - sum(kkvkin(1:3,i))
     enddo
  endif

  call intttvc6(kkvkin, xin,voltet, &
       frhis, nwhis, &
       wtthis )
end subroutine inttetrac6

!---------------------------------------------------------------------------------------------
subroutine intttvc6(kkvkin, v,voltet, &
     frhis, nwhis, &
     wtthis ) !this is accumlating variable
  !- Histgram weights for each bin.
  !r The i-th bin of the Histgrams is [frhis(i) frhis(i+1)].
  !r  wtthis(ihis) += \int_{frhis(ihis)^{frhis(ihis+1)} d \omega  \int d^3k \delta(\omega +v(k) )
  !r  Total sum of the histgram weight is equal to voltet.
  !o  wtthis(ihis): Note this is accumulating variable.

  !r Note wtthis(nwhis,0:3)
  !r if matrix_linear()=T, wtthis is calculated assuming
  !r                       the linear-dependency of matrix elements. takao /dec/2003
  implicit none
  integer::  inihis, iedhis, nwhis,ichk=2313
  integer :: ieaord(1:4),i,isig,idif(3),idf,ix,n,itmp,ihis
  real(8)   ::  v(4), voltet, WW(4),  norm,s1,s2,pvn, &
       integb3p, integb3m , integb2p, integb2m,ww2p,ww3p
  real(8):: frhis(nwhis+1),wtthis(nwhis,4),intega,integb,stot,xxx,wx
  logical ::chkwrt=.false., matrix_linear
  real(8):: kkvkin(4,4),www(1:4)=.25d0,ec,wcg(4),wttt
! #ifdef EXPAND_SORTEA
!   n=4
!   !      isig = 1
!   ! option noparallel
!   do i = 1,n
!      ieaord(i) = i
!   enddo
!   ! option noparallel
!   do ix= 2,n
!      ! option noparallel
!      do i=ix,2,-1
!         if( -v(ieaord(i-1)) > -v(ieaord(i) ) ) then
!            itmp = ieaord(i-1)
!            ieaord(i-1) = ieaord(i)
!            ieaord(i) = itmp
!            !            isig= -isig
!            cycle
!         endif
!         exit
!      enddo
!   enddo
! #else
  call sortea( -v,ieaord,4,isig)
!#endif
  WW(1:4) = -v( ieaord(1:4) )

  !  ww(1)<ww(2)<ww(3)<ww(4)
  if(( .NOT. (WW(1)<=WW(2))) .OR. ( .NOT. (WW(2)<=WW(3))) .OR. ( .NOT. (WW(3)<=WW(4))) ) then
     write(6,"(/,' --- intttvc6: wrong order WW=',4d14.6)") WW
     ! top2rx 2013.08.09 kino        stop 'intttvc6: wrong order of WW'
     call rx( 'intttvc6: wrong order of WW')
  endif

  if(chkwrt) then
     write(ichk,"(/,' --- intttvc6: e=',4d23.16)") WW
  endif

  inihis= -999
  iedhis= -999
  ix=1
  do ihis = 1,nwhis
     if(ix==1 .AND. WW(ix)<frhis(ihis+1)) then
        inihis = ihis
        ix=4
     endif
     if(ix==4 .AND. WW(ix)<frhis(ihis+1)) then
        iedhis = ihis
        exit
     endif
  enddo
  if(iedhis==-999 .OR. inihis==-999) then
     print *,' intttvc6: can not find inihis iedhis'
     ! top2rx 2013.08.09 kino        stop    ' intttvc6: can not find inihis iedhis'
     call rx( ' intttvc6: can not find inihis iedhis')
  endif
  !      if(chkwrt) print *,' inihis iedhis=', inihis,iedhis
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      print *,' ###########integtetn test #########'
  !      ww(1:4)=(/1d0,1.5d0,9.5d0,10d0/)
  !      do ix=1,101
  !        wx = (ix-1)/10d0
  !        call integtetn(WW, Wx, xxx)
  !        write(61,"(i3,2d23.16)")ix,wx,xxx
  !      enddo
  !      stop 'test end integtetn---'
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  call integtetn(WW, WW(4), norm)
  if(chkwrt) write(ichk,"(' norm voltet=',2d24.16)") norm,voltet
  pvn = voltet/norm
  if(chkwrt) write(ichk,"(' norm voltet pvn=',3d14.6)") norm,voltet,pvn
  intega = 0d0
  do ihis = inihis, iedhis
     if( frhis(ihis+1)>ww(4) ) then
        integb = norm
     else
        call integtetn(WW, frhis(ihis+1), integb)
     endif

     if(matrix_linear()) then !-- Weight for each conrner sum(wcg)=1d0
        ec= (frhis(ihis)+frhis(ihis+1))/2d0
        call mkwcg(WW, ec, wcg)
        do ix=1,4
           www(ix)= sum(kkvkin(ix,ieaord(1:4))*wcg(1:4))
        enddo
        !          write(6,"('sum(www)= ',d13.4,' www(1:4)=',4f10.4)") sum(www),www(1:4)
        wttt = pvn*(integb - intega)*4d0
        wtthis(ihis,:) = wtthis(ihis,:) + wttt * www(:)
     else
        wtthis(ihis,1) = wtthis(ihis,1) + pvn*(integb - intega)
     endif

     if(chkwrt) then
        write(ichk,"(' ihis [Init End] wtt=', i5,3f11.6)") &
             ihis, frhis(ihis), frhis(ihis+1), pvn*(integb-intega)
     endif

     intega = integb
  enddo
  !      if(chkwrt) then
  !        dini = 0d0
  !        do ihis=inihis,iedhis
  !          if(ihis>1) dini=wtthis(ihis-1)
  !          write(6,"(' ihis [Init  End]', i4,2f13.6,' wtt=',f13.6)")
  !     &    ihis,dini,frhis(ihis),wtthis(ihis)
  !          dini = frhis(ihis)
  !        enddo
  !      endif
  if(chkwrt) write(ichk,*) ' end of intttvc6'
end subroutine intttvc6


!--------------------------------------------------------------
subroutine integtetn(e, ee, integb)
  !> Calculate primitive integral of integb = 1/pi Imag[\int^ee dE' 1/(E' -e(k))] = \int^ee dE' S[E']
  !! \remark
  !!  S[E] : is area of the cross-section between the omega-constant plane and the tetrahedron.
  !!  [here we assumee e1<e2<e3<e4].
  !!  Normalization is not considered!
  !!  Rath&Freeman Integration of imaginary part on Eq.(17)
  implicit none
  integer:: dummy4doxygen
  real(8) ::  e(1:4), ee, norm, integb,a,b, &
       V1,V2,V3,V4,D1,D2,D3,D4,e1,e2,e3,e4
  integer::i1,i2,i3
  !-------------------------------------------------------------
  !      if((.not.(e(1)<=e(2))).or.(.not.(e(2)<=e(3))).or.(.not.(e(3)<=e(4))) ) then
  !        write(6,"(/,' --- integtetn wrong order of e=',4d14.6)") e
  !        call rx('integtetn: wrong order of e')
  !      endif
  !      write(6,"(' --- integtetn e =',4d23.16)") e
  !      write(6,"(' --- integtetn ee=',d23.16)") ee
  if ( ee<e(1) ) then  !jan2008 ee<=e(1) ---> ee<e(1)
     integb = 0d0
     return
  elseif (ee>e(4) ) then
     call rx( ' integtetn: ee>e(4)')
  endif

  !--- case1. poor numerical accuracy for cases.
  if( .FALSE. ) then
     e1 = e(1)-3d-6
     e2 = e(2)-2d-6
     e3 = e(3)-1d-6
     e4 = e(4)
     V1= ee - e1
     V2= ee - e2
     V3= ee - e3
     V4= ee - e4
     D2 = (V2-V4)*(V2-V3)*(V2-V1) !<0
     D3 = (V3-V4)*(V3-V2)*(V3-V1) !>0
     D1 = (V1-V4)*(V1-V3)*(V1-V2) !>0
     if( e1<=ee ) integb =          V1**3/D1
     if( e2<=ee ) integb = integb + V2**3/D2
     if( e3<=ee ) integb = integb + V3**3/D3
     return
  endif

  !--- case2
  e1 = e(1)-3d-8
  e2 = e(2)-2d-8
  e3 = e(3)-1d-8
  e4 = e(4)
  V1= ee - e1
  V2= ee - e2
  V3= ee - e3
  V4= ee - e4
  !      D1 = (V1-V4)*(V1-V3)*(V1-V2) !>0
  !      D2 = (V2-V4)*(V2-V3)*(V2-V1) !<0
  !      D3 = (V3-V4)*(V3-V2)*(V3-V1) !>0
  !      D4 = (V4-V1)*(V4-V2)*(V4-V3) !<0
  if    ( e1<=ee .AND. ee<e2 ) then
     D1 = (e4-e1)*(e3-e1)*(e2-e1) !>0
     integb =  V1**3/D1
  elseif( e2<=ee .AND. ee<e3 ) then
     D1 = (e4-e1)*(e3-e1)*(e2-e1) !>0
     D2 = (e4-e2)*(e3-e2)*(e1-e2) !<0
     a  =  V1/   D1**(1d0/3d0)
     b  =  V2/(-D2)**(1d0/3d0)
     integb = (a-b) * (a**2+a*b+b**2)
  elseif( e3<=ee .AND. ee<e4 ) then
     D4 = (e1-e4)*(e2-e4)*(e3-e4) !<0
     integb = 1d0 - V4**3/D4
  elseif( ee==e4 ) then
     integb = 1d0
  endif

  !-----------------------------------------------------------
  !      D2 = (V2-V4)*(V2-V3)*(V2-V1) !>0
  !      D3 = (V3-V4)*(V3-V2)*(V3-V1) !<0
  !      integb = 0d0
  !      if( e(1)<=ee ) integb =          V1**3* D2*D3
  !c      write(6,*) ' integb1=',integb
  !      D1 = (V1-V4)*(V1-V3)*(V1-V2) !<0
  !      if( e(2)<=ee ) integb = integb + V2**3* D3*D1
  !c      write(6,*) ' integb2=',integb
  !      if( e(3)<=ee ) integb = integb + V3**3* D1*D2
  !c      write(6,*) ' integb3=',integb
  !--
  !      if(V4==0d0) then
  !       write(6,*)
  !       write(6,"(' V  =',4d18.10)") V1, V2, V3, V4
  !       write(6,"(' int=',3d13.5,16x,d13.5)") V1**3/D1,V2**3/D2,V3**3/D3,
  !     &             max(abs(V1**3/D1),abs(V2**3/D2),abs(V3**3/D3))
  !         write(6,*) ' integb=',integb
  !      endif
end subroutine integtetn


!-----------------------------------------------------
subroutine mkwcg(e, ee, wcg)
  !- calculate wweight for each corners.---------
  !i e(1:4),ee
  !o wcg
  !r normarization is sum(scg(1:4))=1d0
  !----------------------------------------------
  implicit none
  real(8) ::  e(1:4), ee, wcg(1:4),e1,e2,e3,e4, &
       w14_1,w14_4,w12_1,w12_2,w13_1,w13_3,w23_2,w23_3,w24_2,w24_4,w34_3,w34_4
  !-------------------------------
  ! cccccccccccccccccccccccccccccccccccccc
  ! This is for original case thru Dec 2003.
  !      wcg=.25d0
  !      return
  ! cccccccccccccccccccccccccccccccccccccc
  e1 = e(1)-3d-8
  e2 = e(2)-2d-8
  e3 = e(3)-1d-8
  e4 = e(4)
  if    ( ee<=e(1) ) then
     wcg(1)  =1d0
     wcg(2:4)=0d0
  elseif( e1<=ee .AND. ee<e2 ) then
     call wab(e1,e2,ee,w12_1,w12_2)
     call wab(e1,e3,ee,w13_1,w13_3)
     call wab(e1,e4,ee,w14_1,w14_4)
     wcg(1) = (w12_1 +w13_1 + w14_1)/3d0
     wcg(2) =  w12_2/3d0
     wcg(3) =  w13_3/3d0
     wcg(4) =  w14_4/3d0
  elseif( e2<=ee .AND. ee<e3 ) then !Is this correct?
     call wab(e1,e3,ee,w13_1,w13_3)
     call wab(e1,e4,ee,w14_1,w14_4)
     call wab(e2,e3,ee,w23_2,w23_3)
     call wab(e2,e4,ee,w24_2,w24_4)
     wcg(1) = (w13_1 + w14_1)/4d0
     wcg(2) = (w23_2 + w24_2)/4d0
     wcg(3) = (w13_3 + w23_3)/4d0
     wcg(4) = (w14_4 + w24_4)/4d0
  elseif( e3<=ee .AND. ee<e4 ) then
     call wab(e1,e4,ee,w14_1,w14_4)
     call wab(e2,e4,ee,w24_2,w24_4)
     call wab(e3,e4,ee,w34_3,w34_4)
     wcg(1) =  w14_1/3d0
     wcg(2) =  w24_2/3d0
     wcg(3) =  w34_3/3d0
     wcg(4) = (w24_4 +w34_4 + w14_4)/3d0
  elseif( ee> e(4) ) then
     wcg(4)  =1d0
     wcg(1:3)=0d0
  endif
  ! ccccccccccccccccccccccc
  !      write(6,"('mkwcg   e=',4d13.4, ' ee=',d13.4)")e,ee
  !      write(6,"('      wcg=',4d13.4)")wcg
  ! cccccccccccccccccccccc
end subroutine mkwcg
!-------------------------------------
subroutine wab(ea,eb,ee,wa,wb)
  implicit none
  real(8)::ea,eb,wa,wb,eet,ee
  eet= eb - ea
  wa= (eb-ee)/eet
  wb= (ee-ea)/eet
end subroutine wab


!---
subroutine addsciss(delta, ef, nnn, eig)
  integer::nnn,i
  real(8):: eig(nnn),ef,delta
  write(6,*)' asssciss delta=', delta
  do i=1,nnn
     if(eig(i)>ef) eig(i)= eig(i)+delta
  enddo
end subroutine addsciss

!-----------------------------------
subroutine midk3(kk,ee,xx,yy,i,j,   kout,xout,yout)
  !- Calculate x and k(3) at the Fermi energy on the like k(i)---k(j).
  implicit none
  integer:: i,j
  real(8) ::  kk(3,1:4),xx(1:4),yy(1:4),ee(1:4), kout(3) &
       ,xout,yout,ratio
  ratio     = ee(i)/(ee(i)-ee(j))
  xout      = xx(i)     + ratio * (xx(j)-xx(i))
  yout      = yy(i)     + ratio * (yy(j)-yy(i))
  kout(1:3) = kk(1:3,i) + ratio * (kk(1:3,j)-kk(1:3,i))
end subroutine midk3

real(8) function det33(am)
  implicit none
  real(8),intent(in) :: am(3,3)
  det33= am(1,1)*am(2,2)*am(3,3) &
       -am(1,1)*am(3,2)*am(2,3) &
       -am(2,1)*am(1,2)*am(3,3) &
       +am(2,1)*am(3,2)*am(1,3) &
       +am(3,1)*am(1,2)*am(2,3) &
       -am(3,1)*am(2,2)*am(1,3)
END function det33

subroutine chkdgn(ene,ndat,  nrank,ixini,ixend,iof,ipr)
  implicit none
  integer :: ndat,i,ix, ixini(ndat),ixend(ndat),nrank,iof
  real(8)    :: ene(ndat), epsx=1d-4
  logical :: ipr
  if(ipr)  write(6,*) 'chgdgn: ndat=',ndat
  if(ndat<1) then
     nrank =0
     return
  endif

  ixini(1) = 1
  if(ndat==1) then
     ixend(1) = 1
     nrank=1
     return
  endif
  i = 1
  ! option noparallel
  do ix = 2, ndat
     if( abs(ene(ix)-ene(ix-1)) >epsx ) then
        ixend(i) = ix-1
        i = i + 1
        ixini(i) = ixend(i-1)+1
        if(ix==ndat) then
           ixend(i)=ix
        endif
     elseif(ix==ndat) then
        ixend(i) = ndat
     endif
  enddo

  nrank = i
  ! option noparallel
  do i =1,nrank
     ixini(i) = ixini(i)+iof
     ixend(i) = ixend(i)+iof
  enddo
  ! check write
  if(ipr) then
     write(6,*)' nrank=',nrank
     do i = 1, ndat
        write(6,"(' i ',i3,' ene=',d15.7)") i,ene(i)
     enddo
     print *
     do i = 1, nrank
        write(6,"(' i ',2i3,' e=',d15.7)") &
             ixini(i),ixend(i),ene(ixini(i)-iof)
     enddo
  endif
end subroutine chkdgn

subroutine midk(kk,ee,xx,i,j,   kout,xout)
  !     - Calculate x and k(3) at the Fermi energy on the like k(i)---k(j).
  integer:: i,j
  real(8) :: ratio, kk(3,1:4),xx(1:4),ee(1:4), kout(3),xout
  ratio     = ee(i)/(ee(i)-ee(j))
  xout      = xx(i)     + ratio * (xx(j)-xx(i))
  kout(1:3) = kk(1:3,i) + ratio * (kk(1:3,j)-kk(1:3,i))
end subroutine midk

subroutine rsvwwk00_4(jpm,iwgt, nqbz,nband,nctot,ncc,nbnbx,  n1b,n2b,noccxv,nbnb)
  !- get (n1b n2b) corresponding to non-zero wgt.
  implicit none
  integer :: jpm,nband, nctot,  ncc, nqbz,nbnbx ,ib,jb, kx,ix
  integer :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),noccxv,nbnb(nqbz)
  logical :: iwgt(nband+nctot,nband+ncc,nqbz)
  noccxv = 0
  do kx  = 1, nqbz
     ix  = 0
     do ib  = 1, nband + nctot
        do jb  = 1, nband + ncc
           if( iwgt(ib,jb,kx)) then
              ix          = ix+1
              n1b(ix, kx) = ib
              n2b(ix, kx) = jb
              if(jpm==1 .AND. ib<=nband .AND. ib>noccxv ) noccxv =ib
              if(jpm==2 .AND. jb<=nband .AND. jb>noccxv ) noccxv =jb
           endif
        enddo
     enddo
     nbnb(kx) = ix
  enddo
  write(6,*) ' rsvwwk: kx nbnbmax=',kx, maxval(nbnb)
end subroutine rsvwwk00_4

!! ---------------------------------
SUBROUTINE TETFBZF_notused(qbas,N1,N2,N3,rk,nqbz, skipgammacell,half,nadd, &
     IDTET,qbzw,ib1bz,ntetf)
  use m_keyvalue,only: getkeyvalue
  !!  Finds tetrahedra in all 1st BZ. takao mod. from tertirr
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i  qb,n1,n2,n3,ipq, output from BZMESH;
  !i  nq, no. of irreducible k-points;
  !o Outputs:
  !o  ntet, No. of different tetrahedra
  !o  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
  !o  idtet(0,i), no. of tetrahedra of the i'th kind
  !m Memory:
  !m  No large internal storage; heap not accessed.
  !r    This require subroutine CCUTUP (lmto-3).
  ! ----------------------------------------------------------------------
  implicit none
  integer:: indexkw(0:n1,0:n2,0:n3),kount, &
       indexk(0:n1,0:n2,0:n3), &
       i,i1,i2,i3,j1,j2,j3, IPQ(N1,N2,N3), &
       IBTR(3,3),KCUT(3,4,6),IMC(0:1,0:1,0:1), &
       idtet(4, 6*n1*n2*n3),iq(4), &
       ntet,k1,k2,k3,itet,ic &
       ,ib1bz((n1+1)*(n2+1)*(n3+1)) !icase=1 fixed
  real(8) :: QB(3,3),QB1(3,3), qbas(3,3), qbzx(3) &
       , qbzw(3,(n1+1)*(n2+1)*(n3+1)),half(3)
  real(8),parameter :: epss = 1d-12
  logical,save ::chk=.true.
  logical:: skipgammacell,qbzreg

  real(8),intent(in) :: rk(3,nqbz)
  integer,intent(in):: n1,n2,n3,nqbz,nadd
  integer,intent(out)::ntetf

  integer:: ngcell,i1x,i2x,i3x,ndx
  real(8):: qqqx(3),rb(3,3),rrr(3),xvec(3),qh(3)
  !------------------------------------------------------------------
  !      hf=0d0
  !      if(icase==2) hf=0.5d0
  QB(1:3,1) = QBAS(1:3,1)/N1
  QB(1:3,2) = QBAS(1:3,2)/N2
  QB(1:3,3) = QBAS(1:3,3)/N3
  qh=matmul(qbas,half)
  call minv33(qb,rb)
  !- index for k in 1st BZ. See genqbz in BZ.FOR.
  kount      = 0
  do      i1 = 1,n1+nadd
     do      i2 = 1,n2+nadd
        do      i3 = 1,n3+nadd
           kount    = kount + 1
           indexk(i1-1,i2-1,i3-1) = kount
           xvec = (/I1-1,I2-1,I3-1/)
           if(chk) then
              !  qbzx(1:3)= qb(1:3,1)*(i1-1+hf) +qb(1:3,2)*(i2-1+hf) +qb(1:3,3)*(i3-1+hf)
              qbzx(1:3)= matmul(qb,xvec) + qh
              ! ccccccccccccccccccc
              !      write(6,"(3d26.18)") rk(1:3,kount)
              !      write(6,"(3d26.18)") qbzx(1:3)
              !      write(6,*)
              ! ccccccccccccccccccc
              if( sum(abs(rk(1:3,kount)-qbzx(1:3))) > epss ) call rx( 'tetfbzf: rk /= qbzx')
           endif
        end do
     end do
  end do
  chk=.false.
  if (kount /= (n1+nadd)*(n2+nadd)*(n3+nadd)) call rx( ' kount: wrong no. k-points 111')
  if (nqbz  /= kount ) call rx( ' kount: wrong no. k-points 222')
  ! cccccccccccccccccccccccccc
  !      DO  j1 = 0,n1-1
  !      DO  j2 = 0,n2-1
  !      DO  j3 = 0,n3-1
  !        write(6,"(' j1j2j3=',3i4,' ix=',i6)") J1,J2,J3,indexk(J1,J2,J3)
  !      enddo
  !      enddo
  !      enddo
  !      write(6,*)
  ! cccccccccccccccccccccccccc
  ! ccccccccccccc
  if(nadd==0) then
     kount      = 0
     do      i1 = 1,n1+1
        do    i2 = 1,n2+1
           do  i3 = 1,n3+1
              kount    = kount + 1
              indexkw(i1-1,i2-1,i3-1) = kount
              qbzw(1:3,kount) = &
                   qb(1:3,1)*(i1-1+half(1)) + qb(1:3,2)*(i2-1+half(2)) + qb(1:3,3)*(i3-1+half(3))
              ib1bz(kount) = indexk(mod(i1-1,n1), mod(i2-1,n2), mod(i3-1,n3))
           end do
        end do
     end do
  elseif(nadd==1) then
     !        qbzw(:,1:nqbz)= rk(:,1:nqbz) !NOTE: nqbz=kount=(n1+nadd)*(n2+nadd)*(n3+nadd)
     kount      = 0
     do     i1 = 1,n1+1
        do   i2 = 1,n2+1
           do i3 = 1,n3+1
              kount    = kount + 1
              indexkw(i1-1,i2-1,i3-1) = kount
              qbzw(1:3,kount) = &
                   qb(1:3,1)*(i1-1+half(1)) + qb(1:3,2)*(i2-1+half(2)) + qb(1:3,3)*(i3-1+half(3))
              ib1bz(kount) = indexk(i1-1, i2-1, i3-1)
              ! cccccccccccccccccccccc
              !              qqqx = qbzw(1:3,kount) -
              !     &         (qb(1:3,1)*(i1-1+hf) + qb(1:3,2)*(i2-1+hf) + qb(1:3,3)*(i3-1+hf))
              !              if(sum(abs(qqqx))>1d-6) then
              !                write(6,*) 'qbzw=',qbzw(1:3,kount)
              !                write(6,*) 'qbx=',(qb(1:3,1)*(i1-1+hf) + qb(1:3,2)*(i2-1+hf) + qb(1:3,3)*(i3-1+hf))
              !                call rx('qbzw error')
              !              endif
              ! cccccccccccccccccccccc
           enddo
        enddo
     enddo
  else
     call rx('tetfbzf:wrong nadd')
  endif

  CALL CCUTUP(QB,QB1,IBTR,KCUT) !This is from LMTO-3
  ntet = 0
  ! ----- START LOOPING OVER MICROCELLS ---------
  call getkeyvalue("GWinput","ngcell",ngcell,default=1)

  DO 20  I3 = 1, N3
     DO 21  I2 = 1, N2
        DO 22  I1 = 1, N1
           if(skipgammacell) then
              i1x = i1-n1
              i2x = i2-n2
              i3x = i3-n3
              if(i1<n1/2) i1x=i1
              if(i2<n2/2) i2x=i2
              if(i3<n3/2) i3x=i3
              ndx=0
              if(qbzreg()) ndx=1
              if( i1x<ngcell .AND. -ngcell< i1x-ndx) then
                 if( i2x<ngcell .AND. -ngcell< i2x-ndx) then
                    if( i3x<ngcell .AND. -ngcell< i3x-ndx) then
                       kount= indexkw(i1-1,i2-1,i3-1)
                       rrr=matmul(rb,qbzw(1:3,kount))
                       write(6,"('qqqqx qbzw=',3i3, 3f9.5)") i1x,i2x,i3x, rrr
                       ! cccccccccccccccccccccccc
                       !             goto 3010
                       ! cccccccccccccccccccccccc
                       cycle
                    endif
                 endif
              endif
           endif

           ! ccccccccccccccccccc
           !        cycle
           ! 3010   continue
           ! cccccccccccccccccccc


           !         write(6,"(' k1k2k3=',3i4,' ix=',i6,' ii=',2i6)")
           !     &    indexkw(J1,J2,J3)
           !     &    ,indexkw(J1,J2,J3), ib1bz(indexkw(J1,J2,J3))
           !     &    ,indexk(mod(j1,n1), mod(j2,n2), mod(j3,n3))
           ! cccccccccccc
           ! ----- SET UP IDENTIFIERS AT 8 CORNERS OF MICROCELL ------
           !        write(6,*)
           !        print *,'xxxxxxxx qqqq xxxxxxxx ',i1,i2,i3
           do K1 = 0, 1
              J1 = I1 -1 + K1
              do K2 = 0, 1
                 J2 = I2 -1 + K2
                 do K3 = 0, 1
                    J3 = I3 -1 + K3
                    IMC(K1,K2,K3)  = indexkw(J1,J2,J3)
                    ! cccccccccccc
                    !         if(k1==0.and.k2==0.and.k3==0) then
                    !           kount= indexkw(j1,j2,j3)
                    !           write(6,"('qqqq qbzw=',3i3,3f9.5)") j1,j2,j3,matmul(rb,qbzw(1:3,kount))
                    !         endif
                    ! cccccccccccc
                 enddo
              enddo
           enddo

           ! ccccccccccccccccccccccccc
           !         if(i1==1.and.i2==1.and.i3==1) then
           !           goto 4010
           !         endif

           !         if(i1==4.and.i2==4.and.i3==4) then
           !           goto 4010
           !         endif
           ! ccccccccc
           !        write(6,"('ffff =',3i3)") i1,i2,i3

           !         if(i1==1.and.i2==1.and.i3==n3) then
           !          goto 4010
           !        endif

           !         if(i1==4.and.i2==4.and.i3==1) then
           !           goto 4010
           !         endif

           !         cycle
           ! 4010    continue
           ! cccccccccccccccccccccccc

           ! ----- LOOP OVER TETRAHEDRA --------------
           ! ccccccccccccccccccccccccc
           do 10 ITET = 1, 6
              !        do 10 ITET = 3,3
              ! ccccccccccccccccccccccccc
              do  IC = 1, 4
                 K1 = KCUT(1,IC,ITET)
                 K2 = KCUT(2,IC,ITET)
                 K3 = KCUT(3,IC,ITET)
                 IQ(IC) = IMC(K1,K2,K3)
              enddo
              ntet=ntet+1
              do  i = 1, 4
                 idtet(i,ntet) = iq(i)
              enddo
10         enddo
           write(6,"('ntet=',3i3,2x,i6)")i1,i2,i3,ntet
22      enddo
21   enddo
20 enddo
  write(6, "(1x,'TETFBZF: ',2i8 )")  ntet, 6*n1*n2*n3
  ntetf=ntet
end SUBROUTINE TETFBZF_notused