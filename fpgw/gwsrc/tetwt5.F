      subroutine  tetwt5x_dtet4(npm, ncc,
     i    q, eband1,eband2,
     i    qbas,ginv,efermi,   
     i    ntetf, nqbzw, nband,nqbz,
     i    nctot,ecore,  idtetf,qbzw,ib1bz,
     i   job,
     o   iwgt,nbnb,           !job=0
     o   demin, demax,       !job=0   Note these are real(4)
     i          frhis,nwhis,
     i          nbnbx,ibjb,nhwtot, ihw,nhw,jhw, !job=1
     o          whw,                            !job=1
     i   iq,isp1,isp2,nqibz,eibzmode,nwgt) ! new input for mtett
      use m_keyvalue,only: getkeyvalue
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

C-----------------------------------
Cr ek   = eigenvalues at k-points in the 1st BZ
Cr ekq  = eigenvalues at k+q, k   in the 1st BZ
Cr
Cr See J.Rath&A.J.Freeman PRB11(6) p.2109(1975).
Cr
Cr subroutine lindtet6 is the main part to calculate
Cr the microcell integral of wtthis(ihis) = \int_lower(ihis)^upper(ihis) d\omega \int d^3k f(e(k)) (1-f(e(q+k))) \delta (omg- e(q+k) + e(k) )
Cr where e(q+k) is unoccupied, e(k) is occupied,
Cr and f(E) denote the Fermi distribution funciton.
Cr
Cr
Crr The numbering of band index is not unique due to degeneracy.
Crr It affects how to choose tetrahedron.  As a result, it affect on
Crr the integration weight, and it might break the crystal symmetry.
Crr So I add symmetrization at the last in this routine so as to recover the symmetry.
Cr
c-----------------------------c
c takao kotani Aug 2000.
c Sergey Faleev 2002
c takao kotaniMar 2002
c takao nov2003; mtet mode.
c     now only for mtet=(/1,1,1/) normal mode, or
c              for mtet=(/2,2,2/) mode (each tetrahedron is devided into 2x2x2 tetrahedrons.
c
c takao dec2003; matrix_linear() mode
      implicit none
      integer:: dummy4doxygen
      integer(4):: npm,jpm,ibxmx,jbxmx,jbx,nrankc1,nrankc2,nnn1,nnn2,nnni,nnnj,ncc

c---in out -------------------------------
      integer(4) :: nband,nqbz,nctot,ntetf,nqbzw,
     &   idtetf(0:3,ntetf),ib1bz(nqbzw), nbnb(nqbz,npm)
      real(8) :: q(3),
     &   eband1(nband,nqbz), eband2(nband,nqbz),
     &   qbas(3,3),ginv(3,3),efermi,
     &   ecore(nctot), qbzw(3,nqbzw) !(n1+1)*(n2+1)*(n3+1))
      complex(8):: omg
      logical :: iwgt(nband+nctot,nband+ncc,nqbz,npm)
c----------------------------------
      integer(4):: itet,ic, ib,jb
      real(8) :: !ek(nband,0:3), ekq(nband,0:3),
     &    qk(3),qkm(3),qbz(3)  !,rk(3,nqbz)
C- For tetrahedra
      integer(4):: kk(0:3),kq(0:3),kr(0:3),i ,j
      real(8)   ::  det33 , 
     &   kvec(3,0:3), ea(0:3), eb(0:3) ,x(0:3),am(3,3)
      complex(8):: wtt(0:3,3)
c
      integer(4):: noccx_kxx, noccx_k, noccx_kq, noccx1
c      real(8),parameter:: eps=0d0 !1d-12 ! cutoff check to determine cancellation.
      logical ::prt=.false.
      real(8),target :: ek_(nband+nctot,0:3), ekq_(nband+nctot,0:3)
c
      integer(4):: irnk1, nrank1, irnk2, nrank2,nibib,kx,
     &  ires((nband+nctot)**2), iof1,iof2,
     &  ini1(nband+nctot),ied1(nband+nctot), ixi1,ixi2,ixe1,ixe2,
     &  ini2(nband+nctot),ied2(nband+nctot)
      real(8):: ekqxx(nband+nctot),summ,voltot,volt
      complex(8)::  wmean,a1,a2,wgt2
      integer(4):: ik,ibx,nbnc,nb1 !kqxx(nqbz),
      logical :: ipr=.false.
      real(8),parameter:: pii=3.1415926535897932d0
      real(8),pointer:: eocc(:),eunocc(:)
      integer:: idim,ivec
c      complex(8),allocatable:: wgt1(:,:)

CCCCCCCCCCCCCCC faleev Jan 2002
      real(8)::  x_(0:3),cut

ctakao-Feb/2002
      integer(4):: job
c job=0
c      logical :: trans(0:3)


      real(4) :: demax(nband+nctot,nband+ncc,nqbz,npm),demax_,demaxx
     &          ,demin(nband+nctot,nband+ncc,nqbz,npm),demin_,deminn
c      real(8) :: demax(nband+nctot,nband,nqbz),demax_,demaxx
c     &          ,demin(nband+nctot,nband,nqbz),demin_,deminn

      integer(4)::ixx
      integer(4)::tmp1,tmp2,tmp3,tmp4
c job=1
      integer(4):: nbnbx,nhwtot,nwhis,inihis,ihis,ikx,ibib, ikk
     &            ,ioff,isum,ini,ied,jini,iini,nnn,    ntetmx,ntetmin
      integer(4):: ihw(nbnbx,nqbz,npm), ! omega pointer
     &             nhw(nbnbx,nqbz,npm), ! number of data
     &             jhw(nbnbx,nqbz,npm), ! histo-weight pointer for whw(*)
     &             ibjb(nctot+nband,nband+ncc,nqbz,npm)
      real(8) :: whw(nhwtot)   ! histo-weight
     &   , frhis(nwhis+1), wtthis2(nwhis,0:3), wtthis(nwhis),piofvoltot
      logical ::chkwrt=.false.,wxx
      real(8),parameter:: pi=3.1415926535897932d0
c
c-devided tet
      integer(4):: nmtet,nqbzwm,nqbzm,ntetfm
      integer(4),allocatable:: idtetfm(:,:,:),ib1bzm(:) !,index_qbzm(:,:,:)
      real(8),allocatable:: qbzm(:,:),qbzwm(:,:)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

      integer(4)::kkm(0:3),kvecm(1:3, 0:3),ifeig,ifmtet,nqnumm,kqxxm,kqxxm2
     &  ,n,nsp,nband_x,kp,im !,n_index_qbzm
      real(8),allocatable:: ekzz1(:,:),ekzz2(:,:),eigtet(:,:,:),wtet(:,:,:)
      real(8)::qdummy(3)

      integer(4)::isp1,isp2,iq,nqibz,mtet(3),  iqx,ispx, ix1,iy1,iz1,verbose,iqindx
      logical ::mtett !,readgwinput

      real(8):: kkv(3,0:3)
      logical :: matrix_linear
      real(8):: scissors_x0, ddw

      integer :: nbandmx_kq,nbandmx_k
c      real(8):: fermi_dS(3,3) !not yet used...

      logical:: eibzmode
      integer:: nwgt(nqbz)

c---------------------------------------------------------------------
      !print *, frhis, nwhis
      write(6,"('  tetwt5_dtet3:init job efermi:=',i2,d13.6)") job,efermi
!      print *, "nctot+nband, nband+ncc, nqbz,npm <-- ibjb allocate2"
!      print *, nctot+nband,nband+ncc,nqbz,npm
      !print *, "ibjb", ibjb
      !print *, "ibjb(6,5,27,1)", ibjb(6,5,27,1)
      !do tmp1 =1,npm
        !do tmp2 = 1, nqbz_mtet
          !do tmp3 = 
      !print *, "swjibjb...", ibjb
      !print *, "nband"
      !print *, nband
c      if(verbose()>=150) chkwrt=.true.
      voltot = abs(det33(qbas))
      piofvoltot = pi/4d0/voltot
c      fermi_dS=0d0
      if(job==0) iwgt = .false.

!! multitet=T ==> micro tetrahedron method (divided-tetrahedron). Not used so much now...
      mtet=(/1,1,1/)
      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
      if(sum(abs(mtet))/=3) then
        mtett=.true.
        if(sum(abs(mtet))<3) 
     &  print *, ' we use devided-tetrahedron scheme mtet=',mtet
      else
        mtett=.false.
      endif

!! read mtet or not.
      if(mtett) then
        ifmtet=501
        open (ifmtet, file='mtet',form='unformatted')
        read(ifmtet) nmtet,nqbzwm,nqbzm,ntetfm !,n_index_qbzm
        allocate( 
     &       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
        ! Index for tetrahedron;    qbzmw(idtetfm) gives extended q vector.
     &       ib1bzm(nqbzwm), qbzm(3,nqbzm) !,index_qbzm(n_index_qbzm,n_index_qbzm,n_index_qbzm)
        ! qbzm(1:3,ib1bz(iq)) gives q vector within the 1st bz.
     &       , wtet(0:3,nmtet,ntetf) )
        read(ifmtet) idtetfm,ib1bzm,qbzm,qbzwm,wtet  !,index_qbzm
        close(ifmtet)
        ifeig=501
        open (ifeig,file='eigmtet',form='unformatted')
        read(ifeig) nband_x,nqnumm,nsp
        print *,'readin eigmtet ', nband_x,nqnumm,nsp
        if(nband_x /=nband ) then
Cstop2rx 2013.08.09 kino          stop 'tetwt5: nband_x /=nband'
          call rx( 'tetwt5: nband_x /=nband')
        endif
        allocate( eigtet(nband,nqnumm,nsp) )
        do iqx=1,nqnumm
          do ispx=1,nsp
            read (ifeig) eigtet(1:nband,iqx,ispx)
            if(verbose()>150) write(6,"('iq=',i3,'  eig(1:5)=',5f10.4)")iqx,eigtet(1:5,iqx,ispx)
          enddo
        enddo
        close(ifeig)
      else
        nmtet = 1
        nqbzwm= nqbzw
        nqbzm = nqbz
        allocate( 
     &       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
     &       ib1bzm(nqbzwm),           qbzm(3,nqbzm) , wtet(0:3,nmtet,ntetf) )
        idtetfm(:,1,:)=idtetf
        qbzwm =qbzw
        ib1bzm=ib1bz
        wtet=0.25d0
      endif

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
      if(mtett) then
        print *, ' mtett mode nqbzm nqibz=',nqbzm,nqibz
        do kx = 1, nqbzm
          ekzz1( 1:nband, kx)              = eigtet(1:nband,kx,isp1)
          ekzz1( nband+1: nband+nctot, kx) = ecore(1:nctot)
          ekzz2( nband+1: nband+nctot, kx) = ecore(1:nctot)
        enddo
        if(iq<=nqibz) then
          do kx = 1, nqbzm
            kqxxm = iqindx(q(1:3)+qbzm(1:3,kx),ginv,qbzm,nqbzm)
            ekzz2(1:nband, kx) = eigtet(1:nband,kqxxm,isp2) !ekzz1(1:nband, kqxxm)
          enddo
        else
          do kx = 1, nqbzm
            kp  = nqbzm *(iq - nqibz) + kx
            ekzz2(1:nband, kx)= eigtet(1:nband,kp,isp2)
          enddo
        endif
        deallocate( eigtet)!,index_qbzm)
      else
        ekzz1=ekxx1
        ekzz2=ekxx2
      endif



!      print *, "nctot, nqbzm", nctot, nqbzm
!      print *, "nbnc", nbnc
!      print *, "eband1"
!      print *, eband1
!      print *, "ekzz1"
!      print *, ekzz1
!      print *, "eband2"
!      print *, eband2
!      print *, "ekzz2"
!      print *, ekzz2





!! Add eigenvalue shift. scissors_x0() is defined in switch.F
!!         call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
      if(scissors_x0()/=0d0) then
        call addsciss(scissors_x0(),efermi,(nband+nctot)*nqbzm, 
     o   ekzz1)
        call addsciss(scissors_x0(),efermi,(nband+nctot)*nqbzm, 
     o   ekzz2)
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
c        write(6,"('itet im vol=',2i5,d13.5)") itet,im,abs(det33(kvec(1:3,1:3))/6d0)
        enddo
      enddo
Cstop2rx 2013.08.09 kino      if(abs(volt-voltot)>1d-10) stop ' tetwt: abs(volt-voltot)>1d-10'
      if(abs(volt-voltot)>1d-10) call rx( ' tetwt: abs(volt-voltot)>1d-10')

      if(job==0) then
        demin=  1d10
        demax= - 1d10
      endif
      !print *, nwgt

!! === Loop over tetrahedron ===
      !print *, "swjq", q
      do 1000 itet = 1, ntetf !;
        !print *, "itet, isp1,  isp2", itet, isp1, isp2
        kk (0:3) = ib1bz( idtetf(0:3,itet) )     !  k
        !print *, "swj_itet, kk", itet, kk
!! may2013 for eizbmode 
        !print *, "swj1", eibzmode, kk, nwgt
        if(eibzmode) then
           if(sum(nwgt(kk(0:3)))==0) cycle 
        endif
        kkv(1:3, 0:3) = qbzw (1:3, idtetf(0:3,itet) )
c        if(mod(itet,1000)==1.and.verbose()>=150) write(6,*)' tetwt5: itet out of ntetf=',itet ,ntetf
        do 1100 im = 1,nmtet
          !print *, "swjcheck1"
          kkm (0:3)       = ib1bzm( idtetfm(0:3,im,itet) ) !  k   in micro-tet
          kvec(1:3, 0:3) = qbzwm ( 1:3, idtetfm(0:3,im,itet) )
          ek_ ( 1:nband+nctot, 0:3) = ekzz1( 1:nband+nctot, kkm(0:3)) ! k
          ekq_( 1:nband+nctot, 0:3) = ekzz2( 1:nband+nctot, kkm(0:3)) ! k+q
          noccx_k  = noccx1 ( ek_(1:nband, 0:3),4,nband, efermi)  
                    !the highest number of occupied states (HNOS)
          print *, "ek_(1:nband,0:3)"
          print *, ek_(1:nband,0:3)
          print *, "efermi", efermi
          print *, "noccx_k", noccx_k
          
          noccx_kq = noccx1 (ekq_(1:nband, 0:3),4,nband, efermi)
          print *, "noccx_kq", noccx_kq
          print *, "ekq_(1:nband,0:3)"
          print *, ekq_(1:nband,0:3)
          nbandmx_k   = nband
          nbandmx_kq  = nband
          !print *, "swjcheck2"
          do i=1,nband
            if( maxval(ek_(i, 0:3)) >1d9) then
              nbandmx_k = i-1
              exit
            endif
          enddo
          do i=1,nband
            if( maxval(ekq_(i, 0:3)) >1d9) then
              nbandmx_kq = i-1
              exit
            endif
          enddo
c  write(*,'("2  before n & np cycles",$)')   !   call cputid(0)
          do jpm = 1,npm !feb2006
            if(jpm==1) then
              ibxmx = noccx_k + nctot
              jbxmx = nbandmx_kq !nband  Apr2009takao
            else
              ibxmx = nbandmx_k  !nband  Apr2009takao
              jbxmx = noccx_kq + nctot
            endif
            do ibx  = 1, ibxmx !noccx_k + nctot  !   occupied
              do jbx  = 1, jbxmx !nband             ! unoccupied
                !print *, "swjcheck3"
c          if(chkwrt.and.job==1) write(6,*) '****** ib jb=',ib,jb
                if(ibx<=noccx_k .or. jpm==2  ) then
                  ib = ibx
                else
                  ib = ibx - noccx_k + nband
                endif
                if(jbx<=noccx_kq .or. jpm==1  ) then
                  jb = jbx
                else
                  jb = jbx - noccx_kq + nband
                endif
!! This mechanism treat ek_ and ekq_ as occpied or unoccupied.
                if(jpm==1) then
c           ibxmx = noccx_k + nctot
c           jbxmx = nband
                  eocc   => ek_ (ib,0:3)
                  eunocc => ekq_(jb,0:3)
                else
c           ibxmx = nband
c           jbxmx = noccx_kq + nctot
                  eunocc => ek_ (ib,0:3)
                  eocc   => ekq_(jb,0:3)
                endif

                if( minval(eocc  ) <= efermi .and. 
     &        maxval(eunocc) >= efermi ) then
                  if( maxval(eunocc(:)-eocc(:)) <0 ) cycle ! this makes a bit effective.
                  if(job==0 ) then    !takao
                    iwgt(ib,jb,kk(0:3),jpm)= .true.
                    x(0:3) = .5d0*(eocc-eunocc)     ! + omg !Denominator. unit in Hartree.
                    demax_ =  maxval(-x(0:3))  ! in Hartree
                    demin_ =  minval(-x(0:3))  ! in Hartree
                    do ixx=0,3
                      demax(ib,jb,kk(ixx),jpm) = max(demax_, demax(ib,jb,kk(ixx),jpm))
                      demin(ib,jb,kk(ixx),jpm) = min(demin_, demin(ib,jb,kk(ixx),jpm))
                    enddo
                  else
                    x(0:3) = .5d0*(eocc-eunocc)  ! + omg !Denominator. unit in Hartree.
                    wtthis2 = 0d0
                    if(chkwrt) then
                      write(6,"('### Goto lindtet6: itet ib jb Ef='
     &                    ,i8,2i5,d12.4,' ###')" ) itet,ib,jb,efermi
                      write(6,"('  eocc  - Ef= ',4f10.3)") eocc  -efermi
                      write(6,"('  eunocc- Ef= ',4f10.3)") eunocc-efermi
                      write(6,"('  -x(a.u.)= ',4f10.3)") -x(0:3)
                    endif
                    call lindtet6( 
     i          kkv,kvec,eocc,eunocc,
     i          x, efermi, 
     i          frhis, nwhis,  
     o          wtthis2) ! kvec(1:3, 0:3), ea(0:3), x(0:3) !  wgt(ib,jb,kk(0:3)) = wgt(ib,jb,kk(0:3)) + wttx/4d0/voltot
c     m          fermi_dS) ! may2012
                    if(chkwrt) then
                      write(6,"('  sum nwhis wtthis2=',i5,d13.5)") nwhis,sum(wtthis2)
                    endif

                    if(chkwrt) then
                      write(6,"('  === ihis [range(a.u.)] wtthis2(0:3) ===')")
                      do ihis = 1,nwhis
                        ddw = frhis(ihis+1)-frhis(ihis)
                        if(chkwrt.and.sum(abs(wtthis2(ihis,:)))/=0d0) then
c                   if(frhis(ihis)>0.15d0/2d0.and.frhis(ihis+1)<.20d0/2d0) then
c                  write(6,"(i4,' [',2f8.3,'] ',4d12.4)") ihis,frhis(ihis),frhis(ihis+1),wtthis2(ihis,:)/ddw
                          write(6,"('ttt',3i4,f10.5,4d12.4)") itet,ib,jb,(frhis(ihis)+frhis(ihis+1))/2d0, wtthis2(ihis,:)/ddw
c                   endif
                        endif
                        do ikx =0,3
                          ibib = ibjb(ib,jb,kk(ikx),jpm)
                          if(ihis< ihw(ibib,kk(ikx),jpm).or.
     &               ihis> ihw(ibib,kk(ikx),jpm)- 1 + nhw(ibib,kk(ikx),jpm) ) then
                            if(sum(abs(wtthis2(ihis,0:3)))/=0d0) 
Cstop2rx 2013.08.09 kino     &               stop 'tetwt5:wtthis2/=0 for out of range'
     &               call rx( 'tetwt5:wtthis2/=0 for out of range')
                          endif
                        enddo
                      enddo
                    endif
!
                    do ikx = 0,3
                      !print *, "swj_itet, im, jpm, ibx, jbx, kk", itet, im, jpm, ibx, jbx, kk
                      !print *, "kk(0), kk(1), kk(2), kk(3)", kk(0), kk(1), kk(2), kk(3)
                      !print *, "swjcheck4, ikx,jpm,npm", ikx,jpm,npm
                      !print *, "swjcheck ib, jb kk(ikx), ikx, jpm"
                      !print *, "ib, jb, kk(ikx), ikx, jpm"
                      !print *, ib, jb, kk(ikx), ikx, jpm
                      !print *, "swjcheck_ibjb,ikx-1", ibjb(ib, jb, kk(ikx-1), jpm)
                      !print *, "swjcheck_ibjb,ikx", ibjb(6,5,27,1)
                      !print *, "ibjb(6,5,27,1)", ibjb(6,5,27,1)
                      !print *, "swjcheck_ibjb,ikx", ibjb(ib, jb, kk(ikx), jpm)
                      !print *, "ibjb(ib,jb,kk(0),jpm)", ibjb(ib,jb,kk(0),jpm)
                      !print *, "ibjb(ib,jb,kk(1),jpm)", ibjb(ib,jb,kk(1),jpm)
                      !print *, "ibjb(ib,jb,kk(2),jpm)", ibjb(ib,jb,kk(2),jpm)
                      !print *, "ibjb(ib,jb,kk(3),jpm)", ibjb(ib,jb,kk(3),jpm)
                      ibib = ibjb(ib,jb,kk(ikx),jpm)
                      !print *, "swjcheck5-ibib", ibib
                      jini = jhw(ibib,kk(ikx),jpm)
                      !print *, "swjcheck5-jini", jini
                      iini = ihw(ibib,kk(ikx),jpm)
                      !print *, "swjcheck5-iini", iini
                      nnn =  nhw(ibib,kk(ikx),jpm)
                      !print *, "swjcheck5-nnn", nnn

                      if(matrix_linear()) then
                        whw(jini:jini+nnn-1) = 
     &          whw(jini:jini+nnn-1) + wtthis2(iini:iini+nnn-1,ikx) * piofvoltot
                      else
                        whw(jini:jini+nnn-1) = 
     &          whw(jini:jini+nnn-1) + wtthis2(iini:iini+nnn-1,0) * piofvoltot*4*wtet(ikx,im,itet) ! piofvoltot= pi/voltot/4
                      endif
                    enddo

                  endif
                endif

              enddo
            enddo
          enddo
 1100   continue
        !print *, "swj_dddddd"
 1000 continue
      !print *, "swj_xxxxxx"
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
          call chkdgn(ecore, nctot, 
     &    nrankc1, ini1(nrank1+1), ied1(nrank1+1), nband,ipr)
        endif
        nrankc2 = 0
        if(nctot/=0) then
          call chkdgn(ecore, nctot, 
     &    nrankc2, ini2(nrank2+1), ied2(nrank2+1), nband,ipr)
        endif
c       if(ipr) print *,' kx nrank =',kx,nrank1,nrank2,nrankc
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
                if( count(iwgt(ixi1:ixe1,ixi2:ixe2,kx,jpm))>0 ) wxx=.true.
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
                if(isum/=0 .and. isum /= (ixe1-ixi1+1)*(ixe2-ixi2+1) ) then
Cstop2rx 2013.08.09 kino                  stop 'tetwt5: isum srrange---1'
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
                    wtthis(ini:ied) = 
     &        wtthis(ini:ied) + whw(ioff:ioff+nhw(ibib,kx,jpm)-1)
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
 1120   continue
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
c        if(abs(wgt(i,j,ik) )>0d0 .and. j>nband) stop 'tetwt5: bug ?' ! bug checker
c        if(abs(wgt(i,j,ik) )>eps ) then
                if(iwgt(i,j,ik,jpm)) then
                  nbnb(ik,jpm) = nbnb(ik,jpm)+1
cc          write(56,"(3i4,2d13.6)") i,j,ik, wgt(i,j,ik)
c          if(abs(wgt(i,j,ik) )>0.1d0 )
c     &      write(6,"(' k i ik=',3i4,'  wgt= ',2d13.6)")
c     &      i,j,ik, wgt(i,j,ik)
                endif
              enddo
            enddo !;  write(6,*)' ik=',ik,' num of nonzero wgt=',nbnb(ik)
          enddo
        enddo
c       write(6,*)' max num of nonzero wgt(k)=',maxval(nbnb)
c       write(6,*)' tot num of nonzero wgt   =',sum(nbnb)
c       write(6,*)' sum of wgt   =',sum(wgt)
c       write(6,*)' tetwt4: end '; call cputid  (0)
        write(6,"('  tetwt5_dtet3: maxval(nbnb(1:nqbz,1:npm)) sum(nbnb)=',2i8
     &         ,' count(iwgt) =',i10)") maxval(nbnb),sum(nbnb),count(iwgt) !,sum(wgt)
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
      end
c------------------------------------------------------------------------





c----------------------------------------------------------------
      subroutine hisrange(frhis,nwhis, demin,demax,
     o  ihw,nhw)
c- determine the pointer ihw and number of pointer.
      implicit none
      integer(4):: nwhis,ihw,nhw,ihis
c      real(8):: frhis(nwhis+1), demin,demax
      real(8):: frhis(nwhis+1)
      real(4):: demin,demax
c  Range for each division is from frhis(ihis-1) to frhis(ihis).
c      print *,' hisrange nwhis=',nwhis,demin,demax
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
      if(ihw==-9999.or.nhw==-9999) then
        print *,' error info=',nwhis, demin,demax,ihw,nhw
        do ihis = 1, nwhis
          write(6,*) ihis, frhis(ihis+1)
        enddo
Cstop2rx 2013.08.09 kino        stop ' hisrange: wrong'
        call rx( ' hisrange: wrong')
      endif
      end


c------------------------------------------------------------------------
      subroutine lindtet6(
     i  kkv,kvec, ea, eb, x, efermi,      ! kvec(1:3, 0:3), ea(0:3), x(0:3) 
     i  frhis, nwhis,  
     o  wtthis ) 
c     m  fermi_ds)
C- Calculate the imaginary part of \int dk1 dk2 dk3 f(ea)(1-f(eb))/(\omega+x ), that is,
C  the microcell integral of wtthis(ihis) = \int_lower(ihis)^upper(ihis) d\omega \int d^3k f(ea(k)) (1-f(eb(q+k))) \delta (omg + x(k) )
C f(E) denote the Fermi distribution funciton. Only for T=0.
C Tetrahedon is specified by values at 4 corners; kvec(1:3, 1:4), ea(1:4), eb(1:4), x(1:4)
Cr This code is based on J.Rath&A.J.Freeman PRB11(6) p.2109(1975).
Cr The \sum_ihis wthis(ihis) = the total volume of the microcell = \int d^3 k
ctakao/Feb/2002 -----------------------------------------------------------------
!! ds is the accumulation variable for the fermi surface.
      implicit none
      integer(4) :: ieaord(1:4),i,isig,n,itmp,ix
      real(8)    ::  kvec(1:3, 1:4), x(1:4), ea(1:4),
     &  kk(3,1:8),xx(1:8),ee(1:4), am(3,3), eb(1:4),ebf(1:8), !ebfKx(1:4),Kx(1:3,4),xKx(4),
     &  vcell,efermi,etest ,  kkv(3,4)
      integer(4):: nwhis
      real(8):: frhis(nwhis+1),wtthis(nwhis,4)
      logical:: chkwrt=.false.
c      real(8):: fermi_ds(3,3)
c---------------------------------------------------------------------------------
      if(chkwrt) then
        print *, ' lindtet6: *****************************'
        write(6,"(' i=',i3,' x ea eb =',3f10.5)") 
     &  (i, x(i), ea(i)-efermi, eb(i)-efermi,i=4,1,-1)
      endif

#ifdef EXPAND_SORTEA
      n=4
c      isig = 1
*poption noparallel
      do i = 1,n
        ieaord(i) = i
      enddo
*poption noparallel
      do ix= 2,n
*poption noparallel
        do i=ix,2,-1
          if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
c            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      ieaord(1:4) = ieaord(4:1:-1)
#else
      call sortea( ea,ieaord, 4 ,isig); ieaord(1:4) = ieaord(4:1:-1)
#endif
      ! the order E_4,E_3,E_2,E_1 This suits for Rath&Freeman.
      kk(1:3,1:4) = kvec(1:3,ieaord(1:4))
      ! 4 corners  denoted by kvec(:,1:4), ee(1:4), and xx(1:4)
      ee (1:4)    = ea (ieaord(1:4)) - efermi
      xx (1:4)    = x  (ieaord(1:4))
      ebf(1:4)    = eb (ieaord(1:4)) - efermi
ccccccccccccc
c      write(6,"(' i iea=',2i4)") (i,ieaord(i),i=1,4)
c      write(6,"(' i=',i3,' xx ee ebf =',3f10.5,'  kk=',3f10.5)")
c     &     (i,xx(i),ee(i),ebf(i),kk(1:3,i),i=4,1,-1)
ccccccccccccc
      if( 0d0<=ee(4) ) then
        if(chkwrt) write(6,*) 'lindtet6: fig 000'
!        wttx = (0d0,0d0)
      elseif( ee(4) < 0d0 .and. 0d0<= ee(3) ) then   !!! Fig 1.
        if(chkwrt) write(6,*) 'lindtet6: fig 1 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1 !K1 is on the like k4---k2.
        call midk3(kk,ee,xx,ebf, 4,1,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
        call midk3(kk,ee,xx,ebf, 4,3,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3

!! Corners of occupied states at the Fermi energy are kk(1:3,i+4) i=1,2,3.
!! This makes a Fermi surface
c        if(sum(abs(ea-eb))<1d-6) then
c          call fermids(kk(:,2+4),kk(:,1+4),kk(:,3+4),fermi_dS)
c        endif

!! still developing
        call  inttetra6(kkv,kk,xx,ebf,(/4,1+4,2+4,3+4/),  ! k4,K1,K2,K3
     i         frhis,  nwhis,  
     o         wtthis )
      elseif( ee(3) < 0d0 .and. 0d0<= ee(2) ) then   !!! Fig 2.
        if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1
        call midk3(kk,ee,xx,ebf, 4,1,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
        call midk3(kk,ee,xx,ebf, 3,1,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3
        call midk3(kk,ee,xx,ebf, 3,2,  kk(1,4+4),xx(4+4),ebf(4+4)) !K4

c        if(sum(abs(ea-eb))<1d-6) then
c          call fermids(kk(:,1+4),kk(:,2+4), kk(:,3+4),fermi_dS) !K1-K2-K3
c          call fermids(kk(:,1+4),kk(:,2+4), kk(:,4+4),fermi_dS) !K1-K2-K4
c        endif

        if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx2'
        call  inttetra6(kkv,kk,xx,ebf,(/4, 3,1+4,2+4/),  ! k4,k3,K1,K2
     i         frhis,  nwhis,  
     o         wtthis )
        if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx3'
        call  inttetra6(kkv,kk,xx,ebf,(/3,2+4,3+4,1+4/),  ! k3,K2,K3,K1
     i         frhis,  nwhis,  
     o         wtthis )
        if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx4'
        call  inttetra6(kkv,kk,xx,ebf,(/3,1+4,3+4,4+4/),  ! k3,K1,K3,K4
     i         frhis,  nwhis,  
     o         wtthis )
        if(chkwrt) write(6,*) 'lindtet6: fig 2 xxx5'
c-----------
      elseif( ee(2) < 0d0 .and. 0d0<= ee(1) ) then   !!! Fig 3.
        if(chkwrt) write(6,*) 'lindtet6: fig 3 xxx'
        call midk3(kk,ee,xx,ebf, 1,4,  kk(1,1+4),xx(1+4),ebf(1+4)) !K1
        call midk3(kk,ee,xx,ebf, 1,2,  kk(1,2+4),xx(2+4),ebf(2+4)) !K2
        call midk3(kk,ee,xx,ebf, 1,3,  kk(1,3+4),xx(3+4),ebf(3+4)) !K3

c        if(sum(abs(ea-eb))<1d-6) then
c          call fermids(kk(:,1+4),kk(:,2+4), kk(:,3+4),fermi_dS) !K1-K2-K3
c        endif

        call  inttetra6(kkv,kk,xx,ebf,(/3, 4,3+4, 2/),  ! k3,k4,K3,k2
     i         frhis,  nwhis,  
     o         wtthis )
        call  inttetra6(kkv,kk,xx,ebf,(/4,1+4,2+4,3+4/),  ! k4,K1,K2,K3
     i         frhis,  nwhis,  
     o         wtthis )
        call  inttetra6(kkv,kk,xx,ebf,(/4, 2,2+4,3+4/),  ! k4,k2,K2,K3
     i         frhis,  nwhis,  
     o         wtthis )
      else
        if(chkwrt) write(6,*) 'lindtet6: fig 4 xxx'
        call  inttetra6(kkv,kk,xx,ebf, (/1,2,3,4/),  ! k1,k2,k3,k4
     i         frhis, nwhis,  
     o         wtthis )
      endif
      if(chkwrt) write(6,*) 'end of lindtet6',wtthis
c      stop 'yyyyyyyyyyyyyyyyyyyyyy'
      end


c-----------------------------------------------------
      subroutine inttetra6(kkv,kk_,xx_,ebf,itetx, 
     i         frhis, nwhis,  
     o         wtthis )
C calculate tetrahedron integral Eq.(16).
C the four corners and denoted by itetx.
Ci kk (k1-k4) and xx (value of denominator)
Ci Kx (K1-K4) and xkx
Cr The four corners are selected from 8 points. It is specified by itetx.
Cr wtthis is accumulated
      implicit none
      integer(4):: itetx(4),ix,i,ieaord(4),isig,n,itmp
      real(8) ::  kk_(3,1:8),xx_(1:8),       !Kx_(3,1:4),xKx_(1:4),ebfKx(1:4), 
     &    kk(3,1:8),xx(1:8),ebf(1:4),ebfin(1:8) ,ee(4)  ,kin(3,4), xin(4)
      integer(4):: nwhis
      real(8):: frhis(nwhis+1),wtthis(nwhis,4)
      logical :: chkwrt=.false.

      real(8)   ::  kkv(3,4)
c ---
      if(chkwrt) then
        print *,' inttetra6: ---------------------------'
c        print *,' ebf  =', ebf(1:4)
c        print *,' ebfKx=', ebfKx(1:4)
      endif

c--- kk xin ebfin ---
c      do i = 1,4
c        if( itetx(i) < 10 ) then   ! kk (k1-k4 in Paper)
c          ix = itetx(i)
c          xin(    i) = xx_(    ix)
c          ebfin(i)   = ebf(ix)
c          kin(1:3,i) = kk_(1:3,ix)
c        else                       ! Kx (K1-K4 in Paper)
c          ix = itetx(i)/10
c          xin(    i) = xKx_(    ix)
c          ebfin(i)   = ebfKx(ix)
c          kin(1:3,i) = Kx_ (1:3,ix)
c        endif
c      enddo

c--- kk xin ebfin ---
      xin(    1:4) = xx_(    itetx(1:4))
      kin(1:3,1:4) = kk_(1:3,itetx(1:4))
      ebfin(  1:4) = ebf(    itetx(1:4))
c
      if(chkwrt) then
        write(6,"(' i=',i3,' xin ebfin =',2f10.5,'  kin=',3f10.5)") 
     &  (i,xin(i),ebfin(i),kin(1:3,i),i=4,1,-1)
      endif

c... ee decendent order ee(4)>ee(3)>ee(2)>ee(1)
#ifdef EXPAND_SORTEA
      n = 4
      isig = 1
*poption noparallel
      do i = 1,n
        ieaord(i) = i
      enddo
*poption noparallel
      do ix= 2,n
*poption noparallel
        do i=ix,2,-1
          if( ebfin(ieaord(i-1)) >ebfin(ieaord(i) ) ) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo

#else
      call sortea( ebfin,ieaord, 4 ,isig)  ! the order E_4,E_3,E_2,E_1 This suits for Rath&Freeman.
#endif
      kk(1:3,1:4) = kin  (1:3,ieaord(1:4))   ! 4 corners  denoted by kkv(:,1:4), ee(1:4), and xx(1:4)
      ee(    1:4) = ebfin(    ieaord(1:4))
      xx(    1:4) = xin  (    ieaord(1:4))
c
      if(chkwrt)write(6,"(' i=',i3,' xx ee =',2f10.5,'  kk=',3f10.5)") 
     &  (i,xx(i),ee(i),kk(1:3,i),i=4,1,-1)
c
      if( 0d0>=ee(4) ) then
!        wttx = (0d0,0d0)
      elseif( ee(4) > 0d0 .and. 0d0>= ee(3) ) then   !!! Fig 1.
        if(chkwrt) write(6,*) 'inttetra5: fig 1'
        call midk(kk,ee,xx, 4,2,  kk(1,1+4),xx(1+4)) !K1 -> Kx(:,1), x(K1) -> xkx(1). K1 is on the like k4---k2.
        call midk(kk,ee,xx, 4,1,  kk(1,2+4),xx(2+4)) !K2
        call midk(kk,ee,xx, 4,3,  kk(1,3+4),xx(3+4)) !K3
        call inttetrac6(kkv,kk,xx, (/4,1+4,2+4,3+4/),  ! k4,K1,K2,K3
     i         frhis, nwhis,  
     o         wtthis )
      elseif( ee(3) > 0d0 .and. 0d0>= ee(2) ) then   !!! Fig 2.
        if(chkwrt) write(6,*) 'inttetra5: fig 2'
        call midk(kk,ee,xx, 4,2,  kk(1,1+4),xx(1+4)) !K1
        call midk(kk,ee,xx, 4,1,  kk(1,2+4),xx(2+4)) !K2
        call midk(kk,ee,xx, 3,1,  kk(1,3+4),xx(3+4)) !K3
        call midk(kk,ee,xx, 3,2,  kk(1,4+4),xx(4+4)) !K4
        call inttetrac6(kkv,kk,xx, (/4, 3,1+4,2+4/),  ! k4,k3,K1,K2
     i         frhis, nwhis,  
     o         wtthis )
        call inttetrac6(kkv,kk,xx, (/3,2+4,3+4,1+4/),  ! k3,K2,K3,K1
     i         frhis, nwhis,  
     o         wtthis )
        call inttetrac6(kkv,kk,xx, (/3,1+4,3+4,4+4/),  ! k3,K1,K3,K4
     i         frhis, nwhis,  
     o         wtthis )
      elseif( ee(2) > 0d0 .and. 0d0>= ee(1) ) then   !!! Fig 3.
        if(chkwrt) write(6,*) 'inttetra5: fig 3'
        call midk(kk,ee,xx, 1,4,  kk(1,1+4),xx(1+4)) !K1
        call midk(kk,ee,xx, 1,2,  kk(1,2+4),xx(2+4)) !K2
        call midk(kk,ee,xx, 1,3,  kk(1,3+4),xx(3+4)) !K3
        call inttetrac6(kkv,kk,xx, (/3, 4,3+4, 2/),  ! k3,k4,K3,k2
     i         frhis, nwhis,  
     o         wtthis )
        call inttetrac6(kkv,kk,xx, (/4,1+4,2+4,3+4/),  ! k4,K1,K2,K3
     i         frhis, nwhis,  
     o         wtthis )
        call inttetrac6(kkv,kk,xx, (/4, 2,2+4,3+4/),  ! k4,k2,K2,K3
     i         frhis, nwhis,  
     o         wtthis )
      else
        if(chkwrt) write(6,*) 'inttetra5: fig 4'
        call inttetrac6(kkv,kk,xx, (/1,2,3,4/),  ! k1,k2,k3,k4
     i         frhis, nwhis,
     o         wtthis )
      endif
      end

c-----------------------------------
      subroutine inttetrac6(kkv,kk,xx, itetx,
     i         frhis, nwhis,
     o         wtthis )
C- Calculate tetrahedron integral Eq.(16).
C The four corners and denoted by itetx(1:4).
Ci kk (k1-k4) and xx (value of denominator)
Ci Kx (K1-K4) and xkx
Cr The four corners are selected from 8 points. It is specified by itetx(1:4).
Cr wtthis is accumulating.
      implicit none
      integer(4):: itetx(4),ix,i
      real(8) ::  kk(3,1:8),xx(1:8), am(3,3)
     &   ,kin(3,4), xin(4), det33, intttvv,  intttv,intvv,work !Kx(3,1:4),xKx(1:4),
      integer(4):: nwhis
      real(8):: frhis(nwhis+1),wtthis(nwhis,4),voltet,kkv4bm(3),bm(3,3)
      logical:: chkwrt=.false.,matrix_linear
      real(8) ::  kkv(3,4), kkvkin(4,4)
c
      if(chkwrt) print *, ' inttetrac6: === '
      xin(    1:4) = xx(    itetx(1:4))
      kin(1:3,1:4) = kk(1:3,itetx(1:4))
      do i = 1,3
        am(1:3,i) = kin(1:3,i) - kin(1:3,4)
      enddo
      voltet = abs(det33(am)/6d0) ! \omega (volume of tetrahedra) = abs(det33(am)/6d0) See Eq. (17).

c--- kkvkin: kin is decomplosed into kkv.!!!!!!!!!!
!   e.g. kin (:,1) =  \sum_i kkv(:,i) * kkvkin(i,1)
      if(matrix_linear()) then
        do i = 1,3
          am(1:3,i) = kkv(1:3,i) - kkv(1:3,4)
        enddo
        call dinv33x(am, bm)
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
!
      call intttvc6(kkvkin, xin,voltet,
     i         frhis, nwhis, 
     o         wtthis )
      end

c---------------------------------------------------------------------------------------------
      subroutine intttvc6(kkvkin, v,voltet,
     i         frhis, nwhis, 
     o         wtthis ) !this is accumlating variable
C- Histgram weights for each bin.
Cr The i-th bin of the Histgrams is [frhis(i) frhis(i+1)].
Cr  wtthis(ihis) += \int_{frhis(ihis)^{frhis(ihis+1)} d \omega  \int d^3k \delta(\omega +v(k) )
Cr  Total sum of the histgram weight is equal to voltet.
Co  wtthis(ihis): Note this is accumulating variable.
C
Cr Note wtthis(nwhis,0:3)
Cr if matrix_linear()=T, wtthis is calculated assuming
Cr                       the linear-dependency of matrix elements. takao /dec/2003
      implicit none
      integer(4)::  inihis, iedhis, nwhis,ichk=2313
      integer(4) :: ieaord(1:4),i,isig,idif(3),idf,ix,n,itmp,ihis
      real(8)   ::  v(4), voltet, WW(4),  norm,s1,s2,pvn,
     &       integb3p, integb3m , integb2p, integb2m,ww2p,ww3p
      real(8):: frhis(nwhis+1),wtthis(nwhis,4),intega,integb,stot,xxx,wx
      logical ::chkwrt=.false., matrix_linear
      real(8):: kkvkin(4,4),www(1:4)=.25d0,ec,wcg(4),wttt
#ifdef EXPAND_SORTEA
      n=4
c      isig = 1
*poption noparallel
      do i = 1,n
        ieaord(i) = i
      enddo
*poption noparallel
      do ix= 2,n
*poption noparallel
        do i=ix,2,-1
          if( -v(ieaord(i-1)) > -v(ieaord(i) ) ) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
c            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
#else
      call sortea( -v,ieaord,4,isig)
#endif
      WW(1:4) = -v( ieaord(1:4) )

!  ww(1)<ww(2)<ww(3)<ww(4)
      if((.not.(WW(1)<=WW(2))).or.(.not.(WW(2)<=WW(3))).or.(.not.(WW(3)<=WW(4))) ) then
        write(6,"(/,' --- intttvc6: wrong order WW=',4d14.6)") WW
Cstop2rx 2013.08.09 kino        stop 'intttvc6: wrong order of WW'
        call rx( 'intttvc6: wrong order of WW')
      endif

      if(chkwrt) then
        write(ichk,"(/,' --- intttvc6: e=',4d23.16)") WW
      endif

      inihis= -999
      iedhis= -999
      ix=1
      do ihis = 1,nwhis
        if(ix==1 .and. WW(ix)<frhis(ihis+1)) then
          inihis = ihis
          ix=4
        endif
        if(ix==4 .and. WW(ix)<frhis(ihis+1)) then
          iedhis = ihis
          exit
        endif
      enddo
      if(iedhis==-999.or.inihis==-999) then
        print *,' intttvc6: can not find inihis iedhis'
Cstop2rx 2013.08.09 kino        stop    ' intttvc6: can not find inihis iedhis'
        call rx( ' intttvc6: can not find inihis iedhis')
      endif
c      if(chkwrt) print *,' inihis iedhis=', inihis,iedhis
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' ###########integtetn test #########'
c      ww(1:4)=(/1d0,1.5d0,9.5d0,10d0/)
c      do ix=1,101
c        wx = (ix-1)/10d0
c        call integtetn(WW, Wx, xxx)
c        write(61,"(i3,2d23.16)")ix,wx,xxx
c      enddo
c      stop 'test end integtetn---'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c          write(6,"('sum(www)= ',d13.4,' www(1:4)=',4f10.4)") sum(www),www(1:4)
          wttt = pvn*(integb - intega)*4d0
          wtthis(ihis,:) = wtthis(ihis,:) + wttt * www(:)
        else
          wtthis(ihis,1) = wtthis(ihis,1) + pvn*(integb - intega)
        endif

        if(chkwrt) then
          write(ichk,"(' ihis [Init End] wtt=', i5,3f11.6)") 
     &    ihis, frhis(ihis), frhis(ihis+1), pvn*(integb-intega)
        endif

        intega = integb
      enddo
c      if(chkwrt) then
c        dini = 0d0
c        do ihis=inihis,iedhis
c          if(ihis>1) dini=wtthis(ihis-1)
c          write(6,"(' ihis [Init  End]', i4,2f13.6,' wtt=',f13.6)")
c     &    ihis,dini,frhis(ihis),wtthis(ihis)
c          dini = frhis(ihis)
c        enddo
c      endif
      if(chkwrt) write(ichk,*) ' end of intttvc6'
      end


c--------------------------------------------------------------
      subroutine integtetn(e, ee, integb)
!> Calculate primitive integral of integb = 1/pi Imag[\int^ee dE' 1/(E' -e(k))] = \int^ee dE' S[E']
!! \remark
!!  S[E] : is area of the cross-section between the omega-constant plane and the tetrahedron.
!!  [here we assumee e1<e2<e3<e4].
!!  Normalization is not considered!
!!  Rath&Freeman Integration of imaginary part on Eq.(17)
      implicit none
      integer:: dummy4doxygen
      real(8) ::  e(1:4), ee, norm, integb,a,b,
     & V1,V2,V3,V4,D1,D2,D3,D4,e1,e2,e3,e4
      integer(4)::i1,i2,i3
c-------------------------------------------------------------
c      if((.not.(e(1)<=e(2))).or.(.not.(e(2)<=e(3)))
c     & .or.(.not.(e(3)<=e(4))) ) then
c        write(6,"(/,' --- integtetn wrong order of e=',4d14.6)") e
c        stop 'integtetn: wrong order of e'
c      endif

c      write(6,"(' --- integtetn e =',4d23.16)") e
c      write(6,"(' --- integtetn ee=',d23.16)") ee
      if    ( ee<e(1) ) then  !jan2008 ee<=e(1) ---> ee<e(1)
        integb = 0d0
        return
      elseif (ee>e(4) ) then
Cstop2rx 2013.08.09 kino        stop ' integtetn: ee>e(4)'
        call rx( ' integtetn: ee>e(4)')
      endif

c--- case1. poor numerical accuracy for cases.
      if(.false.) then
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

c--- case2
      e1 = e(1)-3d-8
      e2 = e(2)-2d-8
      e3 = e(3)-1d-8
      e4 = e(4)
      V1= ee - e1
      V2= ee - e2
      V3= ee - e3
      V4= ee - e4
c      D1 = (V1-V4)*(V1-V3)*(V1-V2) !>0
c      D2 = (V2-V4)*(V2-V3)*(V2-V1) !<0
c      D3 = (V3-V4)*(V3-V2)*(V3-V1) !>0
c      D4 = (V4-V1)*(V4-V2)*(V4-V3) !<0
      if    ( e1<=ee .and. ee<e2 ) then
        D1 = (e4-e1)*(e3-e1)*(e2-e1) !>0
        integb =  V1**3/D1
      elseif( e2<=ee .and. ee<e3 ) then
        D1 = (e4-e1)*(e3-e1)*(e2-e1) !>0
        D2 = (e4-e2)*(e3-e2)*(e1-e2) !<0
        a  =  V1/   D1**(1d0/3d0)
        b  =  V2/(-D2)**(1d0/3d0)
        integb = (a-b) * (a**2+a*b+b**2)
      elseif( e3<=ee .and. ee<e4 ) then
        D4 = (e1-e4)*(e2-e4)*(e3-e4) !<0
        integb = 1d0 - V4**3/D4
      elseif( ee==e4 ) then
        integb = 1d0
      endif

c-----------------------------------------------------------
c      D2 = (V2-V4)*(V2-V3)*(V2-V1) !>0
c      D3 = (V3-V4)*(V3-V2)*(V3-V1) !<0
c      integb = 0d0
c      if( e(1)<=ee ) integb =          V1**3* D2*D3
cc      write(6,*) ' integb1=',integb
c      D1 = (V1-V4)*(V1-V3)*(V1-V2) !<0
c      if( e(2)<=ee ) integb = integb + V2**3* D3*D1
cc      write(6,*) ' integb2=',integb
c      if( e(3)<=ee ) integb = integb + V3**3* D1*D2
cc      write(6,*) ' integb3=',integb
c--
c      if(V4==0d0) then
c       write(6,*)
c       write(6,"(' V  =',4d18.10)") V1, V2, V3, V4
c       write(6,"(' int=',3d13.5,16x,d13.5)") V1**3/D1,V2**3/D2,V3**3/D3,
c     &             max(abs(V1**3/D1),abs(V2**3/D2),abs(V3**3/D3))
c         write(6,*) ' integb=',integb
c      endif
      end


c-----------------------------------------------------
      subroutine mkwcg(e, ee, wcg)
C- calculate wweight for each corners.---------
Ci e(1:4),ee
Co wcg
Cr normarization is sum(scg(1:4))=1d0
c----------------------------------------------
      implicit none
      real(8) ::  e(1:4), ee, wcg(1:4),e1,e2,e3,e4,
     &  w14_1,w14_4,w12_1,w12_2,w13_1,w13_3,w23_2,w23_3,w24_2,w24_4,w34_3,w34_4
c-------------------------------
cccccccccccccccccccccccccccccccccccccccc
! This is for original case thru Dec 2003.
c      wcg=.25d0
c      return
cccccccccccccccccccccccccccccccccccccccc
      e1 = e(1)-3d-8
      e2 = e(2)-2d-8
      e3 = e(3)-1d-8
      e4 = e(4)
      if    ( ee<=e(1) ) then
        wcg(1)  =1d0
        wcg(2:4)=0d0
      elseif( e1<=ee .and. ee<e2 ) then
        call wab(e1,e2,ee,w12_1,w12_2)
        call wab(e1,e3,ee,w13_1,w13_3)
        call wab(e1,e4,ee,w14_1,w14_4)
        wcg(1) = (w12_1 +w13_1 + w14_1)/3d0
        wcg(2) =  w12_2/3d0
        wcg(3) =  w13_3/3d0
        wcg(4) =  w14_4/3d0
      elseif( e2<=ee .and. ee<e3 ) then !Is this correct?
        call wab(e1,e3,ee,w13_1,w13_3)
        call wab(e1,e4,ee,w14_1,w14_4)
        call wab(e2,e3,ee,w23_2,w23_3)
        call wab(e2,e4,ee,w24_2,w24_4)
        wcg(1) = (w13_1 + w14_1)/4d0
        wcg(2) = (w23_2 + w24_2)/4d0
        wcg(3) = (w13_3 + w23_3)/4d0
        wcg(4) = (w14_4 + w24_4)/4d0
      elseif( e3<=ee .and. ee<e4 ) then
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
ccccccccccccccccccccccccc
c      write(6,"('mkwcg   e=',4d13.4, ' ee=',d13.4)")e,ee
c      write(6,"('      wcg=',4d13.4)")wcg
cccccccccccccccccccccccc
      end
c-------------------------------------
      subroutine wab(ea,eb,ee,wa,wb)
      implicit none
      real(8)::ea,eb,wa,wb,eet,ee
      eet= eb - ea
      wa= (eb-ee)/eet
      wb= (ee-ea)/eet
      end


c---
      subroutine addsciss(delta, ef, nnn, eig)
      real(8):: eig(nnn),ef,delta
      write(6,*)' asssciss delta=', delta
      do i=1,nnn
        if(eig(i)>ef) eig(i)= eig(i)+delta
      enddo
      end

c-----------------------------------
      subroutine midk3(kk,ee,xx,yy,i,j,   kout,xout,yout)
C- Calculate x and k(3) at the Fermi energy on the like k(i)---k(j).
      implicit none
      integer(4):: i,j
      real(8) ::  kk(3,1:4),xx(1:4),yy(1:4),ee(1:4), kout(3)
     &      ,xout,yout,ratio
      ratio     = ee(i)/(ee(i)-ee(j))
      xout      = xx(i)     + ratio * (xx(j)-xx(i))
      yout      = yy(i)     + ratio * (yy(j)-yy(i))
      kout(1:3) = kk(1:3,i) + ratio * (kk(1:3,j)-kk(1:3,i))
      end

      real(8) function det33(am)
      implicit none
      real(8),intent(in) :: am(3,3)
      det33= am(1,1)*am(2,2)*am(3,3)
     &       -am(1,1)*am(3,2)*am(2,3)
     &       -am(2,1)*am(1,2)*am(3,3)
     &       +am(2,1)*am(3,2)*am(1,3)
     &       +am(3,1)*am(1,2)*am(2,3)
     &       -am(3,1)*am(2,2)*am(1,3)
      end

      subroutine chkdgn(ene,ndat,  nrank,ixini,ixend,iof,ipr)
      implicit none
      integer(4) :: ndat,i,ix, ixini(ndat),ixend(ndat),nrank,iof
      real(8)    :: ene(ndat), epsx=1d-4
      logical ipr
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
*poption noparallel
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
c
      nrank = i
*poption noparallel
      do i =1,nrank
        ixini(i) = ixini(i)+iof
        ixend(i) = ixend(i)+iof
      enddo
c-check write
      if(ipr) then
        write(6,*)' nrank=',nrank
        do i = 1, ndat
          write(6,"(' i ',i3,' ene=',d15.7)") i,ene(i)
        enddo
        print *
        do i = 1, nrank
          write(6,"(' i ',2i3,' e=',d15.7)")
     &  ixini(i),ixend(i),ene(ixini(i)-iof)
        enddo
      endif
      end

      subroutine midk(kk,ee,xx,i,j,   kout,xout)
C- Calculate x and k(3) at the Fermi energy on the like k(i)---k(j).
      real(8) ::  kk(3,1:4),xx(1:4),ee(1:4), kout(3),xout
      ratio     = ee(i)/(ee(i)-ee(j))
      xout      = xx(i)     + ratio * (xx(j)-xx(i))
      kout(1:3) = kk(1:3,i) + ratio * (kk(1:3,j)-kk(1:3,i))
      end

      subroutine rsvwwk00_4(jpm,iwgt, nqbz,nband,nctot,ncc,nbnbx, !feb2006 ncc added
     o   n1b,n2b,noccxv,nbnb)
C- get (n1b n2b) corresponding to non-zero wgt.
      implicit none
      integer(4) :: jpm,nband, nctot,  ncc, nqbz,nbnbx ,ib,jb, kx,ix
      integer(4) :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),noccxv,nbnb(nqbz)
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
              if(jpm==1.and. ib<=nband .and. ib>noccxv ) noccxv =ib
              if(jpm==2.and. jb<=nband .and. jb>noccxv ) noccxv =jb
            endif
          enddo
        enddo
        nbnb(kx) = ix
      enddo
      write(6,*) ' rsvwwk: kx nbnbmax=',kx, maxval(nbnb)
      end

      SUBROUTINE TETFBZF(icase,qbas, N1,N2,N3, rk,nqbz,  IDTET,qbzw,ib1bz)
C-  Finds tetrahedra in all 1st BZ. takao mod. from tertirr
C ----------------------------------------------------------------------
Ci Inputs:
Ci  qb,n1,n2,n3,ipq, output from BZMESH;
Ci  nq, no. of irreducible k-points;
Co Outputs:
Co  ntet, No. of different tetrahedra
Co  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
Co  idtet(0,i), no. of tetrahedra of the i'th kind
Cm Memory:
Cm  No large internal storage; heap not accessed.
cr    This require subroutine CCUTUP (lmto-3).
C ----------------------------------------------------------------------
      IMPLICIT none
      integer(4):: n1,n2,n3, indexkw(0:n1,0:n2,0:n3),nqbz,kount,
     &  indexk(0:n1-1,0:n2-1,0:n3-1),
     &  i,i1,i2,i3,j1,j2,j3, IPQ(N1,N2,N3),
     .  IBTR(3,3),KCUT(3,4,6),IMC(0:1,0:1,0:1),
     .  idtet(4, 6*n1*n2*n3),iq(4),
     &  ntet,k1,k2,k3,itet,ic 
     & ,ib1bz((n1+1)*(n2+1)*(n3+1))  ,icase
      real(8) :: QB(3,3),QB1(3,3), qbas(3,3), rk(3,nqbz) ,qbzx(3)
     & , qbzw(3,(n1+1)*(n2+1)*(n3+1)),hf
      real(8),parameter :: epss = 1d-12
      logical,save ::chk=.true.
c------------------------------------------------------------------
      hf=0d0
      if(icase==2) hf=0.5d0
      QB(1:3,1) = QBAS(1:3,1)/N1
      QB(1:3,2) = QBAS(1:3,2)/N2
      QB(1:3,3) = QBAS(1:3,3)/N3
C- index for k in 1st BZ. See genqbz in BZ.FOR.
      kount      = 0
      do      i1 = 1,n1
        do      i2 = 1,n2
          do      i3 = 1,n3
            kount    = kount + 1
            indexk(i1-1,i2-1,i3-1) = kount

c  for check
            if(chk) then
              qbzx(1:3)= qb(1:3,1)*(i1-1+hf) +qb(1:3,2)*(i2-1+hf) +qb(1:3,3)*(i3-1+hf)
ccccccccccccccccccccc
c      write(6,"(3d26.18)") rk(1:3,kount)
c      write(6,"(3d26.18)") qbzx(1:3)
c      write(6,*)
ccccccccccccccccccccc
              if(  sum(abs(rk(1:3,kount)-qbzx(1:3)) ) > epss ) 
Cstop2rx 2013.08.09 kino     &    stop 'tetfbzf: rk /= qbzx'
     &    call rx( 'tetfbzf: rk /= qbzx')
            endif
          end do
        end do
      end do
      chk=.false.
Cstop2rx 2013.08.09 kino      if (kount /= n1*n2*n3) stop ' kount: wrong no. k-points'
      if (kount /= n1*n2*n3) call rx( ' kount: wrong no. k-points')
Cstop2rx 2013.08.09 kino      if (nqbz  /= n1*n2*n3) stop ' kount: wrong no. k-points'
      if (nqbz  /= n1*n2*n3) call rx( ' kount: wrong no. k-points')
cccccccccccccccccccccccccccc
c      DO  j1 = 0,n1-1
c      DO  j2 = 0,n2-1
c      DO  j3 = 0,n3-1
c        write(6,"(' j1j2j3=',3i4,' ix=',i6)") J1,J2,J3,indexk(J1,J2,J3)
c      enddo
c      enddo
c      enddo
c      write(6,*)
cccccccccccccccccccccccccccc
      kount      = 0
      do      i1 = 1,n1+1
        do      i2 = 1,n2+1
          do      i3 = 1,n3+1
            kount    = kount + 1
            indexkw(i1-1,i2-1,i3-1) = kount
            qbzw(1:3,kount) =
     &    qb(1:3,1)*(i1-1+hf) + qb(1:3,2)*(i2-1+hf) + qb(1:3,3)*(i3-1+hf)
            ib1bz(kount) = indexk(mod(i1-1,n1), mod(i2-1,n2), mod(i3-1,n3))
          end do
        end do
      end do
c
      CALL CCUTUP(QB,QB1,IBTR,KCUT) !This is from LMTO-3
      ntet = 0
C ----- START LOOPING OVER MICROCELLS ---------
      DO  20  I3 = 1, N3
      DO  20  I2 = 1, N2
      DO  20  I1 = 1, N1
C ----- SET UP IDENTIFIERS AT 8 CORNERS OF MICROCELL ------
c        write(6,*)
        DO  1  K1 = 0, 1
          J1 = I1 -1 + K1
        DO  1  K2 = 0, 1
          J2 = I2 -1 + K2
        DO  1  K3 = 0, 1
          J3 = I3 -1 + K3
          IMC(K1,K2,K3)  = indexkw(J1,J2,J3)
cccccccccccccc
c         write(6,"(' k1k2k3=',3i4,' ix=',i6,' ii=',2i6)")J1,J2,J3
c     &    ,indexkw(J1,J2,J3), ib1bz(indexkw(J1,J2,J3))
c     &    ,indexk(mod(j1,n1), mod(j2,n2), mod(j3,n3))
cccccccccccccc
   1    continue
C ----- LOOP OVER TETRAHEDRA --------------
        DO  10  ITET = 1, 6
          DO  2  IC = 1, 4
            K1 = KCUT(1,IC,ITET)
            K2 = KCUT(2,IC,ITET)
            K3 = KCUT(3,IC,ITET)
            IQ(IC) = IMC(K1,K2,K3)
    2     continue
          ntet=ntet+1
          do  i = 1, 4
            idtet(i,ntet) = iq(i)
          enddo
   10   CONTINUE
   20 CONTINUE
      write(6, "(1x,'TETFBZF: ',2i8 )")  ntet, 6*n1*n2*n3
      END
