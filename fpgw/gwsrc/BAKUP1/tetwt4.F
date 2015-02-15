      subroutine tetwt4x_dev(omg,q, eband1,eband2,
     i    qbas,ginv,efermi,   
     i    ntetf, nqbzw, nband,nqbz,
     i    nctot,ecore,  idtetf,qbzw,ib1bz, chipm, 
     o    wgt,nwgt,
     i   iq,isp1,isp2,nqibz) ! new input for mtett
      use keyvalue
C- Obtain weights for Dielectric function by tetrahedron method.
Ci omg     = omega in Hartree unit! complex!  omg might be 1d-8 or so
Ci           to choose the correct log branch for intvvc and so on.
Ci q       = q-vector in x(q,iw)  2*pi*q(1:3)/alat is the true q.
Ci qbas    = base reciprocal lattice vectors
Ci ginv    = inverse of qbas s. indxrk.f
Ci efermi  = Fermi level in Rydberg
Ci eband1 = eigenvalues for   k
Ci eband2 = eigenvalues for q+k
Ci ecore  = eigenvalues for core
Ci      ,where k runs in the Full 1st BZ, not only in IBZ.
Ci idtetf,qbzw,ib1bz = outputs from TETFBZF, which should be called before calling tetwt4.
Ci nctot   = the number of core
Ci nband   = total number of bands
Ci nqbz    = number of k-points in the 1st BZ
Ci n1,n2,n3= divisions along base reciprocal lattice vectors
Co wgt : integration weight for the Lindhard function. complex. wgt(band indexfor k, band indexfor q+k,k)
Co nwgt: the number of non-zero weights for each k.
C-----------------------------------
Cr ek   = eigenvalues at k-points in the 1st BZ
Cr ekq  = eigenvalues at k+q, k   in the 1st BZ
C
Cr
Cr See J.Rath&A.J.Freeman PRB11(6) p.2109(1975).
Cr
Crr The numbering of band index is not unique due to degeneracy.
Crr It affects how to choose tetrahedron.  As a result, it affect on
Crr the integration weight, and it might break the crystal symmetry.
Crr So I add symmetrization at the last in this routine so as to recover the symmetry.
Cr
Cr subroutine lindtet3 is the main part to calculate
Cr the microcell integral, \int dk1 dk2 dk3 f(e(k)) (1-f(e(q+k)))/ (omg- e(q+k) + e(k) ),
Cr where e(q+k) is unoccupied, e(k) is occupied,
Cr and f(E) denote the Fermi distribution funciton. Only for T=0.
Cr omg is complex. omg could be 1d-8 or so to choose the correct log branch.
Cr
Cr Their algorism requires rather tedious classifications due to the input eigenvalues.
Cr So I am afraid that the smoothness might be not fine.
Cr There are some room to improve this routine, mainly in intttvc,
Cr which is the core part of the tetrahedron method.
Cr The parameters eps and so on in the funtion intttvc might not be optimum.
Cr These parameters (and classification) are necessary to avoid to give out of range argumets.
Cr e.g., to avoid denominators in intv4c divergent.
c takao kotani Aug 2000.
c
c  I improved intttvc to intttvc2. It sould be a bit better but
c  I am not completely sure whetehre there are no numerically dengerous case or not.
c takao kotani Apr 2002
      implicit none
c---in out -------------------------------
      integer(4) :: ntetf, nqbzw,nband,nqbz,nctot,
     &   idtetf(0:3,ntetf),ib1bz(nqbzw), nwgt(nqbz)
      real(8) :: q(3),
     &   eband1(nband,nqbz), eband2(nband,nqbz),
     &   qbas(3,3),ginv(3,3),efermi,
     &   ecore(nctot), qbzw(3,nqbzw)
      complex(8):: omg, wgt(nband+nctot,nband,nqbz)
c----------------------------------
      integer(4):: ntet, itet,ic, ib,jb
      real(8) :: ek(nband,0:3), ekq(nband,0:3),
     &    qk(3),qkm(3),qbz(3) !,qbzm(3)  !,rk(3,nqbz)
C- For tetrahedra
      integer(4):: kk(0:3),kq(0:3),kr(0:3),i ,j
      real(8)   ::  det33 , 
     &   kvec(3,0:3), ea(0:3), eb(0:3) ,x(0:3),am(3,3)
      complex(8):: wtt(0:3,3)
c
      integer(4):: noccx_kxx, noccx_k, noccx_kq, noccx1
      real(8),parameter:: eps=0d0 !1d-12 ! cutoff check to determine cancellation.
      logical ::prt=.false.
      real(8) :: ek_(nband+nctot,0:3), ekq_(nband+nctot,0:3),voltot,volt
c
      integer(4):: irnk1, nrank1, irnk2, nrank2,nibib,kx,nrankc,
     &  ires((nband+nctot)**2), iof1,iof2,
     &  ini1(nband+nctot),ied1(nband+nctot), ixi1,ixi2,ixe1,ixe2,
     &  ini2(nband+nctot),ied2(nband+nctot)
      real(8):: ekqxx(nband+nctot),summ
c      real(8):: ekxx1(nband+nctot,nqbz),ekxx2(nband+nctot,nqbz)
      complex(8)::  wmean,wttx,a1,a2,wgtA,wgtB
      integer(4):: kqxx(nqbz), ik,ibx,nbnc,nb1
      logical :: ipr=.false.
      real(8),parameter:: pii=3.1415926535897932d0

      integer:: idim,ivec !,ixc
      logical :: chipm
      complex(8),allocatable:: wgt1(:,:)

c-devided tet
      integer(4):: nmtet,nqbzwm,nqbzm,ntetfm
      integer(4),allocatable:: idtetfm(:,:,:),ib1bzm(:)!,index_qbzm(:,:,:)
      real(8),allocatable:: qbzm(:,:),qbzwm(:,:)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

      integer(4)::kkm(0:3),kvecm(1:3, 0:3),ifeig,ifmtet,nqnumm,kqxxm,kqxxm2
     &  ,n,nsp,nband_x,kp,im !,n_index_qbzm
      real(8),allocatable:: ekzz1(:,:),ekzz2(:,:),eigtet(:,:,:),wtet(:,:,:)
      real(8)::qdummy(3)

      integer(4)::iq,nqibz,mtet(3),  iqx,ispx, ix1,iy1,iz1,verbose,iqindx
      logical ::mtett !,readgwinput

      real(8):: scissors_x0
      logical :: testsemif,chkwrt=.false.,iwarn=.true.
      integer(4):: isp1,isp2
c-----------------------------------------------------------------------
      if(verbose()>=80) chkwrt=.true.
      ipr=chkwrt
      if(ipr) write(6,*) ' tetwt4x_dev: '
      voltot = abs(det33(qbas))
      ntet   = ntetf !6*n1*n2*n3
*poption noparallel
      wgt    = 0d0

c
c      isp1=isp
c      isp2=isp
c      if(ixc==22.or.ixc==23) then
c         if(isp==1) isp2=2
c         if(isp==2) isp2=1
c      endif
c--------------------------------------------------------------------
c      call cputid  (0)
c      write(6,*) ' *** q=',q(1:3)
c      write(6,*) ' qbas ='; write(6,"(3f10.5)")(qbas(1:3,i),i=1,3)
c      write(6,*) ' vol tot    =',voltot
c      write(6,*) ' vol of tet =',det33(qbas)/(n1*n2*n3*6)
c      write(6,*) ' efermi=',efermi
c      write(6,*) ' nctot nband =',nctot,nband
c      write(6,*) ' nqbz  =',nqbz
c      do i=1,nctot
c      write(6,"('  ecore=',i3,f10.5)")i, ecore(i)
c      enddo

c--- devided-tetrahedron method.
      mtet=(/1,1,1/)
c      if(readgwinput())
      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
      if(sum(abs(mtet))/=3) then
        mtett=.true.
        if(sum(abs(mtet))>3) then
          write(6,*) ' we use devided-tetrahedron scheme mtet=',mtet
        elseif(sum(mtet)/=6) then
Cstop2rx 2013.08.09 kino          stop ' now only mtet 2 2 2 is allowed '
          call rx( ' now only mtet 2 2 2 is allowed ')
        endif
      else
        mtett=.false.
      endif

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
        read(ifmtet) idtetfm,ib1bzm,qbzm,qbzwm, wtet !,index_qbzm
cccccccccccccccccccccccccccccccccccccccc
c        do ix1=1,n_index_qbzm
c        do iy1=1,n_index_qbzm
c        do iz1=1,n_index_qbzm
c          write(6,"(' ix iy iz=',3i5,' index=',i5)")ix1,iy1,iz1,index_qbzm(ix1,iy1,iz1)
c        enddo
c        enddo
c        enddo
c        stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxx test xxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccc
        close(ifmtet)
        ifeig=501
        open (ifeig,file='eigmtet',form='unformatted')
        read(ifeig) nband_x,nqnumm,nsp
        write(6,*)'readin eigmtet ', nband_x,nqnumm,nsp
        if(nband_x /=nband ) then
Cstop2rx 2013.08.09 kino          stop 'tetwt5: nband_x /=nband'
          call rx( 'tetwt5: nband_x /=nband')
        endif
        allocate( eigtet(nband,nqnumm,nsp) )
        do iqx=1,nqnumm
          do ispx=1,nsp
            read (ifeig) eigtet(1:nband,iqx,ispx)
            if(verbose()>50) write(6,"('iq=',i3,'  eig(1:5)=',5f10.4)")iqx,eigtet(1:5,iqx,ispx)
          enddo
        enddo
        close(ifeig)
      else
        nmtet = 1
        nqbzwm= nqbzw
        nqbzm = nqbz
        allocate( 
     &       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
     &       ib1bzm(nqbzwm),           qbzm(3,nqbzm), wtet(0:3,nmtet,ntetf) )
        idtetfm(:,1,:)=idtetf
        qbzwm =qbzw
        ib1bzm=ib1bz
        wtet=0.25d0
      endif
c---------------------------------------------------------
      allocate( wgt1(nband+nctot,nband) )
      nbnc = nband+nctot
      nb1  = nband+1
      allocate( ekxx1(nband+nctot,nqbz),ekxx2(nband+nctot,nqbz))
      do kx = 1, nqbz
        ekxx1( 1:nband, kx) = eband1(1:nband,kx)
        ekxx2( 1:nband, kx) = eband2(1:nband,kx)
        ekxx1( nband+1: nband+nctot, kx) = ecore(1:nctot)
        ekxx2( nband+1: nband+nctot, kx) = ecore(1:nctot)
      enddo

! eigenvalues at q and q+k ---------------------------------------
!  ekzz1 for q
!  ekzz2 for q+k.
      allocate( ekzz1(nband+nctot,nqbzm),ekzz2(nband+nctot,nqbzm))
      if(mtett) then
        do kx = 1, nqbzm
          ekzz1( 1:nband, kx)              = eigtet(1:nband,kx,isp1)
          ekzz1( nband+1: nband+nctot, kx) = ecore(1:nctot)
          ekzz2( nband+1: nband+nctot, kx) = ecore(1:nctot)
        enddo
        if(iq<=nqibz) then
          do kx = 1, nqbzm
!          q(1:3)+qbzm(1:3,kx) ----> kqxxm

c      write(6,*)' xxx3 ginv=',ginv
c      write(6,*)' q+qbzm= ', q(1:3) + qbzm(1:3,kx)
c          call fbz2 ( q(1:3) + qbzm(1:3,kx),
c     i      ginv ,index_qbzm,n_index_qbzm,qbzm,nqbzm,
c     o      qdummy,kqxxm)
cc--------
c          call iqindx2(q(1:3)+qbzm(1:3,kx),ginv,qbzm,nqbzm,kqxxm, qdummy)
c         write(6,*) 'kqxxm2 kqxxm=',kqxxm2,kqxxm
            kqxxm = iqindx(q(1:3)+qbzm(1:3,kx),ginv,qbzm,nqbzm)
            ekzz2(1:nband, kx) = eigtet(1:nband,kqxxm,isp2) !ekzz1(1:nband, kqxxm)
          enddo
        else
          do kx = 1, nqbzm
            kp  = nqbzm *(iq - nqibz) + kx
            ekzz2(1:nband, kx)= eigtet(1:nband,kp,isp2)
          enddo
        endif
      else
        ekzz1=ekxx1
        ekzz2=ekxx2
      endif

c--- 2004dec
      if(scissors_x0()/=0d0) then
        call addsciss(scissors_x0(),efermi,(nband+nctot)*nqbzm, 
     o   ekzz1)
        call addsciss(scissors_x0(),efermi,(nband+nctot)*nqbzm, 
     o   ekzz2)
      endif

C- Check
c      volt = 0d0
c      do itet = 1, ntet
c        kvec(1:3,0:3) = qbzw (1:3, idtetf(0:3,itet) )
c        do i = 1,3
c        kvec(1:3,i) = kvec(1:3,i) - kvec(1:3,0)
c        enddo
c        volt = volt + abs(det33(kvec(1:3,1:3))/6d0)
c      enddo
c      if(abs(volt-voltot)>1d-10) stop ' tetwt: abs(volt-voltot)>1d-10'
C- Check
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
c      write(6,*)' xxx4 volt voltot=',volt,voltot
Cstop2rx 2013.08.09 kino      if(abs(volt-voltot)>1d-10) stop ' tetwt: abs(volt-voltot)>1d-10'
      if(abs(volt-voltot)>1d-10) call rx( ' tetwt: abs(volt-voltot)>1d-10')
c      endif

C- kqxx(kx)   ekxx
c      do kx = 1, nqbz
c        qk  (1:3) =  rk(1:3,kx) + q(1:3)
c        qkm (1:3) = -qk(1:3)                     ! qkm = -k-q
c        call fbz (qk,  n1,n2,n3,qbas,ginv,indexk,
c     o     qbz,  kqxx(kx) )
c       write(1212,*) kx, kqxx(kx)
c      write(6,*)' ekkdf chk=',
c     & sum(ekxx1( 1:nband, kqxx(kx)) - ekxx2( 1:nband, kx)),
c     & sum(ekxx1( 1:nband, kx))
c      enddo
ccccccccc
c      do kx = 1, nqbz
c        call rwdd1 (ifev, kx, nband, ekxx(1:nband, kx))
c        ekxx1(nband+1: nband+nctot, kx) = ecore(1:nctot)
c      enddo

C- Loop over tetrahedron
*poption parallel
*poption psum(wgt( nband+nctot,nband,nqbz ) )
*poption tlocal( kk,kvec,kq,ek,ek_,ekq,ekq_ )
*poption tlocal( wttx,x, wgt1 )
*poption indep ( omg ,efermi )

      do 1000 itet = 1, ntet
        kk (0:3)       = ib1bz( idtetf(0:3,itet) ) !  k
        do 1100 im=1,nmtet
          kkm (0:3)      = ib1bzm( idtetfm(0:3,im,itet) ) !  k
          kvec(1:3, 0:3) = qbzwm ( 1:3, idtetfm(0:3,im,itet) )
          ek  ( 1:nband,       0:3) = ekzz1( 1:nband,       kkm(0:3)) ! k
          ek_ ( 1:nband+nctot, 0:3) = ekzz1( 1:nband+nctot, kkm(0:3)) ! k
          ekq ( 1:nband,       0:3) = ekzz2( 1:nband,       kkm(0:3)) ! k+q
          ekq_( 1:nband+nctot, 0:3) = ekzz2( 1:nband+nctot, kkm(0:3)) ! k+q
          noccx_k  = noccx1 (ek  ,4,nband, efermi)  !the highest number of occupied states (HNOS)
          noccx_kq = noccx1 (ekq ,4,nband, efermi)

C*poption parallel
C*poption indep ( wgt,omg,kvec,ek_ ,ekq_)
C*poption tlocal( wttx,x )
c----------
          do ibx = 1, noccx_k + nctot  !   occupied
            do jb = 1, nband         ! unoccupied
c          if(prt) write(6,*) '****** ib jb=',ib,jb
              if(ibx <= noccx_k ) then
                ib = ibx
              else
                ib = ibx - noccx_k + nband
              endif
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$           if(ib==4.and.jb==4.and.itet==2) then
c$$$             goto 1820
c$$$           endif
c$$$           wgt1(ib,jb) =  0d0
c$$$           cycle
c$$$ 1820      continue
c$$$           write(6,*)' test itet ib jb=',itet,ib,jb
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccc
c----------
c        do ib = 1, nband + nctot !   occupied
c        do jb = 1, nband         ! unoccupied
c          if(prt) write(6,*) '****** ib jb=',ib,jb
c          if( .not.(ib <= noccx_k .or. ib >= nband+1)) cycle
c----------
              x(0:3)  = .5d0*(ek_ (ib, 0:3) - ekq_(jb, 0:3))         ! + omg !Denominator. unit in Hartree.
ccccccccccccccccc
c          if(testsemif().and.ibx==1) then
c!             x(0:3) = -0.5d0*3.7d0 !test for Na
c             x(0:3) = -0.5d0*1.14d0 !test for K
c          endif
ccccccccccccccccc
              if( sum(abs(x(0:3)+omg))<1d-8 ) cycle
              if(chkwrt) then
                write(6,"('### Goto lindtet3: itet ib jb Ef='
     &                    ,i8,2i5,d12.4,' ###')" ) itet,ib,jb,efermi
                write(6,"('  ek_ - Ef= ',4f10.3)") ek_  (ib, 0:3)-efermi
                write(6,"('  ekq_- Ef= ',4f10.3)") ekq_ (jb, 0:3)-efermi
                write(6,"('  -x(a.u.)= ',4f10.3)") -x(0:3)
              endif
cccccccccccccccccccc
c            omg= 0.150d0/2d0+ (0d0,1d-8)
ccccccccccccccccccc
c 1818     continue
ccccccccccccccccccccc
c          if(chkwrt) then
c            omg= omg+0.0001d0
c          endif
cccccccccccccccccccc
              if(chipm.and.isp1==2) then
                wgtA=0d0
              else
                call lindtet3(omg,kvec,ek_(ib, 0:3),ekq_(jb, 0:3),
     &      x,efermi,wttx)   ! kvec(1:3, 0:3), ea(0:3), x(0:3)
                wgtA = wttx/4d0/voltot
              endif

              if(chipm.and.isp1==1) then
                wgtB = 0d0
              elseif(abs(omg)<1d-10.and.(.not.chipm)) then  ! in the static case.
                wgtB = wgtA
              else
                call lindtet3(-omg,kvec,ek_(ib, 0:3),ekq_(jb, 0:3),
     &      x,efermi,wttx)
                wgtB = wttx/4d0/voltot
              endif

              wgt1(ib,jb) =  wgtA + wgtB

              if(chkwrt) then
                write(6,"('ttt',3i4,f10.5,d12.4)") itet,ib,jb,dreal(omg),-dimag(wttx/pii)
c            endif
c            if(dreal(omg)>0.20d0/2d0) goto 1819
c            goto 1818
c 1819       continue
c            if(itet==ntetf) then
c             write(6,*)' tetwt5: itet==ntetf'
c             stop 'test end'
c            endif
              endif

c--- check for too large wttx -------------
#ifndef COMMENTOUTfor_PARALLEL
c            write(6,"(' ib jb itet= ',3i5,'  wgt= ',4d13.6)")
c     &      ib,jb,itet, wgt1(ib,jb)-wgt2,wgt2
              if((abs(wgtA)>1d0.or.abs(wgtB)>1d0)
     &       .and.omg/=1d30.and.iwarn) then
                write(6,*)' --- tetwt4: WARN! large tetrahedron weigt.    ',  
     &      ' The denominator omg+x for lindtet3 is rather large.   ',
     &      ' This warn does not mean bug or fail of the calculation. ',
     &      ' But you may have to check the convergence for num of k. ',
     &      ' Or you could try another real-axis energy mesh!'
                write(6,"(' ib jb itet= ',3i5,'  wgt= ',4d13.6)") 
     &      ib,jb,itet, wgtA,wgtB
                write(6,*) ' 1 goto lindtet3  ef= ',efermi
                write(6,*) ' omg=',omg
                write(6,*) (kvec(1:3,i),i=0,3)
                write(6,"('x    =',4d23.16)") x    (0:3)
                write(6,"('ek_  =',4d23.16)") ek_  (ib, 0:3)
                write(6,"('ekq_ =',4d23.16)") ekq_ (jb, 0:3)
                write(6,*)' tetwt4: abs(wttx/4d0/voltot)>1d0 too large ?' !It seems too large
                write(6,"('------- we omit this WARN hereafter! -------')")
                iwarn=.false.
              endif
#endif
c-----------
            enddo
          enddo
          wgt(1:noccx_k,:,kk(0))= wgt(1:noccx_k,:,kk(0))+ wgt1(1:noccx_k,:)* 4*wtet(0,im,itet)
          wgt(nb1:nbnc,: ,kk(0))= wgt(nb1:nbnc,: ,kk(0))+ wgt1(nb1:nbnc, :)* 4*wtet(0,im,itet)
          wgt(1:noccx_k,:,kk(1))= wgt(1:noccx_k,:,kk(1))+ wgt1(1:noccx_k,:)* 4*wtet(1,im,itet)
          wgt(nb1:nbnc,: ,kk(1))= wgt(nb1:nbnc,: ,kk(1))+ wgt1(nb1:nbnc, :)* 4*wtet(1,im,itet)
          wgt(1:noccx_k,:,kk(2))= wgt(1:noccx_k,:,kk(2))+ wgt1(1:noccx_k,:)* 4*wtet(2,im,itet)
          wgt(nb1:nbnc,: ,kk(2))= wgt(nb1:nbnc,: ,kk(2))+ wgt1(nb1:nbnc, :)* 4*wtet(2,im,itet)
          wgt(1:noccx_k,:,kk(3))= wgt(1:noccx_k,:,kk(3))+ wgt1(1:noccx_k,:)* 4*wtet(3,im,itet)
          wgt(nb1:nbnc,: ,kk(3))= wgt(nb1:nbnc,: ,kk(3))+ wgt1(nb1:nbnc, :)* 4*wtet(3,im,itet)
ccccccccccccccccccccccccccc
c        write(6,*)' itet=',itet, sum(abs(wgt1)),' xxx=',sum(abs(wgt))
c        write(6,*)' itet=',itet, noccx_k,nb1,nbnc,kk
cccccccccccccccccccccccccc
 1100   continue
 1000 continue


C- Symmetrization of wgt. ! We need special cares for degenerate cases.
*poption parallel
*poption indep( ekxx1,ekxx2,ekxx1, nband ,ipr, ecore, nctot ,wgt)
*poption tlocal( nrank1, ini1,ied1, nrank2, ini2,ied2, nrankc )
      do kx = 1, nqbz     ! ipr = .false.; if(ipr)  write(6,*)' kx =',kx
        call chkdgn( ekxx1(:,kx), nband, nrank1, ini1,ied1,0 ,ipr)
        call chkdgn( ekxx2(:,kx), nband, nrank2, ini2,ied2,0 ,ipr)
        nrankc = 0
        if(nctot/=0) then
          call chkdgn(ecore, nctot, 
     &    nrankc, ini1(nrank1+1), ied1(nrank1+1), nband,ipr)
        endif
c       if(ipr) write(6,*)' kx nrank =',kx,nrank1,nrank2,nrankc
        do irnk1 = 1, nrank1 + nrankc
          do irnk2 = 1, nrank2
            ixi1  = ini1(irnk1);  ixe1 = ied1(irnk1)
            ixi2  = ini2(irnk2);  ixe2 = ied2(irnk2)
            wmean = sum( wgt(ixi1:ixe1,ixi2:ixe2, kx ) )
     &            / ((ixe1-ixi1+1)*(ixe2-ixi2+1) )
            wgt(ixi1:ixe1,ixi2:ixe2, kx ) = wmean
          enddo
        enddo
      enddo
c
      do ik =1, nqbz
        nwgt(ik) = 0
        do i =1, nband+nctot
          do j =1, nband
Cstop2rx 2013.08.09 kino            if(abs(wgt(i,j,ik) )>0d0 .and. j>nband) stop 'tetwt4: bug ?' ! bug checker
            if(abs(wgt(i,j,ik) )>0d0 .and. j>nband) call rx( 'tetwt4: bug ?')
            if(abs(wgt(i,j,ik) )>eps ) then
              nwgt(ik) = nwgt(ik)+1
cc          write(56,"(3i4,2d13.6)") i,j,ik, wgt(i,j,ik)
c          if(abs(wgt(i,j,ik) )>0.1d0 )
c     &      write(6,"(' k i ik=',3i4,'  wgt= ',2d13.6)")
c     &      i,j,ik, wgt(i,j,ik)
            endif
          enddo
        enddo !;  write(6,*)' ik=',ik,' num of nonzero wgt=',nwgt(ik)
      enddo
c      write(6,*)' max num of nonzero wgt(k)=',maxval(nwgt)
c      write(6,*)' tot num of nonzero wgt   =',sum(nwgt)
c      write(6,*)' sum of wgt   =',sum(wgt)
c      write(6,*)' tetwt4: end '; call cputid  (0)

      write(6,"(' tetwt4_dev: omega maxval(nwgt(1:nqbz)) sum(nwgt)=',2f8.4,2i7
     &  ,' sum(wgt)=',2d13.6)") omg*2d0,maxval(nwgt),sum(nwgt),sum(wgt)
      end




c----------------------------------------------------------------------
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

c-------------------
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

c-------------------------------------------------------------------
      subroutine lindtet3( omg, kvec, ea, eb, x, efermi,  ! input
     o    wttx)                                           ! output
C- Calculate the integral of \int dk1 dk2 dk3 f(ea)(1-f(eb))/x
C f(E) denote the Fermi distribution funciton. Only for T=0.
C Tetrahedon is specified by values at 4 corners; kvec(1:3, 1:4), ea(1:4), eb(1:4), x(1:4)
Cr This code is based on J.Rath&A.J.Freeman PRB11(6) p.2109(1975).
C---------------------------------------------------------------------------------
      implicit none
      integer(4) :: ieaord(1:4),i,isig,n,itmp,ix
      real(8)    ::  kvec(1:3, 1:4), x(1:4), ea(1:4),
     &  kk(3,1:4),xx(1:4),ee(1:4), am(3,3), eb(1:4),ebf(1:4),ebfKx(1:4),
     &  vcell,efermi,etest , Kx(1:3,4),xKx(4)
      complex(8) :: omg, wttx,inttetra3
      integer(4):: verbose
      logical :: debug=.false.
c-------------------------------
      debug=.false.
      if(verbose()>=110) debug=.true.
      if(debug) then
        write(996,*)' lindtet3: *****************************'
        write(996,"(' i=',i3,' x ea eb =',3f10.5)") 
     &  (i, x(i), ea(i)-efermi, eb(i)-efermi,i=4,1,-1)
      endif
ccccccccc
#ifdef EXPAND_SORTEA
      n=4
      isig = 1
*poption noparallel
      do i = 1,n
        ieaord(i) = i
      enddo
*poption noparallel
      do ix= 2,n
*poption noparallel
        do i=ix,2,-1
          if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
C#ifdef EXPAND_ISWAP
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp

C#else
C            call iswap (ieaord(i-1),ieaord(i))
C#endif
            isig= -isig
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
      if(debug) then
        write(996,"(' i iea=',2i4)") (i,ieaord(i),i=1,4)
        write(996,"(' i=',i3,' xx ee ebf =',3f10.5,'  kk=',3f10.5)") 
     &     (i,xx(i),ee(i),ebf(i),kk(1:3,i),i=4,1,-1)
      endif
ccccccccccccc
      if( 0d0<=ee(4) ) then
        wttx = (0d0,0d0)
      elseif( ee(4) < 0d0 .and. 0d0<= ee(3) ) then   !!! Fig 1.
        if(debug) write(996,*) 'lindtet3: fig 1 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  Kx(1,1),xKx(1),ebfKx(1))
        !K1 -> Kx(:,1), x(K1) -> xkx(1). K1 is on the like k4---k2.
        call midk3(kk,ee,xx,ebf, 4,1,  Kx(1,2),xKx(2),ebfKx(2)) !K2
        call midk3(kk,ee,xx,ebf, 4,3,  Kx(1,3),xKx(3),ebfKx(3)) !K3
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4,10,20,30/) )  ! k4,K1,K2,K3
      elseif( ee(3) < 0d0 .and. 0d0<= ee(2) ) then   !!! Fig 2.
        if(debug) write(996,*) 'lindtet3: fig 2 xxx'
        call midk3(kk,ee,xx,ebf, 4,2,  Kx(1,1),xKx(1),ebfKx(1)) !K1
        call midk3(kk,ee,xx,ebf, 4,1,  Kx(1,2),xKx(2),ebfKx(2)) !K2
        call midk3(kk,ee,xx,ebf, 3,1,  Kx(1,3),xKx(3),ebfKx(3)) !K3
        call midk3(kk,ee,xx,ebf, 3,2,  Kx(1,4),xKx(4),ebfKx(4)) !K4
ccccccccccc
c         write(6,*)'xxxx inttetra='
c     &       ,inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4, 3,10,20/))  ! k4,k3,K1,K2
c         write(6,*) omg,kk,xx,ebf
c         write(6,*) ' Kx=',Kx,xKx,ebfKx
c         write(6,*)'xxxx inttetra='
c     &       ,inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,20,30,10/))  ! k3,K2,K3,K1
c         write(6,*) omg,kk,xx,ebf
c         write(6,*) ' Kx=',Kx,xKx,ebfKx
c         write(6,*)'xxxx inttetra='
c     &       ,inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,10,30,40/))  ! k3,K1,K3,K4
c         write(6,*) omg,kk,xx,ebf
c         write(6,*) ' Kx=',Kx,xKx,ebfKx
c     stop ' ***************** end test ttt *****************'
ccccccccccc
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4, 3,10,20/))  ! k4,k3,K1,K2
     &       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,20,30,10/))  ! k3,K2,K3,K1
     &       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3,10,30,40/))  ! k3,K1,K3,K4
      elseif( ee(2) < 0d0 .and. 0d0<= ee(1) ) then   !!! Fig 3.
        if(debug) write(6,*) 'lindtet3: fig 3 xxx'
        call midk3(kk,ee,xx,ebf, 1,4,  Kx(1,1),xKx(1),ebfKx(1)) !K1
        call midk3(kk,ee,xx,ebf, 1,2,  Kx(1,2),xKx(2),ebfKx(2)) !K2
        call midk3(kk,ee,xx,ebf, 1,3,  Kx(1,3),xKx(3),ebfKx(3)) !K3
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/3, 4,30, 2/))  ! k3,k4,K3,k2
     &       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4,10,20,30/))  ! k4,K1,K2,K3
     &       + inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfKx,(/4, 2,20,30/))  ! k4,k2,K2,K3
cccccccccc test fulled cccccccccccccccccccccccccccccccccccc
cc     &        + inttetra(kk,xx,Kx,xKx, (/1, 10,20,30/) )  ! test fulling
cc           write(6,*)' fulled test'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
        if(debug) write(996,*) 'lindtet3: fig 4 xxx'
        wttx = inttetra3(omg,kk,xx,ebf,Kx,xKx,ebfkx,(/1,2,3,4/))  ! k1,k2,k3,k4
      endif
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

c-----------------------------------
      complex(8) function 
     & inttetra3(omg,kk_,xx_,ebf,Kx_,xKx_,ebfKx,itetx)
C calculate tetrahedron integral Eq.(16).
C the four corners and denoted by itetx.
Ci kk (k1-k4) and xx (value of denominator)
Ci Kx (K1-K4) and xkx
Cr The four corners are selected from 8 points. It is specified by itetx.
      implicit none
      integer(4):: itetx(4),ix,i,ieaord(4),isig,n,itmp
      real(8) ::  kk_(3,1:4),xx_(1:4),Kx_(3,1:4),xKx_(1:4),
!!!!!!!!!!!!!!!!!      !BUGFIX ----> kx xkx ___>kx_,xkx_ 1999 mar
     &    kk(3,1:4),xx(1:4),ebf(1:4),ebfin(1:4)
     &   ,Kx(3,1:4),xKx(1:4),ebfKx(1:4),ee(4)
     &   ,kin(3,4), xin(4)
      complex(8) :: omg,wttx,inttetrac
      logical :: debug=.false.
c--- kk xin ebfin ---
      do i = 1,4
        if( itetx(i) < 10 ) then   ! kk (k1-k4 in Paper)
          ix = itetx(i)
          xin(    i) = xx_(    ix)
          ebfin(i)   = ebf(ix)
          kin(1:3,i) = kk_(1:3,ix)
        else                       ! Kx (K1-K4 in Paper)
          ix = itetx(i)/10
          xin(    i) = xKx_(    ix)
          ebfin(i)   = ebfKx(ix)
          kin(1:3,i) = Kx_ (1:3,ix)
        endif
      enddo
c      write(6,*) ' ebf  =', ebf
c      write(6,*) ' ebfkx=', ebfkx
c      write(6,"(' i=',i3,' xin ebfin =',2f10.5,'  kin=',3f10.5)")
c     &  (i,xin(i),ebfin(i),kin(1:3,i),i=4,1,-1)

c ee decendent order ee(4)>ee(3)>ee(2)>ee(1)
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
C#ifdef EXPAND_ISWAP
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp

C#else
C            call iswap (ieaord(i-1),ieaord(i))
C#endif
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo

#else
      call sortea( ebfin,ieaord, 4 ,isig)  ! the order E_4,E_3,E_2,E_1 This suits for Rath&Freeman.
#endif
      kk(1:3,1:4) = kin  (1:3,ieaord(1:4))   ! 4 corners  denoted by kvec(:,1:4), ee(1:4), and xx(1:4)
      ee(    1:4) = ebfin(    ieaord(1:4))
      xx(    1:4) = xin  (    ieaord(1:4))
cccc
      if(debug) then
        write(996,*)' inttetra3: ***** '
        write(996,"(' i=',i3,' xx ee =',2f10.5,'  kk=',3f10.5)") 
     &  (i,xx(i),ee(i),kk(1:3,i),i=4,1,-1)
      endif
cccc
      if( 0d0>=ee(4) ) then
        wttx = (0d0,0d0)
      elseif( ee(4) > 0d0 .and. 0d0>= ee(3) ) then   !!! Fig 1.
        if(debug)  write(996,*)' intterra3: fig1'
        call midk(kk,ee,xx, 4,2,  Kx(1,1),xKx(1)) !K1 -> Kx(:,1), x(K1) -> xkx(1). K1 is on the like k4---k2.
        call midk(kk,ee,xx, 4,1,  Kx(1,2),xKx(2)) !K2
        call midk(kk,ee,xx, 4,3,  Kx(1,3),xKx(3)) !K3
        wttx = inttetrac(omg,kk,xx,Kx,xKx,  (/4,10,20,30/))  ! k4,K1,K2,K3
      elseif( ee(3) > 0d0 .and. 0d0>= ee(2) ) then   !!! Fig 2.
        if(debug)  write(996,*)' intterra3: fig2'
        call midk(kk,ee,xx, 4,2,  Kx(1,1),xKx(1)) !K1
        call midk(kk,ee,xx, 4,1,  Kx(1,2),xKx(2)) !K2
        call midk(kk,ee,xx, 3,1,  Kx(1,3),xKx(3)) !K3
        call midk(kk,ee,xx, 3,2,  Kx(1,4),xKx(4)) !K4
        wttx  = inttetrac(omg,kk,xx,Kx,xKx, (/4, 3,10,20/))  ! k4,k3,K1,K2
     &        + inttetrac(omg,kk,xx,Kx,xKx, (/3,20,30,10/))  ! k3,K2,K3,K1
     &        + inttetrac(omg,kk,xx,Kx,xKx, (/3,10,30,40/))  ! k3,K1,K3,K4
      elseif( ee(2) > 0d0 .and. 0d0>= ee(1) ) then   !!! Fig 3.
        if(debug)  write(996,*)' intterra3: fig3'
        call midk(kk,ee,xx, 1,4,  Kx(1,1),xKx(1)) !K1
        call midk(kk,ee,xx, 1,2,  Kx(1,2),xKx(2)) !K2
        call midk(kk,ee,xx, 1,3,  Kx(1,3),xKx(3)) !K3
        wttx  = inttetrac(omg,kk,xx,Kx,xKx, (/3, 4,30, 2/))  ! k3,k4,K3,k2
     &        + inttetrac(omg,kk,xx,Kx,xKx, (/4,10,20,30/))  ! k4,K1,K2,K3
     &        + inttetrac(omg,kk,xx,Kx,xKx, (/4, 2,20,30/))  ! k4,k2,K2,K3
      else
        wttx  = inttetrac(omg,kk,xx,Kx,xKx, (/1,2,3,4/) )  ! k1,k2,k3,k4
        if(debug)  write(996,*)' intterra3: fig4'
      endif
      inttetra3 = wttx
      end

c-----------------------------------
      complex(8) function inttetrac(omg,kk,xx, Kx,xKx,itetx)
C calculate tetrahedron integral Eq.(16).
C the four corners and denoted by itetx.
Ci kk (k1-k4) and xx (value of denominator)
Ci Kx (K1-K4) and xkx
Cr The four corners are selected from 8 points. It is specified by itetx.
      implicit none
      integer(4):: itetx(4),ix,i
      real(8) ::  kk(3,1:4),xx(1:4), Kx(3,1:4),xKx(1:4), am(3,3)
     &   ,kin(3,4), xin(4), det33, intttvv,  intttv,intvv,work
      complex(8) :: omg,intttvc2,aaa
      logical:: debug=.false.
c
      if(omg==1d30.or. omg==-1d30) then
        inttetrac = 1d0 !this mode is only for cont the non-zero pairs, occu and unocc
        return
      endif
      do i = 1,4
        if( itetx(i) < 10 ) then   ! kk (k1-k4 in Paper)
          ix = itetx(i)
          xin(    i) = xx(    ix)
          kin(1:3,i) = kk(1:3,ix)
        else                       ! Kx (K1-K4 in Paper)
          ix = itetx(i)/10
          xin(    i) = xKx(    ix)
          kin(1:3,i) = Kx (1:3,ix)
        endif
      enddo
      do i = 1,3
        am(1:3,i) = kin(1:3,i) - kin(1:3,4)
      enddo
      inttetrac = intttvc2(omg,xin) * abs(det33(am)/6d0) ! \omega (volume of tetrahedra) = abs(det33(am)/6d0) See Eq. (17).

cccccccccccccccccccccccccc
      if(debug) then
        write(996,"(' omg  abs(det33(am)/6d0) =',2d13.5,3x,4d13.5)") 
     &     omg,abs(det33(am)/6d0)
        write(996,"(' xin       =',4d13.5)") xin
        write(996,"(' inttetrac =',4d13.5)") inttetrac
        write(996,*)
      endif
ccccccccccccccccccccccccccccc
      end

c------------------------------------------------------
      complex(8) function intttvc2(omg,v)
C- Tetrahedron integral. Eq.(16). (except volume factor).
Ci v(1:4) + omg : denominator for integration.
Ci includ imaginary part
      implicit none
      complex(8) :: vvvc,intvvc,intv4c,intv3c,intv2c
      integer(4) :: ieaord(1:4),i,isig,idif(3),idf,ix,n,itmp
      real(8)    :: v(4),VV(4),v1,v2,v3,v4,d1,d2,d3,
     &              vvd,vvas,vvv,diffc
     &  ,aaa2,e1,e2,e3,e4,f0,f1,f3,ee
      logical ::debug=.false.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      real(8),parameter :: eps = 1d-4  ! Judge coincidence.
      real(8),parameter :: eps = 1d-3  ! Judge coincidence. Safer
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      complex(8) omg
#ifdef EXPAND_SORTEA
      n=4
      isig = 1
*poption noparallel
      do i = 1,n
        ieaord(i) = i
      enddo
*poption noparallel
      do ix= 2,n
*poption noparallel
        do i=ix,2,-1
          if( v(ieaord(i-1)) > v(ieaord(i) ) ) then
C#ifdef EXPAND_ISWAP
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
C#else
C            call iswap (ieaord(i-1),ieaord(i))
C#endif
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      VV(1:4) = V(ieaord(1:4)) !vv is ordered

#else
      call sortea( v,ieaord,4,isig);  VV(1:4) = V(ieaord(1:4)) !vv is ordered
#endif

      vvas = sum(abs(vv(1:4)))   ! ;write(6,*)'vvas=',vvas
      diffc=  max(eps,eps*vvas)
      idf  = 0
      do i = 1,3
        if( abs( VV(i)-VV(i+1) ) < diffc ) idf= idf+10**(3-i)
      enddo

      if(debug) then
        write(996,"(' xxx vv=',4d23.15,' idf= ',i3)")
     & vv(1),vv(2),vv(3),vv(4),idf
      endif

c      if(   abs( VV(1)-VV(4)   ) < eps2*vvas+1d-4)  idf = 111
      if(    idf==111 ) then  !2002 apr
        VVVC = (VV(1)+VV(2)+VV(3)+VV(4))/4
        vv(1) = vvvc -.5d0*diffc
        vv(2) = vvvc -.5d0*diffc
        vv(3) = vvvc +.5d0*diffc
        vv(4) = vvvc +.5d0*diffc
        idf = 101
      elseif(idf==011) then
        vv(4)= vv(3)+diffc
        idf = 010
      elseif(idf==110) then
        vv(1)=vv(2) -diffc
        idf = 010
      endif

      if(debug) then
        write(996,"(' yyy vv=',4d23.15,' idf= ',i3)")
     & vv(1),vv(2),vv(3),vv(4),idf
      endif

! Possible values of idf are denoted by ijk, where i,j and k can be 0 or 1.
! 1 means corresponding vv are the same.
      intttvc2 = (0d0,0d0)
      if(    idf==  000 ) then
        intttvc2 = intvvc(omg,vv)
      elseif(idf==  100 ) then
        VVV = (VV(1)+VV(2))/2
        V3  = VV(3)
        V4  = VV(4)
        intttvc2 = intv4c(omg,vvv,v3,v4)
      elseif(idf==  010 ) then
        VVV = (VV(2)+VV(3))/2
        V3  = VV(1)
        V4  = VV(4)
        intttvc2 = intv4c(omg,vvv,v3,v4)
      elseif(idf==  001 ) then
        VVV = (VV(3)+VV(4))/2
        V3  = VV(1)
        V4  = VV(2)
        intttvc2 = intv4c(omg,vvv,v3,v4)
      elseif(idf==  101 ) then
        VVV = (VV(1)+VV(2))/2
        VVd = (VV(3)+VV(4))/2
        intttvc2 = intv3c(omg,vvv,vvd)
      else
Cstop2rx 2013.08.09 kino        stop 'intttvc2: idf wrong...'
        call rx( 'intttvc2: idf wrong...')
c      elseif(idf==  110 ) then
c        VVV = (VV(1)+VV(2)+VV(3))/3
c        V4  = VV(4)
c        if( abs(vvv+omg)>eps3 ) intttvc2 = intv2c(omg,VVV,v4)
c      elseif(idf==  011 ) then
c        VVV = (VV(2)+VV(3)+VV(4))/3
c        V4  = VV(1)
c        if( abs(vvv+omg)>eps3 ) intttvc2 = intv2c(omg,VVV,v4)
c      elseif(idf==  111 ) then
c        VVVC = (VV(1)+VV(2)+VV(3)+VV(4))/4 + omg
c        if( abs(vvvc)>eps3 ) intttvc2 = 1d0/vvvc
      endif
      end

c-----------------------------------
      complex(8) function intvvc(omg,vv)
      implicit none
      complex(8):: omg,v1c,v2c,v3c,v4c,intvvcc
      real(8) :: vv(4),v1,v2,v3,v4,d1,d2,d3
      V1=VV(1); V2=VV(2); V3=VV(3); V4=VV(4)
      if(v1 ==0) v1 =1d-15
      if(v2 ==0) v2 =1d-15
      if(v3 ==0) v3 =1d-15
      if(v4 ==0) v4 =1d-15
      D1 = (V1-V4)*(V1-V2)*(V1-V3)
      D2 = (V2-V4)*(V2-V3)*(V2-V1)
      D3 = (V3-V4)*(V3-V1)*(V3-V2)
      v1c = V1+omg
      v2c = V2+omg
      v3c = V3+omg
      v4c = V4+omg
      intvvcc = V1c**2/D1*log(V1c/V4c)  ! Rath&Freeman Eq.(17)
     &          + V2c**2/D2*log(V2c/V4c)
     &          + V3c**2/D3*log(V3c/V4c)
      intvvc = 3*intvvcc
      end

c------------
      complex(8) function intv4c(omg,v,v3,v4) !case (iv) at page 2113.
      implicit none
      complex(8):: omg,v3c,v4c,vc,intv4cc,xfun,x3,x4,dxfun
      real(8) :: v,v3,v4
      if(v ==0) v =1d-15
      if(v3==0) v3=1d-15
      if(v4==0) v4=1d-15
      v3c = V3+omg
      v4c = V4+omg
      vc  = V +omg
c ---case4
      x3 = ( V -V3 )/V3c
      x4 = ( V -V4 )/V4c
      intv4c = 3d0* ( xfun(x3) - xfun(x4) )/(v3-v4)
c      write(996,"('intv4c: case4 x3',2d13.5)") x3
c      write(996,"('              x4',2d13.5)") x4


c--- case1
cc      write(996,"('intv4c: case1'))")
c      intv4cc= V3c**2/(v3-v)**2/(v3-v4)*log(V3c/Vc)
c     &     +   V4c**2/(v4-v)**2/(v4-v3)*log(V4c/Vc)
c     &     +   Vc/(v3-v)/(v4-v)
c      intv4c = 3*intv4cc

c--- case2
c      write(996,"('intv4c: case2'))")
c      intv4cc= 1d0/(v3-v4) * (
c     &     1d0/(v3-v) * ( V3c**2* log(V3c/Vc)/(v3-v) - Vc )
c     &   - 1d0/(v4-v) * ( V4c**2* log(V4c/Vc)/(v4-v) - Vc )
c     &         )
c      intv4c = 3*intv4cc

c---case3
c      write(996,"('intv4c: case3'))")
c      intv4c= 3d0/(v3-v4) * (
c     &   V3c**2/(v3-v)**2 * (log(V3c/Vc) -1d0 + vc/v3c)
c     & - V4c**2/(v4-v)**2 * (log(V4c/Vc) -1d0 + vc/v4c)
c     &         )
c      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(996,"('intv4c-1:', 2d13.5))")
c     &   1d0/(v3-v4) *( 1d0/(v3-v) * ( V3c**2* log(V3c/Vc)/(v3-v) - Vc))
c      write(996,"('intv4c-2:', 2d13.5))")
c     &   1d0/(v3-v4) *( 1d0/(v4-v) * ( V4c**2* log(V4c/Vc)/(v4-v) - Vc))
c      write(996,"('intv4c:', 2d13.5))")
c     &        V4c**2/(v4-v)**2/(v4-v3)*log(V4c/Vc)
c      write(996,"('intv4c:', 2d13.5))")
c     &        Vc/(v3-v)/(v4-v)
cccccccccccccccccccccccccccccccccccccccccccccdddddddd
      end
c----------------------------
      complex(8) function xfun(x)
      complex(8) :: x,y,a
      if(abs(x)<1d-1) then
        a = -x
        xfun = 0d0
        do n = 20,3,-1
          xfun = xfun*a +1d0/n !honar procedure
        enddo
        xfun = xfun*a
      else
        xfun= ( (1d0-x/2d0)*x - log(1d0+x) )/x**2
      endif
      end

c------------
      complex(8) function intv3c(omg,v,vd)
      implicit none
      complex(8):: omg,vc,vdc,intv3cc,x,xfun2
     &  ,intv3cx,intv3cy,xfun3
      real(8) :: v,vd
c      integer(4),save :: ic=0
      if(v ==0) v =1d-15
      if(vd==0) vd=1d-15
c-case2 numerically better
      vc  = V +omg
      vdc = Vd+omg
      x = vdc/vc -1d0
      intv3c = 3*vc*xfun2(x)/(v-vd)**2

c-case1
c      vc  = V +omg
c      vdc = Vd+omg
c      intv3cc = (2*vc*vdc/(v-vd)*log(vdc/vc) + (vc+vdc))
c     &          /(v-vd)**2
c      intv3c = 3*intv3cc

cccccccccccccccccccccccccccccccccccccccccccccccccc
c      ic=ic+1
c      intv3cy = 3*vc*xfun3(x)/(v-vd)**2
c      write(996,"('3',i7,2d12.4,2f10.5))") ic,omg,v,vd
c      write(996,"('3',i7,2d16.8))") ic,intv3c
c      write(996,"(8x,2d16.8))") intv3cx
c      write(996,"(8x,2d16.8))") intv3cy
c      write(996,*)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      end
c------------
c      complex(8) function xfun3(x)
c      complex(8) x
c      xfun3 =  x + 2 - 2*(1+x)*log(1+x)/ x
c      end

c------------
      complex(8) function xfun2(x)
      complex(8) x,a,xxx,xxx1,xxx2
      if(abs(x)<1d-1) then
        a = -x
        xfun2 = 0d0
        do n = 20,3,-1
          xfun2 = xfun2*a +1d0/n !honar procedure
        enddo
        xfun2 = (-xfun2* 2*(1+x) +1d0 )* x**2
      else
        xfun2 =  x + 2 - 2*(1+x)*log(1+x)/ x
      endif
      end

c------------
      complex(8) function intv2c(omg,v,v4)
      implicit none
      complex(8):: omg,vc,v4c,intv2cc
      real(8) :: v,v4
      if(v ==0) v =1d-15
      if(v4==0) v4=1d-15
      vc  = V +omg
      v4c = V4+omg
      intv2cc = v4c**2/(v-v4)**3* log(vc/v4c)
     &     + (1.5d0*v4c**2 + .5d0*vc**2 - 2*vc*v4c) /(v-v4)**3
      intv2c = 3*intv2cc
      end
c--------------------------------------
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
c--------------------------------------
      subroutine midk(kk,ee,xx,i,j,   kout,xout)
C- Calculate x and k(3) at the Fermi energy on the like k(i)---k(j).
      real(8) ::  kk(3,1:4),xx(1:4),ee(1:4), kout(3),xout
      ratio     = ee(i)/(ee(i)-ee(j))
      xout      = xx(i)     + ratio * (xx(j)-xx(i))
      kout(1:3) = kk(1:3,i) + ratio * (kk(1:3,j)-kk(1:3,i))
      end

c---------------------------------------------------------------
      subroutine rsvwwk(wgt, nqbz,nband,nctot,nbnbx,
     o   n1b,n2b,wwk,noccxv,nbnb)
C- Convert wgt to wwk
      implicit none
      integer(4) :: nband, nctot, nqbz, nbnbx ,ib,jb, kx,ix
      integer(4) :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),noccxv,nbnb(nqbz)
      complex(8) :: wgt(nband+nctot,nband,nqbz), wwk(nbnbx,nqbz)
cccccccccccccccccccccc
c cutoff for tetrahedron weight!
      real(8),parameter :: eps = 0d0
cccccccccccccccccccccc
      wwk    = 0d0
      noccxv = 0
      do kx  = 1, nqbz
        ix  = 0
        do ib  = 1, nband + nctot
          do jb  = 1, nband
            if( abs(wgt(ib,jb,kx)) >eps) then
              ix          = ix+1
              n1b(ix, kx) = ib
              n2b(ix, kx) = jb
              wwk(ix, kx) = wgt(ib,jb,kx)
              if( ib<=nband .and. ib>noccxv ) noccxv =ib
            endif
          enddo
        enddo
        nbnb(kx) = ix
      enddo
      write(6,*) ' rsvwwk: kx nbnbmax=',kx, maxval(nbnb)
      end

c---------------------------------------------------------------
      subroutine rsvwwk2(wgt, nqbz,nband,nctot,nbnbx, n1b0,n2b0,nbnb0,
     o   n1b,n2b,wwk,noccxv,nbnb)
C- Convert wgt to wwk
      implicit none
      integer(4) :: nband, nctot, nqbz, nbnbx ,ib,jb, kx,ix
      integer(4) :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),noccxv,nbnb(nqbz)
      complex(8) :: wgt(nband+nctot,nband,nqbz), wwk(nbnbx,nqbz),wtot
c
      integer(4) :: n1b0(nbnbx,nqbz),n2b0(nbnbx,nqbz),nbnb0(nqbz),ibib
      complex(8),allocatable ::wgt0(:,:,:)
      allocate(wgt0(nband+nctot,nband,nqbz))
      wgt0=wgt
      wwk = 0d0
      noccxv = 0
      do kx   = 1,nqbz
        nbnb(kx) = nbnb0(kx)
        ix=0
        do ibib = 1,nbnb0(kx)
          ix          = ix+1
          ib =  n1b0(ibib,kx)
          jb =  n2b0(ibib,kx)
          n1b(ix, kx) = ib
          n2b(ix, kx) = jb
          wwk(ix, kx) = wgt(ib,jb,kx)
          wgt0(ib,jb,kx) = 0d0
          if( ib<=nband .and. ib>noccxv ) noccxv =ib
        enddo
      enddo
      if( sum(abs(wgt0)) > 1d-10) then
        write(6,*)' rsvwwk2: sum(abs(wgt0))>0 the set n1b n2b at iw==1 is
     &    not enough to reserve wwk. 
     &    If you calculate from omega=0. it could be a problem.
     &    You might start rather large omega. 
     &    This problem is because inttetra3 could give complete zero 
     &    for small omega. Also the symmetrization could generate zero.'
     &    ,sum(abs(wgt0))
Cstop2rx 2013.08.09 kino        stop 'rsvwwk2: sum(abs(wgt0))'
        call rx( 'rsvwwk2: sum(abs(wgt0))')
      endif
      end

c---------------------------------------------------------------
      subroutine rsvwwk0(wgt, nqbz,nband,nctot,nbnbx,
     o   n1b,n2b,noccxv,nbnb)
C- get (n1b n2b) corresponding to non-zero wgt.
      implicit none
      integer(4) :: nband, nctot, nqbz, nbnbx ,ib,jb, kx,ix
      integer(4) :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),noccxv,nbnb(nqbz)
      complex(8) :: wgt(nband+nctot,nband,nqbz)
c
      noccxv = 0
      do kx  = 1, nqbz
        ix  = 0
        do ib  = 1, nband + nctot
          do jb  = 1, nband
            if( abs(wgt(ib,jb,kx)) /=0d0 ) then
              ix          = ix+1
              n1b(ix, kx) = ib
              n2b(ix, kx) = jb
              if( ib<=nband .and. ib>noccxv ) noccxv =ib
            endif
          enddo
        enddo
        nbnb(kx) = ix
      enddo
      write(6,*) ' rsvwwk: kx nbnbmax=',kx, maxval(nbnb)
      end
c---------------------------------------------------------------
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

c---------------------------------------------------------------
      subroutine rsvwwk00(iwgt, nqbz,nband,nctot,nbnbx,
     o   n1b,n2b,noccxv,nbnb)
C- get (n1b n2b) corresponding to non-zero wgt.
      implicit none
      integer(4) :: nband, nctot, nqbz, nbnbx ,ib,jb, kx,ix
      integer(4) :: n1b(nbnbx,nqbz),n2b(nbnbx,nqbz),noccxv,nbnb(nqbz)
      logical :: iwgt(nband+nctot,nband,nqbz)
      noccxv = 0
      do kx  = 1, nqbz
        ix  = 0
        do ib  = 1, nband + nctot
          do jb  = 1, nband
            if( iwgt(ib,jb,kx)) then
              ix          = ix+1
              n1b(ix, kx) = ib
              n2b(ix, kx) = jb
              if( ib<=nband .and. ib>noccxv ) noccxv =ib
            endif
          enddo
        enddo
        nbnb(kx) = ix
      enddo
      write(6,*) ' rsvwwk: kx nbnbmax=',kx, maxval(nbnb)
      end
c---------------------------------------------------------------


c---------------------------------------------------------------
      subroutine rsvwwk4(wgt, nqbz,nband,nctot,nbnbx, n1b0,n2b0,nbnb0,
     o  wwk)
C- Convert wgt to wwk
      implicit none
      integer(4) :: nband, nctot, nqbz, nbnbx ,ib,jb, kx,ix
      complex(8) :: wgt(nband+nctot,nband,nqbz), wwk(nbnbx,nqbz),wtot
      integer(4) :: n1b0(nbnbx,nqbz),n2b0(nbnbx,nqbz),nbnb0(nqbz),ibib
      complex(8),allocatable ::wgt0(:,:,:)
      allocate(wgt0(nband+nctot,nband,nqbz))
      wwk = 0d0
      wgt0= wgt
      do kx   = 1,nqbz
        ix = 0
        do ibib = 1,nbnb0(kx)
          ix = ix+1
          ib =  n1b0(ibib,kx)
          jb =  n2b0(ibib,kx)
          wwk(ix, kx) = wgt(ib,jb,kx)
          wgt0(ib,jb,kx) = 0d0
        enddo
      enddo
c      write(6,*)' rsvwwk4: sum(abs(wgt0))=',sum(abs(wgt0))
cccccccccccccccccccccccccccccccccc
c      write(6,*) 'xxxxxxxxxxxxxx nqbz nbnb0', nqbz,nbnb0
c      do kx   = 1,nqbz
c      do ibib = 1,nbnb0(kx)
c        ib =  n1b0(ibib,kx)
c        jb =  n2b0(ibib,kx)
c        write(6,"(' kx ib jb =',3i4,2d13.5)") kx,ib,jb,wgt0(ib,jb,kx)
c      enddo
c      enddo
ccccccccccccccccccccccccccccccccccc
      if( sum(abs(wgt0))/=0d0 ) then
        write(6,*)' rsvwwk4: sum(abs(wgt0))/=0d0 the set n1b n2b is',
     &          ' not enough to reserve wwk. ',sum(abs(wgt0))
        do kx = 1,nqbz
          do ib=1,nband+nctot
            do jb=1,nband
              if(wgt0(ib,jb,kx) /= 0d0) then
                write(6,"(' ib jb kx=',3i5,2d14.6)") ib,jb,kx,wgt0(ib,jb,kx)
              endif
            enddo
          enddo
        enddo
Cstop2rx 2013.08.09 kino        stop 'rsvwwk4: sum(abs(wgt0))'
        call rx( 'rsvwwk4: sum(abs(wgt0))')
      endif
      deallocate(wgt0)
      end

      subroutine rsvwwk4_w(wgt, nqbz,nband,nctot,nbnbx, n1b0,n2b0,nbnb0,
     i  ifwwk)
C- Convert wgt to wwk
      implicit none
      integer(4) :: nband, nctot, nqbz, nbnbx ,ib,jb, kx,ix ,ifwwk
      complex(8) :: wgt(nband+nctot,nband,nqbz), wtot
      integer(4) :: n1b0(nbnbx,nqbz),n2b0(nbnbx,nqbz),nbnb0(nqbz),ibib
      complex(8),allocatable ::wgt0(:,:,:),wwk(:,:)
c      allocate(wgt0(nband+nctot,nband,nqbz))
      allocate(wwk(nbnbx,nqbz))
      wwk = 0d0
c      wgt0= wgt
      do kx   = 1,nqbz
        ix = 0
        do ibib = 1,nbnb0(kx)
          ix = ix+1
          ib =  n1b0(ibib,kx)
          jb =  n2b0(ibib,kx)
          wwk(ix, kx) = wgt(ib,jb,kx)
c        wgt0(ib,jb,kx) = 0d0
        enddo
      enddo
c      write(6,*)' rsvwwk4: sum(abs(wgt0))=',sum(abs(wgt0))
c      if( sum(abs(wgt0))/=0d0 ) then
c        write(6,*)' rsvwwk4: sum(abs(wgt0))/=0d0 the set n1b n2b is
c     &    not enough to reserve wwk. ',sum(abs(wgt0))
c        stop 'rsvwwk4: sum(abs(wgt0))'
c      endif
c      deallocate(wgt0)
      write(ifwwk) wwk
      deallocate(wwk)
      end



















