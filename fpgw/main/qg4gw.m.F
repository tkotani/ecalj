      program qg4gw
!> Generate required q+G vectors and so on for GW calculations.
!! input file
!!   LATTC: contains these lattice informations;
!!    alat       : lattice constant in a.u.
!!    QpGcut_psi : maxmum of |q+G| in a.u. in the expansion of the eigenfunction.
!!    QpGcut_Cou : maxmum of |q+G| in a.u. in the expansion of the Coulomb matrix.
!!    plat(1:3,1): 1st primitive translation vector in the unit of alat
!!    plat(1:3,2): 2nd primitive translation vector
!!    plat(1:3,3): 3rd primitive translation vector
!!   SYMOPS file : include point group operation. See sample.
!!
!! outtput files:
!!   QGpsi: q and G vector for the eigenfunction
!!   QGcou: q and G vector for the Coulomb matrix
!!   Q0P  : offset Gamma point around \Gamma points
!!   EPSwklm : offset Gamma method.
!! and so on.
!!ccc   Qmtet: q vectors for devided-tetrahedron.
!! --------------------------
!! For exampl,e QGpsi is written in the following manner. See mkqg2 in mkqg.F
!!     open(ifiqg, file='QGpsi',)
!!      write(ifiqg ) nqnum,ngpmx,QpGcut_psi,nqbz,nqi,imx,nqibz
!!      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx) ) !inverse mapping table
!!      ngvecprev=9999
!!      ngveccrev=9999
!!      do iq = 1, nqnum
!!         q = qq(1:3,iq)
!!         write (ifiqg) q, ngp, irr(iq) ! irr=1 for irreducible points
!!         do ig = 1,ngp
!!            nnn3 = ngvecp(1:3, ig) 
!!            ngvecprev( nnn3(1), nnn3(2),nnn3(3)) = ig
!!         enddo
!!         write (ifiqg)  ngvecp,ngvecprev !ngvecprev is added on mar2012takao
!!         do ig = 1,ngc
!!            nnn3 = ngvecc(1:3, ig) 
!!            ngveccrev( nnn3(1), nnn3(2),nnn3(3)) = ig
!!         enddo
!!       enddo  
!!     close(ifiqg)
!! -----------------------------------------------------
!! True q (in a.u. in Cartesian coordinate) is given by
!!    q(1:3)     = 2*pi/alat * q(1:3)
!! True q+G is given by
!!    qplusG(1:3,igp) = 2*pi/alat * (q + matmul(qlat * ngvec(1:3,igp))), for igp=1,ngp
!! ----------------------------------------------------------------------
      use m_keyvalue,only: getkeyvalue
      implicit none
      integer(4) ::n1q,n2q,n3q,ifiqg,ifiqgc,ifigw0,ngrp,ifi,i,ig,iq0pin,idummy
      real(8) ::  alat,QpGcut_psi, QpGcut_Cou,dummy ,plat(3,3)
      real(8) :: volum,q0(3),qlat0(3,3),QpGx2,a1,a2,pi,unit !,QpGx1
      real(8),allocatable :: symops(:,:,:)
      character(len=150):: recrdxxx
      character(len=10) :: keyw1='unit_2pioa',keyw2
      logical ::unit2=.false. !  readgwinput,
      integer(4)::nnn(3),ret
      integer(4):: verbose,q0pchoice,wgtq0p    !,normcheck !version,
      logical:: GaussSmear,KeepEigen,core_orth,ldummy, lnq0iadd=.false. !keepppovl,
      integer(4):: iq0pinxxx ,ifile_handle,n1,n2,n3
      integer:: gammacellctrl=0
      pi= 4d0* atan(1d0)
      call cputid (0)
      write(6,*)' qg4gw: Generate Q0P->1; Readin Q0P->2; band mode->3; SW(chipm)->4'
      write(6,*)'        Generate Q0P->101(old offset Gamma)'
      write(6,*)'        Generate Q0P and Q0P for xyz ->201 '
      read (5,*) iq0pin
      call headver('qg4gw',iq0pin)
      write(6,*) ' mode iq0pin = ',iq0pin
      if(iq0pin==-100.or.iq0pin==1.or.iq0pin==2.or.iq0pin==3.or.iq0pin==101) then
        iq0pinxxx=iq0pin
      elseif(iq0pin==10002) then
         iq0pinxxx=2
         gammacellctrl=1 !Gammacell skip mode
      elseif(iq0pin==20002) then
         iq0pinxxx=2
         gammacellctrl=2 !Gammacell only mode
      elseif(iq0pin==4) then
         iq0pinxxx=2
      elseif(iq0pin==201) then
        iq0pinxxx=1
        lnq0iadd=.true.
      else
        call rx( 'Not allowed iq0pin')
      endif
c this is moved to gwinit.m.F march2016
c!! Generate templeta of GWinput for iq0pin=-100
c      if(iq0pin==-100) then
c         call conv2gwinput()
c         call rx0( ' OK! qg4gw mode=-100 to generate GWinput')
c      endif   
      idummy=q0pchoice()
      write(6,"(' q0pchoice() = ',i4)")  q0pchoice()

      ifi=ifile_handle()
      open (ifi, file='LATTC')
      read(ifi,*) alat
      read(ifi,*) plat(1:3,1)
      read(ifi,*) plat(1:3,2)
      read(ifi,*) plat(1:3,3)
      read(ifi,*) !dummy
      close(ifi)
!! --- readin SYMOPS. point group operations. r'=matmul(symops(:,:),r) for any ig.
      ifi=ifile_handle()
      open (ifi, file='SYMOPS')
      read(ifi,*) ngrp
      write(6,*) ' SYMOPS ngrp=',ngrp
      allocate(symops(3,3,ngrp))
      do ig = 1,ngrp
        read(ifi,*)
        do i=1,3
          read(ifi,*) symops(i,1:3,ig)
        enddo
      enddo
      close(ifi)
!! --- check write
      write(6,*) ' --- primitive vectors ---'
      write(6,"(' unit(a.u.) alat  =',f13.6 )") alat
      write(6,"(' primitive_1 =',3f13.6)") plat(1:3,1)
      write(6,"(' primitive_2 =',3f13.6)") plat(1:3,2)
      write(6,"(' primitive_3 =',3f13.6)") plat(1:3,3)
      write(6,*) ' --- point group operations --- '
      do ig = 1,ngrp
        print *, ' ig=',ig
        do i=1,3
          write(6,"(3f14.6)") symops(i,1:3,ig)
        enddo
      enddo
!! --- Readin GWinput
      call getkeyvalue("GWinput", "n1n2n3", nnn,3)
      n1q=nnn(1); n2q=nnn(2); n3q = nnn(3)
      call getkeyvalue("GWinput", "QpGcut_psi",QpGx2)
      call getkeyvalue("GWinput", "QpGcut_cou",QpGcut_Cou)
      call getkeyvalue("GWinput", "unit_2pioa",unit2)
      if(unit2) then
        unit = 2d0*pi/alat
        QpGx2     = QpGx2      *unit
        QpGcut_cou= QpGcut_cou *unit
      endif
      QpGcut_psi = QpGx2
      write(6,"(' --- k points for GW from GWinput =',3i3)") nnn(1:3) 
      write(6,"(' ---  |k+G| < QpG(psi) QpG(Cou)=',2d13.6)") QpGcut_psi, QpGcut_Cou
      ifiqg  = ifile_handle()
      open(ifiqg ,file='QGpsi',form='unformatted')
      ifiqgc = ifile_handle()
      open(ifiqgc,file='QGcou',form='unformatted')
      if(iq0pin==4) then
         QpGcut_psi=0d0
         QpGcut_Cou=0d0
      endif
!!
      call mkQG2(alat,plat,symops,ngrp,nnn,iq0pinxxx, 
     &     QpGcut_psi, QpGcut_Cou, ifiqg, ifiqgc, gammacellctrl,lnq0iadd)
      write(6,*) ' OK! End of qg4gw '
      if(iq0pin ==1)     call rx0( ' OK! qg4gw mode=1 normal mode')
      if(iq0pin ==2)     call rx0( ' OK! qg4gw mode=2 Readin Q0P mode')
      if(iq0pin ==10002) call rx0( ' OK! qg4gw mode=10002 Readin Q0P. GammaCell skipped.')
      if(iq0pin ==20002) call rx0( ' OK! qg4gw mode=20002 Readin Q0P. GammaCell Only.')
      if(iq0pin ==3)     call rx0( ' OK! qg4gw mode=3 band-plot mode')
      if(iq0pin ==4)     call rx0( ' OK! qg4gw mode=4 Readin Q0P mode. Set ngp=ngc=0')      
      end
