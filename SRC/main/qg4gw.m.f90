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
      use m_mpi,only: MPI__Initialize
      use m_lgunit,only: m_lgunit_init
      implicit none
      integer(4) :: ifiqg,ifiqgc,ifigw0,ngrp,ifi,i,ig,iq0pin
      real(8) :: alat,QpGcut_psi, QpGcut_Cou,dummy ,plat(3,3),volum,q0(3),qlat0(3,3),a1,a2,unit
      real(8),allocatable :: symops(:,:,:)
      character(len=150):: recrdxxx
c      character(len=10) :: keyw1='unit_2pioa',keyw2
      integer(4)::nnn(3),ret,verbose,q0pchoice,wgtq0p,iq0pinxxx ,ifile_handle,n1,n2,n3
      logical:: GaussSmear,KeepEigen,core_orth,ldummy, lnq0iadd=.false. !keepppovl,
      integer:: gammacellctrl=0
      logical:: lmagnon = .false., cmdopt2
      character(20):: outs=''
      call MPI__Initialize()
      call M_lgunit_init()
      call show_programinfo(6)
      call cputid (0)
      write(6,"(a)")'qg4gw: Generate Q0P->1; Readin Q0P->2; SW(chipm)->4'
      write(6,"(a)")'       Generate Q0P and Q0P for xyz ->201 '
      if( cmdopt2('--job=',outs) ) then 
         read(outs,*) iq0pin
      else   
         read (5,*) iq0pin
      endif   
      write(6,*) ' mode iq0pin = ',iq0pin
      if(iq0pin==1.or.iq0pin==2) then
         iq0pinxxx=iq0pin
      elseif(iq0pin==4) then
         iq0pinxxx=2
      elseif(iq0pin==201) then
         iq0pinxxx=1
         lnq0iadd=.true.
      elseif(iq0pin==10002) then
         iq0pinxxx=2
         gammacellctrl=1        !Gammacell skip mode
      elseif(iq0pin==20002) then
         iq0pinxxx=2
         gammacellctrl=2        !Gammacell only mode
      elseif(iq0pin==40001) then
         iq0pinxxx=2
         lmagnon=.true.         !! QforEPSL for magnon
      else
         call rx( 'Not allowed iq0pin')
      endif
      call mkQG2(iq0pinxxx, gammacellctrl,lnq0iadd,lmagnon)
      write(6,*) ' OK! End of qg4gw '
      if(iq0pin ==1)     call rx0( ' OK! qg4gw mode=1 normal mode')
      if(iq0pin ==2)     call rx0( ' OK! qg4gw mode=2 Readin Q0P mode')
      if(iq0pin ==4)     call rx0( ' OK! qg4gw mode=4 Readin Q0P mode. Set ngp=ngc=0')      
      if(iq0pin ==201)   call rx0( ' OK! qg4gw mode=201 Generate Q0P and Q0P for xyz ')
      if(iq0pin ==10002) call rx0( ' OK! qg4gw mode=10002 Readin Q0P. GammaCell skipped.')
      if(iq0pin ==20002) call rx0( ' OK! qg4gw mode=20002 Readin Q0P. GammaCell Only.')
      if(iq0pin ==40001) call rx0( ' OK! qg4gw mode=40001 QforEPSL for magnon')
      call rx('qg4gw: iq0pin wrong?')
      end

      include 'show_programinfo.fpp'
