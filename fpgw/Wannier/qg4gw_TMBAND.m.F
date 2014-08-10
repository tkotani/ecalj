      program qg4gw
C- Generate required q and G vectors for GW calculations.
C ----------------------------------------------------------------------
Cr input file
Cr   GWIN0: it is for GW calculation. See sample.
Cr   LATTC: contains these lattice informations;
Cr    alat       : lattice constant in a.u.
Cr    QpGcut_psi : maxmum of |q+G| in a.u. in the expansion of the eigenfunction. 
Cr    QpGcut_Cou : maxmum of |q+G| in a.u. in the expansion of the Coulomb matrix. 
Cr    plat(1:3,1): 1st primitive translation vector in the unit of alat
Cr    plat(1:3,2): 2nd primitive translation vector  
Cr    plat(1:3,3): 3rd primitive translation vector  
Cr   SYMOPS file : include point group operation. See sample.
Cr   GWIN_V2 is read in mkqg to see wether your GW run requires tetrahedron info or not.
Cr
Cr outtput files 
Cr   QGpsi: q and G vector for the eigenfunction
Cr   QGcou: q and G vector for the Coulomb matrix
Cr   Q0P  : q points near \Gamma points
Cr   QIBZ : q points in the Irreducible BZ.
Cr   Qmtet: q vectors for devided-tetrahedron.
Cr --------------------------
Cr The total number of q points are nqnum=(n1q*n2q*n3q)*(1+ nq0i).
Cr You can read q and G, for example, through
Cr
Cr     open(ifiqg, file='QGpsi',)
Cr     read(ifiqg ) nqnum, ngpmx, QpGcut_psi
Cr     allocate( q(3,nqnum),ngp(nqnum),ngvecp(3,ngpmx,nqnum) )
Cr     do  iq    = 1, nqnum
Cr       read (ifiqg) q(1:3), ngp(iq) 
Cr       read (ifiqg) ngvecp(1:3,1:ngp(iq),iq)
Cr     enddo
Cr     enddo
Cr     close(ifiqg)
Cr --------------------------
Cr True q is given by
Cr    True_q(1:3)     = 2*pi/alat * q(1:3)
Cr  True G is given by 
Cr    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvec(1:3,igp)) ,igp=1,ngp
Cr
Cr QpGcut_psi are given in both 
Cr  GWIN0 and LATTC. One of them is used and
Cr  other one is the dummy. The dummy should 
Cr  have negative value.
C ----------------------------------------------------------------------
      use keyvalue
      implicit none
      integer(4) ::n1q,n2q,n3q,ifiqg,ifiqgc,ifigw0,ngrp,ifi,i,ig,iq0pin,idummy
      real(8) :: plat(3,3), alat,QpGcut_psi, QpGcut_Cou,dummy
      real(8) :: volum,q0(3),qlat0(3,3),QpGx1,QpGx2,a1,a2,pi,unit
      real(8),allocatable :: symops(:,:,:)

      character(len=150):: recrdxxx
      character(len=10) :: keyw1='unit_2pioa',keyw2
      logical ::unit2=.false. !  readgwinput,
      integer(4)::nnn(3),ret
      integer(4):: version,verbose,q0pchoice,bzcase,wgtq0p,normcheck
      logical:: keepppovl,GaussSmear,KeepEigen,core_orth,ldummy

      integer(4):: iii
ccccccccccccccccccccccccccccccccc
c NaN test See NaNdble().
c      real(8) :: NaNdble,zzz
c      complex(8):: ccc,Nancmpx
c      real(4):: qqq
c      integer(4):: ixx
c      zzz = NaNdble()
c      qqq =NaNdble()
c      ccc = NaNcmpx()
c      ixx = zzz
c      write(6,*) zzz,zzz+1d0,ixx,qqq,ccc
c      stop 'xxxxxxxxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccc

      pi= 4d0* atan(1d0)
c --- mode
       write(6,*)
     & ' qg4gw: Generate Q0P->1; Readin Q0P ->2; band mode ->3;'
       read (5,*) iq0pin
      call headver('qg4gw',iq0pin)
       write(6,*) ' mode iq0pin = ',iq0pin

c--- If GWinput is not exist, GWIN0 GWIN_V2 and QPNT are converted to GWinput.
c      if(readgwinput()) then
         call conv2gwinput(iq0pin)
         if(iq0pin==-100) stop ' OK! qg4gw mode=-100 for GWinput'
c      else
c        if(iq0pin==-100) stop
c      endif
c-----------------------------------------------------------------

c ... for bug fix of ifc6. ifc6 does not allow 
c      write(6,"(' keepppovl() = ',L)")  KeepPPOVL() 
c     as first call is not allowed by ifc6.
c     init mode of KeepPPOVL() and so... Recursive procedure of read-write (i/0) seems to be not allowd by ifc6.
      ldummy=KeepPPOVL()
      ldummy=KeepEigen()
      ldummy=core_orth()
      idummy=verbose()
      ldummy=GaussSmear()
      idummy=bzcase()
      idummy=q0pchoice()
      idummy=normcheck()

c --- list switches
      print *,' --- switches defined in gwsrc/switch.f ---'
      write(6,"(' version()   = ',I3)")  version()
      write(6,"(' keepppovl() = ',L)")  KeepPPOVL()
      write(6,"(' keepeigen() = ',L)")  KeepEigen()
      write(6,"(' core_orth()  = ',L)")  core_orth()
      write(6,"(' verbose ()  = ',i4)")  verbose()
      write(6,"(' GaussSmear()= ',L4)")  GaussSmear()
      write(6,"(' q0pchoice() = ',i4)")  q0pchoice()
      write(6,"(' bzcase()    = ',i4)")  bzcase()
      write(6,"(' normcheck() = ',i4)")  normcheck()



C --- readin LATTC
c---------------------------------- 
c A sample LATTICE file for a Si case, 
c      10.26d0        ! alat        = lattice constant in a.u.
c       0d0 .5d0 .5d0 ! plat(1:3,1) = 1st primitive translation vector in the unit of alat
c      .5d0 .0d0 .5d0 ! plat(1:3,2)  = 2nd ...
c      .5d0 .5d0 .0d0 ! plat(1:3,3)  = 3rd ...
c      2.7d0          ! QpGcut_psi  = maxmum of |q+G| in a.u.
c---------------------------------- 
      ifi = 6661
      open (ifi, file='LATTC')
      read(ifi,*) alat
      read(ifi,*) plat(1:3,1)
      read(ifi,*) plat(1:3,2)
      read(ifi,*) plat(1:3,3)
      read(ifi,*) QpGx1 !If QpGcut_psi is given by the code (LAPW case)
      close(ifi)

C --- readin SYMOPS.
      ifi = 6661
      open (ifi, file='SYMOPS')
      read(ifi,*) ngrp
      write(6,*) ' SYMOPS ngrp=',ngrp
      allocate(symops(3,3,ngrp))
c      write(667,*) ngrp
      do ig = 1,ngrp
        read(ifi,*)  
c        write(667,*) ig
        do i=1,3
          read(ifi,*) symops(i,1:3,ig)
c          write(667,"(3d24.16)") symops(i,1:3,ig)
        enddo
      enddo
      close(ifi)

C --- check write
      write(6,"(' alat      =',f13.6 )") alat
      write(6,"(' plat a1   =',3f13.6)") plat(1:3,1)
      write(6,"(' plat a2   =',3f13.6)") plat(1:3,2)
      write(6,"(' plat a3   =',3f13.6)") plat(1:3,3)
      write(6,*) ' --- point group operation --- '
      do ig = 1,ngrp
        print *, ' ig=',ig
        do i=1,3
          write(6,"(3f14.6)") symops(i,1:3,ig) 
        enddo
      enddo

C --- Readin GWIN0. nx ny nz for k-points for GW calculations.
c      if(readgwinput()) then
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
c      else !old part.  from GWIN0
c       ifigw0 = 101
c       open(ifigw0,file='GWIN0',status="old",err=1010)
c       write(6,*) 'Reading GWIN0'
c       read(ifigw0,*)
c       read(ifigw0,*) n1q,n2q,n3q !GWIN0
c       read(ifigw0,*)
c
cc ...Jun2002
c       read(ifigw0,"(150a)") recrdxxx
c       read(recrdxxx,*) keyw2
c       if(keyw1==keyw2) then
c        read(recrdxxx,*) keyw2,a1,a2
c        unit = 2d0*pi/alat 
c        write(6,*) ' unit_2pioa find: 2pi/alat = ',unit
c        QpGx2     = unit*a1
c        QpGcut_Cou= unit*a2               ! cut off for G. |k+G|< QpGcut in (2pi/alat)
c       else
c        read(recrdxxx,*) QpGx2, QpGcut_Cou ! cut off for G. |k+G|< QpGcut in  (a.u.)
c       endif

cccccccccc
c      write(6,*) QpGx2, QpGcut_Cou  ! cut off for G. |k+G|< QpGcut in  (a.u.)
c      stop 'testend'
ccccccccccc
c       read(ifigw0,*)
c       read(ifigw0,*) alp
c       close(ifigw0)
c      endif
      if     ( QpGx1>=0d0 .and. QpGx2<0d0 ) then
        QpGcut_psi = QpGx1
        write(6,*) 'QpGcut_psi is given by GWIN0 =', QpGx1
       elseif(  QpGx1< 0d0 .and. QpGx2>=0d0) then
        QpGcut_psi = QpGx2
        write(6,*) 'QpGcut_psi is given by GWIN0 =', QpGx2
       else
        stop ' qg4gw: QpGcut_psi are given in both 
     &    GWIN0 and LATTC. One of them is used and
     &    other one is the dummy. The dummy should 
     &    have negative value like -1.'
      endif

      write(6,*) ' *** k points for GW from GWIN0 =',n1q,n2q,n3q
      write(6,"(' *** |k+G| < QpG(psi) QpG(Cou)=',2d13.6)")
     &             QpGcut_psi, QpGcut_Cou
      goto 1011

 1010 print *," --- a GWIN0 sample for the case of si ---" 
      write(6,"(a)")' n1 n2 n3. the number of k-points for GW '
      write(6,"(a)")' 2 2 2 '
      write(6,"(a)")' Cutoff of |q+G| for Coulomb. a.u.'
      write(6,"(a)")' 0 2.7 !QpGcut_psi(=0 -> dummy) QpGcut_Cou' 
      write(6,"(a)")
     &' Parameter to choose the offsetted gamma points Q0P.'
      write(6,"(a)")' 1d0 '
      write(6,"(a)")' number of bands.'
      write(6,"(a)")' 99999 '
      print *," --- end GWIN0 sample --- " 
      print *,
     & " ERROR! Not find GWIN0. You have to write GWIN0 by hand! 
     &          See above sample!"
      stop " ERROR! Not find GWIN0. You have to write GWIN0 by hand! "
 1011 continue

C --- Open output files QGpsi QGcou  
      ifiqg  = 401
      ifiqgc = 402
      open(ifiqg ,file='QGpsi',form='unformatted')
      open(ifiqgc,file='QGcou',form='unformatted')
!      call mkQG2(alat,plat,symops,ngrp,n1q,n2q,n3q,iq0pin,
!     &  QpGcut_psi, QpGcut_Cou, ifiqg, ifiqgc) 
      call mkQG2_TMBAND(alat,plat,symops,ngrp,n1q,n2q,n3q,iq0pin,
     &                   QpGcut_psi, QpGcut_Cou, ifiqg, ifiqgc)
      write(6,*) ' OK! End of qg4gw '
      if(iq0pin ==1) stop ' OK! qg4gw mode=1 normal mode'
      if(iq0pin ==2) stop ' OK! qg4gw mode=2 Readin Q0P mode'
      if(iq0pin ==3) stop ' OK! qg4gw mode=3 band-plot mode'
      end

