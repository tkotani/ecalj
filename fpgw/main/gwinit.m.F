      program gwinit_v2
C- Generate required q and G vectors for GW calculations.
C ----------------------------------------------------------------------
Ci input file
Ci  (GWIN0): it is for GW calculation. See sample. If no GWIN0, GWIN0.tmp is generated.
Ci           Then you have to edit it and save it as GWIN0, and do gwinit once more.
Ci  LATTC: contains these lattice informations;
Ci    alat       : lattice constant in a.u.
Ci    QpGcut_psi : maxmum of |q+G| in a.u. in the expansion of the eigenfunction.
Ci    QpGcut_Cou : maxmum of |q+G| in a.u. in the expansion of the Coulomb matrix.
Ci    plat(1:3,1): 1st primitive translation vector in the unit of alat
Ci    plat(1:3,2): 2nd primitive translation vector
Ci    plat(1:3,3): 3rd primitive translation vector
Ci    QpGcut_psi  = maxmum of |q+G| in a.u.
Ci    nbas nclass
Ci    iclass(1:nbas)
Ci    lmax(1:nclass),konf(0:lmax,1:nclass)
Ci
Ci  SYMOPS : includes point group operation. See sample.
Ci
Ci  (SYML) : Symmetry line for BZ. See bandplot routines.
Ci           It it exists, it gives QPNTforSYML.tmp, which is suitable to plot bands for SYML.
Ci
Ci  NLAindx: This specify the order of cphi(1:mnla,iband). See rdata4gw_v2 cphi stored in CPHI.
Co output files
Co  GWIN_V2.tmp   : Template for GWIN_V2
Co  QPNT.tmp   : Template for QPNT, which specifies for which k-points
Co                      we will calculate the self-enery.
Co  (QPNTforSYML.tmp) : It is suitable to plot bands for SYML.
Co  KPNTin1BZ: k point in 1st BZ, it is a check write.
C
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
Cr takao kotani Apr 2002
C ----------------------------------------------------------------------
      implicit none
      integer(4) ::n1q,n2q,n3q,ifiqg,ifiqgc,ifigw0,ngrp,ifi,i,ig
      real(8) :: plat(3,3), alat,alp,QpGcut_psi, QpGcut_Cou,dummy
      real(8) :: volum,q0(3),qlat0(3,3),QpGx1,QpGx2,
     & dw,delta,deltaw,esmr,tolopt,qm(3,3)
      real(8),allocatable :: symops(:,:,:),qibz(:,:)

      integer(4) :: nbas,lmxax,ibas,l,ixxx,lmxa,ibasx,ifigw0t,mxkp,
     &  irs,niw,ic,iclass,
     &  nqibz,ifiqibz,iqibz,ifqpnt,iqall,iaf,iii, ifigwin,lk,
     &  nocc,nunocc,
     &  kkk,noccc,nunoccc,ncinc,ncinc2
      integer(4),allocatable :: lmxaa(:), konf(:,:),nncx(:,:),lcutmx(:)
      logical :: gwin0exist
      integer(4),allocatable:: IPQ(:)
      real(8),allocatable   :: WGT(:)
c
      real(8):: qlat(3,3),qp(3)
      integer(4):: nnn,ifkpt,i1,nlinemax=50,ifsyml,nline,iopenxx,nqs
      integer(4),allocatable:: nqq(:),symlon(:)
      real(8),allocatable:: qbz(:,:),qbzs(:,:),wbz(:),qq1(:,:),qq2(:,:),nstbz(:)
      logical :: extsyml
      integer(4)::checksymlon

      integer(4)::ifinla,izz,izzn,izzz,iatbk
      integer(4),allocatable::iat(:),lindx(:),nindx(:),nnvv(:,:)

      integer(4):: mnla, mnla_,idummy
      integer(4),parameter::       NLAindxmx=1000
      character(len=6)::  caption(NLAindxmx)
      character(len=100)::recordnla
      character(len=104)::recordnlax
      character(len=1)::seg2
      character(len=6)::seg1

      character(len=150):: recrdxxx
      character(len=10) :: keyw1='unit_2pioa',keyw2
      real(8)::a1,a2,unit,pi
      integer(4):: ibzcase

      pi= 4d0* atan(1d0)
C --- read ibzcase
      write(6,*) ' gwinit_v2: ibzcase?: Regular BZ mesh -->1; Without q=0 ->2'
      read (5,*)  ibzcase
      call headver('gwinit_v2',ibzcase)
      write(6,*) ' run mode ibzcase = ',ibzcase
c'

c -----------
      inquire(file='GWIN0',exist= gwin0exist)
      if(.not.gwin0exist) then
        ifigw0t = 102
        open(ifigw0t,file='GWIN0.tmp')
        write(6,*) " --- GWIN0.tmp is generated --- "
        write(ifigw0t,"(a)")' n1 n2 n3. the number of k-points for GW '
        write(ifigw0t,"(a)")'  2 2 2 '
        write(ifigw0t,"(a)")
     &      ' Cutoff of |q+G| for Psi and Coulomb. (a.u.)'
        write(ifigw0t,"(a)")'  3.2  2.7  ! QpGcut_psi QpGcut_Cou '
        write(ifigw0t,"(a)")
     &  ' Parameter to choose the offsetted gamma points Q0P.'
        write(ifigw0t,"(a)")'  1d0 '
        write(ifigw0t,"(a)")' number of bands.'
        write(ifigw0t,"(a)")'  99999 '
        write(6,*) " gwinit_v2: GWIN0.tmp is generated! "
        write(6,*) " You have to edit GWIN0.tmp "//  
     &    " and save it as GWIN0. Then repeat!"
        stop " gwinit_v2: Repeat with GWIN0! "
      endif

C --- Readin LATTC
c----------------------------------
c A sample LATTICE file for a Si case with empty sphere,
c  konf should contains minimum principle quantum number for valence electron.
c  In other words, konf(l)-1 for each l is taken as the maximum principle quantum number for core.
c-----------------------------------
c      10.26d0        ! alat        = lattice constant in a.u.
c       0d0 .5d0 .5d0 ! plat(1:3,1) = 1st primitive translation vector in the unit of alat
c      .5d0 .0d0 .5d0 ! plat(1:3,2)  = 2nd ...
c      .5d0 .5d0 .0d0 ! plat(1:3,3)  = 3rd ...
c      2.7d0          ! QpGcut_psi  = maxmum of |q+G| in a.u.
c      ----------------
c      4  4             ! nbas lmaxx
c      1  4   3 3 3 4 5 ! ibas lmax(ibas) konf(l=0:lmax,ibas)
c      2  4   3 3 3 4 5 !
c      3  4   1 2 3 4 5 !  empty spheres ---> nocore apparently.
c      4  4   1 2 3 4 5 !
cc----------------------------------
      ifi = 6661
      open (ifi, file='LATTC')
      read(ifi,*) alat
      read(ifi,*) plat(1:3,1)
      read(ifi,*) plat(1:3,2)
      read(ifi,*) plat(1:3,3)
      read(ifi,*) QpGx1        !If QpGcut_psi is given by the code (LAPW case)
      read(ifi,*)
      read(ifi,*) nbas,lmxax
      read(ifi,*)
      read(ifi,*)
      allocate(lmxaa(nbas), konf(0:lmxax,nbas),nncx(0:lmxax,nbas) )
      do ibas = 1,nbas
        read(ifi,*) ibasx, lmxaa(ibas), konf(0:lmxaa(ibas),ibas)
        write(6,"(100i4)") ibasx, lmxaa(ibas), konf(0:lmxaa(ibas),ibas)
        if(ibasx/=ibas       ) stop ' gwinit_v2: ibasx/=ibas'
        do l    = 0,lmxaa(ibas)
          nncx(l,ibas) = konf(l,ibas) -1 -l   ! number of cores for each l ibas
        enddo
      enddo
      close(ifi)

C --- Readin SYMOPS.
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

C --- Check write
      write(6,"(' alat      =',f13.6 )") alat
      write(6,"(' plat a1   =',3f13.6)") plat(1:3,1)
      write(6,"(' plat a2   =',3f13.6)") plat(1:3,2)
      write(6,"(' plat a3   =',3f13.6)") plat(1:3,3)
c      write(6,*) ' --- point group operation --- '
c      do ig = 1,ngrp
c        print *, ' ig=',ig
c        do i=1,3
c          write(6,"(3f14.6)") symops(i,1:3,ig)
c        enddo
c      enddo

C --- Readin GWIN0. nx ny nz for k-points for GW calculations.
      ifigw0 = 101
      open(ifigw0,file='GWIN0')
      write(6,*) 'Reading GWIN0'
      read(ifigw0,*)
      read(ifigw0,*) n1q,n2q,n3q !GWIN0
      read(ifigw0,*)

! ... Jun2002 it now recognize a keyword unit_2pioa.
!       read(ifigw0,*) QpGx2, QpGcut_Cou  ! cut off for G. |k+G|< QpGcut in  (a.u.)
!!! this is a false block not used now.
      if(.false.) then
        read(ifigw0,"(150a)") recrdxxx
        read(recrdxxx,*) keyw2
        if(keyw1==keyw2) then
          read(recrdxxx,*) keyw2,a1,a2
          unit = 2d0*pi/alat
          write(6,*) ' unit_2pioa find: 2pi/alat = ',unit
          QpGx2     = unit*a1
          QpGcut_Cou= unit*a2               ! cut off for G. |k+G|< QpGcut in (2pi/alat)
        else
          read(recrdxxx,*) QpGx2, QpGcut_Cou ! cut off for G. |k+G|< QpGcut in  (a.u.)
        endif
!
        read(ifigw0,*)
        read(ifigw0,*) alp
        if     ( QpGx1>=0d0 .and. QpGx2<0d0 ) then
          QpGcut_psi = QpGx1
          write(6,*) 'QpGcut_psi is given by GWIN0 =', QpGx1
        elseif(  QpGx1< 0d0 .and. QpGx2>=0d0) then
          QpGcut_psi = QpGx2
          write(6,*) 'QpGcut_psi is given by GWIN0 =', QpGx2
        else
          stop ' qg4gw: QpGcut_psi are given in both 
     &     GWIN0 and LATTC. One of them is used and
     &    other one is the dummy. The dummy should 
     &    have negative value.'
        endif
        close(ifigw0)
        write(6,"(' *** k points mesh for GW from GWIN0 =',3i4)") 
     &     n1q,n2q,n3q
        write(6,"(' *** |k+G| < QpG(psi) QpG(Cou)=',2d13.6)")
     &             QpGcut_psi, QpGcut_Cou
      endif

C --- Make q-points in IBZ.
      mxkp   = n1q*n2q*n3q
      allocate(qibz(3,mxkp),IPQ(mxkp),WGT(mxkp))
      call bzmesh(ibzcase,plat,qm,n1q,n2q,n3q,symops,ngrp,IPQ,qibz
     &  ,WGT,nqibz,mxkp)
      deallocate(IPQ,WGT)
c      ifiqibz = 6661
c      open (ifiqibz, file='QIBZ') !write q-points in IBZ.
c      write(ifiqibz,"(i10)") nqibz
c      do iqibz = 1,nqibz
c        write(ifiqibz,"(3d24.16)") qibz(1:3,iqibz)
c      enddo
c      close(ifiqibz)
c      write(6,*)
c      write(6,"( '  ngrp  = ',i3)")ngrp
c      write(6,'("  qibz=",i3,3f12.5)')(i,qibz(1:3,i),i=1,nqibz)

c --- readin SYML if it exists.
      inquire(file='SYML',exist=extsyml)
      print *,'--- Readin SYML'
      allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax)
     & ,symlon(mxkp))

      if(extsyml) then
        ifsyml = iopenxx('SYML')
        nline = 0
        do
          nline = nline + 1
          read(ifsyml,*,err=602,end=601) 
     &   nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
          write(6,"(' ndiv+1=',i4, ' qinit= ',3f10.6,' qend=',3f10.6)") 
     &   nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
        enddo
  602   stop ' gwinit_v2: Error readin SYML; You follow correct format?'
  601   continue
        close(ifsyml)
        nline = nline - 1
      endif

C --- q in 1BZ KPT1BZ.
      nnn = n1q*n2q*n3q
      allocate ( qbz(3,nnn),qbzs(3,nnn),wbz(nnn),nstbz(nnn) )
      call dinv33x (plat,qlat) !it was dinv33(plat,1,qlat) by Ferdi
      call genqbz (ibzcase,qlat,n1q,n2q,n3q,
     o  qbz,wbz,nstbz )
c ... write to file KPNTin1BZ
      ifkpt = 335
      nqs=0
      open(ifkpt,file='KPTin1BZ.gwinit.chk')
      do      i1 = 1,nnn
        call shorbz(qbz(1,i1),qp,qlat,plat)
        if(extsyml) then
          symlon(i1) 
     &    = checksymlon(qbz(1,i1),qq1,qq2,nline,qlat,plat,ifkpt)
          write (ifkpt,"(1x,i4,4f10.5,2x,3f10.5,i3)") 
     &     i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp,symlon(i1)
          if(symlon(i1)/=0) then
            nqs=nqs+1
            qbzs(:,nqs)=qbz(:,i1)
          endif
        else
          write (ifkpt,"(1x,i4,4f10.5,2x,3f10.5,i3)") 
     &     i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp
        endif
      end do
      close (ifkpt)
      deallocate(wbz)
      write(6,"(' --- TOTAL num of q is n1*n2*n3=',i10)")nnn

C --- Sample QPNT file ---------------
      ifqpnt = 888
      open (ifqpnt,file='QPNT.tmp')
      write(ifqpnt,"(a,a)") " --- Specify the q and band indeces", 
     & " for which we evaluate the self-energy ---"
      write(ifqpnt,*)
      write(ifqpnt,"(a)")
     & "*** all q -->1, otherwise 0;  up only -->1, otherwise 0"
      iqall = 0;      iaf   = 0
      write(ifqpnt,*)iqall,iaf
      write(ifqpnt,"(a)") 
     &  "*** no. states and band index for calculation."
      iii = 2
      write(ifqpnt,*)  iii ! nband
      write(ifqpnt,"(99i3)") (i,i=4,5)
      write(ifqpnt,"(a,a)") "*** q-points, which shoud be in qbz.", 
     & ",See KPNTin1BZ."
      write(ifqpnt,*) min(nqibz,3)
      write(ifqpnt,'(i3,3f23.16)')(i,qibz(1:3,i),i=1,nqibz)
      close(ifqpnt)
c
      if(extsyml) then
        ifqpnt = 888
        open (ifqpnt,file='QPNT.forSYML.tmp')
        write(ifqpnt,"(a)") " --- Specify the q and band indeces, 
     & for which we evaluate the self-energy ---"
        write(ifqpnt,*)
        write(ifqpnt,"(a)")
     & "*** all q -->1, otherwise 0;  up only -->1, otherwise 0"
        iqall = 0;      iaf   = 0
        write(ifqpnt,*)iqall,iaf
        write(ifqpnt,"(a)") 
     &  "*** no. states and band index for calculation."
        iii = 2
        write(ifqpnt,*)  iii ! nband
        write(ifqpnt,"(99i3)") (i,i=4,5)
        write(ifqpnt,"(a)") "*** q-points, which shoud be in qbz. 
     & See KPNTin1BZ."
        write(ifqpnt,*) nqs
        write(ifqpnt,'(i3,3f23.16)')(i,qbzs(1:3,i),i=1,nqs)
        close(ifqpnt)
      endif

c readin NLAindx.
      ifinla = 889
      open(ifinla,file='NLAindx')
      write(6,
     & "(' --- readin NLAindx---[specify the order of cphi(1:mnla)]')")
      write(6,"(' -- n    l    a -- ')")
      read(ifinla,*)
      read(ifinla,*) mnla
      allocate(nindx(mnla),lindx(mnla),iat(mnla),nnvv(0:lmxax,nbas))
      nnvv  = 0
      mnla_ = 0
      iatbk = 0
      do izz=1, mnla
        read(ifinla,"(a100)",err=1039,end=1039) recordnla
c       write(6,"(a100)")recordnla
        recordnlax = recordnla//' $$$'
        read(recordnlax,*,err=1040,end=1040) 
     &   nindx(izz),lindx(izz),iat(izz),idummy,caption(izz)
        write(6,"(4i5,2x,a)") 
     &   nindx(izz),lindx(izz),iat(izz),idummy,caption(izz)
        izzn = izz
        mnla_= mnla_ + 2*lindx(izz)+1
        if(nnvv(lindx(izz),iat(izz))< nindx(izz)) then
          nnvv(lindx(izz),iat(izz)) = nindx(izz)
        endif
      enddo
 1039 continue
      write(6,*)' --- end of readin NLAindx--- nla mnla=',izzn,mnla
      if(mnla /=mnla_) stop ' gwinit_v2:mnla/= mnla_ '
      goto 1041
 1040 continue
      stop 'gwinit: can not readin recordnla'
 1041 continue

C --- A sample GWIN_V2 file -----------------------------------------------
c You should change some points before GW calculation!
c Especially, FLEQUENIES and PRODUCT BASIS
c
c Required inputs lmxa(1:ibas)       : lmax for argumantation.
c                 konf(lmxa(ic),ic)  : number of cores
c--------------------------
      ifigwin = 888
      dw = 0.005d0  !dw =0.02
      niw = 10      !niw =10
      delta  = -1d-6
      deltaw = 0.02d0
      esmr   = 0.003d0 ! 0.01d0

      open (ifigwin,file='GWIN_V2.tmp')
c ... Frequencies
      write(ifigwin,*) 
     & "FREQUENCIES----------------------------------------------------"
      write(ifigwin,"(a)") 
     & " Mesh size along Re axis: "//
     & ' !  (Num. of mesh is in NW given by hchknw.)'
      write(ifigwin,"(f13.5,a)") dw
     &,'  ! =dw(a.u.)  Try at least dw=0.04 and dw=0.02 '

      write(ifigwin,"(a)") 
     & " Number of frequencies along Im axis. "//
     & "  Used for integration to get SEC* in hsfp0"
      write(ifigwin,"(i10,a)") niw,
     & "     ! =niw  Try at least niw=6,10, or 12 "

      write(ifigwin,"(a)") 
     & ' Broadening of x0. negative means tetrahedron method;'
      write(ifigwin,"(d14.5,a)")delta,' ! =delta(a.u.) ;used by hx0fp0.'
     &  //' You get smeard x0 with abs(delta). '

      write(ifigwin,"(a)") 
     & ' Mesh for numerical derivative to get the Z factor'
      write(ifigwin,"(d14.5,a)")deltaw,
     &  ' ! =deltaw(a.u.); used by hsfp0.'
     & //' Check errors on Z--->See two Z values shown in SEC*'

      write(ifigwin,"(a)") 
     & ' Poles of G^LDA are treated as if they have width '//
     & 'esmr in hsfp0. Change esmr for metals. '//
     & ' See DOSACC*---especailly around Ef.'
      write(ifigwin,"(d14.5,a)") esmr,' GaussSmear ! =esmr(Ry)  used by hsfp0.'//
     &' Keep esmr as it is (smaller than band gap) for insulator'

      write(ifigwin,"(a)") "   ------------ dummy line ----------- "
      write(ifigwin,"(a)") "   ------------ dummy line ----------- "

c ... PRODUCT BASIS -------------------------
      write(ifigwin,*)
      write(ifigwin,"(a)") 
     & "PRODUCT BASIS--------------------------------------------------"
      write(ifigwin,"(a)") 
     & " tolerance to remove products due to poor linear-independency"
      tolopt =1d-2
      write(ifigwin,"(d14.6,a)") tolopt, 
     & " ! =tolopt; larger gives smaller num. of product basis."//
     & " See lbas and lbasC, which are output of hbasfp0."

      write(ifigwin,"(a)") 
     & " lcutmx(atom) = maximum l-cutoff for the product basis. "
     & //" =4 is required for atoms with valence d, like Ni Ga"
      allocate(lcutmx(nbas)); lcutmx=4
      write(ifigwin,"(1000i3)") lcutmx(1:nbas)

      write(ifigwin,"(a)")  "  atom   l  nnvv  nnc "
     &  //"! nnvv: num. of radial functions (valence) on the "//
     &   "augmentation-waves, nnc: num. for core."
      do ibas =1,nbas
        do lk   =0,lmxaa(ibas)
          write(ifigwin,"(4i5)") ibas,lk, nnvv(lk,ibas), nncx(lk,ibas)
        enddo
      enddo
c
      write(ifigwin,"(a)")
     & "  atom   l    n  occ unocc  ! Valence(1=yes,0=no) "
c     & //" (n=1->phi n=2->phidot n=3->LocalOrbital)"
      do ibas= 1, nbas
        do lk   = 0, lmxaa(ibas)
          do izz = 1, izzn
            if(iat(izz)==ibas.and.lk==lindx(izz)) then
              nocc = 1; nunocc = 1
              if(lindx(izz) >2 ) nocc   = 0
              if(lindx(izz) >3 ) nunocc = 0
              if(nindx(izz)==2) then
                nocc   = 0
                nunocc = 0
              endif
              seg1='';if(iat(izz)/=iatbk) seg1='----- '
              iatbk=iat(izz)
              write(ifigwin,"(5i5,3x,a )") iat(izz),lindx(izz),nindx(izz)
     &                        , nocc,nunocc, '! '//caption(izz)//seg1
            endif
          enddo
        enddo
      enddo
      write(ifigwin,"(a)") 
     &'  atom   l    n  occ unocc  '//
     &'ForX0 ForSxc ! Core (1=yes, 0=no)'
      do ibas  = 1,nbas
        do lk    = 0,lmxaa(ibas)
          if(lk==0) seg2='S'
          if(lk==1) seg2='P'
          if(lk==2) seg2='D'
          if(lk==3) seg2='F'
          if(lk==4) seg2='G'
          if(lk==5) seg2='H'
          if(lk==6) seg2='I'
          do kkk   = lk+1,konf(lk,ibas)-1
            noccc=0; nunoccc=0; ncinc=0; ncinc2=0
            seg1='';if(lk==0.and.kkk==lk+1) seg1=' -----'
            write(ifigwin,"(5i5,2x,2i5,a)") 
     &   ibas,lk,kkk-lk,noccc,nunoccc,ncinc,ncinc2
     &    ,'    ! '//char(48+kkk)//seg2//seg1
          enddo
        enddo
      enddo
      if(extsyml) then
        print *, 
     & ' OK! gwinit_v2: generated GWIN_V2.tmp QPNT.tmp QPNTforSYML.tmp.'
      else
        print *, 
     & ' OK! gwinit_v2: generated GWIN_V2.tmp QPNT.tmp. '
      endif
      stop ' Write GWIN_V2 and QPNT based on these template! '
      end





