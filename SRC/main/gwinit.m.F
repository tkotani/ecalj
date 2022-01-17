      program gwinit_v2
!!  Generate GWinput.tmp
C ----------------------------------------------------------------------
Ci input file
Ci  GWIN0: n1,n2,n3
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
Ci  NLAindx: This specify the order of cphi(1:mnla,iband). See rdata4gw_v2 cphi stored in CPHI.
Co Output files
Co  GWIN_V2.tmp : Template for GWIN_V2
Co  QPNT.tmp    : Template for QPNT, which specifies for which k-points
Co                      we will calculate the self-enery.
Co  KPNTin1BZ: k point in 1st BZ, it is for check.
C
Cr --------------------------
Cr The total number of q points are nqnum=(n1q*n2q*n3q)*(1+ nq0i).
Cr You can read q and G, for example, through
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
!! -----------------------------------------------------------
      use m_hamindex,only:   Readhamindex, symops, ngrp
      use m_get_bzdata1,only:  getbzdata1, 
     &  nqbz, nqibz, qbz,wbz,qibz,wibz
      implicit none
      integer ::n1q,n2q,n3q,ifiqg,ifiqgc,ifigw0,ifi,i,ig,nsp
      real(8) :: plat(3,3), alat,alp,QpGcut_psi, QpGcut_Cou,dummy
      real(8) :: volum,q0(3),qlat0(3,3),QpGx1,QpGx2,
     & dw,delta,deltaw,esmr,tolopt,qm(3,3)
c      real(8),allocatable :: symops(:,:,:) !,qibz(:,:)

      integer :: nbas,lmxax,ibas,l,ixxx,lmxa,ibasx,ifigw0t,mxkp,
     &  irs,niw,ic,iclass,
     &  ifiqibz,iqibz,ifqpnt,iqall,iaf,iii, ifigwinv2,lk, !nqibz,
     &  nocc,nunocc,
     &  kkk,noccc,nunoccc,ncinc,ncinc2
      integer,allocatable :: lmxaa(:), konf(:,:),nncx(:,:),lcutmx(:)
      logical :: gwin0exist
      integer,allocatable:: IPQ(:)
      real(8),allocatable   :: WGT(:)
      real(8):: qlat(3,3),qp(3)
      integer:: nnn,ifkpt,i1,nlinemax=50,ifsyml,nline,iopenxx,nqs
      integer,allocatable:: nqq(:),symlon(:)!,nstbz(:)
      real(8),allocatable:: qbzs(:,:),qq1(:,:),qq2(:,:)  !qbz(:,:),wbz(:),
      logical :: extsyml
      integer::checksymlon
      integer::ifinla,izz,izzn,izzz,iatbk
      integer,allocatable::iat(:),lindx(:),nindx(:),nnvv(:,:)
      integer:: mnla, mnla_,idummy
      integer,parameter::       NLAindxmx=1000
      character(len=6)::  caption(NLAindxmx)
      character(len=100)::recordnla
      character(len=104)::recordnlax
      character(len=1)::seg2
      character(len=6)::seg1
      character(len=150):: recrdxxx
      character(len=10) :: keyw1='unit_2pioa',keyw2
      real(8)::a1,a2,unit
      integer:: ibzcase,ifile_handle
!! conv2gwinput
      integer :: natom
      character(8),allocatable::spid(:)
      character(len=550):: pppx
      character(len=10) :: add,esmr_char
      logical :: gwinputexist,qpntexist,anfexist
      integer(4):: ifig,ifigwinp,nband_chi0,nband_sigm
      integer:: ifatomlist,ifiq
      integer(4):: nband_sigm2, iatom, nwin, incwfin
      real,parameter:: pi= 4d0* atan(1d0)
      character*255:: readin
!! Read ibzcase
      write(6,"(a)")' === gwinit: ibzcase?: Regular BZ mesh -->1; Without q=0 ->2 + option(n1q,n2q,n3q)'
      read (5,"(a)")  readin
      read(readin,*) ibzcase
c      call headver('gwinit',ibzcase)
      read(readin,*,end=101) idummy,n1q,n2q,n3q
      goto 102
 101  continue
      n1q=4 !defaultvalues 
      n2q=4
      n3q=4
 102  continue
      write(6,"(a,i3)")  ' Run mode ibzcase   = ',ibzcase
      write(6,"(a,3i5)") ' Readin n1q n2q n3q = ',n1q,n2q,n3q
!! Readin LATTC
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
      open (newunit=ifi, file='LATTC',form='unformatted')
      read(ifi) alat,plat,nbas,lmxax,nsp
      allocate(lmxaa(nbas), konf(0:lmxax,nbas),nncx(0:lmxax,nbas) )
      read(ifi)konf,lmxaa
      close(ifi)
      do ibas = 1,nbas
        write(6,"(100i4)") ibas, lmxaa(ibas), konf(0:lmxaa(ibas),ibas)
        do l= 0,lmxaa(ibas)
          nncx(l,ibas) = konf(l,ibas) -1 -l   ! number of cores for each l ibas
        enddo
      enddo
      close(ifi)
      call readhamindex()
      write(6,"(' alat      =',f13.6 )") alat
      write(6,"(' plat a1   =',3f13.6)") plat(1:3,1)
      write(6,"(' plat a2   =',3f13.6)") plat(1:3,2)
      write(6,"(' plat a3   =',3f13.6)") plat(1:3,3)
C --- Make q-points in IBZ.
      mxkp   = n1q*n2q*n3q
      call minv33tp(plat,qlat) !it was dinv33(plat,1,qlat) by Ferdi
      call getbzdata1(qlat,(/n1q,n2q,n3q/), !n1q,n2q,n3q
     &    symops,ngrp,tetrai=.false.,mtet=(/1,1,1/),gammacellctrl=0) 

!! Write to file KPNTin1BZ
      nnn = n1q*n2q*n3q
      ifkpt = ifile_handle()
      nqs=0
      open(ifkpt,file='KPTin1BZ.gwinit.chk')
      do      i1 = 1,nnn
        call shorbz(qbz(1,i1),qp,qlat,plat)
c        if(extsyml) then
c          symlon(i1) 
c     &    = checksymlon(qbz(1,i1),qq1,qq2,nline,qlat,plat,ifkpt)
c          write (ifkpt,"(1x,i4,4f10.5,2x,3f10.5,i3)") 
c     &     i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp,symlon(i1)
c          if(symlon(i1)/=0) then
c            nqs=nqs+1
c            qbzs(:,nqs)=qbz(:,i1)
c          endif
c        else
          write (ifkpt,"(1x,i4,4f10.5,2x,3f10.5,i3)") 
     &     i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp
c        endif
      end do
      close (ifkpt)
c      deallocate(wbz)
      write(6,"(' --- TOTAL num of q is n1*n2*n3=',i10)")nnn

C --- Sample QPNT file ---------------
      ifqpnt = ifile_handle()
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
      rewind ifqpnt
c      close(ifqpnt)
c
c$$$      if(extsyml) then
c$$$        ifqpnt = 888
c$$$        open (ifqpnt,file='QPNT.forSYML.tmp')
c$$$        write(ifqpnt,"(a)") " --- Specify the q and band indeces, 
c$$$     & for which we evaluate the self-energy ---"
c$$$        write(ifqpnt,*)
c$$$        write(ifqpnt,"(a)")
c$$$     & "*** all q -->1, otherwise 0;  up only -->1, otherwise 0"
c$$$        iqall = 0;      iaf   = 0
c$$$        write(ifqpnt,*)iqall,iaf
c$$$        write(ifqpnt,"(a)") 
c$$$     &  "*** no. states and band index for calculation."
c$$$        iii = 2
c$$$        write(ifqpnt,*)  iii ! nband
c$$$        write(ifqpnt,"(99i3)") (i,i=4,5)
c$$$        write(ifqpnt,"(a)") "*** q-points, which shoud be in qbz. 
c$$$     & See KPNTin1BZ."
c$$$        write(ifqpnt,*) nqs
c$$$        write(ifqpnt,'(i3,3f23.16)')(i,qbzs(1:3,i),i=1,nqs)
c$$$        close(ifqpnt)
c$$$      endif

!! Readin NLAindx.
      ifinla = ifile_handle()
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
      write(6,*)' --- end of readin NLAindx--- nla mnla=',izzn,mnla,mnla_
      if(mnla /=mnla_) call rx( ' gwinit_v2:mnla/= mnla_ ')
      goto 1041
 1040 continue
      call rx( 'gwinit: can not readin recordnla')
 1041 continue
      write(6,*)'gwinit: 1041'
!! Write GWinput.tmp
      ifigwinp= ifile_handle()
      open(ifigwinp,file='GWinput.tmp')
      ifi = ifigwinp
      write(ifi,"(a)")'!!! Starting from ! (or nonkeyword) is comment line !!! '
      write(ifi,"(a)")'!!! Each lines consists of "keyword value(s)"  !!! '
      write(ifi,"(a)")'!!! Each tag section in <...>... </...> has its own format. !!! '
      write(ifi,"(a)")
      write(ifi,"(a)")'!EIBZmode off  !no symmetrization for hx0fp0* (default on);Only affects comp. effort. off may faster.'
      write(ifi,"(a)")'!chi_RegQbz off !Use no Gamma mesh for dielectric function. This automaticall set EIBZmode off.'
      write(ifi,"(a)")'!Verbose    0  ! 0-->default; 100--->debug '
      write(ifi,"(a)")'!LFC@Gamma off !(on is default) if on, eps with Local field correction is used at Gamma point'
      write(ifi,"(a)")'!Q0Pchoice 1 !1(default):qzerolimit(in practice, See generated Q0P), 2:1/q^2 average in Gamma region'
c      write(ifi,"(a)")'!CoreOrth off ! off  --> Not enforce core to be orthogonalized with valence (default)'
c      write(ifi,"(a)")'               ! on   --> Enforce cores orthogonalied to valence (this may give strange core functions!).'
c      write(ifi,"(a)")'!multitet 2 2 2 ! tetrahedron is devided into micro tetrahedron. This probably now work'
      write(ifi,"(a)" )'! ##### From GWIN0 ################ '
      write(ifi,"(a,3i4,a)" )'n1n2n3',n1q,n2q,n3q,' ! for BZ meshing in GW, Wannier function and cRPA'
      write(ifi,"(a)") 'QpGcut_psi 4.0  !(See unit_2pioa for unit) |q+G| cutoff for eigenfunction.'
      write(ifi,"(a)") 'QpGcut_cou 3.0  !(See unit_2pioa for unit) |q+G| cutoff for Coulomb and W.'
      write(ifi,"(a)") 'unit_2pioa off ! off --> a.u.; on--> unit of QpGcut_* are in 2*pi/alat '
      write(ifi,"(a)") 'alpha_OffG 1.0 !(a.u.) Used in auxially function in the offset-Gamma method.'
c      write(ifi,"(a,i8,a)")     'nband_chi0 ', nband_chi0,' !    nband cutoff for chi0  (Optional)'
c      write(ifi,"(a,2i8,a)")     'nband_sigm ', nband_sigm
c     &      ,nband_sigm2,' !    nband cutoff for Sigma  (Optional) (1st:num in sigma; 2nd: num of G used in hsfp0)'
c      write(ifi,"(a,2f10.3,a)")  'emax_sigm  ', emax_sigm ,
c     &       emax_sigm2,'  !(Ry)  (Optional) emax cutoff for Sigma (as in the nband_sigm)'
c      write(ifi,"(a,2i8,a)")     'nband_sigm ', nband_sigm
c     &      ,nband_sigm2,' !    nband cutoff for Sigma  (Optional) (1st:num in sigma; 2nd: num of G used in hsfp0)'
c     &       emax_sigm2,'  !(Ry)  (Optional) emax cutoff for Sigma (as in the nband_sigm)'

      write(ifi,"(a)") '!emax_chi0  999 !(Ry) emax cutoff for chi0  (Optional)'
      write(ifi,"(a)") 'emax_sigm  3.0  !(Ry)  emax cutoff for Sigma'
!!
      write(ifi,*)
      write(ifi,"(a)" ) '! ##### FREQUENCIES from GWIN_V2 ################ '
c      write(ifi,"(a,f10.6,a)") 'dw    ', dw,   ' !(a.u.) energy-mesh (bin width size) along real axis.'
c      write(ifi,"(a,f10.3,a)") 'omg_c ', omg_c,' !(a.u.) energy-mesh is twiced at omg_c'
      write(ifi,"(a)")  'HistBin_dw    2d-3 ! 1d-5 is fine mesh (good for metal?) !(a.u.) BinWidth along real axis at omega=0.'
      write(ifi,"(a)")  'HistBin_ratio 1.08 ! 1.03 maybe safer. frhis(iw)= b*(exp(a*(iw-1))-1), where a=ratio-1.0 and dw=b*a' 
      write(ifi,"(a)")  '                   ! This "ba mesh" is from 9Mar2016'
      write(ifi,"(a)")  '                   ! See fpgw/gwsrc/m_freq.F'
      write(ifi,"(a)")  'iSigMode  3   ! QSGW mode switch for gwsc. use =3.'
      write(ifi,"(a)")  'niw      10   ! Number of frequencies along Im axis. Used for integration to get Sigma_c'
      write(ifi,"(a)")  '              ! To test, try niw=6 and niw=12'
      write(ifi,"(a)")  'delta  -1d-6  !(a.u.)  Broadening of x0. negative means tetrahedron method.'
      write(ifi,"(a)")  '              ! used by hx0fp0. You get smeard x0 witth abs(delta).'
      write(ifi,"(a)")  'deltaw  0.02  !(a.u.) Mesh for numerical derivative to get the Z factor'
      write(ifi,"(a)")  'esmr   0.003  !(Ry) used by hsfp0. Keep esmr smaller than band gap for insulators'
      write(ifi,"(a)")  '              ! Poles of G^LDA are treated as if they have width esmr in hsfp0. '
      write(ifi,"(a)")  '              ! Change esmr for metals.  See DOSACC*---especailly around Ef.'
      write(ifi,"(a)")  'GaussSmear on ! Gaussian or Rectangular smearing for Pole of G^LDA with esmr for hsfp0.'
      write(ifi,"(a)")  '!GaussianFilterX0 0.0001 !(a.u.) Gaussian smearing for the polarization function x0. '
      write(ifi,"(a)")  '                         ! This stabilize convergence for metallic systems'
      write(ifi,"(a)")  '                         ! This can be a default setting in the future'
      write(ifi,*)
      write(ifi,"(a)" ) '! ################################################# '
      write(ifi,"(a)")'<PRODUCT_BASIS> '
      ifigwinv2=ifi
c!! Write GWIN_V2.tmp
c      open(newunit=ifigwinv2,file='GWIN_V2.tmp')
!! PRODUCT BASIS section
      write(ifigwinv2,"(a)") " tolerance to remove products due to poor linear-independency"
      write(ifigwinv2,"(a)") 
     & " 1d-3 ! =tolopt; larger gives smaller num. of product basis."//
     & " See lbas and lbasC, which are output of hbasfp0."
      write(ifigwinv2,"(a)") 
     & " lcutmx(atom) = maximum l-cutoff for the product basis. "
     & //" =4 is required for atoms with valence d, like Ni Ga"
      allocate(lcutmx(nbas)); lcutmx=4
      write(ifigwinv2,"(1000i3)") lcutmx(1:nbas)
      write(ifigwinv2,"(a)")  "  atom   l  nnvv  nnc "
     &  //"! nnvv: num. of radial functions (valence) on the "//
     &   "augmentation-waves, nnc: num. for core."
      do ibas =1,nbas
        do lk   =0,lmxaa(ibas)
          write(ifigwinv2,"(4i5)") ibas,lk, nnvv(lk,ibas), nncx(lk,ibas)
        enddo
      enddo
      write(ifigwinv2,"(a)") "  atom   l    n  occ unocc  ! Valence(1=yes,0=no) "
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
              write(ifigwinv2,"(5i5,3x,a )") iat(izz),lindx(izz),nindx(izz)
     &                        , nocc,nunocc, '! '//caption(izz)//seg1
            endif
          enddo
        enddo
      enddo
      write(ifigwinv2,"(a)") '  atom   l    n  occ unocc  ForX0 ForSxc ! Core (1=yes, 0=no)'
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
            write(ifigwinv2,"(5i5,2x,2i5,a)") 
     &       ibas,lk,kkk-lk,noccc,nunoccc,ncinc,ncinc2,'    ! '//char(48+kkk)//seg2//seg1
          enddo
        enddo
      enddo
      write(ifi,"(a)")'</PRODUCT_BASIS>'
c      rewind ifigwinv2
c      do
c        read(ifigwinv2,"(a)",end=755) pppx
c        write(ifi,"(a)")  trim(pppx)
c      enddo
c 755  write(ifi,"(a)")'</PRODUCT_BASIS>'
c      close(ifigwinv2,status='delete')
!! QPNT section
      write(ifi,*)
      write(ifi,"(a)" )'! ################################################# '
      write(ifi,"(a)")'<QPNT> ! This block is the same as QPNT.'
      do
        read(ifqpnt,"(a)",end=756) pppx
        write(ifi,"(a)") trim(pppx)
      enddo
 756  write(ifi,"(a)")'</QPNT>'
      close(ifqpnt,status='delete')

c ... EPScond
c      inquire(file='EPScond',exist= qpntexist)
c      if(qpntexist) then
c        write(ifi,*)
c        open(ifiq, file='EPScond',status='old')
c        read(ifiq,*) epsrng, dwry !epsrng dw in Ry
c        write(ifi,"(a,f8.3,a,a)") '!EPSrange ', epsrng*2d0,' !(Ry) [0,EPSrange] for dielectric function plot.'
c        write(ifi,"(a,f8.3,a,a)") '!EPSdw    ', dwry,      ' !(Ry) energy mesh  for dielectric function plot.'
c      else
      write(ifi,"(a,f8.3,a,a)") '!EPSrange  1    !(Ry) [0,EPSrange] for dielectric function plot.'
      write(ifi,"(a,f8.3,a,a)") '!EPSdw     0.05 !(Ry) energy mesh  for dielectric function plot.'
      write(ifi,*)
      write(ifi,"(a,f8.3,a,a)") '!QforEPSIBZ on ! Use all q in IBZ for the calculation of eps mode.'
      write(ifi,"(a)") '<QforEPS>'
      write(ifi,"(a)") ' 0d0 0d0 0.01d0'
      write(ifi,"(a)") ' 0d0 0d0 0.02d0'
      write(ifi,"(a)") ' 0d0 0d0 0.04d0'
      write(ifi,"(a)") '</QforEPS>'
      write(ifi,"(a)") '!<QforEPSL>'
      write(ifi,"(a)") '! 0d0 0d0 0d0   1d0   0d0  0d0 8'
      write(ifi,"(a)") '! 0d0 0d0 0d0  .5d0  .5d0  0d0 8'
      write(ifi,"(a)") '!</QforEPSL>'
!!
      write(ifi,*)
      write(ifi,"(a)")'!!! ##### Maximally localized Wannier function ################ '
      write(ifi,"(a)")'!!! For s,p,d,f the indices 1-16 correspond to: '
      write(ifi,"(a)")'!!! index l m polynomial '
      write(ifi,"(a)")'!!! 1 0 0 1 '
      write(ifi,"(a)")'!!! -----------------------------  '
      write(ifi,"(a)")'!!! 2 1 -1 y '
      write(ifi,"(a)")'!!! 3 1 0 z  '
      write(ifi,"(a)")'!!! 4 1 1 x  '
      write(ifi,"(a)")'!!! -----------------------------  '
      write(ifi,"(a)")'!!! 5 2 -2 xy '
      write(ifi,"(a)")'!!! 6 2 -1 yz  '
      write(ifi,"(a)")'!!! 7 2 0 3z^2-1 '
      write(ifi,"(a)")'!!! 8 2 1 xz  '
      write(ifi,"(a)")'!!! 9 2 2 x^2-y^2 ' 
      write(ifi,"(a)")'!!! -----------------------------  '
      write(ifi,"(a)")'!!! 10 3 -3 y(3x^2-y^2) '
      write(ifi,"(a)")'!!! 11 3 -2 xyz '
      write(ifi,"(a)")'!!! 12 3 -1 y(5z^2-1) '
      write(ifi,"(a)")'!!! 13 3 0 z(5z^2-3) '
      write(ifi,"(a)")'!!! 14 3 1 x(5z^2-1) '
      write(ifi,"(a)")'!!! 15 3 2 z(x^2-y^2) '
      write(ifi,"(a)")'!!! ------------------------ '
      write(ifi,"(a)")'!!! higher is lm ordered. See Ylm definition in lmto/fpgw doc.'
      write(ifi,"(a)")

      ifatomlist = ifile_handle()
      open(ifatomlist,file="ATOM_LIST")
      read(ifatomlist,*)natom
      allocate(spid(natom))
      read(ifatomlist,*)spid
      close(ifatomlist)
      write(ifi,"(a)")'<Worb> Site '
      do iatom = 1,natom
        write(ifi,"(1a,i3,1x,a,2x,a)") '!', iatom,trim(spid(iatom)),' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16' 
      enddo
      write(ifi,"(a)")'</Worb> '
      write(ifi,"(a)")
      write(ifi,"(a)")'!wan_out_ewin off'
      write(ifi,"(a)")'!wan_out_bmin 16  !band index for outer window'
      write(ifi,"(a)")'!wan_out_bmax 18  !band index for outer window'
      write(ifi,"(a)")'wan_out_emin  -1.05  !eV relative to Efermi'
      write(ifi,"(a)")'wan_out_emax  2.4  !eV relative to Efermi'
      write(ifi,"(a)")'!wan_in_ewin on '
      write(ifi,"(a)")'!wan_in_emin  -1.0  !eV relative to Efermi'
      write(ifi,"(a)")'!wan_in_emax  -0.3  !eV relative to Efermi'
      write(ifi,"(a)")
      write(ifi,"(a)")'wan_tb_cut 15'
      write(ifi,"(a)")'wan_maxit_1st 300'
      write(ifi,"(a)")'wan_conv_1st 1d-7'
      write(ifi,"(a)")'wan_max_1st   0.1'
      write(ifi,"(a)")'wan_maxit_2nd 1500'
      write(ifi,"(a)")'wan_max_2nd   0.3'
      write(ifi,"(a)")'wan_conv_end  1d-8'
      write(ifi,"(a)")'!wmat_all .true.'
      write(ifi,"(a)")'!wmat_rcut1 8'
      write(ifi,"(a)")'!wmat_rcut2 0.01'
      write(ifi,"(a)")
      write(ifi,"(a)")'!vis_wan_band_n 3'
      write(ifi,"(a)")'!vis_wan_band_id 1 2 3  !integer x vis_wan_band_n, this is index for hmaxloc, as you like.'
      write(ifi,"(a)")'!vis_wan_tvec 0 0 0 !1 1 1   !integer x 3, tlat(R)'
      write(ifi,"(a)")'!vis_wan_mesh 5 5 5          !integer x 3, # of mesh'
      write(ifi,"(a)")'!vis_wan_lbound -1.2  -1.2 -1.2 !real x 3, lower bound in alat unit or abc unit'
      write(ifi,"(a)")'!vis_wan_ubound 1.2  1.2 1.2    !real x 3, upper bound in alat or abc unit'
      write(ifi,"(a)")'!vis_wan_outputformat xsf       ! opendx, cube, xsf , default=xsf'
      write(ifi,"(a)" ) '! ################################################# '
      call rx0s( ' OK! We have generated GWinput.tmp! ')
      end
