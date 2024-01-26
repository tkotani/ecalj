subroutine gwinit_v2() !  Generate GWinput.tmp.
  ! ----------------------------------------------------------------------
  ! nput file        (this doc touched at 2022jan)
  !i    HAMindex0 via readhamindex0

  !    alat       : unit in a.u. to measure plat and so on.
  !    plat(1:3,1): 1st primitive translation vector in the unit of alat
  !    plat(1:3,2): 2nd primitive translation vector
  !    plat(1:3,3): 3rd primitive translation vector
  !    QpGcut_psi : maxmum of |q+G| in a.u. in the expansion of the eigenfunction.
  !    QpGcut_Cou : maxmum of |q+G| in a.u. in the expansion of the Coulomb matrix.

  !    nbas nclass,  iclass(1:nbas)
  !    lmax(1:nclass),konf(0:lmax,1:nclass)

  !  symops : includes point group operation. See sample.
  !  nindx,lindx,ibasindex(NLAindx group): These specify the order of cphi(1:mnla,iband).
  ! utput files
  !o  GWinput.tmp : Template for GWinput
  !o  QPNT.chk    : QPNT lists, among which we specify  q for the self-enery.
  !o  KPNTin1BZ.chk: k point in 1st BZ, it is for check.
  !r True q is given by
  !r    True_q(1:3)     = 2*pi/alat * q(1:3)
  !r  True G is given by
  !r    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvec(1:3,igp)) ,igp=1,ngp
  !! -----------------------------------------------------------
  use m_hamindex0,only: readhamindex0,alat,plat,nbas,lmxax,nsp,konft,lmxaa=>lmxa, &
       nindx,lindx,mnla=>ndima,iat=>ibasindx,caption,spid,symops, ngrp,npqn,qlat,zz,pqn,ndima
  use m_get_bzdata1,only:  getbzdata1, nqbz, nqibz, qbz,wbz,qibz,wibz
  implicit none
  integer ::n1q,n2q,n3q,ifiqg,ifiqgc,ifigw0,ifi,i,ig
  real(8) :: alp,QpGcut_psi, QpGcut_Cou,dummy
  real(8) :: volum,q0(3),qlat0(3,3),QpGx1,QpGx2, &
       dw,delta,deltaw,esmr,tolopt,qm(3,3)
  integer :: ibas,l,ixxx,lmxa,ibasx,ifigw0t,mxkp, &
       irs,niw,ic,iclass, ifiqibz,iqibz,ifqpnt,iqall,iaf,iii, ifigwinv2,lk,  nocc,nunocc, &
       kkk,noccc,nunoccc,ncinc,ncinc2,isp
  integer,allocatable :: nncx(:,:),lcutmx(:)
  logical :: gwin0exist
  integer,allocatable:: IPQ(:),konf(:,:)
  real(8),allocatable   :: WGT(:)
  real(8):: qp(3)
  integer:: nnn,ifkpt,i1,nlinemax=50,ifsyml,nline,nqs
  integer,allocatable:: nqq(:),symlon(:)!,nstbz(:)
  real(8),allocatable:: qbzs(:,:),qq1(:,:),qq2(:,:)  !qbz(:,:),wbz(:),
  logical :: extsyml
  integer::checksymlon
  integer::ifinla,izz,izzn,izzz,iatbk
  integer,allocatable::nnvv(:,:)
  integer:: mnla_,idummy
  character(len=1)::seg2
  character(len=6)::seg1
  character(len=10) :: keyw1='unit_2pioa',keyw2
  real(8)::a1,a2,unit
  integer:: ibzcase 
  character(len=550):: pppx
  character(len=10) :: add,esmr_char
  logical :: gwinputexist,qpntexist,anfexist
  integer:: ifig,ifigwinp,nband_chi0,nband_sigm
  integer:: ifatomlist,ifiq
  integer:: nband_sigm2, iatom, nwin, incwfin,nx
  real,parameter:: pi= 4d0* atan(1d0)
  real(8):: pq,zc
  integer:: lx,nxx
  n1q=4                     !defaultvalues
  n2q=4
  n3q=4
  write(6,"(a,3i5)") ' Default n1q n2q n3q = ',n1q,n2q,n3q
  ! readin HAMindx0
  call readhamindex0()
  allocate(konf(0:lmxax,nbas),nncx(0:lmxax,nbas) )
  isp=1 !konf is not spin-dependent
  konf=konft(:,:,isp)
  do ibas = 1,nbas
     write(6,"(i4,f9.4,100i4)") ibas, zz(ibas), lmxaa(ibas), konf(0:lmxaa(ibas),ibas)
     do l= 0,lmxaa(ibas)
        nncx(l,ibas) = konf(l,ibas) -1 -l   ! number of cores for each l ibas
     enddo
  enddo
  write(6,"(' alat      =',f13.6 )") alat
  write(6,"(' plat a1   =',3f13.6)") plat(1:3,1)
  write(6,"(' plat a2   =',3f13.6)") plat(1:3,2)
  write(6,"(' plat a3   =',3f13.6)") plat(1:3,3)
  ! --- Make q-points in IBZ.
  mxkp   = n1q*n2q*n3q
  call getbzdata1(qlat,(/n1q,n2q,n3q/),symops,ngrp,tetrai=.false.,tetraf=.false.,mtet=(/1,1,1/),gammacellctrl=0)
  !! Write to file KPNTin1BZ
  nnn = n1q*n2q*n3q
  nqs=0
!  open(ifkpt,file='KPTin1BZ.gwinit.chk')
!  do      i1 = 1,nnn
!     call shorbz(qbz(1,i1),qp,qlat,plat)
!     write (ifkpt,"(1x,i4,4f10.5,2x,3f10.5,i3)") &
!          i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp
!  end do
!  close (ifkpt)
  write(6,"(' --- TOTAL num of q is n1*n2*n3=',i10)")nnn
  ! --- Sample QPNT file ---------------
  open (newunit=ifqpnt,file='QPNT.chk')
!  write(ifqpnt,"(a,a)") " --- Specify the q and band indeces for which we evaluate the self-energy ---"
  write(ifqpnt,"(a)")"*** all qibz for gw_lmfh -->1, otherwise 0" !;  up only -->1, otherwise 0"
  iqall = 1 !;      iaf   = 0
  write(ifqpnt,*) iqall !,iaf
  write(ifqpnt,"(a)")"*** # of states and band index for gw_lmfh calculation."
  iii = min(ndima,100) !this is the number of MTO (max of iii is ndima+napw(lowest))
  write(ifqpnt,*)  iii ! nband
  write(ifqpnt,'(*(g0,x))') (i,i=1,iii)
!  write(ifqpnt,"(a)") "*** q-points, which shoud be in qbz. See KPNTin1BZ."
!  write(ifqpnt,*) min(nqibz,3)
!  write(ifqpnt,'(i3,3f23.16)')(i,qibz(1:3,i),i=1,nqibz)
  rewind ifqpnt
  allocate(nnvv(0:lmxax,nbas))
  nnvv = 0
  mnla_=0
  do izz=1, mnla
     !        write(6,"(4i5,2x,'!',a)") izz,nindx(izz),lindx(izz),iat(izz),caption(izz)
     mnla_= mnla_ + 2*lindx(izz)+1
     if(nnvv(lindx(izz),iat(izz)) < nindx(izz)) then
        nnvv(lindx(izz),iat(izz)) = nindx(izz)
     endif
  enddo

  !! Write GWinput.tmp ===============================================
  open(newunit=ifigwinp,file='GWinput.tmp')
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
  write(ifi,"(a)" )'! ##### From GWIN0 ################ '
  write(ifi,"(a,3i4,a)" )'n1n2n3',n1q,n2q,n3q,' ! for BZ meshing in GW, Wannier function and cRPA'
  write(ifi,"(a)") 'QpGcut_psi 4.0  !(See unit_2pioa for unit) |q+G| cutoff for eigenfunction.'
  write(ifi,"(a)") 'QpGcut_cou 3.0  !(See unit_2pioa for unit) |q+G| cutoff for Coulomb and W.'
  write(ifi,"(a)") 'unit_2pioa off ! off --> a.u.; on--> unit of QpGcut_* are in 2*pi/alat '
  write(ifi,"(a)") 'alpha_OffG 1.0 !(a.u.) Used in auxially function in the offset-Gamma method.'
  !      write(ifi,"(a,i8,a)")     'nband_chi0 ', nband_chi0,' !    nband cutoff for chi0  (Optional)'
  !      write(ifi,"(a,2i8,a)")     'nband_sigm ', nband_sigm
  !     &      ,nband_sigm2,' !    nband cutoff for Sigma  (Optional) (1st:num in sigma; 2nd: num of G used in hsfp0)'
  !      write(ifi,"(a,2f10.3,a)")  'emax_sigm  ', emax_sigm ,
  !     &       emax_sigm2,'  !(Ry)  (Optional) emax cutoff for Sigma (as in the nband_sigm)'
  !      write(ifi,"(a,2i8,a)")     'nband_sigm ', nband_sigm
  !     &      ,nband_sigm2,' !    nband cutoff for Sigma  (Optional) (1st:num in sigma; 2nd: num of G used in hsfp0)'
  !     &       emax_sigm2,'  !(Ry)  (Optional) emax cutoff for Sigma (as in the nband_sigm)'

  write(ifi,"(a)") '!emax_chi0  999 !(Ry) emax cutoff for chi0  (Optional)'
  write(ifi,"(a)") 'emax_sigm  3.0  !(Ry)  emax cutoff for Sigma'
  !!
  write(ifi,*)
  write(ifi,"(a)" ) '! ##### FREQUENCIES from GWIN_V2 ################ '
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
!  ifi=ifi
  !! PRODUCT BASIS section
  write(ifi,"(a)") " tolerance to remove products due to poor linear-independency"
  write(ifi,"(a)") &
       " 1d-3 ! =tolopt; larger gives smaller num. of product basis."// &
       " See lbas and lbasC, which are output of hbasfp0."
  write(ifi,"(a)") &
       " lcutmx(atom) = maximum l-cutoff for the product basis. " &
       //" =4 is required for atoms with valence d, like Ni Ga"
  allocate(lcutmx(nbas)); lcutmx=4
  do ibas=1,nbas
     if(zz(ibas)<10.5) lcutmx(ibas)=2
     if(57.001<zz(ibas).and.zz(ibas)<71.001) lcutmx(ibas)=6
     if(89.001<zz(ibas).and.zz(ibas)<71.001) lcutmx(ibas)=6
  enddo
  write(ifi,"(1000i3)") lcutmx(1:nbas)
  write(ifi,"(a)")  "  atom   l  nnvv  nnc " &
       //"! nnvv: num. of radial functions (valence) on the "// &
       "augmentation-waves, nnc: num. for core."
  do ibas =1,nbas
     do lk   =0,lmxaa(ibas)
        write(ifi,"(4i5)") ibas,lk, nnvv(lk,ibas), nncx(lk,ibas)
     enddo
  enddo
  write(ifi,"(a)") "  atom   l    n  occ unocc  ! Valence(1=yes,0=no) "
  iatbk=0
  do ibas= 1, nbas
     do lk   = 0, lmxaa(ibas)
        do nx = 1, npqn
           do izz = 1, mnla
              if(iat(izz)==ibas .AND. lk==lindx(izz) .AND. nx==nindx(izz)) then
                 zc = zz(ibas)
                 pq = pqn(izz)
                 lx = lindx(izz)
                 nxx = nindx(izz)
                 nunocc = 1
                 nocc=0
                 ! s
                 if(lx==0.and.pq==2.and.zc>1.5)  nocc= 1
                 if(lx==0.and.pq==3.and.zc>4.5)  nocc= 1 !sp boundary
                 if(lx==0.and.pq==4.and.zc>12.5) nocc= 1
                 if(lx==0.and.pq==5.and.zc>30.5) nocc= 1
                 if(lx==0.and.pq==6.and.zc>48.5) nocc= 1
                 if(lx==0.and.pq==7.and.zc>80.5) nocc= 1
                 if(lx==0.and.pq==2.and.zc>10.5) nunocc= 0 !Ne+
                 if(lx==0.and.pq==3.and.zc>18.5) nunocc= 0
                 if(lx==0.and.pq==4.and.zc>36.5) nunocc= 0
                 if(lx==0.and.pq==5.and.zc>54.5) nunocc= 0
                 if(lx==0.and.pq==6.and.zc>86.5) nunocc= 0
                 ! p
                 if(lx==1.and.pq==2.and.zc>1.5)  nocc= 1
                 if(lx==1.and.pq==3.and.zc>4.5)  nocc= 1
                 if(lx==1.and.pq==4.and.zc>12.5) nocc= 1
                 if(lx==1.and.pq==5.and.zc>30.5) nocc= 1
                 if(lx==1.and.pq==6.and.zc>48.5) nocc= 1
                 if(lx==1.and.pq==7.and.zc>80.5) nocc= 1
                 if(lx==1.and.pq==2.and.zc>10.5) nunocc= 0
                 if(lx==1.and.pq==3.and.zc>18.5) nunocc= 0
                 if(lx==1.and.pq==4.and.zc>36.5) nunocc= 0
                 if(lx==1.and.pq==5.and.zc>54.5) nunocc= 0
                 if(lx==1.and.pq==6.and.zc>86.5) nunocc= 0
                 ! d occupied
                 if(lx==2 .and. zc >20.5 .and. pq==3 ) nocc=1 !3d sd boundary
                 if(lx==2 .and. zc >38.5 .and. pq==4 ) nocc=1 !4d
                 if(lx==2 .and. zc >56.5 .and. pq==5 ) nocc=1 !5d
                 if(lx==2 .and. zc >30.5 .and. pq==3 ) nunocc=0 !3d completely filled for Zn< zc
                 if(lx==2 .and. zc >48.5 .and. pq==4 ) nunocc=0 !4d 
                 ! f occupied
                 if(lx==3 .and. zc >57.5 .and. pq==4 ) nocc=1 !4f La+
                 if(lx==3 .and. zc >89.5 .and. pq==5 ) nocc=1 !5f
                 if(lx==2 .and. zc >71.5 .and. pq==4 ) nunocc=0 !4f completely filled.
                 ! g or more
                 if(lx >3 ) nocc   = 0 
                 if(lx >3 ) nunocc = 0
                 !phidot skip
                 if(nxx==2)  nocc   = 0 
                 if(nxx==2)  nunocc = 0 
                 seg1=''
                 if(iat(izz)/=iatbk) seg1=' -----'
                 iatbk=iat(izz)
                 write(ifi,"(5i5,3x,a )") iat(izz),lindx(izz),nindx(izz) &
                      , nocc,nunocc, '! '//caption(izz)//trim(seg1)
                 exit
              endif
           enddo
        enddo
     enddo
  enddo
  write(ifi,"(a)") '  atom   l    n  occ unocc  ForX0 ForSxc ! Core (1=yes, 0=no)'
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
           write(ifi,"(5i5,2x,2i5,a)") &
                ibas,lk,kkk-lk,noccc,nunoccc,ncinc,ncinc2,'    ! '//char(48+kkk)//seg2//trim(seg1)
        enddo
     enddo
  enddo
  write(ifi,"(a)")'</PRODUCT_BASIS>'
  !! QPNT section
  write(ifi,*)
  write(ifi,"(a)" )'! ################################################# '
  write(ifi,"(a)")'<QPNT> ! This block is the same as QPNT.'
  do
     read(ifqpnt,"(a)",end=756) pppx
     write(ifi,"(a)") trim(pppx)
  enddo
756 write(ifi,"(a)")'</QPNT>'
  close(ifqpnt,status='delete')

  write(ifi,"(a,f8.3,a,a)") '!EPSrange  1    !(Ry) [0,EPSrange] for dielectric function plot.'
  write(ifi,"(a,f8.3,a,a)") '!EPSdw     0.05 !(Ry) energy mesh  for dielectric function plot.'
  write(ifi,*)
  write(ifi,"(a,f8.3,a,a)") '!QforEPSIBZ on ! Use all q in IBZ for the calculation of eps mode.'
  write(ifi,"(a)") 'QforEPSunita on'
  write(ifi,"(a)") '<QforEPS>'
  write(ifi,"(a)") ' 0 0 0.00001   '
  write(ifi,"(a)") ' 0 0 0.001     '
  write(ifi,"(a)") ' 0 0 0.0014142 '
  write(ifi,"(a)") ' 0 0 0.002     '
  write(ifi,"(a)") ' 0 0 0.0028284 '
  write(ifi,"(a)") ' 0 0 0.004'
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
  write(ifi,"(a)")'<Worb> Site '
  do iatom = 1,nbas
     write(ifi,"(1a,i3,1x,a,2x,a)") &
          '!', iatom,trim(spid(iatom)),' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
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
  stop ' OK! We have generated GWinput.tmp! '
end subroutine gwinit_v2
