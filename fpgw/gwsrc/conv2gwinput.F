      subroutine conv2gwinput()
c      use m_keyvalue,only:getkeyvalue
      implicit none
      integer(4) :: ifigw0,ifi,i,ig
      real(8) :: QpGcut_psi, QpGcut_Cou
      real(8) ::  dw,delta,deltaw,esmr
      logical :: gwin0exist

      character(len=150):: recrdxxx
      character(len=550):: pppx
      character(len=10) :: keyw1='unit_2pioa',keyw2,add,esmr_char,eeg='GaussSmear'
c'
      real(8)::a1,a2,unit,pi

      logical :: gwinputexist,qpntexist,anfexist
      integer(4):: ifig,ifigwinp,n1,n2,n3,nband_chi0,nband_sigm
     & ,isigmode,ifigwin,niw,ifiq,ixx
      real(8)::  alpha_OffG, emax_chi0, emax_sigm ,omg_c,
     &  epsrng, dwry
      integer(4)::nnn(3),n

      real(8):: emax_sigm2
      integer(4):: nband_sigm2
      pi= 4d0* atan(1d0)

c -----------
c      inquire(file='GWinput',exist= gwinputexist)
c      write(*,*) 'gwinputexist=', gwinputexist
c      if(gwinputexist) return

c      if(gwinputexist) then
c        call getkeyvalue("GWinput","n1n2n3",nnn,3)
c        call getkeyvalue("GWinput","QpGcut_psi",QpGcut_psi)
c        call getkeyvalue("GWinput","QpGcut_cou",QpGcut_cou)
c        call getkeyvalue("GWinput","nband_sigm",nband_sigm,default=99999)
c        call getkeyvalue("GWinput","nband_chi0",nband_chi0,default=99999)
c        n1=nnn(1);n2=nnn(2);n3=nnn(3)
c        call wGWIN0(n1,n2,n3,QpGcut_psi, QpGcut_Cou, max(nband_sigm,nband_chi0))
c        return
c      endif

c--- Write GWinput from GWIN0 GWIN_V2 QPNT.
      ifigwinp=3501
      ifigwin =3502
      ifigw0  =3503
      ifig    =3504
      open(ifigw0, file='GWIN0',status='old')
      open(ifig,   file='GWIN_V2',status='old')
      write(6,*) " --- Convert data from GWIN0 to GWinput. --- "
      open(ifigwinp,file='GWinput')

c ... GWIN0
      nband_chi0=9999
      emax_chi0 =9999
      nband_sigm=9999
      emax_sigm = 2.0
      emax_sigm2=9999
      nband_sigm2=9999

      read(ifigw0,*)
      read(ifigw0,*) n1,n2,n3
      read(ifigw0,*)

      read(ifigw0,"(150a)") recrdxxx
      read(recrdxxx,*) keyw2
      if(keyw1==keyw2) then
        read(recrdxxx,*) keyw2, QpGcut_psi, QpGcut_Cou
        write(6,*) ' unit_2pioa find: 2pi/alat  '
        ixx=1
      else
        read(recrdxxx,*) QpGcut_psi, QpGcut_Cou
        ixx=0
      endif

      read(ifigw0,*)
      read(ifigw0,*) alpha_OffG
      read(ifigw0,*,end=556)
      read(ifigw0,*,err=555,end=556) nband_chi0, emax_chi0
 555  read(ifigw0,*,err=556,end=556) nband_sigm, emax_sigm
 556  continue
      close(ifigw0)
c
      ifi = ifigwinp
      write(ifi,"(a)")'!!! Starting from ! (or nonkeyword) is comment line !!! '
      write(ifi,"(a)")'!!! Each lines consists of "keyword value(s)"  !!! '
      write(ifi,"(a)")'!!! Each tag section in <...>... </...> has its own format. !!! '
      write(ifi,"(a)")'!EIBZmode off   !no symmetrization for hx0fp0* (default on);Only affects comp. effort. off may faster.'
      write(ifi,"(a)")'!chi_RegQbz off !Use no Gamma mesh for dielectric function. This automaticall set EIBZmode off.'
      write(ifi,"(a)")'!Verbose    0  ! 0-->default; 100--->debug '
      write(ifi,"(a)")'!LFC@Gamma off !(on is default) if on, eps with Local field correction is used at Gamma point'
      write(ifi,"(a)")'!Q0Pchoice 1 !1(default):qzerolimit(in practice, See generated Q0P), 2:1/q^2 average in Gamma region'
      write(ifi,"(a)")'!CoreOrth  off  ! off  --> Not enforce core to be orthogonalized with valence (default)'
      write(ifi,"(a)")'                ! on   --> Enforce cores orthogonalied to valence (this may give strange core functions!).'
      write(ifi,"(a)")'!multitet 2 2 2 ! tetrahedron is devided into micro tetrahedron. This probably now work'
c
      write(ifi,"(a)" )'! ##### From GWIN0 ################ '
      write(ifi,"(a,3i5,a)" )   'n1n2n3     ', n1,n2,n3,  ' ! for BZ meshing in GW '
      write(ifi,"(a,f8.3,a,a)") 'QpGcut_psi ', QpGcut_psi,' !(See unit_2pioa for unit) |q+G| cutoff for eigenfunction.'
      write(ifi,"(a,f8.3,a,a)") 'QpGcut_cou ', QpGcut_cou,' !(See unit_2pioa for unit) |q+G| cutoff for Coulomb and W.'
      if(ixx==0)write(ifi,"(a)")'unit_2pioa off ! off --> a.u.; on--> unit of QpGcut_* are in 2*pi/alat '
      if(ixx==1)write(ifi,"(a)")'unit_2pioa on  ! on--> unit of QpGcut_* are in 2*pi/alat; off --> a.u.'
      write(ifi,"(a,f8.3,a)")   'alpha_OffG ', alpha_offG,' !(a.u.) Used in auxially function in the offset-Gamma method.'
c      write(ifi,"(a,i8,a)")     'nband_chi0 ', nband_chi0,' !    nband cutoff for chi0  (Optional)'
c      write(ifi,"(a,2i8,a)")     'nband_sigm ', nband_sigm
c     &      ,nband_sigm2,' !    nband cutoff for Sigma  (Optional) (1st:num in sigma; 2nd: num of G used in hsfp0)'
c      write(ifi,"(a,2f10.3,a)")  'emax_sigm  ', emax_sigm ,
c     &       emax_sigm2,'  !(Ry)  (Optional) emax cutoff for Sigma (as in the nband_sigm)'
c      write(ifi,"(a,2i8,a)")     'nband_sigm ', nband_sigm
c     &      ,nband_sigm2,' !    nband cutoff for Sigma  (Optional) (1st:num in sigma; 2nd: num of G used in hsfp0)'
c     &       emax_sigm2,'  !(Ry)  (Optional) emax cutoff for Sigma (as in the nband_sigm)'

      write(ifi,"(a,f10.3,a)") '!emax_chi0  ', emax_chi0,' !(Ry) emax cutoff for chi0  (Optional)'
      write(ifi,"(a,f10.3,a)") 'emax_sigm  ', emax_sigm,' !(Ry)  emax cutoff for Sigma'

c ... GWIN_V2
c      dw = 0.005d0
c      omg_c = 0.04d0
c      iSigMode= 3
c      niw = 10      !niw =10
c      delta  = -1d-6
c      deltaw = 0.02d0
c      esmr   = 0.003d0 ! 0.01d0

      read(ifig,*)
      read(ifig,*)
      read(ifig,*)
      read(ifig,*)
      read(ifig,*) 
      read(ifig,*)
      read(ifig,*) 
      read(ifig,*)
      read(ifig,*) 
      read(ifig,*)
      read(ifig,*) 
      read(ifig,*)
      read(ifig,*)
      read(ifig,*)
      read(ifig,*)
c
      write(ifi,*)
      write(ifi,"(a)" ) '! ##### FREQUENCIES from GWIN_V2 ################ '
c      write(ifi,"(a,f10.6,a)") 'dw    ', dw,   ' !(a.u.) energy-mesh (bin width size) along real axis.'
c      write(ifi,"(a,f10.3,a)") 'omg_c ', omg_c,' !(a.u.) energy-mesh is twiced at omg_c'
      write(ifi,"(a)")  'HistBin_dw    2d-3  ! 1d-5 for metal  !(a.u.) bin width along real axis at omega=0.'
      write(ifi,"(a)")  'HistBin_ratio 1.08  ! 1.03 for check. frhis(iw)= b*(exp(a*(iw-1))-1), where a=ratio-1.0 and dw=b*a' 
      write(ifi,"(a)")  '                    ! This "ba mesh" is from 9Mar2016'
      write(ifi,"(a)")  '                    ! See fpgw/gwsrc/m_freq.F'
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
      write(ifi,*)
      write(ifi,"(a)" ) '! ################################################# '
      write(ifi,"(a)")'<PRODUCT_BASIS> '
      write(ifi,"(a)")" tolerance to remove products due to poor linear-independency"
      do
        read(ifig,"(a)",end=755) pppx
        call cuttail(pppx,n)
        write(ifi,"(a)")  pppx(1:n)
      enddo
 755  write(ifi,"(a)")'</PRODUCT_BASIS>'

c ... QPNT
      inquire(file='QPNT',exist= qpntexist)
      if(qpntexist) then
        ifiq=3505
        open(ifiq,    file='QPNT',status='old')
        write(ifi,*)
        write(ifi,"(a)" )'! ################################################# '
        write(ifi,"(a)")'<QPNT> ! This block is the same as QPNT.'
        do
          read(ifiq,"(a)",end=756) pppx
          call cuttail(pppx,n)
          write(ifi,"(a)")   pppx(1:n)
        enddo
 756    write(ifi,"(a)")'</QPNT>'
      endif

cc ... ANFcond
c      inquire(file='ANFcond',exist= anfexist)
c      if(anfexist) then
c        open(ifiq,    file='ANFcond',status='old')
c        write(ifi,*)
c        write(ifi,"(a)" )'! ################################################# '
c        write(ifi,"(a)")'<ANFcond> ! This block is the same as ANFcond'
c        do
c        read(ifiq,"(a)",end=757)  pppx(1:n)
c        write(ifi,"(a)")  pppx(1:n)
c        enddo
c 757    write(ifi,"(a)")'</ANFcond>'
c      endif

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
c      endif
c      close(ifi)
      return
      end
c
!      subroutine wGWIN0(n1,n2,n3,QpGcut_psi, QpGcut_Cou,nband)
!      real(8):: QpGcut_psi,QpGcut_cou
!      open(ifigw0,  file='GWIN0') !This GWIN0 is only for lmfgw and lmf2gw
!      write(ifigw0,"(a)")' --- These are taken from GWinput by qg4gw. Used in lmfgw and lmf2gw. Edit GWinput! ---'
!      write(ifigw0,"(3i8)") n1,n2,n3
!      write(ifigw0,*) ' --- '
!      write(ifigw0,"(2d13.5)") QpGcut_psi, QpGcut_Cou
!      write(ifigw0,*) ' --- '
!      write(ifigw0,*) ' --- '
!      write(ifigw0,*) ' --- '
!      write(ifigw0,*) nband
!      write(ifigw0,*) ' --- End --- '
!      close(ifigw0)
!      end

      subroutine cuttail(a,n)
      integer(4) :: n
      character*(*) :: a
      do n=550,1,-1
        if(a(n:n)/=' ') exit
      enddo
      end
