real(8) function deltaq_scale()
  ! Q0Pchoice=1: qzerolimit. (not too small because of numerical reason.)
  ! Q0Pchoice=2: =1d0/3.0**.5d0/Q is the mean value of \int_{|q|<Q} d^3 q <1/q^2> for a sphere.
  use m_keyvalue,only: getkeyvalue
  integer,save :: ttt=1
  logical,save:: init=.true.
  if(init) then
     call getkeyvalue("GWinput","Q0Pchoice",ttt,default=1)
     write(6,"('  Q0Pchoice=',i3)") ttt
     init=.false.
  endif
  if(ttt==1) then
     deltaq_scale=0.1d0 !this is essentially q to zero limit.
  elseif(ttt==2) then
     deltaq_scale=1d0/3.0**.5d0
  else
     call rx( 'Use Q0Pchoice = 1 or 2 (=1 is default)')
  endif
END function deltaq_scale
logical function localfieldcorrectionllw()
  use m_keyvalue,only: getkeyvalue
  logical,save:: init=.true.,ttt
  if(init) then
     call getkeyvalue("GWinput","LFC@Gamma",ttt,default=.true.)
     write(6,*)'LFC@Gamma=',ttt
     init=.false.
  endif
  localfieldcorrectionllw=ttt
end function localfieldcorrectionllw
logical function addbasnew()
  addbasnew=.true. !new version of adding polynomial-like product basis in MTs.
  ! If false, use old version (less product basis set is added.)
end function addbasnew
!$$$      logical function newaniso()
!$$$      newaniso=.true. !new GW mode. not the offset Gamma method.
!$$$                      !combines methods by two Christoph.
!$$$      end
real(8) function screenfac()
  ! The Coulomb interaction is given as exp(- screeenfac()*r)/r (see hvccfp0.m.F)
  ! Formally scrennfac=0d0 is correct, however, we can not choose too small screenfac
  ! in out current implementation. For example,
  ! screenfac=-1d-8 gives NaN for GaAs222 test-->This gives negative eigenvalue of Vcoul for q=0
  use m_keyvalue,only: getkeyvalue
  real(8):: ddd
  real(8),save :: tss
  logical,save:: init=.true.
  if(init) then
     call getkeyvalue("GWinput","TFscreen",tss, default=1d-5**.5)
     ! 1d-5**.5 is just given by rough test.
     ! Results should not depend on this value as long as default is small enough.
     write(6,*)'TFscreen=',tss
     init=.false.
  endif
  ! screenfac = - TFscreen**2 = energy (negative) ==> (\nabla^2 + e) v= \delta(r-r')
  screenfac= -tss**2 !-ttt  !note negative sign for exp(-sqrt(e)r)
END function screenfac
!! testmode
logical function testomitq0()
  testomitq0=.false.
end function testomitq0
!      integer function nomatm()
!c In cases (corei7+gfortran and so on in my case),zgemm did not work.
!c now in ppbafp.fal.F
!      nomatm=1 !use matmul instead of zgemm called from a subroutine matm
!      end
!===========================================================
!      subroutine headver(head,id1)
!      character*(*) head
!      write(6,"(a,a, i3)") head,": VerNum= xxx: imode=",id1
!      if(id1==-9999) call rx( '---end of version check!')
!      end
!============================================================
logical function is_mix0vec()
  ! s_mis0vec=.false. is original version. But it caused a problem at BZ bounday.
  is_mix0vec=.true.
end function is_mix0vec
logical function evaltest()
  evaltest=.false.
end function evaltest
!      logical function test_symmetric_W()
!      use m_keyvalue
!      logical,save:: init=.true.,ttt
!      if(init) then
!        call getkeyvalue("GWinput","TestSymmetricW",ttt,default=.false.)
!        init=.false.
!      endif
!      test_symmetric_W= ttt
!      end

!      logical function testtr()
!      testtr=.true.
!      end
!      logical function negative_testtr()
!      negative_testtr=.true.
!      end

logical function TimeReversal()
  use m_keyvalue,only: getkeyvalue
  logical,save:: init=.true.,trevc
  if(init) then
     call getkeyvalue("GWinput","TimeReversal",trevc,default=.true.)
     init=.false.
  endif
  timereversal= trevc
end function TimeReversal

logical function oncew()
  logical,save::init=.true.
  if(init) then
     oncew=.true.
     init=.false.
  else
     oncew=.false.
  endif
end function oncew

logical function onceww(i)
  integer:: i
  logical,save::init(100)=.true.
  if(init(i)) then
     onceww=.true.
     init(i)=.false.
  else
     onceww=.false.
  endif
end function onceww
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! for future use. NaN generator.
real(8) function NaNdble()
  real(8):: d
  d = 1d0-1d0
  NaNdble= (1d0-1d0)/d
END function NaNdble
!      real(8) function NaNdble2()
!      NaNdble2= (1d0-1d0)/(1d0-1d0)
!      end
complex(8) function NaNcmpx()
  real(8):: NaNdble
  NaNcmpx=cmplx(NaNdble(),NaNdble())
END function NaNcmpx
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function rmeshrefine()
  use m_keyvalue,only: getkeyvalue
  call getkeyvalue("GWinput","rmeshrefine",rmeshrefine,default=.true.)
end function rmeshrefine
real(8) function delrset()
  ! dr/dI at rmat. used for rmeshrefin=T case
  use m_keyvalue,only: getkeyvalue
  call getkeyvalue("GWinput","dRdIatRmax",delrset,default=0.003d0)
END function delrset

logical function qbzreg()
  use m_keyvalue,only: getkeyvalue
  logical,save:: init=.true.,ccrq
  if(init) then
     call getkeyvalue("GWinput","chi_RegQbz",ccrq,default=.true.)
     init=.false.
  endif
  qbzreg= ccrq
end function qbzreg
logical function smbasis() !
  use m_keyvalue,only: getkeyvalue
  integer(4),save:: smbasis0
  logical,save:: init=.true.,smbasis00=.false.
  if(init) then
     call getkeyvalue("GWinput","smbasis",smbasis0,default=0)
     init=.false.
     if(smbasis0>0) then
        smbasis00 = .true.
     endif
  endif
  smbasis = smbasis00
end function smbasis
integer(4) function smbasiscut() !
  use m_keyvalue,only: getkeyvalue
  integer(4),save:: smbasis0
  logical,save:: init=.true.
  if(init) then
     call getkeyvalue("GWinput","smbasis",smbasis0,default=0)
     init=.false.
  endif
  smbasiscut = smbasis0
END function smbasiscut
integer function smbasis_case() !
  smbasis_case = 1
END function smbasis_case

!      logical function ngczero()
!! ngczero=T: Use ngc for given q even for iqi>nqibz
!! This is for "regular meshing" epsmode (now only for ix=23 mode).
!!      ngczero=.false. !default
!!---------------------------------
!      ngczero= .true. !false.  ! ! false is for gw_lmfh
!      end                        ! true is now only for eps_lmf_chi

logical function qreduce() !
  ! remove the inequivalent q points (G vector shifts)
  qreduce= .true. ! false is safer for usual mode, gw_lmfh
end function qreduce                !(But I think true is OK---not tested completely).
! rue may reduce the size of eigen function files (Cphi Geig).

!      integer(4) function saveiq()
!! I think saveiq=1 is (maybe a little) better to accelate GW calculation.
!! But saveiq=0 may stop calculation when you do multi-k point eps mode.
!!      ---> Then set sqveiq=0
!      saveiq=1
!      end

! Long-range-only Coulomb interaction
real(8) function eees()
  use m_keyvalue,only: getkeyvalue
  logical,save:: init=.true.
  real(8),save:: eee
  real(8):: r0cs
  if(init) then
     call getkeyvalue("GWinput","removed_r0c",r0cs,default=1d60)
     eee = -1d0/r0cs**2
     if(r0cs>1d10) eee=0d0
  endif
  eees = eee
end function eees

!      logical function testsemif() !test for semicore
!      testsemif=.false.
!      end
real(8) function scissors_x0()
  use m_keyvalue,only: getkeyvalue
  use m_ReadEfermi,only: readefermi,ef,bandgap
  !      integer(4):: iopen
  logical,save:: init=.true.
  real(8),save:: sciss !,bandgap,ef
  if(init) then
     call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
     !        call readefermi() !ef bandgap !this can cause a problem to change ef very remote manner via m_ReadEfermi
     !        ifi  = iopen('EFERMI',1,0,0)
     !        read(ifi,*) ef,bandgap
     !        close(ifi)
     init=.false.
  endif
  scissors_x0 = (sciss-1d0) * bandgap
END function scissors_x0


!---------------------------------------------------
integer(4) function zvztest()
  !---------------------
  zvztest=0
  !----No test:
  !      zvztest=0 !! not zvztest mode
  !----
  !     zvztest=1  ! test1  <psi_i psi_j  M_I ><M_I v M_J><M_J psi_i psi_j >
  !----
  !     zvztest=2  ! test2   |M_1> =  phi_s*phi_s basis case for Li. Set product basis as only
  !                    1    0    3    1    1   ! 1S_l
END function zvztest

!      integer function version()
!      version=0
!      end

!      logical function onlyimagaxis()
!      use m_keyvalue,only: getkeyvalue
!      logical,save ::init=.true.,onlyi
!      integer(4):: ret
!      if(init) then
!        call getkeyvalue("GWinput","OnlyImagAxis",onlyi,default=.false.,status=ret )
!        init=.false.
!      endif
!      onlyimagaxis=onlyi
!      end

!      logical function cphigeig_mode()
!c Whether you get cphi and geig from CphiGeig, or DATA4GW.
!c See cphigeig_mode()=.false. is for older method to store eigenfunctions.
!      cphigeig_mode=.true.
!      end

logical function matrix_linear()
  use m_keyvalue,only: getkeyvalue
  ! Use linear interpolation for matrix elements (numerator) in tetrahdron-weight's calculation.
  ! matrix_linear=T seems to give little improvements.
  logical,save::init=.true.,matrix_linear0
  if(init) then
     call getkeyvalue("GWinput","tetrahedron_matrix_linear",matrix_linear0,default=.false.)
     init=.false.
  endif
  matrix_linear=matrix_linear0
end function matrix_linear

!      logical function ifgeigb()
! See rdpp_v2 x0kf_v2hx. only for epsPP_lmfh mode now.
! This option reduce memory usage.
!      use m_keyvalue,only: getkeyvalue
!      call getkeyvalue("GWinput","UseGeigBFile",ifgeigb,default=.false.)
!      end

!      logical function KeepPPOVL()
!      use m_keyvalue,only: getkeyvalue
!c! Keep data from PPOVL in memory or not; in getppx in rdppovl.f.
!c KeepPPOVL=T : speed up
!c KeepPPOVL=F : efficient memory usage
!      logical,save:: init=.true.,Keepppovl0
!      if(init) then
!        call getkeyvalue("GWinput","KeepPPOVL",KeepPPOVL0,default=.true.)
!        init=.false.
!      endif
!      keepppovl = keepppovl0
!      end

logical function KeepEigen()
  use m_keyvalue,only: getkeyvalue
  !! Keep data from CPHI and GEIG in memory or not; in readeigen
  ! KeepEigen=T : speed up
  ! KeepEigen=F : efficient memory usage
  logical,save::init=.true.,keepeigen0
  if(init) then
     call getkeyvalue("GWinput","KeepEigen",KeepEigen0,default=.true.)
     init=.false.
  endif
  keepeigen = keepeigen0
end function KeepEigen

!      logical function readgwinput()
!c Use GWinput instead of GWIN0, GWIN_V2, QPNT
!      readgwinput=.true.
!      end

!$$$      logical function core_orth()
!$$$      use m_keyvalue,only: getkeyvalue
!$$$      logical,save::init=.true.,core_orthx
!$$$      integer(4):: ret
!$$$      if(init) then
!$$$        call getkeyvalue("GWinput","CoreOrth",core_orthx,default=.false. )
!$$$        init=.false.
!$$$      endif
!$$$      core_orth=core_orthx
!$$$  end

integer(4) function verbose()
  use m_keyvalue,only: getkeyvalue
  logical,save ::init=.true.,ggg
  !      logical:: readgwinput
  integer(4):: ret
  integer(4),save::verbosex
  if(init) then
     inquire(file='GWinput',exist=ggg)
     if(ggg) then
        call getkeyvalue("GWinput","Verbose",verbosex,default=0 )
     else
        verbosex=0
     endif
     !        write(6,*)' verbose=',verbosex
     init=.false.
  endif
  verbose=verbosex
END function verbose

!$$$      logical function GaussSmear()
!$$$      use m_keyvalue,only: getkeyvalue
!$$$C- smergin switch for SEx and SEc.
!$$$c      GaussSmear=.true. ! Gaussian smering.
!$$$c      GaussSmear=.false.! original rectoangular smering.
!$$$c It seems that you might need to use narrower esmer in GWIN_V2
!$$$c when you use GaussSmear=.true.
!$$$      logical ::init=.true. !,readgwinput
!$$$      logical,save :: GaussSmearx=.false.
!$$$      character(len=150):: recrdxxx
!$$$      character(len=130):: recrdxxx0
!$$$      character(len=10)  :: keyw1='GaussSmear',keyw2
!$$$      if(init) then
!$$$c        if(readgwinput()) then
!$$$        call getkeyvalue("GWinput","GaussSmear",GaussSmearx )
!$$$c        else
!$$$c         ifinin = 8087
!$$$c         open(ifinin,file='GWIN_V2')
!$$$c         do i=1,10; read(ifinin,*); enddo
!$$$c         read(ifinin,"(130a)") recrdxxx0
!$$$c         recrdxxx = recrdxxx0//' #'
!$$$c         read(recrdxxx,*) a1, keyw2
!$$$c         if(keyw1==keyw2) GaussSmearx=.true.
!$$$c         close(ifinin)
!$$$c         write(6,*)' GaussSmear=',GaussSmearx
!$$$c        endif
!$$$        init=.false.
!$$$      endif
!$$$      GaussSmear=GaussSmearx
!$$$      end

integer function q0pchoice()
  use m_keyvalue,only: getkeyvalue
  !- Switch whether you use new seeting Q0P (offsetted Gamma).
  ! q0pchoice=0: old---along plat
  ! q0pchoice=1: new---along Ex Ey Ez.

  ! See q0irre.f
  logical,save ::init=.true.
  !      logical:: readgwinput
  integer(4),save:: ret,q0pchoicex
  if(init) then
     !       if(readgwinput()) then
     !         write(6,*)' goto getkeyvalue'
     call getkeyvalue("GWinput","Q0P_Choice",q0pchoicex,default=0) !,status=ret )
     !       endif
     init=.false.
  endif
  q0pchoice=q0pchoicex
end function q0pchoice

logical function tetra_hsfp0()
  ! for tetrahedron method of hsfp0. See hsfp0.m.f or so.
  !     & , tetraex  = .false. ! This switch is only meaningful for mode=1,5,6
  !                            ! If you want to calculate exchange, use tetraex=T .
  !                            ! Note that you have to supply EFERMI by the tetrahedon method.
  tetra_hsfp0=.false.
end function tetra_hsfp0

! bzcase==1 only now
!$$$c--------------------------------------------
!$$$      integer(4) function bzcase()
!$$$      bzcase=1
!$$$      end

!$$$      use m_keyvalue,only: getkeyvalue
!$$$c      bzcase==2 is for regular-mesh in BZ without gamma point.
!$$$      logical,save ::init=.true.
!$$$c      logical:: readgwinput
!$$$      integer(4),save::bzcasex
!$$$      if(init) then
!$$$c        if(readgwinput()) then
!$$$        call getkeyvalue("GWinput","BZmesh",bzcasex,default=1)
!$$$c        endif
!$$$        init=.false.
!$$$      endif
!$$$      bzcase=bzcasex
!$$$      end
real(8) function wgtq0p() !essentially dummy
  wgtq0p=0.01d0
END function wgtq0p
!$$$c  This is effective only for bzcase==2.
!$$$c  Offset gamma has integration weight of 8*wgtq0p*(number of BZmesh)
!$$$c  On the other hand, we subtract weight of wgtq0p*(number of BZmesh)
!$$$c  from the eight nearest mesh points to Gamma.

real(8) function escale()
  use m_keyvalue,only: getkeyvalue   !c--- used q0pchoice<0 mode -------
  call getkeyvalue("GWinput","q0scale",escale,default=0.8d0)
END function escale
integer(4) function normcheck()
  use m_keyvalue,only: getkeyvalue
  ! write normcheck files or not
  ! normcheck=0: not
  ! normcheck=1: only dia
  ! normcheck=2: dia and off
  integer(4),save::nnn
  logical,save ::init=.true.
  if(init) then
     call getkeyvalue("GWinput","NormChk",nnn,default=1)
     init=.false.
     !        write(6,"('NormChk mode=',i3)")nnn
  endif
  normcheck=nnn
END function normcheck
integer(4) function auxfunq0p()
  ! =0: usual auxially function  exp(-alpha |q|^2) /q^2
  ! =1: new auxally function     exp(-alpha |q|) /q^2
  !      use m_keyvalue,only: getkeyvalue
  !      integer(4)::nnn
  !      logical,save ::init=.true.
  !      if(init) then
  !        call getkeyvalue("GWinput","AuxFunQ0P",nnn,default=0)
  !        init=.false.
  !c        write(6,"('AuxFun mode=',i3)")nnn
  !      endif
  !      auxfunq0p=nnn
  auxfunq0p=0
END function auxfunq0p
subroutine wgtscale(q1,w1,w2)
  implicit none
  real(8):: q1,q2,w1,qav2,w2,qav
  q2 = 1d0
  !      qav2= (q1**2 +q2**2+q1*q2)/3d0
  !      w1 =  (1d0/qav2 - 1d0/q2**2)/ (1d0/q1**2 - 1d0/q2**2)
  qav = 3d0/4d0*(q2**4-q1**4)/(q2**3-q1**3)
  w1 =  (qav - q2)/ (q1 - q2)
  w2 =  1d0 - w1
end subroutine wgtscale

