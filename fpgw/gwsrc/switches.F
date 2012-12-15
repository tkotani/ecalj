      logical function localfieldcorrectionllw()
      use keyvalue,only: getkeyvalue
      logical,save:: init=.true.,ttt
      if(init) then
        call getkeyvalue("GWinput","LFC@Gamma",ttt,default=.true.)
        print *,'LFC@Gamma=',ttt
        init=.false.
      endif
      localfieldcorrectionllw=ttt
      end

      logical function eibz4sig()
      eibz4sig=.true.
      end

      logical function mtosigmaonly() !this is also in lm7k/subs/switches.F
      mtosigmaonly=.true.
      end

      logical function eibz4x0()
      use keyvalue,only: getkeyvalue
      logical,save:: init=.true.,eibzmode
      if(init) then
        call getkeyvalue("GWinput","EIBZmode",eibzmode,default=.true.)
        init=.false.
      endif
      eibz4x0=eibzmode
! If T, use EIBZ procedure in PRB125102(2010). Not completed yet... 
!WARN: time-reversal is not included yet. (see hx0fp0.m.F L1195 call eibzgen).
!      eibz. In addition, inefficient method of symmetrization at the bottom of x0kf_v4h
!      call eibzgen generates EIBZ.
!      Probably, it will make things easier to start from "inversion mesh for q point".
      end

      logical function addbasnew()
      addbasnew=.true. !new version of adding polynomial-like product basis in MTs.
      ! If false, use old version (less product basis set is added.)
      end

      logical function newaniso()
      newaniso=.true. !new GW mode. not the offset Gamma method.
                      !combines methods by two Christoph. 
      end

      real(8) function screenfac()
c      srceenfac=0d0
c      screenfac=-1d-8 !NaN for GaAs222 test. This gives negative eigenvalue of Vcoul for q=0
      screenfac=-1d-5 
      end

!! testmode
      logical function testomitq0()
      testomitq0=.false.
      end

c      integer function nomatm()
cc In cases (corei7+gfortran and so on in my case),zgemm did not work.
cc now in ppbafp.fal.F
c      nomatm=1 !use matmul instead of zgemm called from a subroutine matm
c      end
c===========================================================
      subroutine headver(head,id1)
      character*(*) head
      write(6,"(a,a, i3)") head,": VerNum= xxx: imode=",id1
      if(id1==-9999) stop '---end of version check!'
      end
c============================================================
      logical function is_mix0vec()
!is_mis0vec=.false. is original version. But it caused a problem at BZ bounday.
      is_mix0vec=.true.
      end

      logical function evaltest()
      evaltest=.false.
      end

      logical function test_symmetric_W()
      use keyvalue
      logical,save:: init=.true.,ttt
      if(init) then
        call getkeyvalue("GWinput","TestSymmetricW",ttt,default=.false.)
        init=.false.
      endif
      test_symmetric_W= ttt
      end

c      logical function testtr()
c      testtr=.true.
c      end
c      logical function negative_testtr()
c      negative_testtr=.true.
c      end

      logical function TimeReversal()
      use keyvalue
      logical,save:: init=.true.,trevc
      if(init) then
        call getkeyvalue("GWinput","TimeReversal",trevc,default=.true.)
        init=.false.
      endif
      timereversal= trevc
      end

      logical function oncew()
      logical,save::init=.true.
      if(init) then
        oncew=.true.
        init=.false.
      else
        oncew=.false.
      endif
      end

      logical function onceww(i)
      logical,save::init(100)=.true.
      if(init(i)) then
        onceww=.true.
        init(i)=.false.
      else
        onceww=.false.
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for future use. NaN generator.
      real(8) function NaNdble()
      real(8):: d
      d = 1d0-1d0
      NaNdble= (1d0-1d0)/d
      end
c      real(8) function NaNdble2()
c      NaNdble2= (1d0-1d0)/(1d0-1d0)
c      end
      complex(8) function NaNcmpx()
      real(8):: NaNdble
      NaNcmpx=cmplx(NaNdble(),NaNdble())
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical function rmeshrefine()
      use keyvalue
      call getkeyvalue("GWinput","rmeshrefine",rmeshrefine,default=.true.)
      end
      real(8) function delrset()
! dr/dI at rmat. used for rmeshrefin=T case
      use keyvalue
      call getkeyvalue("GWinput","dRdIatRmax",delrset,default=0.003d0)
      end

      logical function qbzreg()
      use keyvalue
      logical,save:: init=.true.,ccrq
      if(init) then
        call getkeyvalue("GWinput","chi_RegQbz",ccrq,default=.true.)
        init=.false.
      endif
      qbzreg= ccrq
      end


      logical function smbasis() !
      use keyvalue
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
      end
      integer(4) function smbasiscut() !
      use keyvalue
      integer(4),save:: smbasis0
      logical,save:: init=.true.
      if(init) then
        call getkeyvalue("GWinput","smbasis",smbasis0,default=0)
        init=.false.
      endif
      smbasiscut = smbasis0
      end
      integer(4) function smbasis_case() !
      smbasis_case = 1
      end

c      logical function ngczero()
c! ngczero=T: Use ngc for given q even for iqi>nqibz
c! This is for "regular meshing" epsmode (now only for ix=23 mode).
c!      ngczero=.false. !default
c!---------------------------------
c      ngczero= .true. !false.  ! ! false is for gw_lmfh
c      end                        ! true is now only for eps_lmf_chi

      logical function qreduce() !
c remove the inequivalent q points (G vector shifts)
      qreduce= .true. ! false is safer for usual mode, gw_lmfh
      end                !(But I think true is OK---not tested completely).
      !true may reduce the size of eigen function files (Cphi Geig).

c      integer(4) function saveiq()
c! I think saveiq=1 is (maybe a little) better to accelate GW calculation.
c! But saveiq=0 may stop calculation when you do multi-k point eps mode.
c!      ---> Then set sqveiq=0
c      saveiq=1
c      end

c Long-range-only Coulomb interaction
      real(8) function eees()
      use keyvalue
      logical,save:: init=.true.
      real(8),save:: eee
      real(8):: r0cs
      if(init) then
        call getkeyvalue("GWinput","removed_r0c",r0cs,default=1d60)
        eee = -1d0/r0cs**2
        if(r0cs>1d10) eee=0d0
      endif
      eees = eee
      end

c      logical function testsemif() !test for semicore
c      testsemif=.false.
c      end
      real(8) function scissors_x0()
      use keyvalue
      integer(4):: iopen
      logical,save:: init=.true.
      real(8),save:: sciss,bandgap,ef
      if(init) then
        call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
        ifi  = iopen('EFERMI',1,0,0)
        read(ifi,*) ef,bandgap
        close(ifi)
        init=.false.
      endif
      scissors_x0 = (sciss-1d0) * bandgap
      end


c---------------------------------------------------
      integer(4) function zvztest()
c---------------------
      zvztest=0
c----No test:
c      zvztest=0 !! not zvztest mode
c----
c     zvztest=1  ! test1  <psi_i psi_j  M_I ><M_I v M_J><M_J psi_i psi_j >
c----
c     zvztest=2  ! test2   |M_1> =  phi_s*phi_s basis case for Li. Set product basis as only
c                    1    0    3    1    1   ! 1S_l
      end

      integer function version()
      version=0
      end

c      logical function onlyimagaxis()
c      use keyvalue
c      logical,save ::init=.true.,onlyi
c      integer(4):: ret
c      if(init) then
c        call getkeyvalue("GWinput","OnlyImagAxis",onlyi,default=.false.,status=ret )
c        init=.false.
c      endif
c      onlyimagaxis=onlyi
c      end

c      logical function cphigeig_mode()
cc Whether you get cphi and geig from CphiGeig, or DATA4GW.
cc See cphigeig_mode()=.false. is for older method to store eigenfunctions.
c      cphigeig_mode=.true.
c      end

      logical function matrix_linear()
      use keyvalue
c Use linear interpolation for matrix elements (numerator) in tetrahdron-weight's calculation.
c matrix_linear=T seems to give little improvements.
      logical,save::init=.true.,matrix_linear0
      if(init) then
        call getkeyvalue("GWinput","tetrahedron_matrix_linear",matrix_linear0,default=.false.)
        init=.false.
      endif
      matrix_linear=matrix_linear0
      end

c      logical function ifgeigb()
c!See rdpp_v2 x0kf_v2hx. only for epsPP_lmfh mode now.
c!This option reduce memory usage.
c      use keyvalue
c      call getkeyvalue("GWinput","UseGeigBFile",ifgeigb,default=.false.)
c      end

      logical function KeepPPOVL()
      use keyvalue
c! Keep data from PPOVL in memory or not; in getppx in rdppovl.f.
c KeepPPOVL=T : speed up
c KeepPPOVL=F : efficient memory usage
      logical,save:: init=.true.,Keepppovl0
      if(init) then
        call getkeyvalue("GWinput","KeepPPOVL",KeepPPOVL0,default=.true.)
        init=.false.
      endif
      keepppovl = keepppovl0
      end

      logical function KeepEigen()
      use keyvalue
c! Keep data from CPHI and GEIG in memory or not; in readeigen
c KeepEigen=T : speed up
c KeepEigen=F : efficient memory usage
      logical,save::init=.true.,keepeigen0
      if(init) then
        call getkeyvalue("GWinput","KeepEigen",KeepEigen0,default=.true.)
        init=.false.
      endif
      keepeigen = keepeigen0
      end

c      logical function readgwinput()
cc Use GWinput instead of GWIN0, GWIN_V2, QPNT
c      readgwinput=.true.
c      end

      logical function core_orth()
      use keyvalue
      logical,save::init=.true.,core_orthx
c      logical:: readgwinput
      integer(4):: ret
      if(init) then
c       if(readgwinput()) then
        call getkeyvalue("GWinput","CoreOrth",core_orthx,default=.false. )
c       endif
        init=.false.
      endif
      core_orth=core_orthx
      end

      integer(4) function verbose()
      use keyvalue
      logical,save ::init=.true.,ggg
c      logical:: readgwinput
      integer(4):: ret
      integer(4),save::verbosex
      if(init) then
        inquire(file='GWinput',exist=ggg)
        if(ggg) then
          call getkeyvalue("GWinput","Verbose",verbosex,default=0 )
        else
          verbosex=0
        endif
c        print *,' verbose=',verbosex
        init=.false.
      endif
      verbose=verbosex
      end

      logical function GaussSmear()
      use keyvalue
C- smergin switch for SEx and SEc.
c      GaussSmear=.true. ! Gaussian smering.
c      GaussSmear=.false.! original rectoangular smering.
c It seems that you might need to use narrower esmer in GWIN_V2
c when you use GaussSmear=.true.
      logical ::init=.true. !,readgwinput
      logical,save :: GaussSmearx=.false.
      character(len=150):: recrdxxx
      character(len=130):: recrdxxx0
      character(len=10)  :: keyw1='GaussSmear',keyw2
      if(init) then
c        if(readgwinput()) then
        call getkeyvalue("GWinput","GaussSmear",GaussSmearx )
c        else
c         ifinin = 8087
c         open(ifinin,file='GWIN_V2')
c         do i=1,10; read(ifinin,*); enddo
c         read(ifinin,"(130a)") recrdxxx0
c         recrdxxx = recrdxxx0//' #'
c         read(recrdxxx,*) a1, keyw2
c         if(keyw1==keyw2) GaussSmearx=.true.
c         close(ifinin)
c         print *,' GaussSmear=',GaussSmearx
c        endif
        init=.false.
      endif
      GaussSmear=GaussSmearx
      end

      integer function q0pchoice()
      use keyvalue
C- Switch whether you use new seeting Q0P (offsetted Gamma).
C q0pchoice=0: old---along plat
C q0pchoice=1: new---along Ex Ey Ez.
C
C See q0irre.f
      logical,save ::init=.true.
c      logical:: readgwinput
      integer(4),save:: ret,q0pchoicex
      if(init) then
c       if(readgwinput()) then
c         print *,' goto getkeyvalue'
        call getkeyvalue("GWinput","Q0P_Choice",q0pchoicex,default=0) !,status=ret )
c       endif
        init=.false.
      endif
      q0pchoice=q0pchoicex
      end

      logical function tetra_hsfp0()
c for tetrahedron method of hsfp0. See hsfp0.m.f or so.
c     & , tetraex  = .false. ! This switch is only meaningful for mode=1,5,6
c                            ! If you want to calculate exchange, use tetraex=T .
c                            ! Note that you have to supply EFERMI by the tetrahedon method.
      tetra_hsfp0=.false.
      end

c--------------------------------------------
      integer(4) function bzcase()
      use keyvalue
c      bzcase==2 is for regular-mesh in BZ without gamma point.
      logical,save ::init=.true.
c      logical:: readgwinput
      integer(4),save::bzcasex
      if(init) then
c        if(readgwinput()) then
        call getkeyvalue("GWinput","BZmesh",bzcasex,default=1)
c        endif
        init=.false.
      endif
      bzcase=bzcasex
      end

c------------------------------------------
      real(8) function wgtq0p()
c  This is effective only for bzcase==2.
c  Offset gamma has integration weight of 8*wgtq0p*(number of BZmesh)
c  On the other hand, we subtract weight of wgtq0p*(number of BZmesh)
c  from the eight nearest mesh points to Gamma.
c------------------------------------------------------
      use keyvalue
      logical,save ::init=.true.
      real(8),save:: www
      real(8):: escale,w1,w2
      integer(4)::  q0pchoice
c      common /wgtq00/www
      logical :: noq0p
      if(init) then
        call getkeyvalue("GWinput","WgtQ0P",www,default=0.01d0)
c        write(6,"('wgtq0p (weight for Q0P for BZmesh=2)=',f13.8)")www
        if(q0pchoice()<0)  then
          www = escale()**3 !1d0-wres
c2          call wgtscale(escale(),w1,w2)
c2           www = 1d0 - w2*( 1d0-escale()**3 )
cc          www = 1d0
        endif
        call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
        if(noq0p) www=0d0
        init =.false.
      endif
      wgtq0p = www
      return
      end

      real(8) function escale()
      use keyvalue
cc--- used q0pchoice<0 mode -------
      call getkeyvalue("GWinput","q0scale",escale,default=0.8d0)
      end

c      subroutine  rewrite_wgtq0p(aaa)
c      real(8):: www,aaa
c      common /wgtq00/www
c      www = aaa
c      print *, ' www=',www
c      return
c      end

c      logical function test_omitq0p()
c      test_omitq0p=.false.
c      end

      integer(4) function normcheck()
      use keyvalue
c write normcheck files or not
c normcheck=0: not
c normcheck=1: only dia
c normcheck=2: dia and off
      integer(4)::nnn
      logical,save ::init=.true.
      if(init) then
        call getkeyvalue("GWinput","NormChk",nnn,default=1)
        init=.false.
c        write(6,"('NormChk mode=',i3)")nnn
      endif
      normcheck=nnn
      end

      integer(4) function auxfunq0p()
c =0: usual auxially function  exp(-alpha |q|^2) /q^2
c =1: new auxally function     exp(-alpha |q|) /q^2
c      use keyvalue
c      integer(4)::nnn
c      logical,save ::init=.true.
c      if(init) then
c        call getkeyvalue("GWinput","AuxFunQ0P",nnn,default=0)
c        init=.false.
cc        write(6,"('AuxFun mode=',i3)")nnn
c      endif
c      auxfunq0p=nnn
      auxfunq0p=0
      end

      subroutine wgtscale(q1,w1,w2)
      implicit none
      real(8):: q1,q2,w1,qav2,w2,qav
      q2 = 1d0
c      qav2= (q1**2 +q2**2+q1*q2)/3d0
c      w1 =  (1d0/qav2 - 1d0/q2**2)/ (1d0/q1**2 - 1d0/q2**2)
      qav = 3d0/4d0*(q2**4-q1**4)/(q2**3-q1**3)
      w1 =  (qav - q2)/ (q1 - q2)
      w2 =  1d0 - w1
      end

