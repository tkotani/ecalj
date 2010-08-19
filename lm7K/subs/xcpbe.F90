! This is taken from ABINIT. http://www.abinit.org/  abinit-6.2.1 on Aug7,2010. 
!
! defs_basis.F, xcpbe.F invcb.F, which are taken from abinit-6.2.1,
! are slightly modified and contained in this file.
! :: call wrtout(std_out,message,'COLL') +  call leave_new('COLL') --->  replaced by subroutine rx(message)
!
!  Since ABINIT is under the GPL, T.Kotani think this file is
!  somehow affeced by the licence.
!  But T.Kotani knows little about what the GPL means.
!

!########################################################################################## 
!### defs_basis.F90
!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_basis
!! NAME
!! defs_basis
!!
!! FUNCTION
!! This module contains definitions for a number of named constants and
!! physical constants.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (HM, XG,XW)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! Of the named constants,
!! by far the most important are those that define the 'kind' types of
!! virtually all the variables used in a (well-written) FORTRAN 90 code
!! the content of this file is derived from 'Numerical Recipes in Fortran 90'
!! W.H. Press et al., volume 2 of 'Fortran Numerical Recipes', Cambridge
!! University Press, Second Edition (1996), p. 937 and 1361
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_basis

 implicit none

!Keyword 'integer' stands for default integer type
!and may be used whenever integer are presumed to be small

!nb of bytes related to an integer subtype n such as -10^(argument) < n < 10^(argument) (this is standard F90)
 integer, parameter :: i1b=selected_int_kind(2)
 integer, parameter :: i2b=selected_int_kind(4)
 integer, parameter :: i4b=selected_int_kind(9)
 integer, parameter :: i8b=selected_int_kind(18)

!nb of bytes related to default simple-precision real/complex subtypes
!(= 4 for many machine architectures, = 8 for e.g. Cray)
 integer, parameter :: sp=kind(1.0)          ! Single precision should not be used
 integer, parameter :: spc=kind((1.0,1.0))

!nb of bytes related to default double-precision real/complex subtypes
!(= 8 for many machine architectures)
 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0_dp,1.0_dp))  ! Complex should not be used presently
                                                  ! except for use of libraries

!nb of bytes related to GW arrays, that can be tuned from sp to dp independently
!of other variables in ABINIT. Presently single precision is the default.
#if defined HAVE_GW_DPC
 integer, parameter :: gwp=kind(1.0d00)
 integer, parameter :: gwpc=kind((1.0_dp,1.0_dp))
#else
 integer, parameter :: gwp=kind(1.0)
 integer, parameter :: gwpc=kind((1.0,1.0))
#endif

!Example:
! integer, parameter :: urp=selected_real_kind((p=)12,(r=)50)
! real((kind=)urp) :: d
! d=5.04876_urp   ! for a real d with 12 significative digits
! and such as 10^-50 < |d| < 10^50

!To modify sp/spc and / or dp/dpc, insert instructions such as 'dp='
! but do not modify the other declarations in this module

!Default logical type
 integer, parameter :: lgt=kind(.true.)

!The default lengths
 integer, parameter :: fnlen=264     ! maximum length of file name variables
 integer, parameter :: strlen=2000000 ! maximum length of input string

!Some constants:
 !integer, parameter :: integer_not_used=0
 !logical, parameter :: logical_not_used=.true.

!UNIX unit numbers : standard input, standard output, ab_out, and a number
!for temporary access to a file.
 integer, parameter :: std_in=5,ab_in=5  ! generally, the number 5 is directly used
 integer, parameter :: std_out=6         ! generally, the number 6 is directly used
 integer, parameter :: ab_out=7
 integer, parameter :: ab_xml_out = 50  ! this unit is used to print output into an XML file
 integer, parameter :: tmp_unit=9,tmp_unit2=10

!The 3x3 identity matrix
!WARNING : this seem not to work ?!
! integer, dimension(3,3), parameter :: &
!& identity3by3=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))

!Real constants
 real(dp), parameter :: zero=0._dp
 real(dp), parameter :: one=1._dp
 real(dp), parameter :: two=2._dp
 real(dp), parameter :: three=3._dp
 real(dp), parameter :: four=4._dp
 real(dp), parameter :: five=5._dp
 real(dp), parameter :: six=6._dp
 real(dp), parameter :: seven=7._dp
 real(dp), parameter :: eight=8._dp
 real(dp), parameter :: nine=9._dp
 real(dp), parameter :: ten=10._dp

!Fractionary real constants
 real(dp), parameter :: half=0.50_dp
 real(dp), parameter :: onehalf=1.50_dp
 real(dp), parameter :: third=one/three
 real(dp), parameter :: quarter=0.25_dp
 real(dp), parameter :: fifth=0.20_dp
 real(dp), parameter :: sixth=one/six
 real(dp), parameter :: seventh=one/seven
 real(dp), parameter :: eighth=0.125_dp
 real(dp), parameter :: ninth=one/nine
 real(dp), parameter :: two_thirds=two*third
 real(dp), parameter :: four_thirds=four*third
 real(dp), parameter :: five_thirds=five*third
 real(dp), parameter :: three_quarters=0.75_dp
 real(dp), parameter :: three_fifth=three/five

!Real constants related to the golden number
 real(dp), parameter :: gold=1.618033988749894848204586834365638117720309179_dp
 real(dp), parameter :: goldenratio=two-gold

!Real constants derived from pi
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: two_pi=two*pi
 real(dp), parameter :: four_pi=four*pi
 real(dp), parameter :: piinv=one/pi
!The following are not used
!real(dp), parameter :: rad_to_deg=180._dp/pi
!real(dp), parameter :: deg_to_rad=one/rad_to_deg
!real(dp), parameter :: half_pi=pi*half
!real(dp), parameter :: third_pi=pi*third
!real(dp), parameter :: quarter_pi=pi*quarter
!real(dp), parameter :: two_thirds_pi=two_thirds*pi


!Real precision
 real(dp), parameter :: greatest_real = huge(one)
 real(dp), parameter :: smallest_real = -greatest_real
 real(dp), parameter :: tol3= 0.001_dp
 real(dp), parameter :: tol4= 0.0001_dp
 real(dp), parameter :: tol5= 0.00001_dp
 real(dp), parameter :: tol6= 0.000001_dp
 real(dp), parameter :: tol7= 0.0000001_dp
 real(dp), parameter :: tol8= 0.00000001_dp
 real(dp), parameter :: tol9= 0.000000001_dp
 real(dp), parameter :: tol10=0.0000000001_dp
 real(dp), parameter :: tol11=0.00000000001_dp
 real(dp), parameter :: tol12=0.000000000001_dp
 real(dp), parameter :: tol13=0.0000000000001_dp
 real(dp), parameter :: tol14=0.00000000000001_dp
 real(dp), parameter :: tol15=0.000000000000001_dp
 real(dp), parameter :: tol16=0.0000000000000001_dp

!real constants derived from sqrt(n.)
 real(dp), parameter :: sqrt2=1.4142135623730950488016887242096939_dp
 real(dp), parameter :: half_sqrt2=0.70710678118654752440084436210484697_dp
 real(dp), parameter :: sqrt3=1.7320508075688772935274463415058739_dp
 real(dp), parameter :: half_sqrt3=0.86602540378443864676372317075293693_dp
 real(dp), parameter :: sqrthalf=0.70710678118654752440084436210484697_dp

!Conversion factors of common use, not directly related to physical quantities.
 real(dp), parameter :: b2Mb=one/1024.0_dp**2  ! conversion factor bytes --> Mbytes
 real(dp), parameter :: b2Gb=b2Mb/1000.0_dp    ! conversion factor bytes --> Gbytes

!Real physical constants
!Revised fundamental constants from http://physics.nist.gov/cuu/Constants/index.html
!(from 2006 least squares adjustment)
 real(dp), parameter :: Bohr_Ang=0.52917720859_dp    ! 1 Bohr, in Angstrom
 real(dp), parameter :: Ha_cmm1=219474.6313705_dp  ! 1 Hartree, in cm^-1
 real(dp), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
 real(dp), parameter :: Ha_K=315774.65_dp ! 1Hartree, in Kelvin
 real(dp), parameter :: Ha_THz=6579.683920722_dp ! 1 Hartree, in THz
 real(dp), parameter :: Ha_J=4.35974394d-18    !1 Hartree, in J
 real(dp), parameter :: e_Cb=1.602176487d-19 ! minus the electron charge, in Coulomb
 real(dp), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
 real(dp), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
!This value is 1Ha/bohr^3 in 1d9 J/m^3
!real(dp), parameter :: HaBohr3_GPa=29421.033_dp ! 1 Ha/Bohr^3, in GPa
 real(dp), parameter :: HaBohr3_GPa=Ha_eV/Bohr_Ang**3*e_Cb*1.0d+21 ! 1 Ha/Bohr^3, in GPa
 real(dp), parameter :: Avogadro=6.02214179d23 ! per mole
!This value is 1 Ohm.cm in atomic units
 real(dp), parameter :: Ohmcm=two*pi*Ha_THz*ninth*ten
!real(dp), parameter :: eps0=8.854187817d-12 ! permittivity of free space in F/m
 real(dp), parameter :: eps0=one/(four_pi*0.0000001_dp*299792458.0_dp**2)
 real(dp), parameter :: AmuBohr2_Cm2=e_Cb*1.0d20/(Bohr_Ang*Bohr_Ang)
 real(dp), parameter :: InvFineStruct=137.035999679_dp  ! Inverse of fine structure constant
 real(dp), parameter :: Sp_Lt=2.99792458d8/2.1876912633d6 ! speed of light in atomic units
 real(dp), parameter :: Time_Sec=2.418884326505D-17 !  Atomic unit of time, in seconds

!Complex constants
 complex(dpc), parameter :: czero=(0._dp,0._dp)
 complex(dpc), parameter :: cone =(1._dp,0._dp)
 complex(dpc) ,parameter :: j_dpc=(0._dp,1.0_dp)

!Character constants
 character(len=1), parameter :: ch10 = char(10)

!Define fake communicator for sequential abinit
 integer, parameter :: abinit_comm_serial = -12345

 ! Error codes used by the bindings.
 integer, parameter, public :: AB6_NO_ERROR                = 0
 integer, parameter, public :: AB6_ERROR_OBJ               = 1
 integer, parameter, public :: AB6_ERROR_ARG               = 2
 integer, parameter, public :: AB6_ERROR_INVARS_ATT        = 3
 integer, parameter, public :: AB6_ERROR_INVARS_ID         = 4
 integer, parameter, public :: AB6_ERROR_INVARS_SIZE       = 5
 integer, parameter, public :: AB6_ERROR_SYM_NOT_PRIMITIVE = 6
 integer, parameter, public :: AB6_ERROR_SYM_BRAVAIS_XRED  = 7

! Values of optdriver corresponding to the different run-levels. 
 integer, parameter, public :: RUNL_GSTATE     = 0
 integer, parameter, public :: RUNL_RESPFN     = 1
 integer, parameter, public :: RUNL_SUSCEP     = 2
 integer, parameter, public :: RUNL_SCREENING  = 3
 integer, parameter, public :: RUNL_SIGMA      = 4
 integer, parameter, public :: RUNL_NONLINEAR  = 5
 integer, parameter, public :: RUNL_RDM        = 7
 integer, parameter, public :: RUNL_SCGW       = 8 
 integer, parameter, public :: RUNL_BSE        = 99 !9

! Flags defining the method used for performing IO (input variable accesswff)
 integer, parameter, public :: IO_MODE_FORTRAN_MASTER = -1
 integer, parameter, public :: IO_MODE_FORTRAN        =  0
 integer, parameter, public :: IO_MODE_MPI            =  1
 integer, parameter, public :: IO_MODE_NETCDF         =  2 ! Only for legacy code, should not be used for new implementations.
 integer, parameter, public :: IO_MODE_ETSF           =  3

CONTAINS  !==============================================================================
!!***

!!****f* defs_basis/get_reclen
!! NAME
!!  get_reclen
!!
!! FUNCTION
!!  Return the record lenght used to store a variable of particular type and kind in processor-dependent units.
!!  The returned value can be used as the value of the lenght to be supplied to a RECL= specifier.
!!
!! INPUTS
!!  rcl=string defining the type and kind of the variable. Possible values are:
!!   i1b, i2b, i4b, dp, gwp, spc, dpc, gwpc (and corresponding capital versions).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_reclen(str) result(rcl)


 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: str
 integer :: rcl

!Local variables ------------------------------
!scalars
 integer(i1b)  :: v_i1b 
 integer(i2b)  :: v_i2b
 integer(i4b)  :: v_i4b 
 !real(sp)    :: v_sp 
 real(dp)      :: v_dp 
 real(gwp)     :: v_gwp
 complex(spc)  :: v_spc 
 complex(dpc)  :: v_dpc 
 complex(gwpc) :: v_gwpc

!************************************************************************

 SELECT CASE (str)

 CASE ("i1b","I1B")
  inquire(iolength=rcl) v_i1b

 CASE ("i2b","I2B")
  inquire(iolength=rcl) v_i2b

 CASE ("i4b","I4B")
  inquire(iolength=rcl) v_i4b

 !CASE ("sp","SP")
 ! inquire(iolength=rcl) v_sp

 CASE ("dp","DP")
  inquire(iolength=rcl) v_dp

 CASE ("gwp","GWP")
  inquire(iolength=rcl) v_gwp

 CASE ("spc","SPC")
  inquire(iolength=rcl) v_spc

 CASE ("dpc","DPC")
  inquire(iolength=rcl) v_dpc

 CASE ("gwpc","GWPC")
  inquire(iolength=rcl) v_gwpc

 CASE DEFAULT
  write(*,*)" Unknown kind: "//TRIM(str)
  STOP
 END SELECT

end function get_reclen
!!***

!----------------------------------------------------------------------

!!****f* defs_basis/print_kinds
!! NAME
!! print_kinds
!!
!! FUNCTION
!! Prints info on the basic data types, e.g. kind, precision...
!!
!! INPUTS
!!   unit = Unit number for output file.
!!
!! OUTPUT
!!   Only printing.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

 subroutine print_kinds(unit)

 !Arguments ------------------------------------

 integer,intent(in) :: unit

 !Local variables-------------------------------
 !character(len=500) :: msg
 ! *********************************************************************

 write(unit,'(a)')' DATA TYPE INFORMATION: '

 write(unit,'(a,/,2(a,i6,/),2(a,e14.8,/),a,e14.8)')&
& ' REAL:      Data type name: REAL(DP) ',&
& '            Kind value: ',KIND(0.0_dp),&
& '            Precision:  ',PRECISION(0.0_dp),&
& '            Smallest nonnegligible quantity relative to 1: ',EPSILON(0.0_dp),&
& '            Smallest positive number:                      ',TINY(0.0_dp),&
& '            Largest representable number:                  ',HUGE(0.0_dp)

 write(unit,'(a,/,2(a,i0,/),a,i0)')&
  ' INTEGER:   Data type name: INTEGER(default) ', &
& '            Kind value: ',KIND(0),              &
& '            Bit size:   ',BIT_SIZE(0),          &
  '            Largest representable number: ',HUGE(0)

 write(unit,'(a,/,a,i0)')&
& ' LOGICAL:   Data type name: LOGICAL ',&
& '            Kind value: ',KIND(.TRUE.)

 write(unit,'(2a,i0)')&
& ' CHARACTER: Data type name: CHARACTER ',&
& '            Kind value: ',KIND('C')

end subroutine print_kinds

!----------------------------------------------------------------------

end module defs_basis
!!***



!########################################################################################## 
!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcpbe
!! NAME
!! xcpbe
!!
!! FUNCTION
!! Treat XC functionals closely linked with the Perdew-Wang 92 LSD
!! and the PBE GGA.
!!
!! For a series of values of the density and, if GGA, the square of the
!! gradient of the density, return the associated Exc energy,
!! potential, and, in case of response-function, functions needed
!! to build the XC kernel.
!!
!! If option==2, Exchange-correlation functional from Perdew-Burke-Ernzerhof,
!! Phys.Rev.Lett. 77, 3866 (1996).
!! If option==1, Reduces to Perdew-Wang LSD , PRB45,13244 (1992).
!! If option==-1 or -2, take only exchange part of PW (-1) or PBE (-2) functionals.
!! If option==-4, C09x exchange functional of V. R. Cooper, PRB 81, 161104(R) (2010).
!! If option==3 or 4, take exchange plus RPA correlation
!!   part of LSD PW (3) or GGA PBE (4) functionals.
!! If option==5, revPBE functional of Zhang and Yang, PRL 80, 890 (1998)
!! If option==6, RPBE functional of Hammer, Hansen and Norskov, PRB 59, 7413 (1999)
!! If option==7, WC functional of Wu and Cohen, PRB 73, 235116 (2006)
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG,MF,LG,CE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  exexch= choice of local exact exchange. Active if exexch=1
!!  npts= number of points to be computed
!!  nspden=1 for unpolarized, 2 for spin-polarized
!!  grho2_updn(npts,ngr2)=square of the gradient of the spin-up,
!!     and, if nspden==2, spin-down, and total density (Hartree/Bohr**2),
!!     only used if gradient corrected functional (option=2,-2,-4 and 4 or beyond)
!!  option= see above
!!  order=its absolute value gives the maximal derivative of Exc to be computed.
!!  rho_updn(npts,nspden)=spin-up and spin-down density (Hartree/bohr**3)
!!  ndvxci= size of dvxci(npts,ndvxci)
!!  ngr2= size of grho2_updn(npts,ngr2)
!!  nd2vxci=size of d2vxci(npts,nd2vxci)
!!
!! OUTPUT
!!
!!  d2vxci=third derivative of the xc energy with respect to the density, only
!!    only if local-density approximation
!!   calculated if order==3
!!   In case of local energy functional (option=1,-1 or 3):
!!    d2vxci(npts,nd2vxc)=              (Hartree*bohr^3)
!!     if(nspden=1): d2vxci(:,1)=-(2/3)*dvxci/d$\rho$ 
!!                                  (dvxci is the second derivative, see below)                      
!!     if(nspden=2): d2vxci(:,1)=3rd order derivative of XC energy with respect to rhouprhouprhoup,
!!                   d2vxci(:,2)=3rd order derivative of XC energy with respect to rhouprhouprhodn
!!                   d2vxci(:,3)=3rd order derivative of XC energy with respect to rhodnrhouprhodn
!!                   d2vxci(:,4)=3rd order derivative of XC energy with respect to rhodnrhodnrhodn
!!  dvxcdgr(npts,3)=partial derivative of the exchange-correlation
!!    energy (exci*$\rho$) with respect to the spin-up (dvxcdgr(:,1)),
!!    spin-down (dvxcdgr(:,2)), or total spin (dvxcdgr(:,3)) gradients of the density
!!    divided by the norm of the gradient (the definition changed in v3.3)
!!
!!  dvxci=partial second derivatives of the xc energy, only if abs(order)>1
!!   In case of local energy functional (option=1,-1 or 3):
!!    dvxci(npts,1+nspden)=              (Hartree*bohr^3)
!!     if(nspden=1 .and. order==2): dvxci(:,1)=dvxc/d$\rho$ , dvxc(:,2) empty
!!     if(nspden=1 .and. order==-2): also compute dvxci(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$
!!     if(nspden=2): dvxci(:,1)=dvxc($\uparrow$)/d$\rho(\uparrow)$,
!!                   dvxci(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!                   dvxci(:,3)=dvxc($\downarrow$)/d$\rho(\downarrow)$
!!   In case of gradient corrected functional (option=2,-2, 4, 5, 6, -4):
!!    dvxci(npts,15)=
!!     dvxci(:,1)= d2Ex/drho_up drho_up
!!     dvxci(:,2)= d2Ex/drho_dn drho_dn
!!     dvxci(:,3)= dEx/d(abs(grad(rho_up))) / abs(grad(rho_up))
!!     dvxci(:,4)= dEx/d(abs(grad(rho_dn))) / abs(grad(rho_dn))
!!     dvxci(:,5)= d2Ex/d(abs(grad(rho_up))) drho_up / abs(grad(rho_up))
!!     dvxci(:,6)= d2Ex/d(abs(grad(rho_dn))) drho_dn / abs(grad(rho_dn))
!!     dvxci(:,7)= 1/abs(grad(rho_up)) * d/drho_up (dEx/d(abs(grad(rho_up))) /abs(grad(rho_up)))
!!     dvxci(:,8)= 1/abs(grad(rho_dn)) * d/drho_dn (dEx/d(abs(grad(rho_dn))) /abs(grad(rho_dn)))
!!     dvxci(:,9)= d2Ec/drho_up drho_up
!!     dvxci(:,10)=d2Ec/drho_up drho_dn
!!     dvxci(:,11)=d2Ec/drho_dn drho_dn
!!     dvxci(:,12)=dEc/d(abs(grad(rho))) / abs(grad(rho))
!!     dvxci(:,13)=d2Ec/d(abs(grad(rho))) drho_up / abs(grad(rho))
!!     dvxci(:,14)=d2Ec/d(abs(grad(rho))) drho_dn / abs(grad(rho))
!!     dvxci(:,15)=1/abs(grad(rho)) * d/drho (dEc/d(abs(grad(rho))) /abs(grad(rho)))
!!
!!  exci(npts)=exchange-correlation energy density (hartree)
!!  vxci(npts,nspden)=partial derivative of the exchange-correlation energy (exci*$\rho$)
!!    with respect to the spin-down (vxci(:,1)) and spin-up (vxci(:,2) densities
!! Normalization: Exc=$\int (exc(r)*\rho (r) d^3 r)$ for $\rho$(r)=electron density.
!!
!! TODO
!!  WARNING: option=4 not yet implemented.
!!
!! NOTES
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!      invcb,leave_new,wrtout
!!
!! SOURCE
!!$#if defined HAVE_CONFIG_H
!!$#include "config.h"
!!$#endif

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,nd2vxci, & !Mandatory Arguments
&                d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                          !Optional Arguments

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!takao comment out
!  use interfaces_14_hidewrite
!  use interfaces_16_hideleave
!  use interfaces_56_xc, except_this_one => xcpbe
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndvxci,nd2vxci,ngr2,npts,nspden,option,order
 integer,intent(in),optional :: exexch
!arrays
 real(dp),intent(in) :: rho_updn(npts,nspden)
 real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
 real(dp),intent(out) :: exci(npts),vxci(npts,nspden)
 real(dp),intent(out),optional :: d2vxci(npts,nd2vxci),dvxcdgr(npts,3)
 real(dp),intent(out),optional :: dvxci(npts,ndvxci)

!Local variables-------------------------------
! The "accurate" value of mu is taken from the PBE code.
! The value of mu_c09 is taken from the paper (see above) in order to recover
! the GEA behaviour of the enhancement factor for small values of s rather than
! the LSD linear response limit used in revPBE. 
!scalars
 integer,save :: initialized=0
 integer :: ipts,ispden
 real(dp),parameter :: alpha_c09=0.0483_dp
 real(dp),parameter :: alpha_zeta2=1.0_dp-1.0e-6_dp,alpha_zeta=1.0_dp-1.0e-6_dp
 real(dp),parameter :: b_wc=0.123456790123_dp,beta=0.066725_dp
 real(dp),parameter :: beta_inv=1.0_dp/beta,c_wc=0.00793746933516_dp
 real(dp),parameter :: fsec_inv=1.0_dp/1.709921_dp,kappa_pbe=0.804_dp
 real(dp),parameter :: kappa_revpbe=1.245_dp,mu=0.2195149727645171_dp
 real(dp),parameter :: kappa_c09=1.245_dp, mu_c09=0.0617_dp 
 real(dp),parameter :: mu_divkappa_pbe=mu/kappa_pbe
 real(dp),parameter :: mu_divkappa_revpbe=mu/kappa_revpbe
 real(dp),parameter :: rsfac=0.6203504908994000_dp,tolgrad=tol10
 real(dp),save :: beta_gamma,coeff_tt,factf_zeta,factfp_zeta,gamma,gamma_inv
 real(dp),save :: sixpi2_1_3,sixpi2m1_3,sq_rsfac,sq_rsfac_inv,threefourth_divpi
 real(dp),save :: twom1_3
 real(dp) :: aa,arg_rr,bb,cc,coeff_aa,alphs2,alphmu,coeff_qq,coeffss,d2aa_drs2,d2aa_drsdzeta
 real(dp) :: d2aa_dzeta2,d2bb_drs2,d2bb_drsdzeta,d2bb_dzeta2,d2cc_dbb2
 real(dp) :: d2cc_drs2,d2cc_drsdzeta,d2cc_dzeta2,d2ecrs0_drs2,d2ecrs1_drs2
 real(dp) :: d2ecrs_drdn2,d2ecrs_drdndrup,d2ecrs_drho2,d2ecrs_drs2
 real(dp) :: d2ecrs_drsdzeta,d2ecrs_drup2,d2ecrs_dzeta2,d2fxdg2,d2fxdn2
 real(dp) :: d2fxdndg,d2fxdss2,d2fzeta4_dzeta2,d2gcrs_drs2,d2hh_drs2
 real(dp) :: d2hh_drsdtt,d2hh_drsdzeta,d2hh_dtt2,d2hh_dttdzeta,d2hh_dzeta2
 real(dp) :: d2macrs_drs2,d2pade_drs2,d2pade_drsdtt,d2pade_drsdzeta,d2pade_dtt2
 real(dp) :: d2pade_dttdzeta,d2pade_dxx2,d2pade_dzeta2,d2qq_drs2,d2qq_drsdtt
 real(dp) :: d2qq_drsdzeta,d2qq_dtt2,d2qq_dttdzeta,d2qq_dzeta2,d2rhohh_drho2
 real(dp) :: d2rhohh_drhodg,d2rr_dqq2,d2rr_drs2,d2rr_drsdtt,d2rr_drsdzeta
 real(dp) :: d2rr_dtt2,d2rr_dttdzeta,d2rr_dzeta2,d2rs_dn2,d2ssdn2,d2ssdndg
 real(dp) :: d2vcrs_drs2,d2xx_drs2,d2xx_drsdtt,d2xx_drsdzeta,d2xx_dttdzeta
 real(dp) :: d2xx_dzeta2,d3ecrs0_drs3,d3ecrs_drup3,d3ecrs_drup2drdn
 real(dp) :: d3ecrs_drupdrdn2,d3ecrs_drdn3,d3ecrs_dzeta3
 real(dp) :: d3ecrs_drs2dzeta,d3ecrs_dzeta2drs,d3ecrs1_drs3,d3gcrs_drs3
 real(dp) :: d3ecrs_drs3,d3macrs_drs3
 real(dp) :: d_wc,daa_drs,daa_dzeta,dbb_drs,dbb_dzeta
 real(dp) :: dcc_dbb,dcc_drs,dcc_dzeta,decrs0_drs,decrs1_drs,decrs_drs
 real(dp) :: decrs_dzeta,dfxdg,dfxdn,dfxdss,dfzeta4_dzeta,dgcrs_drs
 real(dp) :: dhh_drs,dhh_dtt,dhh_dzeta,div_rr,divss,dmacrs_drs,dpade_drs
 real(dp) :: dpade_dtt,dpade_dxx,dpade_dzeta,dqq_drs,dqq_dtt,dqq_dzeta
 real(dp) :: drhohh_drho,drr_dqq,drr_drs,drr_dtt,drr_dzeta,drs_dn,dssdg,dssdn
 real(dp) :: dtt_dg,dvcrs_drs,dxx_drs,dxx_dtt,dxx_dzeta,ec0_a1,ec0_aa,ec0_b1
 real(dp) :: ec0_b2,ec0_b3,ec0_b4,ec0_den,ec0_f1,ec0_f2,ec0_log,ec0_q0,ec0_q1
 real(dp) :: ec0_q1p,ec0_q1pp,ec0_q1ppp,ec1_a1,ec1_aa,ec1_b1,ec1_b2,ec1_b3
 real(dp) :: ec1_b4,ec1_den,ec1_f1,ec1_f2,ec1_log,ec1_q0,ec1_q1,ec1_q1p,ec1_q1pp
 real(dp) :: ec1_q1ppp,ecrs,ecrs0,factfppp_zeta
 real(dp) :: ecrs1,ex_gga,ex_lsd,exc,exp_pbe,expss,f_zeta,factfpp_zeta
 real(dp) :: fp_zeta,fpp_zeta,fppp_zeta,fx,gamphi3inv,gcrs,grrho2,hh,kappa
 real(dp) :: mac_a1,mac_aa,mac_b1,mac_b2,mac_b3,mac_b4,mac_den,mac_f1,mac_f2,mac_log,mac_q0
 real(dp) :: mac_q1,mac_q1ppp
 real(dp) :: mac_q1p,mac_q1pp,macrs,mu_divkappa,p1_wc,p2_wc,pade,pade_den
 real(dp) :: phi3_zeta,phi_logder,phi_zeta,phi_zeta_inv,phip_zeta,phipp_zeta,qq
 real(dp) :: rho,rho_inv,rhomot
 real(dp) :: rhotmo6,rhotmot,rhoto6,rhotot,rhotot_inv,rr,rs,rsm1_2,sqr_rs
 real(dp) :: sqr_sqr_rs,ss,tt,vxcadd,xx,zeta,zeta4,zetm_1_3,zetp_1_3
 real(dp) :: a1fa,a2fa,b1fa,b2fa,c1fa,c2fa,e1fa,e2fa,f1fa,f2fa,g1fa,g2fa,h1fa,h2fa
 real(dp) :: i1fa,i2fa,m1fa,m2fa,n1fa,n2fa
 real(dp) :: sp1_up3,sp1_up2dn,sp1_updn2,sp1_dn3
 real(dp) :: sp2_up3,sp2_up2dn,sp2_updn2,sp2_dn3
 real(dp) :: sp3_up3,sp3_up2dn,sp3_updn2,sp3_dn3
 real(dp) :: d3ecrs_sp0,d3ecrs_sp1,d3ecrs_sp2,d3ecrs_sp3




 character(len=500) :: message
!arrays
 real(dp),allocatable :: rho_updnm1_3(:,:),rhoarr(:),rhom1_3(:),zetm(:)
 real(dp),allocatable :: zetmm1_3(:),zetp(:),zetpm1_3(:)
!no_abirules
!integer :: debug 
!real(dp) :: delta,factor,grr,rho_dn,rho_dnm,rho_dnp,rho_up,rho_upm,rho_upp,zeta_mean
!real(dp), allocatable :: wecrsz(:,:),d1wecrsz(:,:),d2wecrsz(:,:),d3wecrsz(:,:)
!real(dp) :: d3ecrs_drho3,d3ecrs_drhodndrho2,d3ecrs_drhoupdrho2
!real(dp) :: ec1_q0p,mac_q0p,sigma1,sigma2,sigma3

! *************************************************************************

!DEBUG
!write(6,*)' xcpbe : enter'
!ENDDEBUG

 d_wc=mu-b_wc

!DEBUG
!allocate(wecrsz(npts,8),d1wecrsz(npts,8),d2wecrsz(npts,8),d3wecrsz(npts,8))
!ENDDEBUG

 if (option<=-4 .or. option==0 .or. option==4 .or. option>=8 ) then
   write(message, '(a,a,a,a,i12,a)' ) ch10,&
&   ' xcpbe : BUG -',ch10,&
&   '  Option must be 1, 2, 3, 5, 6, 7, 8, -1 or -2 ; argument was ',option,'.'
!  call wrtout(std_out,message,'COLL')
!  call leave_new('COLL')
 end if
!Checks the compatibility between the presence of dvxci and ndvxci

!cccccccccccccccccccccccccccc
! print *,' present=', present(d2vxci),present(dvxcdgr),present(dvxci), &
!&  present(exexch),present(grho2_updn)
! print *,'xxxxxxx', ndvxci, present(dvxci)
!cccccccccccccccccccccccccc

 if(ndvxci /=0 .neqv. present(dvxci))then
   write(message, '(4a)' )ch10,&
&   ' xcpbe : BUG -',ch10,&
&   '  If ndvxci/=0 there must the optional arguments dvxci'
    call rx(message)
!   call wrtout(std_out,message,'COLL')
!   call leave_new('COLL')
 end if
!Checks the compatibility between the inputs and the presence of the optional arguments
 if(ndvxci /= 0 .and. abs(order) <= 1)then
   write(message, '(6a,i8,a)' )ch10,&
&   ' xcpbe : BUG -',ch10,&
&   '  The order does not require the presence of dvxci',ch10,&
&   '  that is allowed when |order|>1, while we have',&
&   order,'.'
!   call wrtout(std_out,message,'COLL')
!   call leave_new('COLL')
    call rx(message)
 end if
 if(ndvxci /= 0 .and. (&
& ((option == 1 .or. option == -1 .or. option == 3) .and. ndvxci /= nspden + 1)&
& .or. (option == -2 .and. ndvxci /= 8)&
& .or. ((option == 2 .or. option == 5 .or. option == 6 .or. option == 7 .or. option == 8) .and. ndvxci /= 15)&
& ))then
   write(message, '(16a,i8,a,i8,a,i8)' )ch10,&
&   ' xcpbe : BUG -',ch10,&
&   '  The option is not consistent with the value of ndvxci',ch10,&
&   '  Allowed values are:',ch10,&
&   '  ndvxci     option',ch10,&
&   ' nspden+1    1,-1,3',ch10,&
&   '    8          -2',ch10,&
&   '    15       2, 5,6,7',ch10,&
&   '  While we have: order=',order,'option=',option,'nspden=',nspden
!   call wrtout(std_out,message,'COLL')
!   call leave_new('COLL')
    call rx(message)
 end if
 if (present(grho2_updn)) then
   if (ngr2/=2*nspden-1 ) then
     write(message, '(4a)' ) ch10,&
&     ' xcpbe : BUG -',ch10,&
&     '  ngr2 must be 2*nspden-1 !'
!     call wrtout(std_out,message,'COLL')
!     call leave_new('COLL')
    call rx(message)
   end if
 end if

 if ((option == 1 .or. option == -1 .or. option ==3) .and.  (present(grho2_updn) .or. present(dvxcdgr))) then
   write(message, '(a,a,a,a,a,a,i6,a)' )ch10,&
&   ' xcpbe : BUG -',ch10,&
&   '  The option chosen does not need the presence',ch10,&
&   '  of the gradient, or of the array dvxcdgr in the input, needed if option/=1,-1,3 , while we have',&
&   option,'.'
!   call wrtout(std_out,message,'COLL')
!   call leave_new('COLL')
    call rx(message)
 end if
 if (order /= 3 .and. present(d2vxci)) then
   write(message, '(a,a,a,a,a,a,i6,a)' )ch10,&
&   ' xcpbe : BUG -',ch10,&
&   '  The order chosen does not need the presence',ch10,&
&   '  of the array d2vxci, needed if order=3 , while we have',&
&   order,'.'
!   call wrtout(std_out,message,'COLL')
!   call leave_new('COLL')
    call rx(message)
 end if

 if(initialized==0)then
   twom1_3=two**(-third)
   sixpi2_1_3=(six*pi**2)**third
   sixpi2m1_3=one/sixpi2_1_3
   threefourth_divpi=three_quarters*piinv
   gamma=(one-log(two))*piinv**2
   gamma_inv=one/gamma
   beta_gamma=beta*gamma_inv
   factf_zeta= one / ( two**(four/three)-two )
   factfp_zeta= four_thirds * factf_zeta * alpha_zeta2
   coeff_tt= one/ (four*four*piinv*(three*pi**2)**third)
!  coeff_tt= two * sqrt(four*piinv*(three*pi**2)**third)
   sq_rsfac=sqrt(rsfac)
   sq_rsfac_inv=one/sq_rsfac
   initialized=1
 end if

!Parameters for the Perdew-Wang 92 LSD as well as LSD-RPA,
!see Table I of Phys.Rev.B 45,13244 (1992)
 ec0_aa=0.031091_dp  ; ec1_aa=0.015545_dp ; mac_aa=0.016887_dp
 if(option/=3 .and. option/=4)then
   ec0_a1=0.21370_dp  ; ec1_a1=0.20548_dp  ; mac_a1=0.11125_dp
   ec0_b1=7.5957_dp   ; ec1_b1=14.1189_dp  ; mac_b1=10.357_dp
   ec0_b2=3.5876_dp   ; ec1_b2=6.1977_dp   ; mac_b2=3.6231_dp
   ec0_b3=1.6382_dp   ; ec1_b3=3.3662_dp   ; mac_b3=0.88026_dp
   ec0_b4=0.49294_dp  ; ec1_b4=0.62517_dp  ; mac_b4=0.49671_dp
 else  ! RPA values
   ec0_a1=0.082477_dp ; ec1_a1=0.035374_dp ; mac_a1=0.028829_dp
   ec0_b1=5.1486_dp   ; ec1_b1=6.4869_dp   ; mac_b1=10.357_dp
   ec0_b2=1.6483_dp   ; ec1_b2=1.3083_dp   ; mac_b2=3.6231_dp
   ec0_b3=0.23647_dp  ; ec1_b3=0.11518_dp  ; mac_b3=0.479_dp
   ec0_b4=0.20614_dp  ; ec1_b4=0.082349_dp ; mac_b4=0.112279_dp
 end if

 if(option/=5 .and. option/=-4)then
   kappa=kappa_pbe
   mu_divkappa=mu_divkappa_pbe
 end if
 if(option==5)then
   kappa=kappa_revpbe
   mu_divkappa=mu_divkappa_revpbe
 end if
 if(option==-4)then
   kappa=kappa_c09  
 end if
!DEBUG
!Finite-difference debugging, do not take away
!Note : here work with collinear gradients. Might be generalized ...
!debug=2  ! Choose 1 (rho grads) or 2 (grho grads)
!if(order==3)debug=1
!factor=1.0_dp
!zeta_mean=0.98_dp
!!zeta_mean=zero
!delta=0.000025*factor
!delta=0.0000125*factor
!if(debug/=0)then
!do ipts=1,npts-4,5
!rho=ipts*0.01_dp*factor
!rho_up=rho*(1.0_dp+zeta_mean)*0.5_dp
!rho_dn=rho*(1.0_dp-zeta_mean)*0.5_dp
!rho_upp=rho_up+delta
!rho_upm=rho_up-delta
!rho_dnp=rho_dn+delta
!rho_dnm=rho_dn-delta
!! Here, vary rho
!if(debug==1)then
!rho_updn(ipts  ,1)=rho_up ; rho_updn(ipts  ,2)=rho_dn
!rho_updn(ipts+1,1)=rho_upp; rho_updn(ipts+1,2)=rho_dn
!rho_updn(ipts+2,1)=rho_upm; rho_updn(ipts+2,2)=rho_dn
!rho_updn(ipts+3,1)=rho_up ; rho_updn(ipts+3,2)=rho_dnp
!rho_updn(ipts+4,1)=rho_up ; rho_updn(ipts+4,2)=rho_dnm
!grho2_updn(ipts:ipts+4,1)=(0.2_dp*factor)**2     ! grad2 of spin up density
!grho2_updn(ipts:ipts+4,2)=(0.2_dp*factor)**2     ! grad2 of spin down density
!grho2_updn(ipts:ipts+4,3)=(0.3_dp*factor)**2     ! grad2 of total density
!else
!!  Here, vary grho (interchange rho and grho)
!grho2_updn(ipts  ,1)=rho_up**2 ; grho2_updn(ipts  ,2)=rho_dn**2
!grho2_updn(ipts+1,1)=rho_upp**2; grho2_updn(ipts+1,2)=rho_dn**2
!grho2_updn(ipts+2,1)=rho_upm**2; grho2_updn(ipts+2,2)=rho_dn**2
!grho2_updn(ipts+3,1)=rho_up**2 ; grho2_updn(ipts+3,2)=rho_dnp**2
!grho2_updn(ipts+4,1)=rho_up**2 ; grho2_updn(ipts+4,2)=rho_dnm**2
!grho2_updn(ipts  ,3)=(ipts*0.01_dp*factor)**2
!grho2_updn(ipts+1,3)=(ipts*0.01_dp*factor+delta)**2
!grho2_updn(ipts+2,3)=(ipts*0.01_dp*factor-delta)**2
!grho2_updn(ipts+3,3)=(ipts*0.01_dp*factor+delta)**2   ! identical to ipts+1
!grho2_updn(ipts+4,3)=(ipts*0.01_dp*factor-delta)**2   ! identical to ipts+2
!rho_updn(ipts:ipts+4,1)=0.2_dp*factor*(1.0_dp+zeta_mean)*0.5_dp    ! spin up density
!rho_updn(ipts:ipts+4,2)=0.2_dp*factor*(1.0_dp-zeta_mean)*0.5_dp    ! spin down density
!end if
!end do
!end if
!Usual option :
!nspden=2 ; order=2
!GGA
!nspden=2 ; order=1
!Might take also, although finite difference later is meaningless
!nspden=1 ; order=-2
!Here, alternative specification, in terms of defined rs and zeta
!do ipts=1,5
!if(ipts==1)then ;rs=0.01_dp ; zeta=0.98_dp ; endif
!if(ipts==2)then ;rs=0.01_dp+delta ; zeta=0.98_dp ; endif
!if(ipts==3)then ;rs=0.01_dp-delta ; zeta=0.98_dp ; endif
!if(ipts==4)then ;rs=0.01_dp ; zeta=0.98_dp+delta ; endif
!if(ipts==5)then ;rs=0.01_dp ; zeta=0.98_dp-delta ; endif
!rho=(rsfac/rs)**3
!rho_up=rho*(1.0_dp+zeta)*0.5_dp
!rho_dn=rho*(1.0_dp-zeta)*0.5_dp
!rho_updn(ipts  ,1)=rho_up ; rho_updn(ipts  ,2)=rho_dn
!enddo
!ENDDEBUG

 if(order**2 >1)then
   factfpp_zeta= third * factfp_zeta * alpha_zeta2
 end if


 allocate(rhoarr(npts),rhom1_3(npts),rho_updnm1_3(npts,2))
 allocate(zetm(npts),zetmm1_3(npts),zetp(npts),zetpm1_3(npts))

 do ispden=1,nspden
   call invcb(rho_updn(:,ispden),rho_updnm1_3(:,ispden),npts)
 end do


 if(nspden==1)then
   rhoarr(:)=two*rho_updn(:,1)
   rhom1_3(:)=twom1_3*rho_updnm1_3(:,1)
   rho_updnm1_3(:,2)=rho_updnm1_3(:,1)
 else
   rhoarr(:)=rho_updn(:,1)+rho_updn(:,2)
   call invcb(rhoarr,rhom1_3,npts)
   do ipts=1,npts
     rhotmot=rhom1_3(ipts)
     rhotot_inv=rhotmot*rhotmot*rhotmot
     zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
     zetp(ipts)=1.0_dp+zeta*alpha_zeta
     zetm(ipts)=1.0_dp-zeta*alpha_zeta
   end do
   call invcb(zetp,zetpm1_3,npts)
   call invcb(zetm,zetmm1_3,npts)
 end if


!fab: eliminate the following restriction

!if (order==3 .and. nspden == 1) d2vxci(:,:)=0._dp


 if (order==3) d2vxci(:,:)=0._dp

!!!Loop unrolling summary
!Completely unrolled for spin non-polarized case
!To be optimized for spin-polarized cases
!The loops are unrolled as follows:
!nspden=1     (line 433)
!order^2<=1  (line 460)
!option=2,5 (line 462)
!option=6,7 (line 630)
!option=-1  (line 825)
!option=-2  (line 853)
!option=1   (line 904)
!option=3   (line 963)
!order=3     (line 1024)
!option=2,5
!option=6,7
!option=-1
!option=-2
!option=1
!option=3
!order=-2    (line 1983)
!option=2,5
!option=6,7
!option=-1
!option=-2
!option=1
!option=3
!order^2>1   (line 2875)
!option=2,5
!option=6,7
!option=-1
!option=-2
!option=1
!option=3
!nspden=2     (line 3750)
!order^2<=1  (line 3754)
!order^2>1 (with if statements inside distinguishing between order=3 or -2)   (line 4000)
!!!End loop unrolling summary

!we separate different cases, depending on nspden
 if (nspden==1) then
!  we separate different cases, depending on order
   if (order**2<=1) then
!    we separate different cases, depending on option
     if(option==2 .or. option==5)then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
         if(present(exexch))then
           if(exexch==1) cycle
         end if
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs


!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta

!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)

!        The GGA correlation energy is added
!ccccccccccccccccccccc
!	 print *,'hh=0 test'
!	 hh=0d0
!ccccccccccccccccccccc
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh -third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

!        End condition of GGA

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do
     else if((option==6) .or. (option==7)) then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss

!        This is RPBE modification
         if (option==6) then
           divss=exp(-mu_divkappa*ss)
           dfxdss= mu*divss
           d2fxdss2=-mu*mu_divkappa*divss

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
!          This is the Wu and Cohen modification
         else
           expss=exp(-ss)
           p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
           p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&           four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
           divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
           dfxdss=p1_wc*divss*divss
           d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
         end if

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
         if(present(exexch))then
           if(exexch==1) cycle
         end if
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs


!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta

!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)

!        The GGA correlation energy is added
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

!        End condition of GGA

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do


     else if (option==-1) then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
       end do

     else if(option==-2) then


       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do

!      !!!!!!!!!!!!!!!!!!!!!!!!  CODING STARTS HERE!!!!!!!!!!!!!!!!!!!
     else if(option==-4) then


       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        VALENTINO R. COOPER C09x GGA, This is an exchange term proposed 
!        to use together with vdw-DF (see above).
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
!        the quarter that is lacking is compensated by the grho2_updn in the
!        next line.
         ss=grho2_updn(ipts,ispden)*coeffss
         alphs2=alpha_c09*ss
         alphmu=alpha_c09*mu_c09
         dfxdss= mu_c09*exp(-alphs2)*(one-alphs2)+kappa*alpha_c09*exp(-alphs2/two)/two
         d2fxdss2=-alphmu*exp(-alphs2)*(two-alphs2)-kappa*(alpha_c09**two)*exp(alphs2/two)/four
         fx    = one+mu_c09*ss*exp(-alphs2)+kappa*(one-exp(-alphs2/two))
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn 
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn) !THIS IS NOT VALID
!        AS IT CORRESPONDS TO THE PBE POTENTIAL
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg !here also
         exc=exc+ex_gga*rho

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do

!      !!!!!!!!!!!!!!!!!!!!!!!!! END OF CODING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
       
     else if(option==1)then


       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

!        Correlation has been added
!        -----------------------------------------------------------------------------

       end do

     else if (option==3) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

!        Correlation has been added
!        -----------------------------------------------------------------------------

       end do


     end if

   else if (order==3) then
!    separate cases with respect to option
     if(option==2 .or. option==5) then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
         ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
         ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
         d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
&         ec0_q0*ec0_q1pp) - &
&         ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
&         ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
&         ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,9)=d2ecrs_drho2
         dvxci(ipts,10)=d2ecrs_drho2
         dvxci(ipts,11)=d2ecrs_drho2

!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
         d2bb_drs2=d2ecrs_drs2*gamphi3inv
         d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
         d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
&         12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta
         d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
         d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
         d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
         d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
         d2aa_drs2=coeff_aa*d2cc_drs2
         d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
         d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
&         coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa
         d2xx_drs2=d2aa_drs2*tt
         d2xx_drsdzeta=d2aa_drsdzeta*tt
         d2xx_drsdtt=daa_drs
         d2xx_dttdzeta=daa_dzeta
         d2xx_dzeta2=d2aa_dzeta2*tt

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta
         d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
         d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
         d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
         d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
         d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
         d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
         d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
         d2qq_drs2=coeff_qq*d2pade_drs2
         d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
         d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
         d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
         d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
&         coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
         d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
&         +six*pade*phi_logder*phi_logder            &
&         -two*pade*phi_zeta_inv*phipp_zeta)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta
         d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
         d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
         d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
         d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
         d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
         d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
         d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
         d2hh_drs2=phi3_zeta*d2rr_drs2
         d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
         d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
         d2hh_dtt2=phi3_zeta*d2rr_dtt2
         d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
         d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
&         six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
&         +three*phi_zeta*phi_zeta*rr*phipp_zeta


!        The GGA correlation energy is added
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)



!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

         d2rhohh_drho2=rhotot_inv*&
&         (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
&         +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
&         +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
&         +two*seven*third*tt*zeta*d2hh_dttdzeta)
         d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
&         -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

!        Component 12 : first derivative with respect to the gradient
!        of the density, div by the grad of the density
         dvxci(ipts,12)=dvxcdgr(ipts,3)
!        Components 9, 10 and 11 : second derivatives with respect to the spin-density
!        Note that there is already a contribution from LSDA
         dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one-two*zeta) &
&         -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
         dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
         dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one+two*zeta) &
&         +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
!        Components 13 and 14 : second derivatives with respect to spin density
!        and gradient, divided by the gradient
         dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
         dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
!        Component 15 : derivative of the (derivative wrt the gradient div by the grad),
!        divided by the grad
         dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


!        End condition of GGA

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do
     else if ((option==6) .or. (option==7)) then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss

         if (option==6) then
           divss=exp(-mu_divkappa*ss)
           dfxdss= mu*divss
           d2fxdss2=-mu*mu_divkappa*divss

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
!          This is the Wu and Cohen modification
         else
           expss=exp(-ss)
           p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
           p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&           four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
           divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
           dfxdss=p1_wc*divss*divss
           d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
         end if

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
         ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
         ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
         d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
&         ec0_q0*ec0_q1pp) - &
&         ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
&         ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
&         ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,9)=d2ecrs_drho2
         dvxci(ipts,10)=d2ecrs_drho2
         dvxci(ipts,11)=d2ecrs_drho2

!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
         d2bb_drs2=d2ecrs_drs2*gamphi3inv
         d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
         d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
&         12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta
         d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
         d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
         d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
         d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
         d2aa_drs2=coeff_aa*d2cc_drs2
         d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
         d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
&         coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa
         d2xx_drs2=d2aa_drs2*tt
         d2xx_drsdzeta=d2aa_drsdzeta*tt
         d2xx_drsdtt=daa_drs
         d2xx_dttdzeta=daa_dzeta
         d2xx_dzeta2=d2aa_dzeta2*tt

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta
         d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
         d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
         d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
         d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
         d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
         d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
         d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
         d2qq_drs2=coeff_qq*d2pade_drs2
         d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
         d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
         d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
         d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
&         coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
         d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
&         +six*pade*phi_logder*phi_logder            &
&         -two*pade*phi_zeta_inv*phipp_zeta)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta
         d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
         d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
         d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
         d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
         d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
         d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
         d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
         d2hh_drs2=phi3_zeta*d2rr_drs2
         d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
         d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
         d2hh_dtt2=phi3_zeta*d2rr_dtt2
         d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
         d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
&         six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
&         +three*phi_zeta*phi_zeta*rr*phipp_zeta


!        The GGA correlation energy is added
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

         d2rhohh_drho2=rhotot_inv*&
&         (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
&         +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
&         +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
&         +two*seven*third*tt*zeta*d2hh_dttdzeta)
         d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
&         -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

!        Component 12 : first derivative with respect to the gradient
!        of the density, div by the grad of the density
         dvxci(ipts,12)=dvxcdgr(ipts,3)
!        Components 9, 10 and 11 : second derivatives with respect to the spin-density
!        Note that there is already a contribution from LSDA
         dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one-two*zeta) &
&         -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
         dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
         dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one+two*zeta) &
&         +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
!        Components 13 and 14 : second derivatives with respect to spin density
!        and gradient, divided by the gradient
         dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
         dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
!        Component 15 : derivative of the (derivative wrt the gradient div by the grad),
!        divided by the grad
         dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


!        End condition of GGA


!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)
       end do

     else if (option==-1) then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho
!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
         if(order==3)then
!          Compute the second derivative of vx
!          vx^(2) = -2*vx^(1)/(3*rhotot)
           d2vxci(ipts,1) = -2._dp*dvxci(ipts,1)/(3._dp*rhotot)
         end if
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

!        Correlation has been added
!        -----------------------------------------------------------------------------

       end do
     else if (option==-2) then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)

!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do
     else if(option==1) then

       do ipts=1,npts


         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
         if(order==3)then
!          Compute the second derivative of vx
!          vx^(2) = -2*vx^(1)/(3*rhotot)
           d2vxci(ipts,1) = -2._dp*dvxci(ipts,1)/(3._dp*rhotot)
         end if
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
         ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
         ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
         d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
&         ec0_q0*ec0_q1pp) - &
&         ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
&         ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
&         ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
         d2vcrs_drs2 = third*(d2ecrs_drs2 - rs*d3ecrs0_drs3)
         drs_dn = -1._dp*four_pi*ninth*rs**4
         d2rs_dn2 = 64._dp*pi*pi*(rs**7)/81._dp
         if(order==3)then
           d2vxci(ipts,1) = d2vxci(ipts,1) + d2vcrs_drs2*drs_dn*drs_dn + &
&           dvcrs_drs*d2rs_dn2





         end if

!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.


!        Correlation has been added
!        -----------------------------------------------------------------------------
       end do
     else if (option==3) then

       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho


!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
         if(order==3)then
!          Compute the second derivative of vx
!          vx^(2) = -2*vx^(1)/(3*rhotot)
           d2vxci(ipts,1) = -2._dp*dvxci(ipts,1)/(3._dp*rhotot)
         end if
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
         ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
         ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
         d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
&         ec0_q0*ec0_q1pp) - &
&         ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
&         ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
&         ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs

         dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
         d2vcrs_drs2 = third*(d2ecrs_drs2 - rs*d3ecrs0_drs3)
         drs_dn = -1._dp*four_pi*ninth*rs**4
         d2rs_dn2 = 64._dp*pi*pi*(rs**7)/81._dp
         if(order==3)then
           d2vxci(ipts,1) = d2vxci(ipts,1) + d2vcrs_drs2*drs_dn*drs_dn + &
&           dvcrs_drs*d2rs_dn2



         end if


       end do

     end if

   else if(order==-2) then
     if(option==2 .or. option==5) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,9)=d2ecrs_drho2
         dvxci(ipts,10)=d2ecrs_drho2
         dvxci(ipts,11)=d2ecrs_drho2

!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
         d2bb_drs2=d2ecrs_drs2*gamphi3inv
         d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
         d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
&         12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta
         d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
         d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
         d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
         d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
         d2aa_drs2=coeff_aa*d2cc_drs2
         d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
         d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
&         coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa
         d2xx_drs2=d2aa_drs2*tt
         d2xx_drsdzeta=d2aa_drsdzeta*tt
         d2xx_drsdtt=daa_drs
         d2xx_dttdzeta=daa_dzeta
         d2xx_dzeta2=d2aa_dzeta2*tt

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta
         d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
         d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
         d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
         d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
         d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
         d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
         d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
         d2qq_drs2=coeff_qq*d2pade_drs2
         d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
         d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
         d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
         d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
&         coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
         d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
&         +six*pade*phi_logder*phi_logder            &
&         -two*pade*phi_zeta_inv*phipp_zeta)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta
         d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
         d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
         d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
         d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
         d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
         d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
         d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
         d2hh_drs2=phi3_zeta*d2rr_drs2
         d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
         d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
         d2hh_dtt2=phi3_zeta*d2rr_dtt2
         d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
         d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
&         six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
&         +three*phi_zeta*phi_zeta*rr*phipp_zeta

!        The GGA correlation energy is added
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

         d2rhohh_drho2=rhotot_inv*&
&         (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
&         +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
&         +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
&         +two*seven*third*tt*zeta*d2hh_dttdzeta)
         d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
&         -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

!        Component 12 : first derivative with respect to the gradient
!        of the density, div by the grad of the density
         dvxci(ipts,12)=dvxcdgr(ipts,3)
!        Components 9, 10 and 11 : second derivatives with respect to the spin-density
!        Note that there is already a contribution from LSDA
         dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one-two*zeta) &
&         -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
         dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
         dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one+two*zeta) &
&         +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
!        Components 13 and 14 : second derivatives with respect to spin density
!        and gradient, divided by the gradient
         dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
         dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
!        Component 15 : derivative of the (derivative wrt the gradient div by the grad),
!        divided by the grad
         dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


!        End condition of GGA


!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do

     else if ((option==6) .or. (option==7)) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss

         if (option==6) then
           divss=exp(-mu_divkappa*ss)
           dfxdss= mu*divss
           d2fxdss2=-mu*mu_divkappa*divss

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
!          This is the Wu and Cohen modification
         else
           expss=exp(-ss)
           p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
           p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&           four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
           divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
           dfxdss=p1_wc*divss*divss
           d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
         end if

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,9)=d2ecrs_drho2
         dvxci(ipts,10)=d2ecrs_drho2
         dvxci(ipts,11)=d2ecrs_drho2

!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
         d2bb_drs2=d2ecrs_drs2*gamphi3inv
         d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
         d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
&         12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta
         d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
         d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
         d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
         d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
         d2aa_drs2=coeff_aa*d2cc_drs2
         d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
         d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
&         coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa
         d2xx_drs2=d2aa_drs2*tt
         d2xx_drsdzeta=d2aa_drsdzeta*tt
         d2xx_drsdtt=daa_drs
         d2xx_dttdzeta=daa_dzeta
         d2xx_dzeta2=d2aa_dzeta2*tt

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta
         d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
         d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
         d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
         d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
         d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
         d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
         d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
         d2qq_drs2=coeff_qq*d2pade_drs2
         d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
         d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
         d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
         d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
&         coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
         d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
&         +six*pade*phi_logder*phi_logder            &
&         -two*pade*phi_zeta_inv*phipp_zeta)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta
         d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
         d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
         d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
         d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
         d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
         d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
         d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
         d2hh_drs2=phi3_zeta*d2rr_drs2
         d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
         d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
         d2hh_dtt2=phi3_zeta*d2rr_dtt2
         d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
         d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
&         six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
&         +three*phi_zeta*phi_zeta*rr*phipp_zeta

!        The GGA correlation energy is added
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

         d2rhohh_drho2=rhotot_inv*&
&         (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
&         +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
&         +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
&         +two*seven*third*tt*zeta*d2hh_dttdzeta)
         d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
&         -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

!        Component 12 : first derivative with respect to the gradient
!        of the density, div by the grad of the density
         dvxci(ipts,12)=dvxcdgr(ipts,3)
!        Components 9, 10 and 11 : second derivatives with respect to the spin-density
!        Note that there is already a contribution from LSDA
         dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one-two*zeta) &
&         -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
         dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
         dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one+two*zeta) &
&         +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
!        Components 13 and 14 : second derivatives with respect to spin density
!        and gradient, divided by the gradient
         dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
         dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
!        Component 15 : derivative of the (derivative wrt the gradient div by the grad),
!        divided by the grad
         dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


!        End condition of GGA


!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do

     else if (option==-1) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
!        Compute the second derivative of vx
!        vx^(2) = -2*vx^(1)/(3*rhotot)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv


!        Correlation has been added
!        -----------------------------------------------------------------------------


       end do

     else if (option==-2) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do

     else if (option==1) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho


!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
!        Compute the second derivative of vx
!        vx^(2) = -2*vx^(1)/(3*rhotot)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)

         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs

         dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
         dvxci(ipts,2)=dvxci(ipts,2)+d2ecrs_drho2-d2ecrs_dzeta2*rhotot_inv


!        Correlation has been added
!        -----------------------------------------------------------------------------

       end do

     else if (option==3) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
!        Compute the second derivative of vx
!        vx^(2) = -2*vx^(1)/(3*rhotot)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp


!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2
         dvxci(ipts,2)=dvxci(ipts,2)+d2ecrs_drho2-d2ecrs_dzeta2*rhotot_inv

!        Correlation has been added
!        -----------------------------------------------------------------------------

       end do

     end if


   else if (order**2>1) then
!    separate cases depending to option
     if(option==2 .or. option==5) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss

         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss

         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,9)=d2ecrs_drho2
         dvxci(ipts,10)=d2ecrs_drho2
         dvxci(ipts,11)=d2ecrs_drho2


!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
         d2bb_drs2=d2ecrs_drs2*gamphi3inv
         d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
         d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
&         12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta
         d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
         d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
         d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
         d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
         d2aa_drs2=coeff_aa*d2cc_drs2
         d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
         d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
&         coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa
         d2xx_drs2=d2aa_drs2*tt
         d2xx_drsdzeta=d2aa_drsdzeta*tt
         d2xx_drsdtt=daa_drs
         d2xx_dttdzeta=daa_dzeta
         d2xx_dzeta2=d2aa_dzeta2*tt

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta
         d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
         d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
         d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
         d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
         d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
         d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
         d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2


!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
         d2qq_drs2=coeff_qq*d2pade_drs2
         d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
         d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
         d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
         d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
&         coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
         d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
&         +six*pade*phi_logder*phi_logder            &
&         -two*pade*phi_zeta_inv*phipp_zeta)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta
         d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
         d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
         d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
         d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
         d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
         d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
         d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
         d2hh_drs2=phi3_zeta*d2rr_drs2
         d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
         d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
         d2hh_dtt2=phi3_zeta*d2rr_dtt2
         d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
         d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
&         six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
&         +three*phi_zeta*phi_zeta*rr*phipp_zeta


!        The GGA correlation energy is added
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)


!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

         d2rhohh_drho2=rhotot_inv*&
&         (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
&         +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
&         +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
&         +two*seven*third*tt*zeta*d2hh_dttdzeta)
         d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
&         -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

!        Component 12 : first derivative with respect to the gradient
!        of the density, div by the grad of the density
         dvxci(ipts,12)=dvxcdgr(ipts,3)
!        Components 9, 10 and 11 : second derivatives with respect to the spin-density
!        Note that there is already a contribution from LSDA
         dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one-two*zeta) &
&         -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
         dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
         dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one+two*zeta) &
&         +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
!        Components 13 and 14 : second derivatives with respect to spin density
!        and gradient, divided by the gradient
         dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
         dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
!        Component 15 : derivative of the (derivative wrt the gradient div by the grad),
!        divided by the grad
         dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


!        End condition of GGA

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do

     else if ((option==6) .or. (option==7)) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss

         if (option==6) then
           divss=exp(-mu_divkappa*ss)
           dfxdss= mu*divss
           d2fxdss2=-mu*mu_divkappa*divss

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
!          This is the Wu and Cohen modification
         else
           expss=exp(-ss)
           p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
           p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&           four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
           divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
           dfxdss=p1_wc*divss*divss
           d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa

           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
         end if

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,9)=d2ecrs_drho2
         dvxci(ipts,10)=d2ecrs_drho2
         dvxci(ipts,11)=d2ecrs_drho2


!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         phi_zeta=1.0_dp
         phip_zeta=0.0_dp
         phi_zeta_inv=1.0_dp
         phi_logder=0.0_dp
         phi3_zeta=1.0_dp
         gamphi3inv=gamma_inv
         phipp_zeta=-two*ninth*alpha_zeta*alpha_zeta

!        From ec to bb
         bb=ecrs*gamphi3inv
         dbb_drs=decrs_drs*gamphi3inv
         dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
         d2bb_drs2=d2ecrs_drs2*gamphi3inv
         d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
         d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
&         12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

!        From bb to cc
         exp_pbe=exp(-bb)
         cc=one/(exp_pbe-one)
         dcc_dbb=cc*cc*exp_pbe
         dcc_drs=dcc_dbb*dbb_drs
         dcc_dzeta=dcc_dbb*dbb_dzeta
         d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
         d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
         d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
         d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

!        From cc to aa
         coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
         aa=coeff_aa*cc
         daa_drs=coeff_aa*dcc_drs
         daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
         d2aa_drs2=coeff_aa*d2cc_drs2
         d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
         d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
&         coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

!        Introduce tt : do not assume that the spin-dependent gradients are collinear
         grrho2=four*grho2_updn(ipts,1)
         dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!        Note that tt is (the t variable of PBE divided by phi) squared
         tt=half*grrho2*dtt_dg

!        Get xx from aa and tt
         xx=aa*tt
         dxx_drs=daa_drs*tt
         dxx_dzeta=daa_dzeta*tt
         dxx_dtt=aa
         d2xx_drs2=d2aa_drs2*tt
         d2xx_drsdzeta=d2aa_drsdzeta*tt
         d2xx_drsdtt=daa_drs
         d2xx_dttdzeta=daa_dzeta
         d2xx_dzeta2=d2aa_dzeta2*tt

!        From xx to pade
         pade_den=one/(one+xx*(one+xx))
         pade=(one+xx)*pade_den
         dpade_dxx=-xx*(two+xx)*pade_den**2
         dpade_drs=dpade_dxx*dxx_drs
         dpade_dtt=dpade_dxx*dxx_dtt
         dpade_dzeta=dpade_dxx*dxx_dzeta
         d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
         d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
         d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
         d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
         d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
         d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
         d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2


!        From pade to qq
         coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
         qq=coeff_qq*pade
         dqq_drs=coeff_qq*dpade_drs
         dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
         dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
         d2qq_drs2=coeff_qq*d2pade_drs2
         d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
         d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
         d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
         d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
&         coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
         d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
&         +six*pade*phi_logder*phi_logder            &
&         -two*pade*phi_zeta_inv*phipp_zeta)

!        From qq to rr
         arg_rr=one+beta*gamma_inv*qq
         div_rr=one/arg_rr
         rr=gamma*log(arg_rr)
         drr_dqq=beta*div_rr
         drr_drs=drr_dqq*dqq_drs
         drr_dtt=drr_dqq*dqq_dtt
         drr_dzeta=drr_dqq*dqq_dzeta
         d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
         d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
         d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
         d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
         d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
         d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
         d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

!        From rr to hh
         hh=phi3_zeta*rr
         dhh_drs=phi3_zeta*drr_drs
         dhh_dtt=phi3_zeta*drr_dtt
         dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
         d2hh_drs2=phi3_zeta*d2rr_drs2
         d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
         d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
         d2hh_dtt2=phi3_zeta*d2rr_dtt2
         d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
         d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
&         six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
&         +three*phi_zeta*phi_zeta*rr*phipp_zeta


!        The GGA correlation energy is added
         exci(ipts)=exci(ipts)+hh

!        Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)


!        From hh to the derivative of the energy wrt the density
         drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
         vxci(ipts,1)=vxci(ipts,1)+drhohh_drho

!        From hh to the derivative of the energy wrt to the gradient of the
!        density, divided by the gradient of the density
!        (The v3.3 definition includes the division by the norm of the gradient)
         dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

         d2rhohh_drho2=rhotot_inv*&
&         (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
&         +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
&         +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
&         +two*seven*third*tt*zeta*d2hh_dttdzeta)
         d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
&         -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

!        Component 12 : first derivative with respect to the gradient
!        of the density, div by the grad of the density
         dvxci(ipts,12)=dvxcdgr(ipts,3)
!        Components 9, 10 and 11 : second derivatives with respect to the spin-density
!        Note that there is already a contribution from LSDA
         dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one-two*zeta) &
&         -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
         dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
         dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
&         ( d2hh_dzeta2*(one+two*zeta) &
&         +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
!        Components 13 and 14 : second derivatives with respect to spin density
!        and gradient, divided by the gradient
         dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
         dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
!        Component 15 : derivative of the (derivative wrt the gradient div by the grad),
!        divided by the grad
         dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg


!        End condition of GGA

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do

     else if (option==-1) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
!        Compute the second derivative of vx
!        vx^(2) = -2*vx^(1)/(3*rhotot)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

       end do

     else if (option==-2) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
         dvxcdgr(ipts,3)=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew-Burke-Ernzerhof GGA, exchange part
         rho_inv=rhomot*rhomot*rhomot
         coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
         ss=grho2_updn(ipts,ispden)*coeffss
         divss=one/(one+mu_divkappa*ss)
         dfxdss= mu*divss*divss
         d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
         fx    = one+kappa*(one-divss)
         ex_gga= ex_lsd*fx
         dssdn=-eight*third*ss*rho_inv
         dfxdn  = dfxdss*dssdn
         vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!        The new definition (v3.3) includes the division by the norm of the gradient
         dssdg =two*coeffss
         dfxdg=dfxdss*dssdg
         dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
         exc=exc+ex_gga*rho

!        Perdew-Burke-Ernzerhof GGA, exchange part
!        Components 3 or 4
         dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!        Components 1 or 2
         d2ssdn2=-11.0_dp*third*dssdn*rho_inv
         d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
         dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&         ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!        Components 5 or 6
         d2ssdndg=-eight*third*dssdg*rho_inv
         d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
         dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!        Components 7 or 8
         d2fxdg2=d2fxdss2*dssdg**2
         dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!        For the time being, treat non-spin-polarized like spin-polarized
         dvxci(ipts,2)=dvxci(ipts,1)
         dvxci(ipts,4)=dvxci(ipts,3)
         dvxci(ipts,6)=dvxci(ipts,5)
         dvxci(ipts,8)=dvxci(ipts,7)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional

!        Correlation has been added
!        -----------------------------------------------------------------------------

!        vxci(ipts,2)=vxci(ipts,1)
         dvxcdgr(ipts,2)=dvxcdgr(ipts,1)

       end do
     else if (option==1) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
!        Compute the second derivative of vx
!        vx^(2) = -2*vx^(1)/(3*rhotot)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs

         dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2

       end do
     else if (option==3) then
       do ipts=1,npts

         rhotot=rhoarr(ipts)
         rhotmot=rhom1_3(ipts)
         rhotot_inv=rhotmot*rhotmot*rhotmot
         rhotmo6=sqrt(rhotmot)
         rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!        -----------------------------------------------------------------------
!        First take care of the exchange part of the functional

         exc=zero
!        loop over the spin
         ispden=1
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
!        Perdew_Wang 91 LSD
         vxci(ipts,ispden)=four_thirds*ex_lsd
         exc=exc+ex_lsd*rho

!        Perdew_Wang 91 LSD
         dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&         threefourth_divpi*sixpi2_1_3*rhomot*rhomot
         dvxci(ipts,2)=zero
!        If non-spin-polarized, first component of dvxci is second
!        derivative with respect to TOTAL density.
         dvxci(ipts,1)=dvxci(ipts,1)*half
!        Compute the second derivative of vx
!        vx^(2) = -2*vx^(1)/(3*rhotot)
!        end of loop over the spin
!        If non spin-polarized, treat spin down contribution now, similar to spin up
         exc=exc*2
         exci(ipts)=exc*rhotot_inv
!        -----------------------------------------------------------------------------
!        Then takes care of the LSD correlation part of the functional


         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
         ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
         ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den
         ecrs=ecrs0
         decrs_drs=decrs0_drs
         decrs_dzeta=0.0_dp
         d2ecrs_drs2=d2ecrs0_drs2
         d2ecrs_dzeta2=alpha_zeta**2*(-macrs)
         d2ecrs_drsdzeta=zero
         zeta=0.0_dp

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxci(ipts,1)=vxci(ipts,1)+ecrs-rs*third*decrs_drs

         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drho2

       end do
     end if

   end if


!  fab: here it starts the spin polarized case

 else if(nspden==2) then

!  we separate different cases depending on order

   if (order**2<=1) then

     do ipts=1,npts

       rhotot=rhoarr(ipts)
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       rhotmo6=sqrt(rhotmot)
       rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!      -----------------------------------------------------------------------
!      First take care of the exchange part of the functional

       exc=zero
       if (present(dvxcdgr)) dvxcdgr(ipts,3)=zero
       do ispden=1,nspden
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
         if(option==1 .or. option==-1 .or. option==3)then
!          Perdew_Wang 91 LSD
           vxci(ipts,ispden)=four_thirds*ex_lsd
           if(present(dvxcdgr)) dvxcdgr(ipts,ispden)=0.0_dp
           exc=exc+ex_lsd*rho
         else          
           rho_inv=rhomot*rhomot*rhomot
           coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
           ss=grho2_updn(ipts,ispden)*coeffss
           if(option==7) then ! This is WC          
             expss=exp(-ss)
             p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
             p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&             four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
             divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
             dfxdss=p1_wc*divss*divss
             d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa
           else
             if(option/=6)then ! This is Perdew-Burke-Ernzerhof GGA, exchange part
               divss=one/(one+mu_divkappa*ss)
               dfxdss= mu*divss*divss
               d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             else  ! This is RPBE modification
               divss=exp(-mu_divkappa*ss)
               dfxdss= mu*divss
               d2fxdss2=-mu*mu_divkappa*divss
             end if
           end if
           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
         end if

       end do
       exci(ipts)=exc*rhotot_inv
       if(present(exexch)) then
         if(exexch==1)cycle
       end if

!      -----------------------------------------------------------------------------
!      Then takes care of the LSD correlation part of the functional

       if(option>0)then

         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         if(option/=3 .and. option/=4)then
           ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
           ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         else
           sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
           ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
           ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
         end if
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
         ec1_q0=-2.0_dp*ec1_aa*(1.0_dp+ec1_a1*rs)
         if(option/=3 .and. option/=4)then
           ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
           ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+4._dp*ec1_b4*rs)
         else
           ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs/sqr_sqr_rs)
           ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+3.5_dp*ec1_b4*rs/sqr_sqr_rs)
         end if
         ec1_den=1.0_dp/(ec1_q1*ec1_q1+ec1_q1)
!        ec1_log=log( 1.0_dp + 1.0_dp / ec1_q1 )
         ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
         ecrs1=ec1_q0*ec1_log
         decrs1_drs= -2.0_dp*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

!        alpha_zeta is introduced in order to remove singularities for fully
!        polarized systems.
         zetp_1_3=(1.0_dp+zeta*alpha_zeta)*zetpm1_3(ipts)**2
         zetm_1_3=(1.0_dp-zeta*alpha_zeta)*zetmm1_3(ipts)**2

         f_zeta=( (1.0_dp+zeta*alpha_zeta2)*zetp_1_3 +                      &
&         (1.0_dp-zeta*alpha_zeta2)*zetm_1_3 - 2.0_dp ) * factf_zeta
         fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
         zeta4=zeta**4

         gcrs=ecrs1-ecrs0+macrs*fsec_inv
!        ecrs=ecrs0+f_zeta*(-macrs*(1.0_dp-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
         ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

         dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
!        decrs_drs=decrs0_drs+f_zeta*&
!        &        (-dmacrs_drs*(1.0_dp-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
         decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
         dfzeta4_dzeta=4.0_dp*zeta**3*f_zeta+fp_zeta*zeta4
         decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
         vxci(ipts,1)=vxci(ipts,1)+vxcadd+decrs_dzeta
         vxci(ipts,2)=vxci(ipts,2)+vxcadd-decrs_dzeta

!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         if(option==2 .or. option==5 .or. option==6 .or. option==7)then
!          The definition of phi has been slightly changed, because
!          the original PBE one gives divergent behaviour for fully
!          polarized points
!          zetpm1_3=(1.0_dp+zeta*alpha_zeta)**(-third)
!          zetmm1_3=(1.0_dp-zeta*alpha_zeta)**(-third)
           phi_zeta=( zetpm1_3(ipts)*(1.0_dp+zeta*alpha_zeta)+ &
&           zetmm1_3(ipts)*(1.0_dp-zeta*alpha_zeta)   )*0.5_dp
           phip_zeta=(zetpm1_3(ipts)-zetmm1_3(ipts))*third*alpha_zeta
           phi_zeta_inv=1.0_dp/phi_zeta
           phi_logder=phip_zeta*phi_zeta_inv
           phi3_zeta=phi_zeta*phi_zeta*phi_zeta
           gamphi3inv=gamma_inv*phi_zeta_inv*phi_zeta_inv*phi_zeta_inv

!          From ec to bb
           bb=ecrs*gamphi3inv
           dbb_drs=decrs_drs*gamphi3inv
           dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
!          From bb to cc
           exp_pbe=exp(-bb)
           cc=one/(exp_pbe-one)
           dcc_dbb=cc*cc*exp_pbe
           dcc_drs=dcc_dbb*dbb_drs
           dcc_dzeta=dcc_dbb*dbb_dzeta

!          From cc to aa
           coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
           aa=coeff_aa*cc
           daa_drs=coeff_aa*dcc_drs
           daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
!          Introduce tt : do not assume that the spin-dependent gradients are collinear
           grrho2=grho2_updn(ipts,3)
           dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!          Note that tt is (the t variable of PBE divided by phi) squared
           tt=half*grrho2*dtt_dg

!          Get xx from aa and tt
           xx=aa*tt
           dxx_drs=daa_drs*tt
           dxx_dzeta=daa_dzeta*tt
           dxx_dtt=aa
!          From xx to pade
           pade_den=one/(one+xx*(one+xx))
           pade=(one+xx)*pade_den
           dpade_dxx=-xx*(two+xx)*pade_den**2
           dpade_drs=dpade_dxx*dxx_drs
           dpade_dtt=dpade_dxx*dxx_dtt
           dpade_dzeta=dpade_dxx*dxx_dzeta

!          From pade to qq
           coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
           qq=coeff_qq*pade
           dqq_drs=coeff_qq*dpade_drs
           dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
           dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)

!          From qq to rr
           arg_rr=one+beta*gamma_inv*qq
           div_rr=one/arg_rr
           rr=gamma*log(arg_rr)
           drr_dqq=beta*div_rr
           drr_drs=drr_dqq*dqq_drs
           drr_dtt=drr_dqq*dqq_dtt
           drr_dzeta=drr_dqq*dqq_dzeta

!          From rr to hh
           hh=phi3_zeta*rr
           dhh_drs=phi3_zeta*drr_drs
           dhh_dtt=phi3_zeta*drr_dtt
           dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)

!          The GGA correlation energy is added
           exci(ipts)=exci(ipts)+hh

!          Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)

!          From hh to the derivative of the energy wrt the density
           drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
           vxci(ipts,1)=vxci(ipts,1)+drhohh_drho+dhh_dzeta
           vxci(ipts,2)=vxci(ipts,2)+drhohh_drho-dhh_dzeta


!          From hh to the derivative of the energy wrt to the gradient of the
!          density, divided by the gradient of the density
!          (The v3.3 definition includes the division by the norm of the gradient)
           dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

!          End condition of GGA
         end if

       else  ! no correlation

!        End condition of including correlation, and not only exchange
       end if

!      Correlation has been added
!      -----------------------------------------------------------------------------

     end do

!    fab: the following is the "else" on order !!!

   else

     do ipts=1,npts

       rhotot=rhoarr(ipts)
       rhotmot=rhom1_3(ipts)
       rhotot_inv=rhotmot*rhotmot*rhotmot
       rhotmo6=sqrt(rhotmot)
       rhoto6=rhotot*rhotmot*rhotmot*rhotmo6
!      -----------------------------------------------------------------------
!      First take care of the exchange part of the functional

       exc=zero
       if (present(dvxcdgr)) dvxcdgr(ipts,3)=zero
       do ispden=1,nspden
         rho   =rho_updn(ipts,ispden)
         rhomot=rho_updnm1_3(ipts,ispden)
         ex_lsd= - threefourth_divpi * sixpi2_1_3*rhomot*rhomot*rho
         if(option==1 .or. option==-1 .or. option==3)then
!          Perdew_Wang 91 LSD
           vxci(ipts,ispden)=four_thirds*ex_lsd
           if(present(dvxcdgr)) dvxcdgr(ipts,ispden)=0.0_dp
           exc=exc+ex_lsd*rho
         else          
           rho_inv=rhomot*rhomot*rhomot
           coeffss=quarter*sixpi2m1_3*sixpi2m1_3*rho_inv*rho_inv*rhomot*rhomot
           ss=grho2_updn(ipts,ispden)*coeffss       
           if(option==7) then ! This is WC          
             expss=exp(-ss)
             p1_wc=b_wc+(mu-b_wc)*(one-ss)*expss+two*c_wc*ss/(one+c_wc*ss*ss)
             p2_wc=d_wc*(ss-two)*expss+two*c_wc/(one+c_wc*ss*ss)-&
&             four*c_wc*c_wc*ss*ss/((one+c_wc*ss*ss)*(one+c_wc*ss*ss))
             divss=one/(one+(b_wc*ss+d_wc*ss*expss+log(one+c_wc*ss*ss))/kappa)
             dfxdss=p1_wc*divss*divss
             d2fxdss2=p2_wc*divss*divss-two*divss*divss*divss*p1_wc*p1_wc/kappa
           else
             if(option/=6)then ! This is Perdew-Burke-Ernzerhof GGA, exchange part
               divss=one/(one+mu_divkappa*ss)
               dfxdss= mu*divss*divss
               d2fxdss2=-mu*two*mu_divkappa*divss*divss*divss
             else  ! This is RPBE modification
               divss=exp(-mu_divkappa*ss)
               dfxdss= mu*divss
               d2fxdss2=-mu*mu_divkappa*divss
             end if
           end if
           fx    = one+kappa*(one-divss)
           ex_gga= ex_lsd*fx
           dssdn=-eight*third*ss*rho_inv
           dfxdn  = dfxdss*dssdn
           vxci(ipts,ispden)=ex_lsd*(four_thirds*fx+rho*dfxdn)
!          The new definition (v3.3) includes the division by the norm of the gradient
           dssdg =two*coeffss
           dfxdg=dfxdss*dssdg
           dvxcdgr(ipts,ispden)=ex_lsd*rho*dfxdg
           exc=exc+ex_gga*rho
         end if

         if(option==1 .or. option==-1 .or. option==3)then

!          Perdew_Wang 91 LSD
           dvxci(ipts,2*ispden-1)=-four_thirds*third*&
&           threefourth_divpi*sixpi2_1_3*rhomot*rhomot
           dvxci(ipts,2)=zero
           if(order==3)then
!            If non-spin-polarized, first component of dvxci is second
!            derivative with respect to TOTAL density.
!            Compute the second derivative of vx
!            vx^(2) = -2*vx^(1)/(3*rhotot)

!            fab: third order derivatives of the exchange part in the spin polarized case

             d2vxci(ipts,3*ispden-2) = -2._dp*dvxci(ipts,2*ispden-1)*(rhomot*rhomot*rhomot)/3._dp

!            mixed thir order derivatives of the exchange energy with respect to rho of 
!            different spin polarization are zero 
             d2vxci(ipts,2)=zero
             d2vxci(ipts,3)=zero

           end if

         else
!          Perdew-Burke-Ernzerhof GGA, exchange part
!          Components 3 or 4
           dvxci(ipts,2+ispden)=dvxcdgr(ipts,ispden)
!          Components 1 or 2
           d2ssdn2=-11.0_dp*third*dssdn*rho_inv
           d2fxdn2=d2fxdss2*dssdn**2+dfxdss*d2ssdn2
           dvxci(ipts,ispden)=third*rho_inv*vxci(ipts,ispden)+&
&           ex_lsd*(seven*third*dfxdn+rho*d2fxdn2)
!          Components 5 or 6
           d2ssdndg=-eight*third*dssdg*rho_inv
           d2fxdndg=d2fxdss2*dssdn*dssdg+dfxdss*d2ssdndg
           dvxci(ipts,4+ispden)=ex_lsd*(four_thirds*dfxdg+rho*d2fxdndg)
!          Components 7 or 8
           d2fxdg2=d2fxdss2*dssdg**2
           dvxci(ipts,6+ispden)=ex_lsd*rho*d2fxdg2
!          For the time being, treat non-spin-polarized like spin-polarized
         end if
       end do
       exci(ipts)=exc*rhotot_inv
!      -----------------------------------------------------------------------------
!      Then takes care of the LSD correlation part of the functional

       if(option>0)then

         rs=rsfac*rhotmot
         sqr_rs=sq_rsfac*rhotmo6
         rsm1_2=sq_rsfac_inv*rhoto6

!        Formulas A6-A8 of PW92LSD
         ec0_q0=-2.0_dp*ec0_aa*(1.0_dp+ec0_a1*rs)
         if(option/=3 .and. option/=4)then
           ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs)
           ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+4._dp*ec0_b4*rs)
         else
           sqr_sqr_rs=max(1.e-15_dp,sqrt(sqr_rs))
           ec0_q1=2.0_dp*ec0_aa*(ec0_b1*sqr_rs+ec0_b2*rs+ec0_b3*rs*sqr_rs+ec0_b4*rs*rs/sqr_sqr_rs)
           ec0_q1p=ec0_aa*(ec0_b1*rsm1_2+2._dp*ec0_b2+3._dp*ec0_b3*sqr_rs+3.5_dp*ec0_b4*rs/sqr_sqr_rs)
         end if
         ec0_den=1.0_dp/(ec0_q1*ec0_q1+ec0_q1)
!        ec0_log=log( 1.0_dp + 1.0_dp / ec0_q1 )
         ec0_log=-log( ec0_q1*ec0_q1*ec0_den )
         ecrs0=ec0_q0*ec0_log
         decrs0_drs= -2.0_dp*ec0_aa*ec0_a1*ec0_log - ec0_q0*ec0_q1p *ec0_den
         ec0_q1pp=0.5_dp*ec0_aa*(-ec0_b1*rsm1_2**3+3._dp*ec0_b3*rsm1_2+8._dp*ec0_b4)
         d2ecrs0_drs2= 4.0_dp*ec0_aa*ec0_a1*ec0_q1p*ec0_den            &
&         -ec0_q0*ec0_q1pp*ec0_den                        &
&         +ec0_q0*ec0_q1p**2*ec0_den**2*(2._dp*ec0_q1+1.0_dp)
         if (order==3) then
           ec0_q1ppp = 0.75_dp*ec0_aa*(rsm1_2**5)*(ec0_b1-ec0_b3*rs)
           ec0_f1 = 1._dp/(ec0_q1*ec0_q1*(1._dp + ec0_q1))
           ec0_f2 = 1._dp/(ec0_q1*(1+ec0_q1))
           d3ecrs0_drs3 = 6._dp*ec0_q1p*ec0_f1*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + &
&           ec0_q0*ec0_q1pp) - &
&           ec0_f2*(-6._dp*ec0_aa*ec0_a1*ec0_q1pp + ec0_q0*ec0_q1ppp + &
&           ec0_f2*(3._dp*ec0_q1p*(-2._dp*ec0_aa*ec0_a1*ec0_q1p + ec0_q0*ec0_q1pp) + &
&           ec0_f2*2._dp*ec0_q0*(ec0_q1p**3)*(1._dp + 3._dp*ec0_q1*(1._dp + ec0_q1))))
         end if

         mac_q0=-2.0_dp*mac_aa*(1.0_dp+mac_a1*rs)
         mac_q1=2.0_dp*mac_aa*(mac_b1*sqr_rs+mac_b2*rs+mac_b3*rs*sqr_rs+mac_b4*rs*rs)
         mac_q1p=mac_aa*(mac_b1*rsm1_2+2._dp*mac_b2+3._dp*mac_b3*sqr_rs+4._dp*mac_b4*rs)
         mac_den=1.0_dp/(mac_q1*mac_q1+mac_q1)
         mac_log=-log( mac_q1*mac_q1*mac_den )
         macrs=mac_q0*mac_log
         dmacrs_drs= -2.0_dp*mac_aa*mac_a1*mac_log - mac_q0*mac_q1p*mac_den

         zeta=(rho_updn(ipts,1)-rho_updn(ipts,2))*rhotot_inv
         ec1_q0=-2.0_dp*ec1_aa*(1.0_dp+ec1_a1*rs)
         if(option/=3 .and. option/=4)then
           ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs)
           ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+4._dp*ec1_b4*rs)
         else
           ec1_q1=2.0_dp*ec1_aa*(ec1_b1*sqr_rs+ec1_b2*rs+ec1_b3*rs*sqr_rs+ec1_b4*rs*rs/sqr_sqr_rs)
           ec1_q1p=ec1_aa*(ec1_b1*rsm1_2+2._dp*ec1_b2+3._dp*ec1_b3*sqr_rs+3.5_dp*ec1_b4*rs/sqr_sqr_rs)
         end if
         ec1_den=1.0_dp/(ec1_q1*ec1_q1+ec1_q1)
!        ec1_log=log( 1.0_dp + 1.0_dp / ec1_q1 )
         ec1_log=-log( ec1_q1*ec1_q1*ec1_den )
         ecrs1=ec1_q0*ec1_log
         decrs1_drs= -2.0_dp*ec1_aa*ec1_a1*ec1_log - ec1_q0*ec1_q1p *ec1_den

!        alpha_zeta is introduced in order to remove singularities for fully
!        polarized systems.
         zetp_1_3=(1.0_dp+zeta*alpha_zeta)*zetpm1_3(ipts)**2
         zetm_1_3=(1.0_dp-zeta*alpha_zeta)*zetmm1_3(ipts)**2

         f_zeta=( (1.0_dp+zeta*alpha_zeta2)*zetp_1_3 +                      &
&         (1.0_dp-zeta*alpha_zeta2)*zetm_1_3 - 2.0_dp ) * factf_zeta
         fp_zeta=( zetp_1_3 - zetm_1_3 ) * factfp_zeta
         zeta4=zeta**4

         gcrs=ecrs1-ecrs0+macrs*fsec_inv
!        ecrs=ecrs0+f_zeta*(-macrs*(1.0_dp-zeta4)*fsec_inv+(ecrs1-ecrs0)*zeta4)
         ecrs=ecrs0+f_zeta*(zeta4*gcrs-macrs*fsec_inv)

         dgcrs_drs=decrs1_drs-decrs0_drs+dmacrs_drs*fsec_inv
!        decrs_drs=decrs0_drs+f_zeta*&
!        &        (-dmacrs_drs*(1.0_dp-zeta4)*fsec_inv+(decrs1_drs-decrs0_drs)*zeta4)
         decrs_drs=decrs0_drs+f_zeta*(zeta4*dgcrs_drs-dmacrs_drs*fsec_inv)
         dfzeta4_dzeta=4.0_dp*zeta**3*f_zeta+fp_zeta*zeta4
         decrs_dzeta=dfzeta4_dzeta*gcrs-fp_zeta*macrs*fsec_inv

         ec1_q1pp=0.5_dp*ec1_aa*(-ec1_b1*rsm1_2**3+3._dp*ec1_b3*rsm1_2+8._dp*ec1_b4)

         d2ecrs1_drs2= 4.0_dp*ec1_aa*ec1_a1*ec1_q1p*ec1_den            &
&         -ec1_q0*ec1_q1pp*ec1_den                        &
&         +ec1_q0*ec1_q1p**2*ec1_den**2*(2._dp*ec1_q1+1.0_dp)


         mac_q1pp=0.5_dp*mac_aa*(-mac_b1*rsm1_2**3+3._dp*mac_b3*rsm1_2+8._dp*mac_b4)
         d2macrs_drs2= 4.0_dp*mac_aa*mac_a1*mac_q1p*mac_den            &
&         -mac_q0*mac_q1pp*mac_den                        &
&         +mac_q0*mac_q1p**2*mac_den**2*(2._dp*mac_q1+1.0_dp)

         d2gcrs_drs2=d2ecrs1_drs2-d2ecrs0_drs2+d2macrs_drs2*fsec_inv
         fpp_zeta=(zetpm1_3(ipts)**2+zetmm1_3(ipts)**2) * factfpp_zeta
         d2fzeta4_dzeta2=12.0_dp*zeta**2*f_zeta  &
&         + 8.0_dp*zeta**3*fp_zeta &
&         +       zeta4  *fpp_zeta

         d2ecrs_drs2=d2ecrs0_drs2+&
&         f_zeta*(zeta4*d2gcrs_drs2-d2macrs_drs2*fsec_inv)
         d2ecrs_drsdzeta=dfzeta4_dzeta*dgcrs_drs-fp_zeta*dmacrs_drs*fsec_inv
         d2ecrs_dzeta2=d2fzeta4_dzeta2*gcrs-fpp_zeta*macrs*fsec_inv

!        End condition of abs(order)>1
!        Add LSD correlation functional to GGA exchange functional
         exci(ipts)=exci(ipts)+ecrs
         vxcadd=ecrs-rs*third*decrs_drs-zeta*decrs_dzeta
!        decrs_drup=vxcadd+decrs_dzeta
!        decrs_drdn=vxcadd-decrs_dzeta
         vxci(ipts,1)=vxci(ipts,1)+vxcadd+decrs_dzeta
         vxci(ipts,2)=vxci(ipts,2)+vxcadd-decrs_dzeta



         dvcrs_drs=third*(2._dp*decrs_drs-rs*d2ecrs_drs2)
!        And d(vxc)/d(rho)=(-rs/(3*rho))*d(vxc)/d(rs)
         d2ecrs_drho2= -rs**4*(four_pi*third)*third*dvcrs_drs
         d2ecrs_drup2=d2ecrs_drho2+&
&         two*(-third*rs*d2ecrs_drsdzeta)*(1._dp-zeta)*rhotot_inv+ &
&         d2ecrs_dzeta2*(1._dp-zeta)**2*rhotot_inv
         d2ecrs_drdndrup=d2ecrs_drho2+&
&         2.0_dp*(-third*rs*d2ecrs_drsdzeta)*(-zeta)*rhotot_inv+ &
&         d2ecrs_dzeta2*(1._dp-zeta)*(-1._dp-zeta)*rhotot_inv
         d2ecrs_drdn2=d2ecrs_drho2+&
&         2.0_dp*(-third*rs*d2ecrs_drsdzeta)*(-1._dp-zeta)*rhotot_inv+ &
&         d2ecrs_dzeta2*(-1._dp-zeta)**2*rhotot_inv




         if (order==3) then


!          fab : INGREDIENTS NEEDED
           
           a1fa=-third*(threefourth_divpi**(third))*((rhotot)**(-4._dp/3._dp))
           a2fa=(1._dp-zeta)/rhotot
           b2fa=(-2._dp/3._dp)*((threefourth_divpi)**(third))*((7._dp/3._dp)*(-1._dp+zeta)/((rhotot)**(7._dp/3._dp)))
           b1fa=a2fa
           c2fa=((1._dp-zeta)**2)*(-3._dp*(1._dp/((rhotot)**2)))
           c1fa=((1._dp-zeta)*(1._dp-zeta))/rhotot
           e2fa=(2._dp/3._dp)*((threefourth_divpi)**(third))*((1._dp-(7._dp/3._dp)*zeta)/((rhotot)**(7._dp/3._dp)))
           e1fa=-zeta/rhotot
           f2fa=(2._dp*zeta)*(1._dp/((rhotot)**2))-(3._dp*zeta*zeta)*(1._dp/((rhotot)**2))+1._dp/(((rhotot)**2))
           f1fa=(zeta*zeta-1._dp)/rhotot
           g1fa=a1fa
           g2fa=(-1._dp-zeta)/rhotot
           h2fa=(2._dp/3._dp)*((threefourth_divpi)**(third))*((-1._dp-(7._dp/3._dp)*zeta)/((rhotot)**(7._dp/3._dp)))
           h1fa=e1fa
           i2fa=((-2._dp*zeta)-(3*zeta*zeta)+1)/((rhotot)**2)
           i1fa=f1fa
           m2fa=(-2._dp/3._dp)*((threefourth_divpi)**(third))*(((7._dp/3._dp)*(zeta+1._dp))/((rhotot)**(7._dp/3._dp)))
           m1fa=g2fa
           n2fa=(-3._dp*(1._dp+zeta)*(1._dp+zeta))/(rhotot*rhotot)
           n1fa=((-1._dp-zeta)**2)/rhotot


!          TERMS APPEARING IN THE THIRD ORDER DERIVATIVES
!          terms appearing in the third order derivatives of the spin polarized 
!          correlation energy with respect to spin densities


!          ec1_q0p=-2.0_dp*ec1_aa*ec1_a1
!          ec1_q1ppp=(3._dp/4._dp)*ec1_aa*(ec1_b1*(rsm1_2**5)-ec1_b3*(rsm1_2**3))
!          This must be erroneous ...
!          d3ecrs1_drs3=(ec1_q1pp*(4._dp*ec1_aa*ec1_a1-ec1_q0p)-ec1_q0*ec1_q1ppp)*ec1_den+ &
!          &           ((-ec1_q1p**2)*(4._dp*ec1_aa*ec1_a1)+ec1_q0*ec1_q1pp*ec1_q1+ec1_q0p*(ec1_q1p**2)+ &
!          &           2._dp*ec1_q0*ec1_q1p*ec1_q1pp)*(ec1_den**2)*(2._dp*ec1_q1+1._dp)- &
!          &           (2._dp*ec1_q0*(ec1_q1p**3)*((2._dp*ec1_q1+1._dp)**2))*(ec1_den**3)+  &
!          &           (2._dp*ec1_q0*(ec1_q1p**3))*(ec1_den**2)      

           ec1_q1ppp = 0.75_dp*ec1_aa*(rsm1_2**5)*(ec1_b1-ec1_b3*rs)
           ec1_f1 = 1._dp/(ec1_q1*ec1_q1*(1._dp + ec1_q1))
           ec1_f2 = 1._dp/(ec1_q1*(1+ec1_q1))
           d3ecrs1_drs3 = 6._dp*ec1_q1p*ec1_f1*(-2._dp*ec1_aa*ec1_a1*ec1_q1p + &
&           ec1_q0*ec1_q1pp) - &
&           ec1_f2*(-6._dp*ec1_aa*ec1_a1*ec1_q1pp + ec1_q0*ec1_q1ppp + &
&           ec1_f2*(3._dp*ec1_q1p*(-2._dp*ec1_aa*ec1_a1*ec1_q1p + ec1_q0*ec1_q1pp) + &
&           ec1_f2*2._dp*ec1_q0*(ec1_q1p**3)*(1._dp + 3._dp*ec1_q1*(1._dp + ec1_q1))))

           
!          mac_q0p=-2.0_dp*mac_aa*mac_a1
!          mac_q1ppp=(3._dp/4._dp)*mac_aa*(mac_b1*((rsm1_2)**5)-mac_b3*((rsm1_2)**3))     
!          This must be erroneous ...
!          d3macrs_drs3=(mac_q1pp*(4._dp*mac_aa*mac_a1-mac_q0p)-mac_q0*mac_q1ppp)*mac_den+ &
!          &           ((-mac_q1p**2)*(4._dp*mac_aa*mac_a1)+mac_q0*mac_q1pp*mac_q1+mac_q0p*(mac_q1p**2)+ &
!          &           2._dp*mac_q0*mac_q1p*mac_q1pp)*(mac_den**2)*(2._dp*mac_q1+1._dp)- &
!          &           (2._dp*mac_q0*(mac_q1p**3)*((2._dp*mac_q1+1._dp)**2))*(mac_den**3)+  &
!          &           (2._dp*mac_q0*(mac_q1p**3))*(mac_den**2)         

           mac_q1ppp = 0.75_dp*mac_aa*(rsm1_2**5)*(mac_b1-mac_b3*rs)
           mac_f1 = 1._dp/(mac_q1*mac_q1*(1._dp + mac_q1))
           mac_f2 = 1._dp/(mac_q1*(1+mac_q1))
           d3macrs_drs3 = 6._dp*mac_q1p*mac_f1*(-2._dp*mac_aa*mac_a1*mac_q1p + &
&           mac_q0*mac_q1pp) - &
&           mac_f2*(-6._dp*mac_aa*mac_a1*mac_q1pp + mac_q0*mac_q1ppp + &
&           mac_f2*(3._dp*mac_q1p*(-2._dp*mac_aa*mac_a1*mac_q1p + mac_q0*mac_q1pp) + &
&           mac_f2*2._dp*mac_q0*(mac_q1p**3)*(1._dp + 3._dp*mac_q1*(1._dp + mac_q1))))

           d3gcrs_drs3=d3ecrs1_drs3-d3ecrs0_drs3+d3macrs_drs3*fsec_inv      
           d3ecrs_drs3=d3ecrs0_drs3+f_zeta*(zeta4*d3gcrs_drs3-d3macrs_drs3*fsec_inv)   
           factfppp_zeta=-two*third*factfpp_zeta*alpha_zeta2                   
           fppp_zeta=factfppp_zeta*(((zetpm1_3(ipts))**5)-((zetmm1_3(ipts))**5))   

           d3ecrs_dzeta3=(24._dp*zeta*f_zeta+36._dp*(zeta**2)*fp_zeta+  &
&           12._dp*(zeta**3)*fpp_zeta+(zeta**4)*fppp_zeta)*gcrs+   &
&           fppp_zeta*(-macrs)*fsec_inv                                   

           d3ecrs_drs2dzeta=dfzeta4_dzeta*(d2gcrs_drs2)+   &
&           fp_zeta*(-d2macrs_drs2)*fsec_inv                        

           d3ecrs_dzeta2drs=d2fzeta4_dzeta2*dgcrs_drs+  &
&           fpp_zeta*(-dmacrs_drs)*fsec_inv                               




!          ***************** all this part has been commented following the suggestion by xavier

!          The following is the calculations of only a part of the third order derivatives.
!          the term d3ecrs_drho3 is the part which remains in the
!          non spin polarized limit
!          THIS CODING IS CORRECT, but the alternative one is also correct ...

!          d3ecrs_drho3=(128._dp*pi*pi/243._dp)*(rs**7)*(decrs_drs)- &
!          &           (48._dp*pi*pi/243._dp)*(rs**8)*d2ecrs_drs2- &
!          &           (16._dp*pi*pi/243._dp)*(rs**9)*d3ecrs_drs3

!          d3ecrs_drhoupdrho2=d3ecrs_drho3+ &
!          &           a2fa*((-8._dp*pi/27._dp)*(rs**4)*d2ecrs_drsdzeta+ &
!          &           (4._dp*pi/27._dp)*(rs**5)*d3ecrs_drs2dzeta)

           
!          d3ecrs_drhodndrho2=d3ecrs_drho3+ &
!          &           g2fa*((-8._dp*pi/27._dp)*(rs**4)*d2ecrs_drsdzeta+ &
!          &           (4._dp*pi/27._dp)*(rs**5)*d3ecrs_drs2dzeta)    


!          third order derivatives of the exchange-correlation part

!          d3ecrs_drup3=d3ecrs_drhoupdrho2+b2fa*d2ecrs_drsdzeta-  &
!          &           (2._dp/3._dp)*rs*b1fa*a1fa*d3ecrs_drs2dzeta- &
!          &           (2._dp/3._dp)*rs*b1fa*b1fa*d3ecrs_dzeta2drs+ &
!          &           c2fa*d2ecrs_dzeta2+c1fa*(a1fa*d3ecrs_dzeta2drs+a2fa*d3ecrs_dzeta3)


!          d3ecrs_drup2drdn=d3ecrs_drhoupdrho2+e2fa*d2ecrs_drsdzeta-  &
!          &           (2._dp/3._dp)*rs*e1fa*a1fa*d3ecrs_drs2dzeta- &
!          &           (2._dp/3._dp)*rs*e1fa*b1fa*d3ecrs_dzeta2drs+ &
!          &           f2fa*d2ecrs_dzeta2+f1fa*(a1fa*d3ecrs_dzeta2drs+a2fa*d3ecrs_dzeta3)

!          d3ecrs_drupdrdn2=d3ecrs_drhodndrho2+h2fa*d2ecrs_drsdzeta-  &
!          &           (2._dp/3._dp)*rs*h1fa*g1fa*d3ecrs_drs2dzeta- &
!          &           (2._dp/3._dp)*rs*h1fa*g2fa*d3ecrs_dzeta2drs+ &
!          &           i2fa*d2ecrs_dzeta2+i1fa*(g1fa*d3ecrs_dzeta2drs+g2fa*d3ecrs_dzeta3)


!          d3ecrs_drdn3=d3ecrs_drhodndrho2+m2fa*d2ecrs_drsdzeta-  &
!          &           (2._dp/3._dp)*rs*m1fa*g1fa*d3ecrs_drs2dzeta- &
!          &           (2._dp/3._dp)*rs*m1fa*g2fa*d3ecrs_dzeta2drs+ &
!          &           n2fa*d2ecrs_dzeta2+n1fa*(g1fa*d3ecrs_dzeta2drs+g2fa*d3ecrs_dzeta3)


!          ********* suggested by xavier  (now corrected, XG100524)

           sp1_up3=three-three*zeta
           sp1_up2dn=one-three*zeta
           sp1_updn2=-one-three*zeta
           sp1_dn3=-three-three*zeta

           sp2_up3=three-six*zeta+three*zeta*zeta
           sp2_up2dn=-one-two*zeta+three*zeta*zeta
           sp2_updn2=-one+two*zeta+three*zeta*zeta
           sp2_dn3=three+six*zeta+three*zeta*zeta

           sp3_up3=(one-zeta)**3
           sp3_up2dn=-one+zeta+zeta**2-zeta**3
           sp3_updn2=one+zeta-zeta**2-zeta**3
           sp3_dn3=-(one+zeta)**3

           d3ecrs_sp0=(eight*rs*decrs_drs-three*rs*rs*d2ecrs_drs2-rs*rs*rs*d3ecrs_drs3)/(rhotot*rhotot*27.0_dp)
           d3ecrs_sp1=(four*rs*d2ecrs_drsdzeta+rs*rs*d3ecrs_drs2dzeta)/(rhotot*rhotot*nine)
           d3ecrs_sp2=(-three*d2ecrs_dzeta2-rs*d3ecrs_dzeta2drs)/(rhotot*rhotot*three)
           d3ecrs_sp3=d3ecrs_dzeta3/(rhotot*rhotot)

           d3ecrs_drup3=d3ecrs_sp0+d3ecrs_sp1*sp1_up3+d3ecrs_sp2*sp2_up3+d3ecrs_sp3*sp3_up3
           d3ecrs_drup2drdn=d3ecrs_sp0+d3ecrs_sp1*sp1_up2dn+d3ecrs_sp2*sp2_up2dn+d3ecrs_sp3*sp3_up2dn
           d3ecrs_drupdrdn2=d3ecrs_sp0+d3ecrs_sp1*sp1_updn2+d3ecrs_sp2*sp2_updn2+d3ecrs_sp3*sp3_updn2
           d3ecrs_drdn3=d3ecrs_sp0+d3ecrs_sp1*sp1_dn3+d3ecrs_sp2*sp2_dn3+d3ecrs_sp3*sp3_dn3

!          **** end of the section suggested by xavier...


!          fab: this is the end of the if over order==3


         end if

!        fab: I think the following lines are wrong..indeed we are in the case option>0, so option cannot be -1
!        I comment them and I put the if only for option 1 and 3

!        if(option==1 .or. option==-1 .or. option==3)then
!        dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drup2
!        dvxci(ipts,2)=dvxci(ipts,2)+d2ecrs_drdndrup
!        dvxci(ipts,3)=dvxci(ipts,3)+d2ecrs_drdn2


!        fab: however, here the thing seems a bit strange...option=3 doesn't seem to be completely implemented 
!        (the second derivatives of ec1_q and ec0_q are the derived only in correspondance of the first derivative in the case option !=3 and !=4)
!        so..I think that here the case "or option==3 should be cancelled

         if(option==1 .or. option==3)then
           dvxci(ipts,1)=dvxci(ipts,1)+d2ecrs_drup2
           dvxci(ipts,2)=dvxci(ipts,2)+d2ecrs_drdndrup
           dvxci(ipts,3)=dvxci(ipts,3)+d2ecrs_drdn2

           if(order==3) then
!            third order derivatives of the spin polarized exchange+correlation energy
             d2vxci(ipts,1)=d2vxci(ipts,1)+d3ecrs_drup3
             d2vxci(ipts,2)=d2vxci(ipts,2)+d3ecrs_drup2drdn
             d2vxci(ipts,3)=d2vxci(ipts,3)+d3ecrs_drupdrdn2
             d2vxci(ipts,4)=d2vxci(ipts,4)+d3ecrs_drdn3
!            DEBUG
!            wecrsz(ipts,1)=ecrs
!            wecrsz(ipts,1)=ecrs*rhotot
!            wecrsz(ipts,2)=rs
!            wecrsz(ipts,2)=rho_updn(ipts,1)
!            wecrsz(ipts,3)=zeta
!            wecrsz(ipts,3)=rho_updn(ipts,2)
!            wecrsz(ipts,5)=ecrs0
!            wecrsz(ipts,6)=gcrs
!            wecrsz(ipts,7)=macrs
!            wecrsz(ipts,8)=ecrs1
!            d1wecrsz(ipts,1)=decrs_drs
!            d1wecrsz(ipts,1)=decrs_drup
!            d1wecrsz(ipts,2)=decrs_dzeta
!            d1wecrsz(ipts,2)=decrs_drdn
!            d1wecrsz(ipts,5)=decrs0_drs
!            d1wecrsz(ipts,6)=dgcrs_drs
!            d1wecrsz(ipts,7)=dmacrs_drs
!            d1wecrsz(ipts,8)=decrs1_drs
!            d2wecrsz(ipts,1)=d2ecrs_drs2
!            d2wecrsz(ipts,1)=d2ecrs_drup2
!            d2wecrsz(ipts,2)=d2ecrs_drsdzeta
!            d2wecrsz(ipts,2)=d2ecrs_drdndrup
!            d2wecrsz(ipts,3)=d2ecrs_dzeta2
!            d2wecrsz(ipts,3)=d2ecrs_drdn2
!            d2wecrsz(ipts,5)=d2ecrs0_drs2
!            d2wecrsz(ipts,6)=d2gcrs_drs2
!            d2wecrsz(ipts,7)=d2macrs_drs2
!            d2wecrsz(ipts,8)=d2ecrs1_drs2
!            d3wecrsz(ipts,1)=d3ecrs_drs3
!            d3wecrsz(ipts,1)=d3ecrs_drup3
!            d3wecrsz(ipts,2)=d3ecrs_drs2dzeta
!            d3wecrsz(ipts,2)=d3ecrs_drup2drdn
!            d3wecrsz(ipts,3)=d3ecrs_dzeta2drs
!            d3wecrsz(ipts,3)=d3ecrs_drupdrdn2
!            d3wecrsz(ipts,4)=d3ecrs_dzeta3
!            d3wecrsz(ipts,4)=d3ecrs_drdn3
!            d3wecrsz(ipts,5)=d3ecrs0_drs3
!            d3wecrsz(ipts,6)=d3gcrs_drs3
!            d3wecrsz(ipts,7)=d3macrs_drs3
!            d3wecrsz(ipts,8)=d3ecrs1_drs3
!            ENDDEBUG
           end if


         else
           dvxci(ipts,9)=d2ecrs_drup2
           dvxci(ipts,10)=d2ecrs_drdndrup
           dvxci(ipts,11)=d2ecrs_drdn2
         end if




!        -----------------------------------------------------------------------------
!        Eventually add the GGA correlation part of the PBE functional
!        Note : the computation of the potential in the spin-unpolarized
!        case could be optimized much further. Other optimizations are left to do.

         if(option==2 .or. option==5 .or. option==6 .or. option==7)then
!          The definition of phi has been slightly changed, because
!          the original PBE one gives divergent behaviour for fully
!          polarized points
!          zetpm1_3=(1.0_dp+zeta*alpha_zeta)**(-third)
!          zetmm1_3=(1.0_dp-zeta*alpha_zeta)**(-third)
           phi_zeta=( zetpm1_3(ipts)*(1.0_dp+zeta*alpha_zeta)+ &
&           zetmm1_3(ipts)*(1.0_dp-zeta*alpha_zeta)   )*0.5_dp
           phip_zeta=(zetpm1_3(ipts)-zetmm1_3(ipts))*third*alpha_zeta
           phi_zeta_inv=1.0_dp/phi_zeta
           phi_logder=phip_zeta*phi_zeta_inv
           phi3_zeta=phi_zeta*phi_zeta*phi_zeta
           gamphi3inv=gamma_inv*phi_zeta_inv*phi_zeta_inv*phi_zeta_inv
           phipp_zeta=-alpha_zeta*alpha_zeta*ninth*&
&           (zetpm1_3(ipts)*zetpm1_3(ipts)*zetpm1_3(ipts)*zetpm1_3(ipts) + &
&           zetmm1_3(ipts)*zetmm1_3(ipts)*zetmm1_3(ipts)*zetmm1_3(ipts)  )

!          From ec to bb
           bb=ecrs*gamphi3inv
           dbb_drs=decrs_drs*gamphi3inv
           dbb_dzeta=gamphi3inv*(decrs_dzeta-three*ecrs*phi_logder)
           d2bb_drs2=d2ecrs_drs2*gamphi3inv
           d2bb_drsdzeta=gamphi3inv*(d2ecrs_drsdzeta-three*decrs_drs*phi_logder)
           d2bb_dzeta2=gamphi3inv*(d2ecrs_dzeta2-six*decrs_dzeta*phi_logder+&
&           12.0_dp*ecrs*phi_logder*phi_logder-three*ecrs*phi_zeta_inv*phipp_zeta)

!          From bb to cc
           exp_pbe=exp(-bb)
           cc=one/(exp_pbe-one)
           dcc_dbb=cc*cc*exp_pbe
           dcc_drs=dcc_dbb*dbb_drs
           dcc_dzeta=dcc_dbb*dbb_dzeta
           d2cc_dbb2=cc*cc*exp_pbe*(two*cc*exp_pbe-one)
           d2cc_drs2=d2cc_dbb2*dbb_drs*dbb_drs+dcc_dbb*d2bb_drs2
           d2cc_drsdzeta=d2cc_dbb2*dbb_drs*dbb_dzeta+dcc_dbb*d2bb_drsdzeta
           d2cc_dzeta2=d2cc_dbb2*dbb_dzeta*dbb_dzeta+dcc_dbb*d2bb_dzeta2

!          From cc to aa
           coeff_aa=beta*gamma_inv*phi_zeta_inv*phi_zeta_inv
           aa=coeff_aa*cc
           daa_drs=coeff_aa*dcc_drs
           daa_dzeta=-two*aa*phi_logder+coeff_aa*dcc_dzeta
           d2aa_drs2=coeff_aa*d2cc_drs2
           d2aa_drsdzeta=-two*daa_drs*phi_logder+coeff_aa*d2cc_drsdzeta
           d2aa_dzeta2=aa*(-two*phi_zeta_inv*phipp_zeta+six*phi_logder*phi_logder)+&
&           coeff_aa*(-four*dcc_dzeta*phi_logder+d2cc_dzeta2)

!          Introduce tt : do not assume that the spin-dependent gradients are collinear
           grrho2=grho2_updn(ipts,3)
           dtt_dg=two*rhotot_inv*rhotot_inv*rhotmot*coeff_tt
!          Note that tt is (the t variable of PBE divided by phi) squared
           tt=half*grrho2*dtt_dg

!          Get xx from aa and tt
           xx=aa*tt
           dxx_drs=daa_drs*tt
           dxx_dzeta=daa_dzeta*tt
           dxx_dtt=aa
           d2xx_drs2=d2aa_drs2*tt
           d2xx_drsdzeta=d2aa_drsdzeta*tt
           d2xx_drsdtt=daa_drs
           d2xx_dttdzeta=daa_dzeta
           d2xx_dzeta2=d2aa_dzeta2*tt

!          From xx to pade
           pade_den=one/(one+xx*(one+xx))
           pade=(one+xx)*pade_den
           dpade_dxx=-xx*(two+xx)*pade_den**2
           dpade_drs=dpade_dxx*dxx_drs
           dpade_dtt=dpade_dxx*dxx_dtt
           dpade_dzeta=dpade_dxx*dxx_dzeta
           d2pade_dxx2=two*(-one+xx*xx*(three+xx))*pade_den*pade_den*pade_den
           d2pade_drs2=d2pade_dxx2*dxx_drs*dxx_drs+dpade_dxx*d2xx_drs2
           d2pade_drsdtt=d2pade_dxx2*dxx_drs*dxx_dtt+dpade_dxx*d2xx_drsdtt
           d2pade_drsdzeta=d2pade_dxx2*dxx_drs*dxx_dzeta+dpade_dxx*d2xx_drsdzeta
           d2pade_dtt2=d2pade_dxx2*dxx_dtt*dxx_dtt
           d2pade_dttdzeta=d2pade_dxx2*dxx_dtt*dxx_dzeta+dpade_dxx*d2xx_dttdzeta
           d2pade_dzeta2=d2pade_dxx2*dxx_dzeta*dxx_dzeta+dpade_dxx*d2xx_dzeta2

!          From pade to qq
           coeff_qq=tt*phi_zeta_inv*phi_zeta_inv
           qq=coeff_qq*pade
           dqq_drs=coeff_qq*dpade_drs
           dqq_dtt=pade*phi_zeta_inv*phi_zeta_inv+coeff_qq*dpade_dtt
           dqq_dzeta=coeff_qq*(dpade_dzeta-two*pade*phi_logder)
           d2qq_drs2=coeff_qq*d2pade_drs2
           d2qq_drsdtt=phi_zeta_inv*phi_zeta_inv*(dpade_drs+tt*d2pade_drsdtt)
           d2qq_drsdzeta=coeff_qq*(d2pade_drsdzeta-two*dpade_drs*phi_logder)
           d2qq_dtt2=phi_zeta_inv*phi_zeta_inv*(two*dpade_dtt+tt*d2pade_dtt2)
           d2qq_dttdzeta=phi_zeta_inv*phi_zeta_inv*(dpade_dzeta-two*pade*phi_logder)+&
&           coeff_qq*(d2pade_dttdzeta-two*dpade_dtt*phi_logder)
           d2qq_dzeta2=coeff_qq*( d2pade_dzeta2-four*dpade_dzeta*phi_logder &
&           +six*pade*phi_logder*phi_logder            &
&           -two*pade*phi_zeta_inv*phipp_zeta)

!          From qq to rr
           arg_rr=one+beta*gamma_inv*qq
           div_rr=one/arg_rr
           rr=gamma*log(arg_rr)
           drr_dqq=beta*div_rr
           drr_drs=drr_dqq*dqq_drs
           drr_dtt=drr_dqq*dqq_dtt
           drr_dzeta=drr_dqq*dqq_dzeta
           d2rr_dqq2=-div_rr**2*beta*beta*gamma_inv
           d2rr_drs2=d2rr_dqq2*dqq_drs*dqq_drs+drr_dqq*d2qq_drs2
           d2rr_drsdtt=d2rr_dqq2*dqq_drs*dqq_dtt+drr_dqq*d2qq_drsdtt
           d2rr_drsdzeta=d2rr_dqq2*dqq_drs*dqq_dzeta+drr_dqq*d2qq_drsdzeta
           d2rr_dtt2=d2rr_dqq2*dqq_dtt*dqq_dtt+drr_dqq*d2qq_dtt2
           d2rr_dttdzeta=d2rr_dqq2*dqq_dtt*dqq_dzeta+drr_dqq*d2qq_dttdzeta
           d2rr_dzeta2=d2rr_dqq2*dqq_dzeta*dqq_dzeta+drr_dqq*d2qq_dzeta2

!          From rr to hh
           hh=phi3_zeta*rr
           dhh_drs=phi3_zeta*drr_drs
           dhh_dtt=phi3_zeta*drr_dtt
           dhh_dzeta=phi3_zeta*(drr_dzeta+three*rr*phi_logder)
           d2hh_drs2=phi3_zeta*d2rr_drs2
           d2hh_drsdtt=phi3_zeta*d2rr_drsdtt
           d2hh_drsdzeta=phi3_zeta*(d2rr_drsdzeta+three*drr_drs*phi_logder)
           d2hh_dtt2=phi3_zeta*d2rr_dtt2
           d2hh_dttdzeta=phi3_zeta*(d2rr_dttdzeta+three*drr_dtt*phi_logder)
           d2hh_dzeta2=phi3_zeta*(six*rr*phi_logder*phi_logder+&
&           six*phi_logder*drr_dzeta+d2rr_dzeta2)  &
&           +three*phi_zeta*phi_zeta*rr*phipp_zeta

!          The GGA correlation energy is added
           exci(ipts)=exci(ipts)+hh

!          Change of variables : from (rs,zeta,tt) to (rhoup,rhodn,grrho)


!          From hh to the derivative of the energy wrt the density
           drhohh_drho=hh-third*rs*dhh_drs-zeta*dhh_dzeta-seven*third*tt*dhh_dtt
           vxci(ipts,1)=vxci(ipts,1)+drhohh_drho+dhh_dzeta
           vxci(ipts,2)=vxci(ipts,2)+drhohh_drho-dhh_dzeta


!          From hh to the derivative of the energy wrt to the gradient of the
!          density, divided by the gradient of the density
!          (The v3.3 definition includes the division by the norm of the gradient)
           dvxcdgr(ipts,3)=rhotot*dtt_dg*dhh_dtt

           d2rhohh_drho2=rhotot_inv*&
&           (-two*ninth*rs*dhh_drs +seven*four*ninth*tt*dhh_dtt &
&           +ninth*rs*rs*d2hh_drs2+zeta*zeta*d2hh_dzeta2+(seven*third*tt)**2*d2hh_dtt2 &
&           +two*third*rs*zeta*d2hh_drsdzeta+two*seven*ninth*rs*tt*d2hh_drsdtt &
&           +two*seven*third*tt*zeta*d2hh_dttdzeta)
           d2rhohh_drhodg=dtt_dg*(-four*third*dhh_dtt-third*rs*d2hh_drsdtt &
&           -zeta*d2hh_dttdzeta-seven*third*tt*d2hh_dtt2)

!          Component 12 : first derivative with respect to the gradient
!          of the density, div by the grad of the density
           dvxci(ipts,12)=dvxcdgr(ipts,3)
!          Components 9, 10 and 11 : second derivatives with respect to the spin-density
!          Note that there is already a contribution from LSDA
           dvxci(ipts,9)=dvxci(ipts,9)+d2rhohh_drho2+rhotot_inv*           &
&           ( d2hh_dzeta2*(one-two*zeta) &
&           -two*third*rs*d2hh_drsdzeta-14.0_dp*third*tt*d2hh_dttdzeta)
           dvxci(ipts,10)=dvxci(ipts,10)+d2rhohh_drho2-rhotot_inv*d2hh_dzeta2
           dvxci(ipts,11)=dvxci(ipts,11)+d2rhohh_drho2+rhotot_inv*           &
&           ( d2hh_dzeta2*(one+two*zeta) &
&           +two*third*rs*d2hh_drsdzeta+14.0_dp*third*tt*d2hh_dttdzeta)
!          Components 13 and 14 : second derivatives with respect to spin density
!          and gradient, divided by the gradient
           dvxci(ipts,13)=d2rhohh_drhodg+dtt_dg*d2hh_dttdzeta
           dvxci(ipts,14)=d2rhohh_drhodg-dtt_dg*d2hh_dttdzeta
!          Component 15 : derivative of the (derivative wrt the gradient div by the grad),
!          divided by the grad
           dvxci(ipts,15)=rhotot*d2hh_dtt2*dtt_dg*dtt_dg

!          End condition of GGA
         end if

       else  ! no correlation

         if(ndvxci > 8)then
!          Must zero the correlation part of the xc kernel
           dvxci(:,9:15)=zero
         end if

!        End condition of including correlation, and not only exchange
       end if

!      Correlation has been added
!      -----------------------------------------------------------------------------

     end do
   end if 

!  fab: this should be the else on nspden
 else
!  Disallowed value for nspden
   write(message, '(a,a,a,a,a,a,i12,a)' ) ch10,&
&   ' xcpbe: BUG -',ch10,&
&   '  Argument nspden must be 1 or 2; ',ch10,&
&   '  Value provided as argument was ',nspden,'.'
!   call wrtout(std_out,message,'COLL')
!   call leave_new('COLL')
    call rx(message)
 end if

!DEBUG
!Finite-difference debugging, do not take away
!if(debug/=0)then
!do ipts=1,5,5

!rho=rho_updn(ipts,1)+rho_updn(ipts,2)
!write(6, '(a,i5,a,2es16.8)' ) ' Point number',ipts,' with rhoup,rhodn=',rho_updn(ipts,1),rho_updn(ipts,2)
!write(6, '(a,i5,a,2es16.8)' ) ' Point number',ipts+1,' with rhoup,rhodn=',rho_updn(ipts+1,1),rho_updn(ipts+1,2)
!write(6, '(a,i5,a,2es16.8)' ) ' Point number',ipts+2,' with rhoup,rhodn=',rho_updn(ipts+2,1),rho_updn(ipts+2,2)
!write(6, '(a,i5,a,2es16.8)' ) ' Point number',ipts+3,' with rhoup,rhodn=',rho_updn(ipts+3,1),rho_updn(ipts+3,2)
!write(6, '(a,i5,a,2es16.8)' ) ' Point number',ipts+4,' with rhoup,rhodn=',rho_updn(ipts+4,1),rho_updn(ipts+4,2)

!! For rho
!if(debug==1)then
!write(6, '(a)' )' Direct values :'
!write(6, '(3es16.8)' )exci(ipts)*rho,vxci(ipts,1),vxci(ipts,2)
!else
!!  For grho2
!write(6, '(4es16.8)' )exci(ipts)*rho,dvxcdgr(ipts,1),&
!&  dvxcdgr(ipts,2),dvxcdgr(ipts,3)
!end if

!write(6, '(4es16.8)' )dvxci(ipts,1:3)  ! For LDA
!write(6, '(a)' )'     3rd-order '
!write(6, '(4es16.8)' )d2vxci(ipts,1:4)  ! For LDA

!write(6, '(4es16.8)' )dvxci(ipts,1:4)  ! For exchange
!write(6, '(4es16.8)' )dvxci(ipts,5:8)  ! For exchange
!write(6, '(4es16.8)' )dvxci(ipts,9:12) ! For correlation
!write(6, '(4es16.8)' )dvxci(ipts,13:15) ! For correlation

!if(debug==1)then
!!  For rho
!write(6, '(a)' )' Finite-difference values :'
!write(6, '(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*(rho+delta) - exci(ipts+2)*(rho-delta) )/2._dp/delta,&
!&      ( exci(ipts+3)*(rho+delta) - exci(ipts+4)*(rho-delta) )/2._dp/delta
!write(6, '(3es16.8)' )&
!&    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2._dp/delta,&
!&    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2._dp/delta,&
!&    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2._dp/delta
!!This is for order 3
!write(6, '(a)' )'     3rd-order by two methods, giving components 1, 2, 3, then on the next line 2, 3, 4 '
!write(6, '(3es16.8)' )&
!&    ( dvxci(ipts+1,1) - dvxci(ipts+2,1) )/2._dp/delta,&
!&    ( dvxci(ipts+1,2) - dvxci(ipts+2,2) )/2._dp/delta,&
!&    ( dvxci(ipts+1,3) - dvxci(ipts+2,3) )/2._dp/delta
!write(6, '(3es16.8)' )&
!&    ( dvxci(ipts+3,1) - dvxci(ipts+4,1) )/2._dp/delta,&
!&    ( dvxci(ipts+3,2) - dvxci(ipts+4,2) )/2._dp/delta,&
!&    ( dvxci(ipts+3,3) - dvxci(ipts+4,3) )/2._dp/delta

!write(6,*)
!write(6,*)' Now for ecrs and derivatives '
!write(6, '(a)' )' ecrs, rs, zeta ='
!write(6, '(3es16.8)' )wecrsz(ipts,1:3)
!write(6, '(3es16.8)' )wecrsz(ipts+1,1:3)
!write(6, '(3es16.8)' )wecrsz(ipts+2,1:3)
!write(6, '(3es16.8)' )wecrsz(ipts+3,1:3)
!write(6, '(3es16.8)' )wecrsz(ipts+4,1:3)
!write(6, '(a)' )' ecrs derivatives :'
!write(6, '(3es16.8)' )d1wecrsz(ipts,1:2)
!write(6, '(3es16.8)' )d2wecrsz(ipts,1:3)
!write(6, '(4es16.8)' )d3wecrsz(ipts,1:4)
!write(6, '(a)' )' Finite-differences :'
!write(6, '(3es16.8)' )&
!&    ( wecrsz(ipts+1,1) - wecrsz(ipts+2,1) )/( wecrsz(ipts+1,2) - wecrsz(ipts+2,2) ),&
!&    ( wecrsz(ipts+3,1) - wecrsz(ipts+4,1) )/( wecrsz(ipts+3,3) - wecrsz(ipts+4,3) )
!write(6, '(3es16.8)' )&
!&    ( d1wecrsz(ipts+1,1) - d1wecrsz(ipts+2,1) )/( wecrsz(ipts+1,2) - wecrsz(ipts+2,2) ),&
!&    ( d1wecrsz(ipts+1,2) - d1wecrsz(ipts+2,2) )/( wecrsz(ipts+1,2) - wecrsz(ipts+2,2) )
!write(6, '(3es16.8)' )&
!&    ( d1wecrsz(ipts+3,1) - d1wecrsz(ipts+4,1) )/( wecrsz(ipts+3,3) - wecrsz(ipts+4,3) ),&
!&    ( d1wecrsz(ipts+3,2) - d1wecrsz(ipts+4,2) )/( wecrsz(ipts+3,3) - wecrsz(ipts+4,3) )
!write(6, '(a)' )' Finite-differences, 3rd order :'
!write(6, '(3es16.8)' )&
!&    ( d2wecrsz(ipts+1,1) - d2wecrsz(ipts+2,1) )/( wecrsz(ipts+1,2) - wecrsz(ipts+2,2) ),&
!&    ( d2wecrsz(ipts+1,2) - d2wecrsz(ipts+2,2) )/( wecrsz(ipts+1,2) - wecrsz(ipts+2,2) ),&
!&    ( d2wecrsz(ipts+1,3) - d2wecrsz(ipts+2,3) )/( wecrsz(ipts+1,2) - wecrsz(ipts+2,2) )
!write(6, '(3es16.8)' )&
!&    ( d2wecrsz(ipts+3,1) - d2wecrsz(ipts+4,1) )/( wecrsz(ipts+3,3) - wecrsz(ipts+4,3) ),&
!&    ( d2wecrsz(ipts+3,2) - d2wecrsz(ipts+4,2) )/( wecrsz(ipts+3,3) - wecrsz(ipts+4,3) ),&
!&    ( d2wecrsz(ipts+3,3) - d2wecrsz(ipts+4,3) )/( wecrsz(ipts+3,3) - wecrsz(ipts+4,3) )


!!This is for GGA

!write(6, '(4es16.8)' )&
!&    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2._dp/delta,&
!&    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2._dp/delta,&
!&    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2._dp/delta,&
!&    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2._dp/delta
!else
!!  For grho2  (should distinguish exchange and correlation ...)
!grr=sqrt(grho2_updn(ipts,1)) ! Analysis of exchange
!grr=sqrt(grho2_updn(ipts,3)) ! Analysis of correlation
!write(6, '(3es16.8)' )exci(ipts)*rho,&
!&      ( exci(ipts+1)*rho - exci(ipts+2)*rho )/2._dp/delta/grr,&
!&      ( exci(ipts+3)*rho - exci(ipts+4)*rho )/2._dp/delta/grr
!write(6, '(3es16.8)' )&
!&    ( vxci(ipts+1,1) - vxci(ipts+2,1) )/2._dp/delta/grr,&
!&    ( vxci(ipts+3,1) - vxci(ipts+4,1) )/2._dp/delta/grr,&
!&    ( vxci(ipts+3,2) - vxci(ipts+4,2) )/2._dp/delta/grr
!write(6, '(4es16.8)' )&
!&    ( dvxcdgr(ipts+1,1) - dvxcdgr(ipts+2,1) )/2._dp/delta/grr,&
!&    ( dvxcdgr(ipts+3,2) - dvxcdgr(ipts+4,2) )/2._dp/delta/grr,&
!&    ( dvxcdgr(ipts+1,3) - dvxcdgr(ipts+2,3) )/2._dp/delta/grr,&
!&    ( dvxcdgr(ipts+3,3) - dvxcdgr(ipts+4,3) )/2._dp/delta/grr
!end if
!end do
!stop
!end if
!ENDDEBUG

 deallocate(rhoarr,rhom1_3,rho_updnm1_3)
 deallocate(zetm,zetmm1_3,zetp,zetpm1_3)

!DEBUG
!deallocate(wecrsz,d1wecrsz,d2wecrsz,d3wecrsz)
!ENDDEBUG 

!DEBUG
!write(6,*)' xcpbe : exit'
!write(6,*)' nspden=',nspden
!if(order==2)stop
!ENDDEBUG

end subroutine xcpbe
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/invcb
!! NAME
!! invcb
!!
!! FUNCTION
!! Compute a set of inverse cubic roots as fast as possible :
!! rspts(:)=rhoarr(:)$^\frac{-1}{3}$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npts=number of real space points on which density is provided
!!  rhoarr(npts)=input data
!!
!! OUTPUT
!!  rspts(npts)=inverse cubic root of rhoarr
!!
!! PARENTS
!!      drivexc,gammapositron,xchcth,xcpbe,xcpositron
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE


!### invcb.F90
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine invcb(rhoarr,rspts,npts)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
!arrays
 real(dp),intent(in) :: rhoarr(npts)
 real(dp),intent(out) :: rspts(npts)

!Local variables-------------------------------
!scalars
 integer :: ii,ipts
 real(dp),parameter :: c2_27=2.0e0_dp/27.0e0_dp,c5_9=5.0e0_dp/9.0e0_dp
 real(dp),parameter :: c8_9=8.0e0_dp/9.0e0_dp,m1thrd=-third
 real(dp) :: del,prod,rho,rhom1,rhomtrd
 logical :: test
!character(len=500) :: message

! *************************************************************************

!Loop over points : here, brute force algorithm
!do ipts=1,npts
!rspts(ipts)=sign( (abs(rhoarr(ipts)))**m1thrd,rhoarr(ipts))
!end do
!

!write(6,*)' invcb : rhoarr, rspts'

 rhomtrd=sign( (abs(rhoarr(1)))**m1thrd, rhoarr(1) )
 rhom1=one/rhoarr(1)
 rspts(1)=rhomtrd
 do ipts=2,npts
!  write(6,*)
!  write(6,*)rhoarr(ipts),rspts(ipts)
   rho=rhoarr(ipts)
   prod=rho*rhom1
!  If the previous point is too far ...
   if(prod < 0.01_dp .or. prod > 10._dp )then
     rhomtrd=sign( (abs(rho))**m1thrd , rho )
     rhom1=one/rho
   else
     del=prod-one
     do ii=1,5
!      Choose one of the two next lines, the last one is more accurate
!      rhomtrd=((one+third*del)/(one+two_thirds*del))*rhomtrd
       rhomtrd=((one+c5_9*del)/(one+del*(c8_9+c2_27*del)))*rhomtrd
       rhom1=rhomtrd*rhomtrd*rhomtrd
       del=rho*rhom1-one
!      write(6,*)rhomtrd,del
       test = del*del < 1.0e-24_dp
       if(test) exit
     end do
     if( .not. test) then
       rhomtrd=sign( (abs(rho))**m1thrd , rho )
     end if
   end if
   rspts(ipts)=rhomtrd
 end do

 end subroutine invcb
!!***





