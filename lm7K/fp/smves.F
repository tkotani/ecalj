      subroutine smves(mode,nbas,ssite,sspec,slat,k1,k2,k3,qmom,gpot0,
     .vval,hpot0,sgp0,smrho,smpot,vconst,smq,qsmc,f,rhvsm,zvnsm,zsum,
     .vrmt,qbg)

      use m_struc_def  !Cgetarg
      use m_globalvariables
C- Electrostatic potential of the smooth density.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 use input vconst
Ci         :1 generate vconst as - average v(RMT)
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Ci   slat  :struct containing information about the lattice
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   smrho :smooth density on real-space mesh
Ci   qbg   : back ground charge
Cio Inputs/Outputs
Cio  vconst:constant potential to be added to total
Cio        :On input  vconst is set to a default value
Cio        :On output vconst may be set to the average estat
Cio        :          at the MT boundary.
Co Outputs (see also Remarks)
Co   gpot0 :integrals of compensating gaussians g_RL * phi0~
Co         :For accuracy, integral is split into
Co         :g_RL phi0 (vesgcm) + g_RL phi [n0~-n0] (ugcomp)
Co         :vesgcm projects g_RL to the mesh to do the integral
Co         :ugcomp does its integrals analytically (structure constants)
Co         :NB: There is a local analog of gpot0 generated in locpt2.
Co   vval  :coffs to YL expansion of es potential at MT boundary
Co   hpot0 :integrals of semicore smooth Hankels * phi0~
Co   sgp0  :sgp0 = sum_RL integral qmom_RL g_RL phi0~
Co         :     = integral [n0~-n0 ] phi0~
Co   smpot :smooth potential phi0~ (includes compensating gaussians)
Co   smq   :integral of smooth density n0
Co   qsmc  :pseudocore charge
Co   f     :electrostatic contribution to force.
Co   rhvsm :integral n0~ [phi0~ + vconst]
Co         :(electrostatic energy of sm. density n0~) + vconst*smq
Co   zvnsm :integral (qcorg-z + rhoc) phi0~
Co   vrmt  :electrostatic potential at rmt, with G=0 term in smpot=0
Cl Local variables
Cl   u00   :integral n0 phi[n0] = n0 phi0
Cl   u0g   :integral n0 [phi0~-phi0]
Cl   ugg   :integral [n0~-n0] [phi0~-phi0]
Cl         :Note: ugg is not used.
Cr Remarks
Cr  The total density is a sum of three terms,
Cr
Cr    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
Cr
Cr  The first term is the smooth density on a mesh of points; the
Cr  second is the true density and is defined on a radial mesh for each
Cr  sphere; the last is the 1-center expansion of the smooth density on
Cr  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
Cr  identical to the one-center expansion of n0(mesh).  The sum of the
Cr  three terms converges rapidly with l because errors in n_R(r) are
Cr  mostly canceled by errors in n0_R(r).)
Cr
Cr  We add and subtract a set of compensating gaussian orbitals
Cr
Cr    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
Cr
Cr  which render the integral of the local part (the last 3 terms)
Cr  zero in each RL channel.  The g_RL must be localized enough that
Cr  their spillout beyond the MT radius is negligible.
Cr
Cr  We define
Cr
Cr    n0~ = n0 + compensating gaussians sum_RL Q_RL g_RL
Cr
Cr  In the interstitial, the electrostatic potential of n0~ is the true
Cr  estat potential.  The potential of n0 is called phi0 and the
Cr  potential of n0~ is called phi0~.  The total electrostatic energy
Cr  is computed as
Cr
Cr    the electrostatic energy of  n0~ + integral n0*vconst +
Cr    the electrostatic energy of (neutral) local parts
Cr
Cr  vconst may either be passed as an input (mode=0) or it is
Cr  generated here as the average ves(RMT).
Cr  This routine computes the estat potential and energy from the
Cr  first two terms.  Some variables used in smves and its subroutines:
Cr    Let n0  = smooth density without the compensating sum_RL Q_RL g_RL
Cr        n0~ = n0 + sum_RL Q_RL g_RL
Cr      phi0  = ves[n0]
Cr      phi0~ = ves[n0~]
Cr      g_RL  = gaussian in RL channel
Cr      h_R   = l=0 sm hankel in RL channel, (represents core densities)
Cr    qmom_RL = multipole moment in RL channel of (n_R(r) - n0_R(r))
Cr              so that int n_RL(r)-n0_RL(r) = qmom_RL * g_RL(r)
Cr      gpot0 = vector of integrals g_RL * phi0~
Cr            =  integral g_RL * (phi0 = phi[n0])
Cr              +integral g_RL * (phi0~-phi0 = phi[n0~-n0])
Cr               The integral is partitioned to minimize mesh errors.
Cr               The first part is done by projecting g_RL to a mesh
Cr               and integrating the product g_RL*phi0 on the mesh
Cr               The second is done analytically by structure constants
Cr      hpot0 = integrals h_R * phi0~ (contributions from core)
Cr            = integrals h_R * (phi0 = phi[n0])
Cr             +integrals h_R * (phi0~-phi0 = phi[n0~-n0])
Cr       u00   :integral n0 phi[n0] = integral n0 phi0
Cr       u0g   :integral n0 [phi0~-phi0]
Cr       sgp0  :integral [n0~-n0] phi0~
Cr   Therefore :u00 + u0g + sgp0 = integral n0~ phi0~
Cr       smq   :integral n0
Cr       vconst:constant potential to be added to total.
Cr             :It is computed from average (v(RMT))
Cr       rhvsm :u00 + u0g + sgp0 + vconst*smq
Cr             := integral n0~ phi0~ + vconst*smq
Cr       zvnsm :integral core density * phi0~
Cr
Cr  Subroutines called by smves:
Cr    vesft    computes the electrostatic potential of n0 = phi0
Cr             (i.e. without the compensating gaussians).  This
Cr             is pretty trivial, since nabla^2 -> G^2 in G-space
Cr
Cr    vesgcm   1. makes the first term in gpot0
Cr                = integral g_RL * (phi0 = phi[n0])
Cr             2. makes the first term in hpot0
Cr             3. adds ves[n0~-n0] to the mesh estat potential
Cr
Cr    ugcomp   1. makes the second term in gpot0
Cr             2. makes the second term in hpot0
Cu Updates
Cb Bugs
Cb   Possible to make vval(l=0) for sites with lmxl=-1, which tells
Cb   value of ves at point.  However, vval doesn't have the
Cb   space allocated.  So skip for now.
Cu Updates
Cu   01 Jul 05 handle sites with lmxl=-1
Cu   19 Sep 02 (WRL) Added background term
Cu   24 Aug 01 Extended to calc vval.  Altered argument list.
Cu   20 Apr 01 Generates vrmt
Cu   21 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp ves_smooth.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,nbas,mode
      double precision qsmc,smq,rhvsm,sgp0,vconst,zsum,zvnsm,qbg
      real(8):: qmom(1) , f(3,nbas) , gpot0(1) , vval(1) , hpot0(nbas) 
     ., vrmt(nbas)
      type(s_lat)::slat
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      double complex smrho(k1,k2,k3,2),smpot(k1,k2,k3,2)
C ... Local parameters
         integer:: ib , igetss , ilm , ipr , iprint , is , iv0 , lfoc 
     .   , lgunit , lmxl , m , n1 , n2 , n3 , ng , ngabc(3) , nglob , 
     .   nlm , nsp , stdo , j1 , j2 , j3 
ckino Dec.12.2011:           integer,pointer :: iv_p_okv(:) =>NULL()

ckino Dec.12.2011:          real(8),pointer :: rv_p_ogv(:) =>NULL()

ckino Dec.12.2011:         real(8),pointer :: rv_p_ocy(:) =>NULL()

      complex(8) ,allocatable :: cg1_zv(:)
      complex(8) ,allocatable :: cgsum_zv(:)
      complex(8) ,allocatable :: cv_zv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision ceh,cofg,cofh,dgetss,hsum,pi,qcorg,qcorh,qsc,
     .rfoc,rmt,s1,s2,sbar,srfpi,sum,sum1,sum2,u00,u0g,ugg,usm,vbar,
     .vcnsto,vol,y0,z,R,eint
C ... Heap
       integer ::iwdummy 


C ... Setup
      call tcn('smves')
      ipr   = iprint()
      stdo  = lgunit(1)
Changenglob      nsp   = nglob('nsp')
      nsp   = globalvariables%nsp
      pi    = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0    = 1d0/srfpi

      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      ng=slat%ng
ckino Dec.12.2011:         rv_p_ogv => slat%rv_p_ogv 

ckino Dec.12.2011:         iv_p_okv => slat%iv_p_okv 

      vol=slat%vol


ckino Dec.12.2011:         rv_p_ocy => slat%rv_p_ocy 



C     Electrostatics depend only on total spin density
      if (nsp .eq. 2) then
        call daxpy(k1*k2*k3*2,1d0,smrho(1,1,1,2),1,smrho,1)
      endif

      allocate(cv_zv(ng))

      call dpzero(f, 3*nbas)

C ... FT of smooth density to reciprocal space
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,-1)

C     debugging ... one-center expansion of smrho
C     call msh21c(1,ssite,sspec,slat,ng,w(ogv),w(okv),k1,k2,k3,smrho)

C --- Estatic potential of smooth density without gaussians ---
ckino Dec.12.2011:          call vesft ( slat , ng , rv_p_ogv , iv_p_okv , cv_zv , k1 , k2 
ckino Dec.12.2011:          call vesft ( slat , ng , rv_p_ogv , slat%iv_p_okv , cv_zv , k1 , k2 
        call vesft ( slat , ng , slat%rv_p_ogv , slat%iv_p_okv , cv_zv , k1 , k2 
     .  , k3 , smrho , smpot , u00 ) 




C ... Add estatic potential of compensating gaussians to smpot
      allocate(cg1_zv(ng))

      allocate(cgsum_zv(ng))

ckino Dec.12.2011:           call vesgcm ( nbas , ssite , sspec , slat , rv_p_ocy , qmom , 
         call vesgcm ( nbas , ssite , sspec , slat , slat%rv_p_ocy , qmom , 
ckino Dec.12.2011:       .   ng , rv_p_ogv , iv_p_okv , cv_zv , cg1_zv , cgsum_zv , k1 , k2 
ckino Dec.12.2011:       .   ng , rv_p_ogv , slat%iv_p_okv , cv_zv , cg1_zv , cgsum_zv , k1 , k2 
     .   ng , slat%rv_p_ogv , slat%iv_p_okv , cv_zv , cg1_zv , cgsum_zv , k1 , k2 
     .   , k3 , smpot , f , gpot0 , hpot0 , qsmc , zsum , vrmt ) 




      if (ipr .ge. 40) write (stdo,230) (ib,(f(m,ib),m=1,3),ib=1,nbas)
  230 format(/' after vesgcomp: forces are:'/(i4,3f12.6))

C ... Compute e.s. potential at MT boundary
ckino Dec.12.2011:          call mshvmt ( nbas , ssite , sspec , slat , ng , rv_p_ogv , iv_p_okv 
ckino Dec.12.2011:          call mshvmt ( nbas , ssite , sspec , slat , ng , rv_p_ogv , slat%iv_p_okv 
        call mshvmt ( nbas , ssite , sspec , slat , ng , slat%rv_p_ogv , slat%iv_p_okv 
     .  , cv_zv , k1 , k2 , k3 , smpot , vval ) 



      call symvvl(nbas,ssite,sspec,slat,vval,vrmt)
      if (allocated(cgsum_zv)) deallocate(cgsum_zv)
      if (allocated(cg1_zv)) deallocate(cg1_zv)

C     call zprm3('smpot',0,smpot,n1,n2,n3)

C --- Make vbar = avg v(RMT) and optionally assign to vconst ---
      vbar = 0
      sbar = 0
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)

        rmt = (sspec(is)%rmt)

        vbar = vbar + rmt**2 * vrmt(ib)
        sbar = sbar + rmt**2
      enddo
      vbar = vbar/sbar
      vcnsto = vconst
      if (mode .ne. 0) vconst = -vbar

      if (ipr .ge. 20) write (stdo,232) vbar,vcnsto,vconst
  232 format(/' avg es pot at rmt=',f9.6,
     .'  avg sphere pot=',f9.6,'  vconst=',f9.6)

C ... Adjust vbar, vval, gpot0 by vconst
      iv0 = 0
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)

        lmxl = int(sspec(is)%lmxl)

        if (lmxl .gt. -1) then
          nlm = (lmxl+1)**2
          vrmt(ib) = vrmt(ib) + vconst
          vval(1+iv0) = vval(1+iv0) + vconst/y0
          gpot0(1+iv0) = gpot0(1+iv0) + vconst/y0
          iv0 = iv0 + nlm
        endif
      enddo
      if (ipr .ge. 40) then
        write (stdo,233)
  233   format(' average electrostatic potential at MT boundaries',
     .  ' after shift')
         call arrprt ( ' Site    ves' , '%,4i%:-3,6;6d' , 'Id' , nbas 
     .   , 0 , 4 , 0 , '  | ' , iwdummy , vrmt , iwdummy , iwdummy , iwdummy 
     .   , iwdummy , iwdummy , iwdummy ) 

      endif

C ... Back transform of density and potential to real-space mesh
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1)
      call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)

C ... Add background to smrho
      do j1=1,k1
        do j2=1,k2
          do j3=1,k3
            smrho(j1,j2,j3,1)=smrho(j1,j2,j3,1)+qbg/vol
          enddo
        enddo
      enddo

      if (qbg .ne. 0) then
        R = (3d0/pi/4d0*vol)**(1d0/3d0)
        eint = qbg*2*9d0/10d0/R
        call info(30,0,0,' cell interaction energy from homogeneous'//
     .  ' background (q=%d) is %;6,6d',qbg,eint)
      endif

C     Integral n0
      call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,sum1,sum2)
      smq = sum1
C     call mshint(vol,1,n1,n2,n3,k1,k2,k3,smpot,sum1,sum2)
C     Integral n0 phi0~
      call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smrho,smpot,s1,s2)
      u0g = s1 - u00

      call ugcomp(nbas,ssite,sspec,slat,qmom,gpot0,hpot0,ugg,f)

      if (ipr .ge. 50) write (stdo,231) (ib,(f(m,ib),m=1,3),ib=1,nbas)
  231 format(/' after ugcomp: forces are'/(i4,3f12.6))
      if (ipr .ge. 50) write(stdo,926) u00,u0g,ugg
  926 format(' u00,u0g,ugg=',3f14.6)

C --- Collect energy terms; make zvnuc for smooth problem ---
      zvnsm = 0d0
      rhvsm = u00 + u0g + vconst*smq
      sum = 0d0
      iv0 = 0
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)

        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        lmxl = int(sspec(is)%lmxl)

        if (lmxl .gt. -1) then
          nlm = (lmxl+1)**2
C       hsum = integral of charge in sm. Hankel
          hsum = -srfpi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
          hpot0(ib) = hpot0(ib) + vconst*hsum
          zvnsm = zvnsm + (qcorg-z)*y0*gpot0(iv0+1) + cofh*hpot0(ib)
          do  ilm = 1, nlm
            rhvsm = rhvsm + qmom(iv0+ilm)*gpot0(iv0+ilm)
            sum = sum + qmom(iv0+ilm)*gpot0(iv0+ilm)
          enddo
          iv0 = iv0+nlm
        endif
      enddo
      sgp0 = sum

c|      write(stdo,991) zvnsm,rhvsm,sum
c|  991 format(' zvnsm=',f12.6,'   rhvsm=',f12.6,
c|     .   /' sum over gpot0*qmom',f12.6)

      usm = 0.5d0*(rhvsm+zvnsm)

      if (ipr .ge. 30) write (stdo,500) usm,smq
  500 format(/' smooth rhoves',f14.6,'   charge',f13.6)

      if (allocated(cv_zv)) deallocate(cv_zv)

      call tcx('smves')

C ... subtract background
      do j1=1,k1
        do j2=1,k2
          do j3=1,k3
            smrho(j1,j2,j3,1)=smrho(j1,j2,j3,1)-qbg/vol
          enddo
        enddo
      enddo
      smq=smq-qbg

C     Restore spin 1 density, copy potential to second spin channel
      if (nsp .eq. 2) then
        call daxpy(k1*k2*k3*2,-1d0,smrho(1,1,1,2),1,smrho,1)
        call dcopy(k1*k2*k3*2,smpot,1,smpot(1,1,1,2),1)
      endif

      end subroutine smves


