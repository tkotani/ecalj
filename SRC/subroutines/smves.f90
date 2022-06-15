subroutine smves(nbas,ssite,sspec,k1,k2,k3,qmom,gpot0, & 
     vval,hpot0,sgp0,smrho,smpot,vconst,smq,qsmc,f,rhvsm,zvnsm,zsum, &
     vrmt,qbg)
  use m_supot,only: iv_a_okv, rv_a_ogv
  use m_struc_def  !Cgetarg
  use m_lmfinit,only: rv_a_ocy,nsp,stdo
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc
  use m_supot,only: lat_ng
  use m_MPItk,only: master_mpi
  use m_ext,only: sname
  !- Electrostatic potential of the smooth density.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode=0 use input vconst
  !ixxxx     (removed  mode=1 generate vconst as - average v(RMT))
  !i   nbas  :size of basis
  !i   ssite :struct containing site-specific information
  !i   sspec :struct containing species-specific information
  !i   slat  :struct containing information about the lattice
  !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
  !i   qmom  :multipole moments of on-site densities (rhomom.f)
  !i   smrho :smooth density on real-space mesh
  !i   qbg   : back ground charge
  ! o Inputs/Outputs
  ! o  vconst:constant potential to be added to total
  ! o        :On input  vconst is set to a default value
  ! o        :On output vconst may be set to the average estat
  ! o        :          at the MT boundary.
  !o Outputs (see also Remarks)
  !o   gpot0 :integrals of compensating gaussians g_RL * phi0~
  !o         :For accuracy, integral is split into
  !o         :g_RL phi0 (vesgcm) + g_RL phi [n0~-n0] (ugcomp)
  !o         :vesgcm projects g_RL to the mesh to do the integral
  !o         :ugcomp does its integrals analytically (structure constants)
  !o         :NB: There is a local analog of gpot0 generated in locpt2.
  !o   vval  :coffs to YL expansion of es potential at MT boundary
  !o   hpot0 :integrals of semicore smooth Hankels * phi0~
  !o   sgp0  :sgp0 = sum_RL integral qmom_RL g_RL phi0~
  !o         :     = integral [n0~-n0 ] phi0~
  !o   smpot :smooth potential phi0~ (includes compensating gaussians)
  !o   smq   :integral of smooth density n0
  !o   qsmc  :pseudocore charge
  !o   f     :electrostatic contribution to force.
  !o   rhvsm :integral n0~ [phi0~ + vconst]
  !o         :(electrostatic energy of sm. density n0~) + vconst*smq
  !o   zvnsm :integral (qcorg-z + rhoc) phi0~
  !o   vrmt  :electrostatic potential at rmt, with G=0 term in smpot=0
  !l Local variables
  !l   u00   :integral n0 phi[n0] = n0 phi0
  !l   u0g   :integral n0 [phi0~-phi0]
  !l   ugg   :integral [n0~-n0] [phi0~-phi0]
  !l         :Note: ugg is not used.
  !r Remarks
  !r  The total density is a sum of three terms,
  !r
  !r    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
  !r
  !r  The first term is the smooth density on a mesh of points; the
  !r  second is the true density and is defined on a radial mesh for each
  !r  sphere; the last is the 1-center expansion of the smooth density on
  !r  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
  !r  identical to the one-center expansion of n0(mesh).  The sum of the
  !r  three terms converges rapidly with l because errors in n_R(r) are
  !r  mostly canceled by errors in n0_R(r).)
  !r
  !r  We add and subtract a set of compensating gaussian orbitals
  !r
  !r    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
  !r
  !r  which render the integral of the local part (the last 3 terms)
  !r  zero in each RL channel.  The g_RL must be localized enough that
  !r  their spillout beyond the MT radius is negligible.
  !r
  !r  We define
  !r
  !r    n0~ = n0 + compensating gaussians sum_RL Q_RL g_RL
  !r
  !r  In the interstitial, the electrostatic potential of n0~ is the true
  !r  estat potential.  The potential of n0 is called phi0 and the
  !r  potential of n0~ is called phi0~.  The total electrostatic energy
  !r  is computed as
  !r
  !r    the electrostatic energy of  n0~ + integral n0*vconst +
  !r    the electrostatic energy of (neutral) local parts
  !r
  !r  vconst may either be passed as an input (mode=0) or it is
  !r  generated here as the average ves(RMT).
  !r  This routine computes the estat potential and energy from the
  !r  first two terms.  Some variables used in smves and its subroutines:
  !r    Let n0  = smooth density without the compensating sum_RL Q_RL g_RL
  !r        n0~ = n0 + sum_RL Q_RL g_RL
  !r      phi0  = ves[n0]
  !r      phi0~ = ves[n0~]
  !r      g_RL  = gaussian in RL channel
  !r      h_R   = l=0 sm hankel in RL channel, (represents core densities)
  !r    qmom_RL = multipole moment in RL channel of (n_R(r) - n0_R(r))
  !r              so that int n_RL(r)-n0_RL(r) = qmom_RL * g_RL(r)
  !r      gpot0 = vector of integrals g_RL * phi0~
  !r            =  integral g_RL * (phi0 = phi[n0])
  !r              +integral g_RL * (phi0~-phi0 = phi[n0~-n0])
  !r               The integral is partitioned to minimize mesh errors.
  !r               The first part is done by projecting g_RL to a mesh
  !r               and integrating the product g_RL*phi0 on the mesh
  !r               The second is done analytically by structure constants
  !r      hpot0 = integrals h_R * phi0~ (contributions from core)
  !r            = integrals h_R * (phi0 = phi[n0])
  !r             +integrals h_R * (phi0~-phi0 = phi[n0~-n0])
  !r       u00   :integral n0 phi[n0] = integral n0 phi0
  !r       u0g   :integral n0 [phi0~-phi0]
  !r       sgp0  :integral [n0~-n0] phi0~
  !r   Therefore :u00 + u0g + sgp0 = integral n0~ phi0~
  !r       smq   :integral n0
  !r       vconst:constant potential to be added to total.
  !r             :It is computed from average (v(RMT))
  !r       rhvsm :u00 + u0g + sgp0 + vconst*smq
  !r             := integral n0~ phi0~ + vconst*smq
  !r       zvnsm :integral core density * phi0~
  !r
  !r  Subroutines called by smves:
  !r    vesft    computes the electrostatic potential of n0 = phi0
  !r             (i.e. without the compensating gaussians).  This
  !r             is pretty trivial, since nabla^2 -> G^2 in G-space
  !r
  !r    vesgcm   1. makes the first term in gpot0
  !r                = integral g_RL * (phi0 = phi[n0])
  !r             2. makes the first term in hpot0
  !r             3. adds ves[n0~-n0] to the mesh estat potential
  !r
  !r    ugcomp   1. makes the second term in gpot0
  !r             2. makes the second term in hpot0
  !u Updates
  !b Bugs
  !b   Possible to make vval(l=0) for sites with lmxl=-1, which tells
  !b   value of ves at point.  However, vval doesn't have the
  !b   space allocated.  So skip for now.
  !u Updates
  !u   01 Jul 05 handle sites with lmxl=-1
  !u   19 Sep 02 (WRL) Added background term
  !u   24 Aug 01 Extended to calc vval.  Altered argument list.
  !u   20 Apr 01 Generates vrmt
  !u   21 Jun 00 spin polarized
  !u   22 Apr 00 Adapted from nfp ves_smooth.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: k1,k2,k3,nbas,i_copy_size
  real(8):: qmom(1) , f(3,nbas) , gpot0(1) , vval(1) , hpot0(nbas) &
       , vrmt(nbas),qsmc,smq,rhvsm,sgp0,vconst,zsum,zvnsm,qbg
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  complex(8):: smrho(k1,k2,k3,2),smpot(k1,k2,k3,2)
  integer:: ib , igetss , ilm , ipr , iprint , is , iv0 , lfoc &
       , lgunit , lmxl , m , n1 , n2 , n3 , ng , ngabc(3) , nglob , nlm , j1 , j2 , j3
  complex(8) ,allocatable :: cg1_zv(:)
  complex(8) ,allocatable :: cgsum_zv(:)
  complex(8) ,allocatable :: cv_zv(:)
  double precision :: ceh,cofg,cofh,dgetss,hsum,pi,qcorg,qcorh,qsc, &
       rfoc,rmt,s1,s2,sbar,srfpi,sumx,sum1,sum2,u00,u0g,ugg,usm,vbar, &
       vcnsto,vol,y0,z,R,eint
  integer ::iwdummy, ifivsmconst
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  ! ... Setup
  call tcn('smves')
  ipr   = iprint()
  !      stdo  = lgunit(1)
  !      nsp   = globalvariables%nsp
  pi    = 4d0*datan(1d0)
  srfpi = dsqrt(4d0*pi)
  y0    = 1d0/srfpi
  !      i_copy_size=size(lat_nabc)
  !     call icopy(i_copy_size,lat_nabc,1,ngabc,1)
  ngabc=lat_nabc
  ng=lat_ng
  vol=lat_vol

  !     Electrostatics depend only on total spin density
  if (nsp == 2) then
     call daxpy(k1*k2*k3*2,1d0,smrho(1,1,1,2),1,smrho,1)
  endif
  allocate(cv_zv(ng))
  call dpzero(f, 3*nbas)
  ! ... FT of smooth density to reciprocal space
  call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,-1)

  !     debugging ... one-center expansion of smrho
  !     call msh21c(1,ssite,sspec,slat,ng,w(ogv),w(okv),k1,k2,k3,smrho)

  ! --- Estatic potential of smooth density without gaussians ---
  call vesft ( ng , rv_a_ogv , iv_a_okv , cv_zv , k1 , k2  &
       , k3 , smrho , smpot , u00 )

  ! ... Add estatic potential of compensating gaussians to smpot
  allocate(cg1_zv(ng))
  allocate(cgsum_zv(ng))
  call vesgcm ( nbas , ssite , sspec ,  rv_a_ocy , qmom , &
       ng , rv_a_ogv , iv_a_okv , cv_zv , cg1_zv , cgsum_zv , k1 , k2 &
       , k3 , smpot , f , gpot0 , hpot0 , qsmc , zsum , vrmt )
  if (ipr >= 40) write (stdo,230) (ib,(f(m,ib),m=1,3),ib=1,nbas)
230 format(/' after vesgcomp: forces are:'/(i4,3f12.6))

  ! --- M. OBATA check
  call esmsmves( nbas , ssite , sspec ,  rv_a_ocy , qmom , &
       ng , rv_a_ogv , iv_a_okv , cv_zv , cg1_zv , cgsum_zv , k1 , k2 &
       , k3 , smrho, qbg, smpot , f , gpot0 , hpot0 , qsmc , zsum , vrmt )
  ! ... Compute e.s. potential at MT boundary
  call mshvmt ( nbas , ssite , sspec ,  ng , rv_a_ogv , iv_a_okv &
       , cv_zv , k1 , k2 , k3 , smpot , vval )
  call symvvl(nbas,ssite,sspec,vval,vrmt) 
  if (allocated(cgsum_zv)) deallocate(cgsum_zv)
  if (allocated(cg1_zv)) deallocate(cg1_zv)
  !     call zprm3('smpot',0,smpot,n1,n2,n3)

  ! --- Make vbar = avg v(RMT) and optionally assign to vconst ---
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
  vconst = -vbar !  if (mode /= 0) vconst = -vbar
  if (ipr >= 20) write (stdo,232) vbar,vcnsto,vconst
232 format(' smves:: avg es pot at rmt=',f9.6,'  avg sphere pot=',f9.6,'  vconst=',f9.6)
!  if(mode==0)call rx(' vsmconst is not vconst due to estatic. need to implement something!')
  if(master_mpi) then
     open(newunit=ifivsmconst,file='vessm.'//trim(sname))
     write(ifivsmconst,"(d23.15)") vconst
     close(ifivsmconst)
  endif

  ! ... Adjust vbar, vval, gpot0 by vconst
  iv0 = 0
  do  ib = 1, nbas
     is = int(ssite(ib)%spec)
     lmxl = int(sspec(is)%lmxl)
     if (lmxl > -1) then
        nlm = (lmxl+1)**2
        vrmt(ib) = vrmt(ib) + vconst
        vval(1+iv0) = vval(1+iv0) + vconst/y0
        gpot0(1+iv0) = gpot0(1+iv0) + vconst/y0
        iv0 = iv0 + nlm
     endif
  enddo
  if (ipr >= 40) then
     write (stdo,233)
233  format(' average electrostatic potential at MT boundaries after shift')
     write(stdo, "(a)") ' Site    ves'
     do ib=1,nbas
        write(stdo,"(i4, f12.6)")ib,vrmt(ib)
     enddo
  endif
  ! ... Back transform of density and potential to real-space mesh
  call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1)
  call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)

  ! ... Add background to smrho
  do j1=1,k1
     do j2=1,k2
        do j3=1,k3
           smrho(j1,j2,j3,1)=smrho(j1,j2,j3,1)+qbg/vol
        enddo
     enddo
  enddo

  if (qbg /= 0) then
     R = (3d0/pi/4d0*vol)**(1d0/3d0)
     eint = qbg*2*9d0/10d0/R
     call info(30,0,0,' cell interaction energy from homogeneous'// &
          ' background (q=%d) is %;6,6d',qbg,eint)
  endif

  !     Integral n0
  call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,sum1,sum2)
  smq = sum1
  !     call mshint(vol,1,n1,n2,n3,k1,k2,k3,smpot,sum1,sum2)
  !     Integral n0 phi0~
  call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smrho,smpot,s1,s2)
  u0g = s1 - u00

  call ugcomp(nbas,ssite,sspec,qmom,gpot0,hpot0,ugg,f) 
  !      write(6,*)'sumcheck=',sum(smrho(:,:,:,1)),sum(smpot(:,:,:,1))

  if (ipr >= 50) write (stdo,231) (ib,(f(m,ib),m=1,3),ib=1,nbas)
231 format(/' after ugcomp: forces are'/(i4,3f12.6))
  if (ipr >= 50) write(stdo,926) u00,u0g,ugg
926 format(' u00,u0g,ugg=',3f14.6)

  ! --- Collect energy terms; make zvnuc for smooth problem ---
  zvnsm = 0d0
  rhvsm = u00 + u0g + vconst*smq
  sumx = 0d0
  iv0 = 0
  do  ib = 1, nbas
     is = int(ssite(ib)%spec)

     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     lmxl = int(sspec(is)%lmxl)

     if (lmxl > -1) then
        nlm = (lmxl+1)**2
        !       hsum = integral of charge in sm. Hankel
        hsum = -srfpi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
        hpot0(ib) = hpot0(ib) + vconst*hsum
        zvnsm = zvnsm + (qcorg-z)*y0*gpot0(iv0+1) + cofh*hpot0(ib)
        do  ilm = 1, nlm
           rhvsm = rhvsm + qmom(iv0+ilm)*gpot0(iv0+ilm)
           sumx = sumx + qmom(iv0+ilm)*gpot0(iv0+ilm)
        enddo
        iv0 = iv0+nlm
     endif
  enddo
  sgp0 = sumx

  !|      write(stdo,991) zvnsm,rhvsm,sum
  !|  991 format(' zvnsm=',f12.6,'   rhvsm=',f12.6,
  !|     .   /' sum over gpot0*qmom',f12.6)

  usm = 0.5d0*(rhvsm+zvnsm)

  if (ipr >= 30) write (stdo,500) usm,smq
500 format(/' smooth rhoves',f14.6,'   charge',f13.6)

  if (allocated(cv_zv)) deallocate(cv_zv)

  call tcx('smves')

  ! ... subtract background
  do j1=1,k1
     do j2=1,k2
        do j3=1,k3
           smrho(j1,j2,j3,1)=smrho(j1,j2,j3,1)-qbg/vol
        enddo
     enddo
  enddo
  smq=smq-qbg

  !     Restore spin 1 density, copy potential to second spin channel
  if (nsp == 2) then
     call daxpy(k1*k2*k3*2,-1d0,smrho(1,1,1,2),1,smrho,1)
     call dcopy(k1*k2*k3*2,smpot,1,smpot(1,1,1,2),1)
  endif

end subroutine smves


