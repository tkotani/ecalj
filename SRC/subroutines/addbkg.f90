subroutine addbkgsm(smrho,k1,k2,k3,nsp,qbg,vol,fac)
  !- add uniform background to smooth charge density
  !i Inputs: smrho: smooth density
  !i         qbg: background charge
  !i         k1,k2,k3, dimensions of smrho mesh
  !i         nsp: number of spins
  !i         vol: volume per cell
  !i         fac: fac*qbg/vol//nsp is added to
  !-------------------------------------------
  implicit none
  integer :: j1,j2,j3,k1,k2,k3,nsp,is
  double precision :: qbg,rhobg,vol,fac
  double complex smrho(k1,k2,k3,nsp)
  rhobg=qbg/vol
  smrho(:,:,:,:)=smrho(:,:,:,:)+rhobg*fac/nsp
end subroutine addbkgsm

subroutine adbkql( sv_p_orhoat , nbas , nsp , qbg , vol , fac )
  use m_struc_def
  use m_lmfinit,only: ispec, sspec=>v_sspec
  !- Add uniform bkg charge density to local smooth rho
  !i orhoat: pointers to local density in spheres
  !i nbas: number of atoms in basis
  !i qbg: background charge
  !i sspec: species structure
  !i nsp: spins
  !i vol: vol of cell
  !i fac: fac * backg density is added
  !u Updates
  !u   01 Jul 05 Zero-radius sites skipped over
  !----------------------------------------
  implicit none
  integer :: nrmx,nlmx,nlml,lmxl,nbas
  parameter (nrmx=1501,nlmx=64)
  integer:: nsp
  type(s_rv1) :: sv_p_orhoat(3,nbas)
  real(8):: qbg , fac
  integer :: ib,nr,is
  double precision :: rhobkg,vol,a,rmt,rofi(nrmx)
  rhobkg = fac*qbg/vol
  do  ib = 1, nbas
     is=ispec(ib) !ssite(ib)%spec !specie(ib)
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     lmxl=sspec(is)%lmxl
     if (lmxl == -1) goto 10
     nlml=(lmxl+1)**2
     call rxx(nr .gt. nrmx,  'addbkgloc: increase nrmx')
     call rxx(nlml .gt. nlmx,'addbkgloc: increase nlmx')
     call radmsh(rmt,a,nr,rofi)
     call addbkgl(sv_p_orhoat(1,ib )%v,sv_p_orhoat(2,ib)%v, rhobkg , nr , nsp , rofi , nlml )
10   continue
  enddo
end subroutine adbkql

subroutine addbkgl(rho1,rho2,rhobkg,nr,nsp,rofi,nlml)
  ! adds uniform background to local smooth density at this site for l=0 component (ilm=1) 
  implicit none
  integer :: nsp,is,nr,nlml,i
  real(8):: rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rofi(nr),rhobkg
  real(8),parameter:: pi = 4d0*datan(1d0),srfpi = dsqrt(4*pi)
  do is = 1, nsp
     rho1(:,1,is) = rho1(:,1,is)+srfpi*rofi(:)**2*rhobkg/nsp
     rho2(:,1,is) = rho2(:,1,is)+srfpi*rofi(:)**2*rhobkg/nsp
  enddo
end subroutine addbkgl




