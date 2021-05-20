      subroutine pvsms2 ( ssite , sspec , rotm , nbas , nsp , sv_p_orhoat )
      use m_struc_def  
      use m_lmfinit,only: stdo


C- Rotate local densities by specified rotation
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: nr lmxl
Ci     Stored:    name
Ci     Passed to: spacks
Ci   rotm  :3x3 cartesian rotation matrix
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   orhoat:vector of offsets containing site density
Co Outputs
Co   orhoat:On output the different m-channels of rhoat(1) and rhoat(2)
Co         :are mixed by the rotation
Cl Local variables
Cr Remarks
Cr   For a rotation matrix R, The density is stored in the 1-center form
Cr      rho_l(r) YL(rhat)
Cr   Given a rotation matrix R, this it transforms as
Cr      rho_l(r) YL(R rhat) = rho_l(r) rYL(rhat)
Cr   where rYL is made by ylmrtg
Cr
Cb Bugs
Cb   No ability is supplied when the Yl are true instead of real
Cb   spherical harmonics
Cu Updates
Cu   21 Dec 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer:: nbas , nsp
      type(s_rv1) :: sv_p_orhoat(3,1)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      double precision rotm(3,3)
      integer ib,i,j,is,lmxl,igetss,nr,nlml,ipr,nlx,nl2,nglob
      parameter (nlx=9, nl2=nlx*nlx)
      double precision rYL(nl2,nl2)
      character spid*8
      call getpr(ipr)
C ... Rotation matrix for real spherical harmonics
      call ylmrtg(nl2,rotm,rYL)
C --- For each site and l, rotate the m-components ---
      if (ipr .ge. 20) then
        call info0(20,0,0,' Rotate local densities using R=')
        write (stdo,350) ((rotm(i,j),j=1,3),i=1,3)
  350   format(3f11.6)
      endif
      do  10  ib = 1, nbas
        is = int(ssite(ib)%spec)
        spid=sspec(is)%name
        nr=sspec(is)%nr
        lmxl=sspec(is)%lmxl
        if (lmxl .eq. -1) goto 10
        nlml = (lmxl+1)**2
        if (nlml .gt. nl2) call rx('increase nl2 in pvsms2')
        call pvsms3 ( nr , nr , nlml , nsp , ryl , nl2 , sv_p_orhoat( 1 , ib )%v )
        call pvsms3 ( nr , nr , nlml , nsp , ryl , nl2 , sv_p_orhoat( 2 , ib )%v )
   10 continue
      end subroutine pvsms2

      subroutine pvsms3(nrx,nr,nlml,nsp,rYL,nl2,rho)
C- Rotation of an l-dependent density
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrx   :leading dimension of rho
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cl   rYL   :rotation matrix that rotates Y_lm
Ci   nl2   :leading dimension of rYL
Co Outputs
Co   rho   :On output the different m-channels of rho are
Co         :mixed by rYL
C ----------------------------------------------------------------------
      implicit none
      integer nrx,nr,nlml,nsp,nl2
      double precision rho(nrx,nlml,nsp)
      double precision rYL(nl2,nl2)
      integer isp
      double precision rwk(nrx,nlml)
      if (nlml .eq. 0) return
      do  isp = 1, nsp
        call dgemm('N','T',nr,nlml,nlml,1d0,rho(1,1,isp),nrx,
     .  rYL,nl2,0d0,rwk,nrx)
        call dcopy(nrx*nlml,rwk,1,rho(1,1,isp),1)
      enddo
      end subroutine pvsms3


