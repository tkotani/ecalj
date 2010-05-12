Cgetarg...info...           structure ['mkorbm', 'site', 'ssite']
Cgetarg...info...           structure ['mkorbm', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['mkorbm', 'm_struc_def']
      subroutine mkorbm(ssite,sspec,isp,nsp,nspc,nlmax,ndham,nev,
     .  wtkp,iq,nbas,ppnl,aus,nl,nkp,orbtm)

      use m_struc_def  !Cgetarg

C- Decomposition of norm from projection of w.f into MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa rmt
Ci     Stored:
Ci     Passed to:
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nlmax :leading dimension of aus
Ci   ndham :dimensions aus,wtkp
Ci   nev   :number of eigenvectors to accumulate orbital moment
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   iq    :current k-point
Ci   nbas  :size of basis
Ci   ppnl  :NMTO potential parameters; see eg potpus.f
Ci   aus   :values of (phi,phidot) MT sphere boundary; see makusq
Ci   nl    :(global maximum l) + 1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Co Outputs
Co   orbtm :orbital moments accumulated for this qp
Cl Local variables
Cl   ispc  :the current spin index in the coupled spins case.
Cl         :Some quantities have no separate address space for each
Cl         :spin in the indepedent-spins case (evec,evl,ewgt) but do
Cl         :in the coupled-spins case.  A separate loop ispc=1..nspc
Cl         :must be added for the latter case
Cl         :ispc is the appropriate index for objects which distinguish
Cl         :spins in the spin-coupled case only
Cl   isp   :isp  is the appropriate index for objects which distinguish
Cl         :spins in the spin-uncoupled case only
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Cl         :ksp is appropriate spin index for quantities that have
Cl         :separate address space for each spin in every case
Cl         :(potential- and density-like objects).
Cu Updates
Cu   25 Apr 05 (A. Chantis) extended to local orbitals
Cu   24 Dec 04 Extended to spin-coupled case
Cu   30 Aug 04 (A. Chantis) first written, adapted from mkpdos
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer isp,nsp,nspc,nlmax,ndham,nbas,nev
      integer n0,nppn,nl,nkp,nab,iq
      parameter (n0=10,nppn=12,nab=9)
Cgetarg       double precision ssite(1),sspec(1),wtkp(ndham,nsp,nkp)
      real(8):: wtkp(ndham,nsp,nkp)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      double precision ppnl(nppn,n0,nsp,*),orbtm(nl,nsp,*)
      double complex aus(nlmax,ndham*nspc,3,nsp,*)
C ... Local parameters
      integer lmxa,lmxax,lmdim,ichan,ib,is,igetss,iv,ilm,l,m,nlma,
     . ll,lc,em,ispc,ksp
      double precision suml(11),s11,s22,s12,s33,s31,s32,s13,s23,
     .  suma,rmt,sab(nab,n0,2)
      double complex au,as,az,iot

      lmxax = ll(nlmax)
      iot = dcmplx(0d0,1d0)
      ichan = 0
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
        is = int(ssite(ib)%spec)

Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec lmxa rmt',sspec,is,lmxa,rmt,0,0)

        lmxa=sspec(is)%lmxa
        rmt=sspec(is)%rmt

        lmxa = min(lmxa,lmxax)
        if (lmxa .eq. -1) goto 10

        nlma = (lmxa+1)**2
        lmdim = nlma

        call phvsfp(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)

C       In noncollinear case, isp=1 always => need internal ispc=1..2
C       ksp is the current spin index in both cases:
C       ksp = isp  in the collinear case
C           = ispc in the noncollinear case
C       ispc=1 for independent spins, and spin index when nspc=2
        do  ispc = 1, nspc
          ksp = max(ispc,isp)
          do  iv = 1, nev
            call dpzero(suml,1+lmxa)
            suma = 0
            ilm = 0

C  ....  Rotate from real to spherical harmonics (order assumed: m,...,-m).
C        |Psi>_l = \Sum_{m}(A_l,m * u_l + B_l,m * s_l)*R_l,m --->
C        |Psi>_l = \Sum_{m}(C_l,m * u_l + D_l,m * s_l)*Y_l,m
C        R_l,m and Y_l,m are the real and spherical harmonics respectively.
C
C              | (-1)^m/sqrt(2)*A_l,-m + i*(-1)^m/sqrt(2)*A_l,m , m>0
C        C_l,m=|  A_l,m                                         , m=0
C              |  1/sqrt(2)*A_l,-m -  i*1/sqrt(2)*A_l,m         , m<0
C
C       Same relationships are valid between D and B.

            do  l = 0, lmxa
              lc = (l+1)**2 - l
              do  m = -l, l
                em = abs(m)
                ilm = ilm+1
C    ...    m,...,-m order
C            if (m .lt. 0) then
C            au = 1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) -
C     .           iot/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
C            as = 1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) -
C     .           iot/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
C            az = 1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) -
C     .           iot/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
C            else if (m .gt. 0) then
C            au = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) +
C     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
C            as = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) +
C     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
C            az = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) +
C     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
C            else
C    ...   -m,...,m order
                if (m .lt. 0) then
                  au = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) + 
     .           1d0/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
                  as = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) + 
     .           1d0/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
                  az = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) + 
     .           1d0/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
                else if (m .gt. 0) then
                  au = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) + 
     .           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
                  as = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) + 
     .           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
                  az = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) + 
     .           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
                else
                  au = aus(ilm,iv,1,ksp,ib)
                  as = aus(ilm,iv,2,ksp,ib)
                  az = aus(ilm,iv,3,ksp,ib)
                endif

C           If (au,as) are coefficients to (u,s), use this
                s11 = dconjg(au)*au*sab(1,l+1,ksp)
                s12 = 2*dconjg(au)*as*sab(2,l+1,ksp)
                s22 = dconjg(as)*as*sab(4,l+1,ksp)
                s33 = dconjg(az)*az*sab(5,l+1,ksp)
                s31 = dconjg(az)*au*sab(6,l+1,ksp)
                s32 = dconjg(az)*as*sab(7,l+1,ksp)
                s13 = dconjg(au)*az*sab(6,l+1,ksp)
                s23 = dconjg(as)*az*sab(7,l+1,ksp)

                orbtm(l+1,ksp,ib)=orbtm(l+1,ksp,ib)+m*(s11+s12+s22+
     .     s33+s32+s23+s31+s13)*wtkp(iv,isp,iq)

              enddo
            enddo

          enddo
        enddo
C        print*, isp, ib, 'ORB.MOMNT=',orbtm(isp,ib)
   10   continue
      enddo
Cgetarg       end
      end subroutine mkorbm











