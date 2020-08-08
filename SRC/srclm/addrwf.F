      module m_addrwf
      contains
      subroutine addrwf(mode,z,l,v,ndg,n1,nr,rofi,rwgt,eadd,ev,fac,gadd,
     .g,s)
C- Add constant * radial wave function to another radial wave function
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :0 use both large and small components of radial w.f.
Ci         :1 use both large component of radial w.f. only
Ci   z     :nuclear charge
Ci         :(used only to compute overlap s, mode=0)
Ci   l     :l-quantum number
Ci   v     :spherical potential, without nuclear part
Ci         :(used only to compute overlap s, mode=0)
Ci   ndg   :leading dimension of g and gadd
Ci   n1    :if 0<n1<=nr, rwgt(n1) is scaled by 2
Ci         :(see vxtrap.f)
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights for numerical integration
Ci   ev    :eigenvalue of input wave function g
Ci         :(used only to compute overlap s, mode=0)
Ci   eadd  :eigenvalue of wave function gadd
Ci         :(used only to compute overlap s, mode=0)
Ci   fac   :Add fac*gadd into g
Ci   gadd  :See fac
Co Inputs/Outputs
Co   g     :g is overwritten by g + fac*g
Co   s     :overlap between gadd and new g
Cr Remarks
Cr   Input g and gadd are assumed to be solutions of the Schrodinger
Cr   equation with eigenvalues ev and eadd.  (For the scalar
Cr   relativistic case, the inner product depends slightly
Cr   on z,v, and eigenvalue)
Cu Updates
Cu   04 Sep 04 Adapted to extended local orbitals
Cu   12 Jul 04 ndg,n1 arguments (altered argument list)
Cu   14 Feb 02 New routine
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,ndg,n1,nr,mode
#define USEOPTIONAL
#ifdef USEOPTIONAL
      real(8):: fac,rofi(nr),rwgt(nr),
     .gadd(ndg,2),g(ndg,2),s
      real(8),optional:: z,v(nr),ev,eadd
#else
      double precision fac,z,rofi(nr),rwgt(nr),v(nr),ev,eadd,
     .gadd(ndg,2),g(ndg,2),s
#endif
C ... Local parameters
      integer ir
      double precision cc,vi,fllp1,gf11,gf22,gf12,r,tmc
      common /cc/ cc

      fllp1 = l*(l+1)
      s = 0
      if (n1 .gt. 0 .and. n1 .lt. nr) rwgt(n1) = 2*rwgt(n1)

      if (mode .eq. 0) then

        do  ir = 2, nr
          r = rofi(ir)

          if (fac .ne. 0) then
            g(ir,1) = g(ir,1) + fac*gadd(ir,1)
            g(ir,2) = g(ir,2) + fac*gadd(ir,2)
          endif

C       Rest of the loop computes overlap between new g and gadd
          vi = v(ir) - 2d0*z/r
          tmc = cc - (vi-ev)/cc
          gf11 = 1d0 + fllp1/(tmc*r)**2
          tmc = cc - (vi-eadd)/cc
          gf22 = 1d0 + fllp1/(tmc*r)**2
          gf12 = (gf11 + gf22)/2

          s = s + rwgt(ir)*(gf12*g(ir,1)*gadd(ir,1) + g(ir,2)*gadd(ir,2))
        enddo

      else
        do  ir = 2, nr
          r = rofi(ir)
          if (fac .ne. 0) then
            g(ir,1) = g(ir,1) + fac*gadd(ir,1)
          endif
          s = s + rwgt(ir)*g(ir,1)*gadd(ir,1)
        enddo
      endif

      if (n1 .gt. 0 .and. n1 .lt. nr) rwgt(n1) = rwgt(n1)/2
      end subroutine addrwf
      end module m_addrwf

      subroutine wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz,
     .phzp,dphzp,pnz,rsml,ehl,g0,g1,gz)
C- Add a linear combination of two w.f. to a 3rd to make a local orbital
C ----------------------------------------------------------------------
Ci Inputs
Ci   l     :l quantum number
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   phi   :1st wave function at rmax, i.e. r*g0
Ci   dphi  :radial derivative of r*g0
Ci   phip  :2nd wave function at rmax, i.e. r*g1
Ci   dphip :radial derivative of r*g1
Ci   g0    :1st radial w.f. to which to orthogonalize gz
Ci   g1    :2nd radial w.f. to which to orthogonalize gz
Ci   ehl   :energy of smoothed Hankel tail for extended local orbital
Ci   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
Cio Inputs/Outputs
Cio  phz   :3rd wave function at rmax, i.e. r*gz
Cio        :Input for standard local orbital; scaled for extended loc. orb
Cio  dphz  :radial derivative of r*gz
Cio        :Input for standard local orbital; scaled for extended loc. orb
Cio  phzp  :energy derivative of phz
Cio        :scaled for extended local orbital
Cio  dphzp :energy derivative of dphz
Cio        :scaled for extended local orbital
Cio  gz    :on input, radial w.f.
Cio        :on output, gz is overwritten by:
Cio        :(gz - alpha g0 - beta g1) so that value and slope
Cio        :of the result are zero at rmt, standard local orbital
Cio        :scaled to match sm. hankel at rmt, extended loc. orb.
Cu Updates
Cu   04 Sep 04 Adapted to extended local orbitals
Cu   06 Mar 02 New routine
C ----------------------------------------------------------------------
C     implicit none
      use m_addrwf
C ... Passed parameters
      integer l,nr
      double precision a,rofi(nr),rwgt(nr),rsml(*),ehl(*)
      double precision phi,dphip,dphi,phip,phz,dphz,phzp,dphzp,pnz(*)
      double precision g0(nr,2),g1(nr,2),gz(nr,2)
C ... Local parameters
      integer k,lpzi,nrmx
      parameter (nrmx=1501)
      double precision det,au,bu,as,bs,fac,x,xx
      double precision gzbig(nrmx*2)

      k = l+1
      lpzi = 0
      if (pnz(k) .gt.  0) lpzi = 1
      if (pnz(k) .ge. 10) lpzi = 2
      if (pnz(k) .ge. 20) lpzi = 3

      if (lpzi .eq. 0) then
        return

      elseif (lpzi .eq. 1) then
        det = phi*dphip - dphi*phip
        au = dphip/det
        bu = -dphi/det
        as = -phip/det
        bs = phi/det
        fac = phz*au + dphz*as
        x = 0
#ifdef USEOPTIONAL
        call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g0(1,2),g=gz(1,2),
     .  s=xx)
        call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g0,g=gz,s=xx)
#else
        call addrwf(1,x,l,x,nr,0,nr,rofi,rwgt,x,x,-fac,g0(1,2),gz(1,2),
     .  xx)
        call addrwf(1,x,l,x,nr,0,nr,rofi,rwgt,x,x,-fac,g0,gz,xx)
#endif

        fac = phz*bu + dphz*bs
#ifdef USEOPTIONAL
        call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g1(1,2),g=gz(1,2),
     .  s=xx)
        call addrwf(mode=1,l=l,ndg=nr,n1=0,nr=nr,rofi=rofi,rwgt=rwgt,fac=-fac,gadd=g1,g=gz,s=xx)
#else
        call addrwf(1,x,l,x,nr,0,nr,rofi,rwgt,x,x,-fac,g1(1,2),gz(1,2),
     .  xx)
        call addrwf(1,x,l,x,nr,0,nr,rofi,rwgt,x,x,-fac,g1,gz,xx)
#endif

      elseif (lpzi .eq. 2 .or. lpzi .eq. 3) then
        call dcopy(nr,gz,1,gzbig,1)
        call dcopy(nr,gz(1,2),1,gzbig(1+nr),1)
        call rwftai(5,rofi(nr),a,nr,nr,rofi,phz,dphz,xx,l,
     .  ehl(k),rsml(k),gzbig)
C       If rwftai scales gzbig, rescale phz,gz
        if (gzbig(nr) .ne. gz(nr,1)) then
          xx = gzbig(nr)/gz(nr,1)
          phz   = phz*xx
          dphz  = dphz*xx
          phzp  = phzp*xx
          dphzp = dphzp*xx
          call dscal(nr,xx,gz(1,1),1)
          call dscal(nr,xx,gz(1,2),1)
        endif
      endif

C     call prrmsh('local orbital',rofi,gz,nr,nr,2)
      end

      subroutine ortrwf(mode,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,ez,g0,g1,gz,
     .D)
C- Orthogonalize a radial wave function gz to a pair of other functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do not change gz, but return scaling factor D that would
Ci         :  normalize gz after orthogonalization to (g0,g1)
Ci         :1 orthonormalize gz.
Ci         :  NB: this routine assumes g0 and g1 are orthogonal
Ci         :2 orthonormalize g0,g1; do not change gz or compute D
Ci         :3 orthogonalize g1 to g0; do not normalize
Ci         :4 orthogonalize g0 and g1 to gz
Ci         :10s digit
Ci         :0 use both large and small components of radial w.f.
Ci         :1 use large component of radial w.f. only.
Ci            In this case, z,v,e0,e1,ez are not used
Ci   z     :nuclear charge
Ci   l     :l quantum number
Ci   v     :spherical potential (atomsr.f)
Ci   ng    :leading dimension of g and gadd
Ci   n1    :if 0<n1<=nr, rwgt(n1) is scaled by 2
Ci         :(see vxtrap.f)
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   e0    :energy eigenvalue of g0
Ci   e1    :energy eigenvalue of g1
Ci   ez    :energy eigenvalue of gz
Ci   g0    :1st radial w.f. to which to orthogonalize gz
Ci   g1    :2nd radial w.f. to which to orthogonalize gz
Ci   gz    :radial w.f. to orthogonalize
Co Outputs
Co   D     :scaling factor that normalizes the orthogonalized gz
Cl Local variables
Cl         :
Cr Remarks
Cr
Cb Bug
Cb   for 1s digit mode=1, this routine assumes g0 and g1 are orthogonal
Cu Updates
Cu   12 Jul 04 Add option 3 to mode.  New argument list
Cu   06 Mar 02 New routine
C ----------------------------------------------------------------------
C     implicit none
      use m_addrwf
C ... Passed parameters
      integer mode,l,ng,n1,nr
      double precision z,v(nr),rofi(nr),rwgt(nr),e0,e1,ez,D
      double precision g0(ng,2),g1(ng,2),gz(ng,2)
C ... Local parameters
      integer mode0,mode1
      double precision s00,s01,s11,s0z,s1z,szz,x,s01hat,s11hat,s1zhat

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

C --- mode 2 : orthonormalize g0 and g1 ---
      if (mode0 .eq. 2) then
C       <g0 g0>
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
        call dscal(nr,1/sqrt(s00),g0(1,1),1)
        call dscal(nr,1/sqrt(s00),g0(1,2),1)
C       <g0 g1>
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
        call daxpy(nr,-s01,g0(1,1),1,g1(1,1),1)
        call daxpy(nr,-s01,g0(1,2),1,g1(1,2),1)
C       <g1 g1>
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
        call dscal(nr,1/sqrt(s11),g1(1,1),1)
        call dscal(nr,1/sqrt(s11),g1(1,2),1)

C       Check
C        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
C        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
C        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
        return

C --- mode 3 : orthogonalize g1 to g0 ---
      elseif (mode0 .eq. 3) then
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
        call daxpy(nr,-s01/s00,g0(1,1),1,g1(1,1),1)
        call daxpy(nr,-s01/s00,g0(1,2),1,g1(1,2),1)
C       Check
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
        return
      endif

      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e0,0d0,g0,g0,s00)
      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,e1,0d0,g1,g1,s11)
      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,e1,0d0,g0,g1,s01)
      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,ez,ez,0d0,gz,gz,szz)
      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,0d0,g0,gz,s0z)
      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,ez,0d0,g1,gz,s1z)

C --- mode 4 : orthogonalize g0 and g1 to gz ---
      if (mode0 .eq. 4) then
        call daxpy(nr,-s0z/szz,gz(1,1),1,g0(1,1),1)
        call daxpy(nr,-s0z/szz,gz(1,2),1,g0(1,2),1)
        call daxpy(nr,-s1z/szz,gz(1,1),1,g1(1,1),1)
        call daxpy(nr,-s1z/szz,gz(1,2),1,g1(1,2),1)

C       Check
C        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,ez,0d0,g0,gz,s0z)
C        call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e1,ez,0d0,g1,gz,s1z)
        return
      endif


      s01hat = s01/sqrt(s00)
      s11hat = s11 - s01hat**2
      s1zhat = s1z - s01*s0z/s00

C     Scaling factor that normalizes the orthogonalized gz
C     D = sqrt(szz - s0z**2/s00 - s1z**2/s11)
      D = sqrt(szz - s0z**2/s00 - s1zhat**2/s11hat)

      if (mode0 .eq. 0) return

C     Orthogonalize
      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,-s0z/s00,g0,gz,x)
      call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,-s1z/s11,g1,gz,x)
C     Normalize
      call dscal(nr,1/D,gz(1,1),1)
      call dscal(nr,1/D,gz(1,2),1)

C     check
C     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,0d0,g0,gz,x)
C     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,e0,ez,0d0,g1,gz,x)
C     call addrwf(mode1,z,l,v,ng,n1,nr,rofi,rwgt,ez,ez,0d0,gz,gz,x)
C     print '('' Normalization of gz is now 1 + '',1pe10.3)', x-1
C     pause

C     check normalization from file
C     mc -qr out.ext -e2 x1 'x2*x2' -av:nr,1 rmax -int 0 rmax
C     call prrmsh('gz',rofi,gz,nr,nr,1)

      end

