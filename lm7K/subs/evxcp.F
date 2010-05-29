      subroutine evxcp(r1,r2,n,nsp,lxcf,ex,ec,exc,vxup,vxdn,vcup,vcdn,
     .vxcup,vxcdn)
C-Local part xc energy and potential for PW91 or PBE (in LDA, PBE=PW91)
C ----------------------------------------------------------------------
Ci Inputs: r1, r2 charge density in a mesh of n points. If nsp=1, r1
Ci         is total charge; if nsp=2 r1 and r2 are up and 
Ci         down charge.
Ci         lxcf=3,4 for PBE (AKA PW91)
Co Outputs:
Co   ex    :exchange-only part of exc
Co   ec    :correlation-only part of exc
Co   exc   :local exchange energy density for the n points
Co   vxup  :exchange-only part of vxc, spin1
Co   vcup  :correlation-only part of vxc, spin1
Co   vxcup :local exchange potential for the n points, spin1
Co   vxdn  :exchange-only part of vxc, spin2
Co   vcdn  :correlation-only part of vxc, spin2
Co   vxcdn :local exchange potential for the n points, spin2
Co         If nsp=1, only spin numbers are returned
Cr Remarks: 
Cr         If nsp=1 r2 points to an arbitrary address. NB evxcp calls
Cr         easypbe written by K. Burke, which uses Hartree atomic units.
Cr         This routine adapted from pbevxc in FP-LMTO and uses names
Cr         longer than six characters.
Cu Updates
Cu   21 Apr 09 Returns exchange and correlation parts separately
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nsp,n,lxcf
      double precision r1(n),r2(n),ex(n),ec(n),exc(n)
      double precision vxup(n),vxdn(n),vcup(n),vcdn(n),vxcup(n),vxcdn(n)
C Local Variables
      integer i,iprint
      double precision up, dn, exlsd, vxuplsd, vxdnlsd, eclsd,
     .vcuplsd, vcdnlsd, expw91, vxuppw91, vxdnpw91,
     .ecpw91, vcuppw91, vcdnpw91, expbe, vxuppbe,
     .vxdnpbe, ecpbe, vcuppbe, vcdnpbe

      if (lxcf .ne. 3 .and. lxcf .ne. 4) then
        if (iprint() .lt. 10) call pshpr(10)
        call rxi('evxcp cannot handle lxcf =',lxcf)
      endif

      do  i = 1, n
        if (nsp .eq. 1) then
          up = r1(i) / 2
          dn = up
        else
          up = r1(i)
          dn = r2(i)
        endif
        call easypbe(up,0d0,0d0,0d0,dn,0d0,0d0,0d0,0d0,0d0,1,1,0,
     .  exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     .  expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     .  expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe) 
C   ... Times two to convert Ha to Ry
        ex(i) = exlsd * 2d0
        ec(i) = eclsd * 2d0
        exc(i) = (exlsd + eclsd) * 2d0
        vxup(i)  = vxuplsd * 2d0
        vcup(i)  = vcuplsd * 2d0
        vxcup(i) = (vxuplsd + vcuplsd) * 2d0
        if (nsp .eq. 2) then
          vxdn(i)  = vxdnlsd * 2d0
          vcdn(i)  = vcdnlsd * 2d0
          vxcdn(i) = (vxdnlsd + vcdnlsd) * 2d0
        endif
      enddo
      end

