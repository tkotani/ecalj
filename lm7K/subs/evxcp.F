      subroutine evxcp(r1,r2,n,nsp,lxcf,exc,vxcup,vxcdn)
C-Local xc energy and potential for PW91 or PBE
C ----------------------------------------------------------------------
Ci Inputs: r1, r2 charge density in a mesh of n points. If nsp=1, r1
Ci         is total charge; if nsp=2 r1 and r2 are up and
Ci         down charge.
Ci         lxcf=3 for PBE
Co Outputs:
Co         exc, vxcup and vxcdn local xc energy and potential.
Co         If nsp=1, vxcup is vxc
Cr Remarks:
Cr         If nsp=1 r2 points to an arbitrary address. NB evxcp calls
Cr         easypbe written by K. Burke, which uses Hartree atomic units.
Cr         This routine adapted from pbevxc in FP-LMTO and uses names
Cr         longer than six characters.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nsp,n,lxcf
      double precision r1(n),r2(n),exc(n),vxcup(n),vxcdn(n)
C Local Variables
      integer i
      double precision up, dn, exlsd, vxuplsd, vxdnlsd, eclsd,
     .vcuplsd, vcdnlsd, expw91, vxuppw91, vxdnpw91,
     .ecpw91, vcuppw91, vcdnpw91, expbe, vxuppbe,
     .vxdnpbe, ecpbe, vcuppbe, vcdnpbe

      if (lxcf .ne. 3 .and. lxcf .ne. 4) call rx('evxcp: bad lxcf')

      do 1 i = 1, n
        if (nsp .eq. 1) then
          up = r1(i) / 2
          dn = up
        else
          up = r1(i)
          dn = r2(i)
        endif
        call easypbe(up,0d0,0d0,0d0,dn,0d0,0d0,0d0,0d0,0d0,1,1,
     .  exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     .  expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     .  expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
C --- times two to convert Ha to Ry
        exc(i) = (exlsd + eclsd) * 2d0
        vxcup(i) = (vxuplsd + vcuplsd) * 2d0
        if (nsp .eq. 2) vxcdn(i) = (vxdnlsd + vcdnlsd) * 2d0
    1 continue
      end

