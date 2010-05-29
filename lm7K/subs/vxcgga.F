      subroutine vxcgga(lxcg,n,nsp,rhop,rhom,grhop,grhom,ggrhop,ggrhom,
     .grho,grpgrm,grggr,grggrp,grggrm,vxc1,vxc2,exc)
C- PW91 and PBE gradient corrections to exc and vxc for a set of points
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rhop  :spin up density if nsp=2; otherwise total density
Ci   rhom  :spin down density if nsp=2; otherwise not used
Ci   grhop :|grad rhop| or |grad rho| if nsp=1
Ci   grhom :|grad rhom| (nsp=2)
Ci   ggrhop:Laplacian of rhop
Ci   ggrhom:Laplacian of rhom
Ci   grho  :|grad total rho| if nsp=2; otherwise not used
Ci   grpgrm:grad rho+ . grad rho- (not used here)
Ci   grggr :(grad rho).(grad |grad rho|)
Ci   grggrp:(grad up).(grad |grad up|) if nsp=2; otherwise ggrgr
Ci   grggrm:(grad dn).(grad |grad dn|) if nsp=2; otherwise not used
Co Outputs
Co   vxc1  :GGA contr. to vxc+ if nsp=2; otherwise GGA contr. to vxc
Co   vxc2  :GGA contr. to vxc- if nsp=2; otherwise GGA contr. to vxc
Co    exc  :GGA contr. to exc
Cl Local variables
Cl   symbols match easypbe
Cb Bugs
Cr Remarks
Cu Updates
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nsp,n,lxcg
      double precision rhop(1),rhom(1),grhop(1),grhom(1),grho(1),
     .ggrhop(1),ggrhom(1),grggr(1),grggrp(1),grggrm(1),grpgrm(1),
     .vxc1(1),vxc2(1),exc(1)
C Local Variables
      integer i,isw
      double precision up,dn,agrup,delgrup,uplap,agrdn,delgrdn,dnlap,
     .agr,delgr
      double precision exlsd,vxuplsd,vxdnlsd,eclsd,
     .vcuplsd,vcdnlsd,expw91,vxuppw91,vxdnpw91,
     .ecpw91,vcuppw91,vcdnpw91,expbe,vxuppbe,
     .vxdnpbe,ecpbe,vcuppbe,vcdnpbe

      if (lxcg .ne. 3 .and. lxcg .ne. 4) call rx('evxcp: bad lxcg')

      do  i = 1, n
        if (nsp .eq. 1) then
          up = rhop(i) / 2
          dn = up
          agrup = grhop(i) / 2
          agrdn = agrup
          uplap = ggrhop(i) / 2
          dnlap = uplap
          delgrup = grggrp(i) / 4
          delgrdn = delgrup
          agr = grhop(i)
          delgr = grggrp(i)
        else
          up = rhop(i)
          dn = rhom(i)
          agrup = grhop(i)
          agrdn = grhom(i)
          uplap = ggrhop(i)
          dnlap = ggrhom(i)
          delgrup = grggrp(i)
          delgrdn = grggrm(i)
          agr = grho(i)
          delgr = grggr(i)
        endif
        call easypbe(up,agrup,delgrup,uplap,
     .  dn,agrdn,delgrdn,dnlap,agr,delgr,1,1,
     .  isw(lxcg.eq.4),
     .  exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd,
     .  expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91,
     .  expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
        if (lxcg .eq. 2) then
          exc(i) = exc(i) + (expw91 + ecpw91 - exlsd - eclsd) * 2d0
          vxc1(i) = vxc1(i) +
     .    (vxuppw91 + vcuppw91 - vxuplsd - vcuplsd) * 2d0
          if (nsp .eq. 2) vxc2(i) = vxc2(i) +
     .    (vxdnpw91 + vcdnpw91 - vxdnlsd - vcdnlsd) * 2d0
        else
          exc(i) = exc(i) + (expbe + ecpbe - exlsd - eclsd) * 2d0
          vxc1(i) = vxc1(i) +
     .    (vxuppbe + vcuppbe - vxuplsd - vcuplsd) * 2d0
          if (nsp .eq. 2) vxc2(i) = vxc2(i) +
     .    (vxdnpbe + vcdnpbe - vxdnlsd - vcdnlsd) * 2d0
        endif
      enddo
      end


