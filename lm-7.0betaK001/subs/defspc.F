      subroutine defspc(sspec)
C     implicit none
      double precision sspec(1)
      double precision dgetss,rg,rsma,rfoca,rsmfa,rmt
      integer nspec,is,nglob

      nspec = nglob('nspec')

      do  is = 1, nspec

        rmt = dgetss('spec rmt',is,sspec)
        call upack2('spec rmt',sspec,is,rmt)
        call upack('spec rg rsma rfoca rsmfa',sspec,is,rg,rsma,rfoca,
     .    rsmfa)

        if (rg    .eq. 0) rg    = -1
        if (rg    .lt. 0) rg    = -rg*0.25d0*rmt
        if (rsma  .eq. 0) rsma  = -1
        if (rsma  .lt. 0) rsma  = -rsma*0.4d0*rmt
        if (rfoca .eq. 0) rfoca = -1
        if (rfoca .lt. 0) rfoca = -rfoca*0.4d0*rmt
        if (rsmfa .eq. 0) rsmfa = -1
        if (rsmfa .lt. 0) rsmfa = -rsmfa*0.5d0*rmt

        call pack5('spec rg rsma rfoca rsmfa',sspec,is,rg,rsma,rfoca,
     .    rsmfa)
        call pack2('spec rsmv',sspec,is,rmt*.5d0)

C        print 333, is,rg,rsma,rfoca,rsmfa
C  333   format(i3,4f12.6)
      enddo
C     stop
      end

