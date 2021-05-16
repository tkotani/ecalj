      subroutine defspc(sspec)
      use m_struc_def,only:s_spec
      use m_lmfinit,only: nspec
      type(s_spec)::sspec(*)
      double precision dgetss,rg,rsma,rfoca,rsmfa,rmt
      integer is,nglob
c      nspec = globalvariables%nspec
      do  is = 1, nspec
        rmt = (sspec(is)%rmt)
        rmt=sspec(is)%rmt
        rg=sspec(is)%rg
        rsma=sspec(is)%rsma
        rfoca=sspec(is)%rfoca
        rsmfa=sspec(is)%rsmfa
        if (rg    .eq. 0) rg    = -1
        if (rg    .lt. 0) rg    = -rg*0.25d0*rmt
        if (rsma  .eq. 0) rsma  = -1
        if (rsma  .lt. 0) rsma  = -rsma*0.4d0*rmt
        if (rfoca .eq. 0) rfoca = -1
        if (rfoca .lt. 0) rfoca = -rfoca*0.4d0*rmt
        if (rsmfa .eq. 0) rsmfa = -1
        if (rsmfa .lt. 0) rsmfa = -rsmfa*0.5d0*rmt
        sspec(is)%rg=rg
        sspec(is)%rsma=rsma
        sspec(is)%rfoca=rfoca
        sspec(is)%rsmfa=rsmfa
        sspec(is)%rsmv=rmt*.5d0
      enddo
      end subroutine defspc


