Cgetarg...info...           structure ['defspc', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['defspc', 'm_struc_def']
      subroutine defspc(sspec)
       
       use m_struc_def  !Cgetarg

C     implicit none
Cgetarg       double precision sspec(1)
       type(s_spec)::sspec(*)

      double precision dgetss,rg,rsma,rfoca,rsmfa,rmt
      integer nspec,is,nglob

      nspec = nglob('nspec')

      do  is = 1, nspec

Cgetarg         rmt = dgetss('spec rmt',is,sspec)
         rmt = (sspec(is)%rmt) 

Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack2('spec rmt',sspec,is,rmt)
         
         rmt=sspec(is)%rmt

Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg...info...           real(8) :: rsma , real(8) :: sspec%rsma 1
Cgetarg...info...           real(8) :: rfoca , real(8) :: sspec%rfoca 1
Cgetarg...info...           real(8) :: rsmfa , real(8) :: sspec%rsmfa 1
Cgetarg         call upack('spec rg rsma rfoca rsmfa',sspec,is,rg,rsma,rfoca,
Cgetarg      .    rsmfa)
         
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

Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg...info...           real(8) :: rsma , real(8) :: sspec%rsma 1
Cgetarg...info...           real(8) :: rfoca , real(8) :: sspec%rfoca 1
Cgetarg...info...           real(8) :: rsmfa , real(8) :: sspec%rsmfa 1
Cgetarg         call pack5('spec rg rsma rfoca rsmfa',sspec,is,rg,rsma,rfoca,
Cgetarg      .    rsmfa)
         
         sspec(is)%rg=rg 
         sspec(is)%rsma=rsma 
         sspec(is)%rfoca=rfoca 
         sspec(is)%rsmfa=rsmfa 

Cgetarg...info...           NG:: rmt*.5d0 , real(8) :: sspec%rsmv 1
Cgetarg         call pack2('spec rsmv',sspec,is,rmt*.5d0)
         
         sspec(is)%rsmv=rmt*.5d0 


C        print 333, is,rg,rsma,rfoca,rsmfa
C  333   format(i3,4f12.6)
      enddo
C     stop
Cgetarg       end
       end subroutine defspc 


