      subroutine praugm(sspec,is)
      use m_struc_def  !Cgetarg
      use m_lmfinit,only: nspec,stdo
C-  Print species information
C ----------------------------------------------------------------------
Ci Inputs
Ci   sspec :struct containing species-specific information
Ci   is    :species index (use 0 to print info for all species)
Co Outputs
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer is
      type(s_spec)::sspec(*)
      integer is1,is2,js,kmxt,lmxa,lgunit
      integer kmxv,lmxl,lfoca
      double precision rmt,rsma,rfoca,rg,rsmv
      character spid*8
      is1 = is
      is2 = is
      if (is .le. 0) then
        is1 = 1
        is2 = nspec !globalvariables%nspec
      endif
      write (stdo,501)
  501 format(/' species data:  augmentation',27x,'density'/
     .' spec       rmt   rsma lmxa kmxa',5x,' lmxl     rg   rsmv  kmxv foca   rfoca')
      do js = is1,is2
        spid=sspec(js)%name
        rmt=sspec(js)%rmt
        rsma=sspec(js)%rsma
        lmxa=sspec(js)%lmxa
        kmxt=sspec(js)%kmxt
        lmxl=sspec(js)%lmxl
        rg=sspec(js)%rg
        rsmv=sspec(js)%rsmv
        kmxv=sspec(js)%kmxv
        lfoca=sspec(js)%lfoca
        rfoca=sspec(js)%rfoca
        write (stdo,500) spid,rmt,rsma,lmxa,kmxt, lmxl,rg,rsmv,kmxv,lfoca,rfoca
 500    format(1x,a,f6.3,f7.3,2i5,6x,i4,2f7.3,i6,i5,f8.3)
      enddo
      end subroutine praugm



