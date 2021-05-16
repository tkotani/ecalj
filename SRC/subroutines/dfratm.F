      subroutine dfratm ( nbas , ssite , sspec , sv_p_orhoat )
      use m_struc_def  !Cgetarg
      use m_lmfinit,only: nsp
C- Allocate arrays for local atomic densities.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Co Outputs
Ci   orhoat:vector of offsets containing site density allocated
Cr Remarks
Cr   rhoat(1,ib):  true local density
Cr   rhoat(2,ib):  smooth local density
Cr   rhoat(3,ib):  core density
Cu Updates
Cu   18 Jun 98 adapted from nfp df_rhoat.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer:: nbas
      type(s_rv1) :: sv_p_orhoat(3,nbas)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
C ... Local parameters
      integer ib,is,igetss,lmxl,nlml,nr
c      nsp = globalvariables%nsp
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)
        nr=sspec(is)%nr
        lmxl=sspec(is)%lmxl
        nlml = (lmxl+1)**2
        if (lmxl .gt. -1) then
          allocate(sv_p_orhoat(1,ib)%v(abs(nr*nlml*nsp)))
          allocate(sv_p_orhoat(2,ib)%v(abs(nr*nlml*nsp)))
          allocate(sv_p_orhoat(3,ib)%v(abs(nr*nsp)))
        endif
      enddo
      end subroutine dfratm


