Cdelw1 undel [['call', 'array', 'wref'], 'orhoat', 'dfratm']
Cgetarg...info...           structure ['dfratm', 'site', 'ssite']
Cgetarg...info...           structure ['dfratm', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['dfratm', 'm_struc_def']
      subroutine dfratm(nbas,ssite,sspec,orhoat)
       
       use m_struc_def  !Cgetarg

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
Cdelw1 do not change  [['call', 'array', 'wref'], 'orhoat', 'dfratm']
      integer nbas,orhoat(3,nbas)
Cgetarg       double precision ssite(1),sspec(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

C ... Local parameters
      integer ib,is,igetss,lmxl,nlml,nr,nsp,nglob

      nsp = nglob('nsp')

      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg         call upack('spec nr lmxl',sspec,is,nr,lmxl,0,0)
         
         nr=sspec(is)%nr
         lmxl=sspec(is)%lmxl

        nlml = (lmxl+1)**2
        if (lmxl .gt. -1) then
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
        call defrr (orhoat(1,ib),  nr*nlml*nsp)
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
        call defrr (orhoat(2,ib),  nr*nlml*nsp)
Cdelw1 do not change orhoat because of ['call', 'array', 'wref']
        call defrr (orhoat(3,ib),  nr*nsp)
        endif
      enddo
Cgetarg       end
Cdelw1 w_varlist remains: orhoat 
Cdelw1 w_varlistundel: [orhoat]
Cdelw1 w_varlist (undel), remains: [orhoat]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate orhoat because of [call array wref]
       end subroutine dfratm 


