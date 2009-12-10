      subroutine dfratm(nbas,ssite,sspec,orhoat)
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
      integer nbas,orhoat(3,nbas)
      double precision ssite(1),sspec(1)
C ... Local parameters
      integer ib,is,igetss,lmxl,nlml,nr,nsp,nglob

      nsp = nglob('nsp')

      do  ib = 1, nbas
        is = igetss('site spec',ib,ssite)
        call upack('spec nr lmxl',sspec,is,nr,lmxl,0,0)
        nlml = (lmxl+1)**2
        if (lmxl .gt. -1) then
        call defrr (orhoat(1,ib),  nr*nlml*nsp)
        call defrr (orhoat(2,ib),  nr*nlml*nsp)
        call defrr (orhoat(3,ib),  nr*nsp)
        endif
      enddo
      end

