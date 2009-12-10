      subroutine getidu(nbas,sspec,ssite,idu)
C- Get idu for all sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idu
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Co Outputs
Co   idu   :idu(l+1,ib)=1 => this l and site has a nonlocal U matrix
Cr Remarks
Cr
Cu Updates
Cu   Mar 29, 07 (Jialei) first created
C------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,idu(4,nbas)
      double precision sspec(1), ssite(1)
C ... Local parameters
      integer igetss,is,ib

      call iinit(idu,4*nbas)
      do  ib = 1, nbas
        is = igetss('site spec',ib,ssite)
        call upack('spec idu',sspec,is,idu(1,ib),0,0,0)
      enddo

      end




