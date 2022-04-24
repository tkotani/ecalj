      integer function iprint()
      implicit none
      integer nstack,istk,iprt,mpipid
      integer:: verbose_in,setprint,ix,set0,setprint0,vb
      integer,save:: verbose0,verbose
      include "mpif.h"
      iprint = verbose
      if(mpipid(1)>0) iprint=0 !write only at master node
      return
      entry setprint0(verbose_in)!base
      verbose0=verbose_in
      verbose=verbose0
      return
      entry set0()
      verbose=verbose0
      return
      entry setprint(ix)
      verbose=ix
      return
      end
      
      subroutine pshpr(vb)
      integer:: setprint,vb,i
      i=setprint(vb)
      end
      
      subroutine poppr()
      integer:: set0,i
      i=set0()
      end

      subroutine getpr(ix)
      integer:: iprint,ix
      ix=iprint()
      end
