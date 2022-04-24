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
      
c$$$C --- pshprt, pshpr, setpr, getpr, popprt, poppr --
c$$$      subroutine pshprt(vb)
c$$$C     implicit none
c$$$      entry pshpr(vb)
c$$$      integer vb,istk,j
c$$$      integer nstack
c$$$      parameter (nstack=7)
c$$$      integer vstack(0:nstack-1),stackp
c$$$      common /iprnt/ vstack,stackp
c$$$      stackp = mod(stackp+(nstack-1),nstack)
c$$$
c$$$      entry setpr(vb)
c$$$      vstack(stackp) = vb
c$$$      return
c$$$
c$$$      entry getpr(vb)
c$$$      vb = vstack(stackp)
c$$$      return
c$$$
c$$$      entry popprt
c$$$      entry poppr
c$$$      stackp = mod(stackp+1,nstack)
c$$$      return
c$$$
c$$$      entry sprt(istk,vb)
c$$$      vstack(mod(stackp+istk,nstack)) = vb
c$$$      return
c$$$
c$$$      entry shopr
c$$$      print 333, (vstack(mod(stackp+j,nstack)), j=0,nstack-1)
c$$$  333 format(' vb stack',10i4)
c$$$      end

c$$$C --- togprt ---
c$$$      subroutine togprt
c$$$C     implicit none
c$$$      integer nstack
c$$$      parameter (nstack=7)
c$$$      integer vstack(0:nstack-1),stackp
c$$$      common /iprnt/ vstack,stackp
c$$$      integer itmp,jtmp
c$$$
c$$$      entry togpr
c$$$      itmp = vstack(stackp)
c$$$      call popprt
c$$$      jtmp = vstack(stackp)
c$$$      call popprt
c$$$      call pshprt(itmp)
c$$$      call pshprt(jtmp)
c$$$      end
c$$$C --- Defaults ---
c$$$      block data iprntb
c$$$C     implicit none
c$$$      integer nstack
c$$$      parameter (nstack=7)
c$$$      integer vstack(0:nstack-1),stackp
c$$$      common /iprnt/ vstack,stackp
c$$$      data vstack /30,30,30,30,30,30,30/, stackp /0/
c$$$      end
c$$$
