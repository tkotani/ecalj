C --- Verbosity stack functions ---
C Functions are:
C iprint(): return value from top of stack
C pshprt(vb): push value onto top of stack
C popprt:     pop verbosity stack
C togprt:     exhanges position of top two values on stack
C setpr(vb):  replace top of stack with vb
C getpr(vb):  set vb to value at top of stack
C iprt(istk): return verbosity istk levels down
C sprt(istk,vb): sets verbosity istk levels down
C --- iprint, iprt ---
      integer function iprint()
C- get last integer off print priority stack
Cu Updates
Cu   28 Aug 04 Increase stack size to 7
C     implicit none
      integer nstack,istk,iprt,mpipid
      parameter (nstack=7)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
#if MPI | MPIK
      include "mpif.h"
#endif
      iprint = vstack(stackp)
      if(mpipid(1)>0) iprint=0 !write only at master node
      return

      entry iprt(istk)
      iprt = vstack(mod(stackp+istk,nstack))
c      print*, istk,iprt,stackp,vstack
      end
C --- pshprt, pshpr, setpr, getpr, popprt, poppr --
      subroutine pshprt(vb)
C     implicit none
      entry pshpr(vb)
      integer vb,istk,j
      integer nstack
      parameter (nstack=7)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      stackp = mod(stackp+(nstack-1),nstack)

      entry setpr(vb)
      vstack(stackp) = vb
      return

      entry getpr(vb)
      vb = vstack(stackp)
      return

      entry popprt
      entry poppr
      stackp = mod(stackp+1,nstack)
      return

      entry sprt(istk,vb)
      vstack(mod(stackp+istk,nstack)) = vb
      return

      entry shopr
      print 333, (vstack(mod(stackp+j,nstack)), j=0,nstack-1)
  333 format(' vb stack',10i4)
      end

C --- togprt ---
      subroutine togprt
C     implicit none
      integer nstack
      parameter (nstack=7)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      integer itmp,jtmp

      entry togpr
      itmp = vstack(stackp)
      call popprt
      jtmp = vstack(stackp)
      call popprt
      call pshprt(itmp)
      call pshprt(jtmp)
      end
C --- Defaults ---
      block data iprntb
C     implicit none
      integer nstack
      parameter (nstack=7)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      data vstack /30,30,30,30,30,30,30/, stackp /0/
      end

