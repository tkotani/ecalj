c-taken from iolib of lmto3 --------------------------------------
C
C Routine iprint is intended as a function that returns an integer used
C to control how verbose the output is.  There is a stack of
C levels; iprint returns the top of the stack.  A new value is pushed
C onto the stack by calling pshprt; the previous value is popped by
C calling popprt.  The top two levels can be toggled by calling togprt.
      integer function iprint()
C- get last integer off print priority stack
C     implicit none
      integer nstack
      parameter (nstack=5)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      iprint = vstack(stackp)
      return
      end
      subroutine pshprt(vb)
C     implicit none
      integer vb
      integer nstack
      parameter (nstack=5)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      stackp = mod(stackp+1,nstack)
      vstack(stackp) = vb
      return
      end
      subroutine popprt
C     implicit none
      integer nstack
      parameter (nstack=5)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      stackp = mod(stackp+(nstack-1),nstack)
      return
      end
      subroutine togprt
C     implicit none
      integer nstack
      parameter (nstack=5)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      integer itmp,jtmp
      itmp = vstack(stackp)
      call popprt
      jtmp = vstack(stackp)
      vstack(stackp) = itmp
      call pshprt(jtmp)
      return
      end

C --- Default initialization ---
      block data data_init_iprint
C     implicit none
C for iprint...
      integer nstack
      parameter (nstack=5)
      integer vstack(0:nstack-1),stackp
      common /iprnt/ vstack,stackp
      data vstack /30,30,30,30,30/, stackp /0/
      end

