      double precision function stkpsh(stack,nstack,stackp,x)
C- pushes variable onto stack
Ci   stack,nstack,stackp:  stack, stack size and current stack pointer
Ci   (stackp should be initially zero)
      integer nstack,stackp
      double precision stack(0:nstack-1),x
      stackp = mod(stackp+1,nstack)
      stack(stackp) = x
      stkpsh = x
      end
      double precision function stkpop(stack,nstack,stackp)
C- Pops variable off top of stack
      integer nstack,stackp
      double precision stack(0:nstack-1)
      stkpop = stack(stackp)
      stackp = mod(stackp+(nstack-1),nstack)
      end

