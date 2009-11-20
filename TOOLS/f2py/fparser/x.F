      subroutine bndfp(alpha)
      call xxxbfp(alpha,beta, aaa) ! Patch to avoid PGI compiler bug on AMD processor
      end

