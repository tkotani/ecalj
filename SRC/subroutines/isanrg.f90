      logical function isanrg(i,i1,i2,t1,t2,lreqd)
      logical lreqd
      integer i,i1,i2,lgunit,iprint,k1,k2,it1
      character*(*) t1,t2
      character strn*80,strn2*80,t3*30
      if (i>=i1 .and. i<=i2) return
      call rx(trim(t1)//' '//trim(t2))
      end

      subroutine fsanrg(f,f1,f2,tol,t1,t2,lreqd)
      logical lreqd
      double precision f,f1,f2,tol
      character*(*) t1,t2
      character strn*100,strn2*100,t3*30
      if (f>=f1 .and. f<=f2) return
      if (f1==f2 .and. f>=f1-tol/2d0 .and. f<=f2+tol/2d0) return
      call rx(trim(t1)//' '//trim(t2))
      end

