      subroutine qdiv(tr,ti,dr,di,t1,t2)
C- Complex divide (t1,t2) = (tr,ti) / (dr,di), quad precision
Cr Remarks
Cr   Adapted from eispack.
Cr   It is permissible for (t1,t2) to occupy the same address space
Cr   as either (tr,ti) or (dr,di)
C      double precision tr,ti,dr,di,t1,t2
C      double precision rr,d,tmp
#if NOQUAD
      double precision tr,ti,dr,di,t1,t2
      double precision rr,d,tmp
#else
      real*16 tr,ti,dr,di,t1,t2
      real*16 rr,d,tmp
#endif
      if (abs(di) .gt. abs(dr)) then
        rr = dr / di
        d = dr * rr + di
        tmp = (tr * rr + ti) / d
        t2 = (ti * rr - tr) / d
        t1 = tmp
      else
        rr = di / dr
        d = dr + di * rr
        tmp = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        t1 = tmp
      endif
      end

