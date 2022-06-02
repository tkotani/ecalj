      subroutine pgvmat2 (ng1,gv1, ng2,gv2,kv2 )
      integer:: ng1,ng2,kv2(1)
      real(8):: gv1(1),gv2(1)
      call pgvmat ( ng1,gv1, ng2,gv2,kv2 )
      end
!
      subroutine pgvmat(ngs,gvs,ngb,gvb,kvb)
      implicit none
      integer ngs,ngb
      integer kvb(ngb,3),kk(ngs),iwk(ngs)
      double precision gvs(ngs,3),gvb(ngb,3),gg(ngs)
      integer ig,low,high,jg,mm
      double precision xx,tol,tol2
      parameter (tol=1d-6,tol2=1d-9)
C ... Generate length of g and ensure length matches for both
      do  10  ig = 1, ngs
        gg(ig) = gvs(ig,1)**2 + gvs(ig,2)**2 + gvs(ig,3)**2
   10 continue
C --- For each G in small, find matching G in big list ---
      xx = -1d0
      low = 0
      high = 0
      do  30  ig = 1, ngs
        iwk(ig) = -1
C   ... Find first and last g-vector list with same length
        if (abs(xx-gg(ig)) .gt. tol) then
          call huntx(gg,ngs,gg(ig)+tol,0,high)
          low = ig
          high = min(high,ngs)
          xx = gg(ig)
        endif
        do  32  jg = low, high
          do  34  mm = 1, 3
            if (abs(gvb(jg,mm)-gvs(ig,mm)) .gt. tol2) goto 32
   34     continue
C     ... Found a match
          iwk(ig) = jg
   32   continue
C   ... Sanity check
        if (iwk(ig) .eq. -1) call rxi('bug in gvmatch, ig=',ig)
   30 continue
C ... Rearrange gvb, kvb
      do  40  mm = 1, 3
        do  42  ig = 1, ngs
          jg = iwk(ig)
          gg(ig) = gvb(jg,mm)
          kk(ig) = kvb(jg,mm)
   42   continue
        do  44  ig = 1, ngs
          gvb(ig,mm) = gg(ig)
          kvb(ig,mm) = kk(ig)
   44   continue
   40 continue
      end


