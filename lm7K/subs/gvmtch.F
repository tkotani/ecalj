cki      subroutine gvmtch(ngs,ogvs,okvs,ngb,ogvb,okvb)
ckiC-  Ensure lower block of big list matches lower block of small list
ckiC     implicit none
cki      integer ngs,ngb,ogvs,okvs,ogvb,okvb
cki      real(8) ,allocatable :: gg_rv(:)
cki      real(8) ,allocatable :: kk_rv(:)
cki      real(8) ,allocatable :: iwk_rv(:)
cki
cki      integer w(1)
cki      common /w/ w
cki
cki      allocate(gg_rv(ngs))
ckiCtakao_ZeroClear_NotRequiered        if (ngs<0) gg_rv(:)=0.0d0
cki
cki      allocate(kk_rv(ngs))
ckiCtakao_ZeroClear_NotRequiered        if (ngs<0) kk_rv(:)=0.0d0
cki
cki      allocate(iwk_rv(ngs))
ckiCtakao_ZeroClear_NotRequiered        if (ngs<0) iwk_rv(:)=0.0d0
cki
cki      call pgvmat ( gg_rv , kk_rv , iwk_rv , ngs , w ( ogvs 
cki     .) , ngb , w ( ogvb ) , w ( okvb ) )
cki
cki      if (allocated(iwk_rv)) deallocate(iwk_rv)
cki      if (allocated(kk_rv)) deallocate(kk_rv)
cki      if (allocated(gg_rv)) deallocate(gg_rv)
cki
cki
cki      end
      subroutine pgvmat2 (ng1,gv1, ng2,gv2,kv2 )
      integer:: ng1,ng2,kv2(1)
      real(8):: gv1(1),gv2(1)
      real(8) ,allocatable :: gg_rv(:)
      real(8) ,allocatable :: kk_rv(:)
      real(8) ,allocatable :: iwk_rv(:)
      allocate( gg_rv(ng1),kk_rv(ng1),iwk_rv(ng1) )
      call pgvmat ( gg_rv , kk_rv , iwk_rv , ng1,gv1, ng2,gv2,kv2 )
      deallocate(iwk_rv,kk_rv,gg_rv)
      end

      subroutine pgvmat(gg,kk,iwk,ngs,gvs,ngb,gvb,kvb)
C     implicit none
      integer ngs,ngb
      integer kvb(ngb,3),kk(ngs),iwk(ngs)
      double precision gvs(ngs,3),gvb(ngb,3),gg(ngs)
      integer ig,low,high,jg,mm
      double precision xx,tol,tol2
      parameter (tol=1d-6,tol2=1d-9)

C      call prm('gvs',gvs,ngs,ngs,3)
C      call prm('gvb',gvb,ngb,ngb,3)

C ... Generate length of g and ensure length matches for both
      do  10  ig = 1, ngs
        gg(ig) = gvs(ig,1)**2 + gvs(ig,2)**2 + gvs(ig,3)**2
   10 continue

C     call prm('gg',gg,ngs,ngs,1)


C ... Sanity check.  Make sure |G| is same for each vector in both lists
C      do  20  ig = 1, ngs
C        xx = gvb(ig,1)**2 + gvb(ig,2)**2 + gvb(ig,3)**2
C        if (abs(xx-gg(ig)) .gt. tol) then
C          print *, xx,gg(ig)
C          call rxi('bug in gvmatch, ig=',ig)
C        endif
C   20 continue

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

c      call prm('gvb',gvb,ngb,ngb,3)

      end


